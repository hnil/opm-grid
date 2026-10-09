/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"

#define BOOST_TEST_MODULE CoarsenInverseOfLgrTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>
#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

struct Fixture
{
    Fixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};
BOOST_GLOBAL_FIXTURE(Fixture);

namespace
{
using Opm::Coarsening::CoarsenRequest;
using Opm::Coarsening::Grdecl;

Grdecl uniformGrid(int n, double dx, double dz)
{
    Grdecl g;
    g.dims = {n, n, n};
    for (int j = 0; j <= n; ++j) {
        for (int i = 0; i <= n; ++i) {
            g.coord.insert(g.coord.end(), {i*dx, j*dx, 0.0, i*dx, j*dx, n*dz});
        }
    }
    g.zcorn.resize(8ull*n*n*n);
    for (int k = 0; k < 2*n; ++k) {
        for (int j = 0; j < 2*n; ++j) {
            for (int i = 0; i < 2*n; ++i) {
                g.zcorn[i + 2ull*n*(j + 2ull*n*k)] = ((k + 1)/2)*dz;
            }
        }
    }
    g.actnum.assign(1ull*n*n*n, 1);
    return g;
}

grdecl view(const Grdecl& g)
{
    grdecl in{};
    for (int d = 0; d < 3; ++d) {
        in.dims[d] = g.dims[d];
    }
    in.coord = g.coord.data();
    in.zcorn = g.zcorn.data();
    in.actnum = g.actnum.data();
    return in;
}

using Key = std::array<long,3>;
Key key(const Dune::FieldVector<double,3>& x)
{
    return {std::lround(x[0]*1e3), std::lround(x[1]*1e3), std::lround(x[2]*1e3)};
}

// Per cell (by centroid): volume, number of faces, number of distinct nodes.
struct Shape
{
    std::map<Key, std::tuple<double,int,int>> cells;
    std::set<Key> vertices;
};

Shape shape(const Dune::CpGrid& g)
{
    Shape s;
    const auto& gv = g.leafGridView();
    for (const auto& v : vertices(gv)) {
        s.vertices.insert(key(v.geometry().center()));
    }
    for (int c = 0; c < gv.size(0); ++c) {
        std::set<int> nodes;
        for (int f = 0; f < g.numCellFaces(c); ++f) {
            const int face = g.cellFace(c, f);
            for (int v = 0; v < g.numFaceVertices(face); ++v) {
                nodes.insert(g.faceVertex(face, v));
            }
        }
        s.cells[key(g.cellCentroid(c))] = {g.cellVolume(c), g.numCellFaces(c),
                                           static_cast<int>(nodes.size())};
    }
    return s;
}

void checkSame(const Shape& a, const Shape& b)
{
    BOOST_REQUIRE_EQUAL(a.cells.size(), b.cells.size());
    BOOST_CHECK(a.vertices == b.vertices);
    for (const auto& [k, va] : a.cells) {
        const auto it = b.cells.find(k);
        BOOST_REQUIRE(it != b.cells.end());
        BOOST_CHECK_CLOSE(std::get<0>(va), std::get<0>(it->second), 1e-10);
        BOOST_CHECK_EQUAL(std::get<1>(va), std::get<1>(it->second));
        BOOST_CHECK_EQUAL(std::get<2>(va), std::get<2>(it->second));
    }
}

Dune::CpGrid merged(const Grdecl& fine, const std::vector<CoarsenRequest>& requests)
{
    const auto layout = Opm::Coarsening::blockLayout(fine.dims, requests);
    Dune::CpGrid grid;
    grid.processEclipseFormatCoarsened(view(fine), layout.blockOfCartesian, layout.boxes,
                                       /*edge_conformal*/ true, /*collapse_coarse_faces*/ true);
    return grid;
}

// The coarse grid with the given boxes refined, edge-conformal.
Dune::CpGrid refined(const Grdecl& fine, int factor,
                     const std::vector<std::pair<std::array<int,3>, std::array<int,3>>>& boxes)
{
    CoarsenRequest all;
    all.endIJK = fine.dims;
    all.cellsPerDim = {fine.dims[0]/factor, fine.dims[1]/factor, fine.dims[2]/factor};
    const auto coarse = Opm::Coarsening::coarsenCornerPoint(fine, {all}).grid;
    Dune::CpGrid grid;
    grid.processEclipseFormat(view(coarse), false, false, true);
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
        coarse.dims, coarse.coord, coarse.zcorn, coarse.actnum, /*edgeConformal*/ true));
    std::vector<std::array<int,3>> perDim, lo, hi;
    std::vector<std::string> names;
    for (const auto& [l, h] : boxes) {
        perDim.push_back({factor, factor, factor});
        lo.push_back(l);
        hi.push_back(h);
        names.push_back("LGR" + std::to_string(names.size() + 1));
    }
    grid.addLgrsUpdateLeafView(perDim, lo, hi, names);
    return grid;
}
} // namespace

// 6^3 fine cells, factor 2: the coarse grid is 3^3 and the box its centre cell.
BOOST_AUTO_TEST_CASE(CoarsenedOutsideIsTheCoarseGridWithTheBoxRefined)
{
    const auto fine = uniformGrid(6, 100.0, 10.0);
    std::vector<CoarsenRequest> requests;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                if (i == 1 && j == 1 && k == 1) {
                    continue;
                }
                CoarsenRequest r;
                r.startIJK = {2*i, 2*j, 2*k};
                r.endIJK = {2*i + 2, 2*j + 2, 2*k + 2};
                requests.push_back(r);
            }
        }
    }
    checkSame(shape(merged(fine, requests)), shape(refined(fine, 2, {{{1,1,1}, {2,2,2}}})));
}

BOOST_AUTO_TEST_CASE(CoarsenedBoxIsTheCoarseGridRefinedAroundIt)
{
    const auto fine = uniformGrid(6, 100.0, 10.0);
    CoarsenRequest box;
    box.startIJK = {2, 2, 2};
    box.endIJK = {4, 4, 4};
    const std::vector<std::pair<std::array<int,3>, std::array<int,3>>> slabs{
        {{0,0,0}, {3,3,1}}, {{0,0,2}, {3,3,3}},
        {{0,0,1}, {3,1,2}}, {{0,2,1}, {3,3,2}},
        {{0,1,1}, {1,2,2}}, {{2,1,1}, {3,2,2}}};
    checkSame(shape(merged(fine, {box})), shape(refined(fine, 2, slabs)));
}

// Neighbouring boxes with the same factor give the grid of the box they make up.
BOOST_AUTO_TEST_CASE(AdjacentBoxesAreOneBox)
{
    const auto fine = uniformGrid(8, 100.0, 10.0);
    CoarsenRequest whole;
    whole.startIJK = {0, 0, 0};
    whole.endIJK = {8, 8, 4};
    whole.cellsPerDim = {4, 4, 2};
    std::vector<CoarsenRequest> quarters;
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
            CoarsenRequest r;
            r.startIJK = {4*i, 4*j, 0};
            r.endIJK = {4*i + 4, 4*j + 4, 4};
            r.cellsPerDim = {2, 2, 2};
            quarters.push_back(r);
        }
    }
    checkSame(shape(merged(fine, {whole})), shape(merged(fine, quarters)));
}

// Where the coarsening is a corner-point grid, the merge gives that grid, also
// where neighbouring columns group their layers differently.
BOOST_AUTO_TEST_CASE(CollapseIsTheCornerPointCoarsening)
{
    const auto fine = uniformGrid(6, 100.0, 10.0);
    CoarsenRequest south, north;
    south.startIJK = {0, 0, 0};
    south.endIJK = {6, 3, 6};
    south.cellsPerDim = {3, 3, 2};
    north.startIJK = {0, 3, 0};
    north.endIJK = {6, 6, 6};
    north.cellsPerDim = {3, 3, 6};
    const auto coarse = Opm::Coarsening::coarsenCornerPoint(fine, {south, north}).grid;
    Dune::CpGrid grid;
    grid.processEclipseFormat(view(coarse), false, false, true);
    checkSame(shape(merged(fine, {south, north})), shape(grid));
}
