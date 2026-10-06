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

// Layers that dip and bend: z = k*dz + a dome in x and y.
Grdecl curvedGrid(int n, double dx, double dz)
{
    auto g = uniformGrid(n, dx, dz);
    const double L = n*dx;
    for (int k = 0; k < 2*n; ++k) {
        for (int j = 0; j < 2*n; ++j) {
            for (int i = 0; i < 2*n; ++i) {
                const double x = ((i + 1)/2)*dx, y = ((j + 1)/2)*dx;
                const double bump = 15.0*std::sin(M_PI*x/L)*std::sin(M_PI*y/L) + 0.02*x;
                g.zcorn[i + 2ull*n*(j + 2ull*n*k)] += bump*(1.0 + 0.1*((k + 1)/2));
            }
        }
    }
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

using Point = Dune::FieldVector<double,3>;
struct CellShape { Point centroid; double volume; int faces; int nodes; };
struct FaceShape { Point centroid; Point normal; double area; int nodes; };

struct Shape
{
    std::map<Key, CellShape> cells;
    std::map<Key, FaceShape> faces;     // unique by centroid
    std::map<Key, Point> vertices;
};

Shape shape(const Dune::CpGrid& g)
{
    Shape s;
    const auto& gv = g.leafGridView();
    for (const auto& v : vertices(gv)) {
        s.vertices[key(v.geometry().center())] = v.geometry().center();
    }
    for (int c = 0; c < gv.size(0); ++c) {
        std::set<int> nodes;
        for (int f = 0; f < g.numCellFaces(c); ++f) {
            const int face = g.cellFace(c, f);
            for (int v = 0; v < g.numFaceVertices(face); ++v) {
                nodes.insert(g.faceVertex(face, v));
            }
            // orient the normal by the lower centroid so both grids agree
            auto n = g.faceNormal(face);
            if (n[0] + 1e-3*n[1] + 1e-6*n[2] < 0) {
                n *= -1.0;
            }
            s.faces[key(g.faceCentroid(face))] = {g.faceCentroid(face), n, g.faceArea(face),
                                                  g.numFaceVertices(face)};
        }
        s.cells[key(g.cellCentroid(c))] = {g.cellCentroid(c), g.cellVolume(c),
                                           g.numCellFaces(c), static_cast<int>(nodes.size())};
    }
    return s;
}

void checkClose(const Point& a, const Point& b, double tol)
{
    for (int d = 0; d < 3; ++d) {
        BOOST_CHECK_SMALL(a[d] - b[d], tol);
    }
}

// extraEdgeNodes: a may list nodes on its face edges that b leaves hanging.
void checkSame(const Shape& a, const Shape& b, const bool extraEdgeNodes = false)
{
    BOOST_REQUIRE_EQUAL(a.cells.size(), b.cells.size());
    BOOST_REQUIRE_EQUAL(a.faces.size(), b.faces.size());
    BOOST_REQUIRE_EQUAL(a.vertices.size(), b.vertices.size());
    for (const auto& [k, p] : a.vertices) {
        const auto it = b.vertices.find(k);
        BOOST_REQUIRE(it != b.vertices.end());
        checkClose(p, it->second, 1e-9);
    }
    for (const auto& [k, ca] : a.cells) {
        const auto it = b.cells.find(k);
        BOOST_REQUIRE(it != b.cells.end());
        const auto& cb = it->second;
        BOOST_CHECK_CLOSE(ca.volume, cb.volume, 1e-10);
        checkClose(ca.centroid, cb.centroid, 1e-9);
        BOOST_CHECK_EQUAL(ca.faces, cb.faces);
        if (extraEdgeNodes) {
            BOOST_CHECK_GE(ca.nodes, cb.nodes);  // b may leave a node hanging on every face
        } else {
            BOOST_CHECK_EQUAL(ca.nodes, cb.nodes);
        }
    }
    for (const auto& [k, fa] : a.faces) {
        const auto it = b.faces.find(k);
        BOOST_REQUIRE(it != b.faces.end());
        const auto& fb = it->second;
        BOOST_CHECK_CLOSE(fa.area, fb.area, 1e-10);
        checkClose(fa.centroid, fb.centroid, 1e-9);
        checkClose(fa.normal, fb.normal, 1e-9*std::max(1.0, fa.area));
        if (extraEdgeNodes) {
            BOOST_CHECK_GE(fa.nodes, fb.nodes);
        } else {
            BOOST_CHECK_EQUAL(fa.nodes, fb.nodes);
        }
    }
}

// Vertices lying inside an edge of a face that does not list them.
int hangingNodes(const Dune::CpGrid& g)
{
    std::vector<Point> xyz;
    for (const auto& v : vertices(g.leafGridView())) {
        xyz.push_back(v.geometry().center());
    }
    int count = 0;
    const int numFaces = g.numFaces();
    for (int f = 0; f < numFaces; ++f) {
        const int n = g.numFaceVertices(f);
        for (int e = 0; e < n; ++e) {
            const int a = g.faceVertex(f, e), b = g.faceVertex(f, (e + 1) % n);
            const Point d = xyz[b] - xyz[a];
            const double len2 = d.two_norm2();
            for (std::size_t v = 0; v < xyz.size(); ++v) {
                if (static_cast<int>(v) == a || static_cast<int>(v) == b) {
                    continue;
                }
                const Point r = xyz[v] - xyz[a];
                const double t = (r*d)/len2;
                if (t > 1e-9 && t < 1 - 1e-9 && (r - t*d).two_norm() < 1e-6) {
                    ++count;
                }
            }
        }
    }
    return count;
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
    // Corner-point processing leaves the thin cells' corners hanging on the thick
    // cells' edges; the merge lists them, which is the only difference.
    const auto collapsed = merged(fine, {south, north});
    BOOST_CHECK_EQUAL(hangingNodes(collapsed), 0);
    BOOST_CHECK_GT(hangingNodes(grid), 0);
    checkSame(shape(collapsed), shape(grid), /*extraEdgeNodes*/ true);
}

// The same on layers that dip and bend: the coarse faces are corner-point faces.
BOOST_AUTO_TEST_CASE(CollapseIsTheCornerPointCoarseningOnCurvedLayers)
{
    const auto fine = curvedGrid(6, 100.0, 10.0);
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
    // Corner-point processing leaves the thin cells' corners hanging on the thick
    // cells' edges; the merge lists them, which is the only difference.
    const auto collapsed = merged(fine, {south, north});
    BOOST_CHECK_EQUAL(hangingNodes(collapsed), 0);
    BOOST_CHECK_GT(hangingNodes(grid), 0);
    checkSame(shape(collapsed), shape(grid), /*extraEdgeNodes*/ true);
}

// A layer of zero thickness in some columns, as edge-conformal MINPV leaves it
// (its volume moved into the cell below), inside the blocks: the merged body has
// no crack there, and is the corner-point coarsening.
BOOST_AUTO_TEST_CASE(ZeroThicknessLayerInsideBlocks)
{
    auto fine = uniformGrid(6, 100.0, 10.0);
    const int n = 6;
    for (int j = 2; j < 4; ++j) {
        for (int i = 2; i < 4; ++i) {
            fine.actnum[i + n*(j + n*3)] = 0;
            for (int dj = 0; dj < 2; ++dj) {
                for (int di = 0; di < 2; ++di) {
                    const auto c = [&](int k) { return (2*i + di) + 2ull*n*((2*j + dj) + 2ull*n*k); };
                    fine.zcorn[c(2*3 + 1)] = fine.zcorn[c(2*3)];   // removed cell: top == bottom
                    fine.zcorn[c(2*4)] = fine.zcorn[c(2*3)];       // the cell below takes its volume
                }
            }
        }
    }
    CoarsenRequest layers;
    layers.startIJK = {0, 0, 2};
    layers.endIJK = {6, 6, 4};
    layers.cellsPerDim = {6, 6, 1};

    Opm::Coarsening::Options options;
    options.activity = Opm::Coarsening::Activity::FillHoles;
    const auto coarse = Opm::Coarsening::coarsenCornerPoint(fine, {layers}, options).grid;
    Dune::CpGrid grid;
    grid.processEclipseFormat(view(coarse), false, false, true);
    const auto collapsed = merged(fine, {layers});

    const auto boundaryFaces = [](const Dune::CpGrid& g) {
        int count = 0;
        for (int f = 0; f < g.numFaces(); ++f) {
            count += (g.faceCell(f, 0) < 0 || g.faceCell(f, 1) < 0) ? 1 : 0;
        }
        return count;
    };
    BOOST_CHECK_EQUAL(boundaryFaces(collapsed), boundaryFaces(grid));
    // The throw's nodes on faces between two uncoarsened cells hang in both: that
    // is the edge-conformal processing of the input, not the merge.
    BOOST_TEST_MESSAGE("hanging nodes: merged " << hangingNodes(collapsed)
                       << ", corner-point " << hangingNodes(grid));
    BOOST_CHECK_LT(hangingNodes(collapsed), hangingNodes(grid));
    checkSame(shape(collapsed), shape(grid), /*extraEdgeNodes*/ true);
}
