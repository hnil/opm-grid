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

#define BOOST_TEST_MODULE PillarBendTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>
#include <opm/grid/cpgrid/GeometricCheck.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <cmath>
#include <limits>
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
constexpr int nx = 2, ny = 1, nz = 3;
constexpr double dx = 100.0, dz = 10.0, sx = 0.5, sy = 0.2;  // pillar slant, m per m of depth

struct Grid
{
    std::vector<double> coord, zcorn;
    std::vector<int> actnum;

    // Column i = 1 is shifted down by `throw_`, a fault along the middle pillars.
    explicit Grid(double throw_ = 0.0)
    {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                coord.insert(coord.end(), {i*dx, j*dx, 0.0, i*dx + 100*sx, j*dx + 100*sy, 100.0});
            }
        }
        zcorn.resize(8*nx*ny*nz);
        for (int k = 0; k < 2*nz; ++k) {
            for (int j = 0; j < 2*ny; ++j) {
                for (int i = 0; i < 2*nx; ++i) {
                    zcorn[i + 2*nx*(j + 2*ny*k)] = ((k + 1)/2)*dz + (i >= 2 ? throw_ : 0.0);
                }
            }
        }
        actnum.assign(nx*ny*nz, 1);
    }

    grdecl view() const
    {
        grdecl g{};
        g.dims[0] = nx; g.dims[1] = ny; g.dims[2] = nz;
        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.data();
        return g;
    }
};

std::vector<double> bends(double top, double bottom)
{
    std::vector<double> b;
    for (int p = 0; p < (nx + 1)*(ny + 1); ++p) {
        b.insert(b.end(), {top, bottom});
    }
    return b;
}

std::vector<Dune::FieldVector<double,3>> vertices(const Dune::CpGrid& grid)
{
    std::vector<Dune::FieldVector<double,3>> v;
    for (const auto& vertex : vertices(grid.leafGridView())) {
        v.push_back(vertex.geometry().center());
    }
    return v;
}

// Pillar (i, j) from Grid, at depth z.
bool onPillar(const Dune::FieldVector<double,3>& v, double zAt)
{
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            if (std::hypot(v[0] - (i*dx + sx*zAt), v[1] - (j*dx + sy*zAt)) < 1e-9) {
                return true;
            }
        }
    }
    return false;
}
} // namespace

BOOST_AUTO_TEST_CASE(PillarsAreVerticalOutsideTheBends)
{
    const Grid g;
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.view(), bends(10.0, 20.0), /*edge_conformal*/ true);
    const auto v = vertices(grid);
    BOOST_CHECK_EQUAL(v.size(), std::size_t((nx + 1)*(ny + 1)*(nz + 1)));
    for (const auto& x : v) {
        BOOST_CHECK_MESSAGE(onPillar(x, std::clamp(x[2], 10.0, 20.0)),
                            "vertex " << x << " is not on its bent pillar");
    }
    BOOST_CHECK(Dune::cpgrid::checkGeometric(grid).ok());
}

BOOST_AUTO_TEST_CASE(NoBendIsTheStraightGrid)
{
    const Grid g(4.0);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Dune::CpGrid straight, unbent;
    straight.processEclipseFormat(g.view(), false, false, /*edge_conformal*/ true);
    unbent.processEclipseFormat(g.view(), bends(nan, nan), /*edge_conformal*/ true);
    const auto a = vertices(straight), b = vertices(unbent);
    BOOST_REQUIRE_EQUAL(a.size(), b.size());
    for (std::size_t n = 0; n < a.size(); ++n) {
        BOOST_CHECK(a[n] == b[n]);
    }
}

BOOST_AUTO_TEST_CASE(FaultedBentGridIsGeometric)
{
    // The bend at 10 m cuts through the shifted column's top cell (4-14 m).
    const Grid g(4.0);
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.view(), bends(10.0, 24.0), /*edge_conformal*/ true);
    const auto check = Dune::cpgrid::checkGeometric(grid);
    BOOST_CHECK_MESSAGE(check.ok(), check.summary());
    double volume = 0.0;
    for (const auto& cell : elements(grid.leafGridView())) {
        volume += cell.geometry().volume();
    }
    // Fault nodes sample the bent pillar, so cells next to a bend are not prisms.
    BOOST_REQUIRE(!check.boundaries.empty());
    BOOST_CHECK_CLOSE(volume, check.boundaries[0].volume, 1e-9);
    BOOST_CHECK_CLOSE(volume, nx*ny*dx*dx*nz*dz, 0.5);
}

BOOST_AUTO_TEST_CASE(WrongBendSizeThrows)
{
    const Grid g;
    Dune::CpGrid grid;
    BOOST_CHECK_THROW(grid.processEclipseFormat(g.view(), std::vector<double>(3, 0.0), true),
                      std::invalid_argument);
}
