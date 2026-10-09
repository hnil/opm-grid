/*
  Copyright 2026 SINTEF Digital.

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
//! \file
//! \brief Micro-benchmark for AdaptiveCpGrid: time the coarse build, the
//!        (static-equivalent) refinement, and iteration/lookup over the refined
//!        leaf. A CpGrid-only reference for now (OPM's ALUGrid has no in-sim LGR
//!        refinement to compare against). Prints timings; always passes.
#include <config.h>

#define BOOST_TEST_MODULE AdaptiveCpGridBench
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/AdaptiveCpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <memory>

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
double msOf(Clock::time_point a, Clock::time_point b)
{ return std::chrono::duration<double, std::milli>(b - a).count(); }

struct MPIFixture {
    MPIFixture() {
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        Dune::MPIHelper::instance(argc, argv);
    }
};

// Uniform Cartesian macro grid (unit cells), expressed as COORD/ZCORN.
void makeCartesian(int nx, int ny, int nz,
                   std::vector<double>& coord, std::vector<double>& zcorn)
{
    coord.resize(6 * (nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            double* p = &coord[6 * (static_cast<std::size_t>(j) * (nx + 1) + i)];
            p[0] = i; p[1] = j; p[2] = 0.0;
            p[3] = i; p[4] = j; p[5] = nz;
        }
    zcorn.resize(8 * static_cast<std::size_t>(nx) * ny * nz);
    for (int k = 0; k < 2 * nz; ++k)
        for (int j = 0; j < 2 * ny; ++j)
            for (int i = 0; i < 2 * nx; ++i)
                zcorn[static_cast<std::size_t>(i) + 2 * static_cast<std::size_t>(nx) * j
                      + 4 * static_cast<std::size_t>(nx) * ny * k] = 0.5 * (k + 1);
}

// Measure access over a refined leaf: iterate+geometry and element->index lookup.
template<class Grid>
std::pair<double,double> measureAccess(const Grid& grid, double& volOut, long& chkOut)
{
    const auto leaf = grid.leafGridView();
    const long leafCells = grid.size(0);
    Dune::MultipleCodimMultipleGeomTypeMapper<std::decay_t<decltype(leaf)>>
        mapper(leaf, Dune::mcmgElementLayout());

    double bestIter = 1e30, bestLook = 1e30, vol = 0.0;
    long chk = 0;
    for (int rep = 0; rep < 3; ++rep) {       // best-of-3 to damp cache/timing noise
        auto a = Clock::now();
        vol = 0.0;
        for (const auto& e : Dune::elements(leaf)) vol += e.geometry().volume();
        auto b = Clock::now();
        chk = 0;
        for (const auto& e : Dune::elements(leaf)) chk += mapper.index(e);
        auto c = Clock::now();
        bestIter = std::min(bestIter, msOf(a, b));
        bestLook = std::min(bestLook, msOf(b, c));
    }
    volOut = vol; chkOut = chk;
    return { leafCells / (bestIter * 1e3), leafCells / (bestLook * 1e3) };  // M cells/s
}

void runBench(int nx, int ny, int nz, int fac)
{
    std::vector<double> coord, zcorn;
    makeCartesian(nx, ny, nz, coord, zcorn);
    const long coarse = static_cast<long>(nx) * ny * nz;
    double vol; long chk;

    // --- static first-level LGR: coarse + addLgrsUpdateLeafView (flow's path) ---
    Dune::CpGrid sgrid;
    grdecl g{};
    g.dims[0] = nx; g.dims[1] = ny; g.dims[2] = nz;
    g.coord = coord.data(); g.zcorn = zcorn.data(); g.actnum = nullptr;
    sgrid.processEclipseFormat(g, false);
    auto sa = Clock::now();
    sgrid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    sgrid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
        std::array<int,3>{nx, ny, nz}, coord, zcorn, std::vector<int>{}));
    sgrid.addLgrsUpdateLeafView({{fac,fac,fac}}, {{0,0,0}}, {{nx,ny,nz}}, {"LGR1"});
    auto sb = Clock::now();
    const double staticMs = msOf(sa, sb);
    auto [sIter, sLook] = measureAccess(sgrid, vol, chk);

    // --- adaptive: AdaptiveCpGrid (coarse in ctor) + adapt() ---
    Opm::AdaptiveCpGrid grid({nx, ny, nz}, coord, zcorn, {});
    grid.markBox({0,0,0}, {nx,ny,nz}, {fac,fac,fac});
    auto aa = Clock::now();
    grid.adapt();
    auto ab = Clock::now();
    const double adaptMs = msOf(aa, ab);
    const long leafCells = grid.grid().size(0);
    auto [aIter, aLook] = measureAccess(grid.grid(), vol, chk);

    std::cout << "\n=== bench " << nx << "x" << ny << "x" << nz
              << " refine whole grid x" << fac << "/dir  (coarse " << coarse
              << " -> leaf " << leafCells << ") ===\n"
              << "  REFINE   static LGR : " << staticMs << " ms ("
              << leafCells / (staticMs * 1e3) << " M cells/s)\n"
              << "  REFINE   adaptive   : " << adaptMs << " ms ("
              << leafCells / (adaptMs * 1e3) << " M cells/s)   [adaptive/static = "
              << adaptMs / staticMs << "x]\n"
              << "  ACCESS iterate+geom : static " << sIter << "  adaptive " << aIter
              << " M cells/s\n"
              << "  ACCESS index lookup : static " << sLook << "  adaptive " << aLook
              << " M cells/s\n";
    BOOST_CHECK_GT(leafCells, coarse);
    BOOST_CHECK_EQUAL(sgrid.size(0), grid.grid().size(0));  // same leaf both ways
}

} // namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(benchmark)
{
    runBench(40, 40, 10, 2);    // 16k -> ~128k
    runBench(60, 60, 20, 2);    // 72k -> ~576k
    runBench(40, 40, 10, 3);    // 16k -> ~432k (x27/cell)
}
