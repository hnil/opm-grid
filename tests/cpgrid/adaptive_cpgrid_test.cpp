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
#include <config.h>

#define BOOST_TEST_MODULE AdaptiveCpGridTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/AdaptiveCpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace
{

struct MPIFixture
{
    MPIFixture()
    {
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        Dune::MPIHelper::instance(argc, argv);
    }
};

// Macro corner-point description plus the raw grdecl view used to build the
// coarse Dune::CpGrid. Mirrors the helper in conforming_builder_test.cpp.
struct TestGrdecl
{
    std::array<int,3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;

    grdecl raw() const
    {
        grdecl g{};
        g.dims[0] = dims[0]; g.dims[1] = dims[1]; g.dims[2] = dims[2];
        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.empty() ? nullptr : actnum.data();
        return g;
    }
};

int cellOf(int doubled) { return doubled / 2; }
int sideOf(int doubled) { return doubled % 2; }

// A distorted vertical-pillar grid: clean hexes (so the trilinear refinement is
// exact), but non-planar tops -- a faithful stand-in for a real corner-point
// CARFIN deck at the grid level.
TestGrdecl makeVerticalPillarGrid(const std::array<int,3>& dims,
                                  const std::function<double(int,int,int)>& depth)
{
    TestGrdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;

    g.coord.resize(6 * (nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6 * (static_cast<std::size_t>(j) * (nx + 1) + i)];
            p[0] = i; p[1] = j; p[2] = 0.0;
            p[3] = i; p[4] = j; p[5] = 100.0;
        }
    }
    g.zcorn.resize(8 * static_cast<std::size_t>(nx) * ny * nz);
    for (int k_ = 0; k_ < 2 * nz; ++k_) {
        for (int j_ = 0; j_ < 2 * ny; ++j_) {
            for (int i_ = 0; i_ < 2 * nx; ++i_) {
                g.zcorn[static_cast<std::size_t>(i_) + 2 * static_cast<std::size_t>(nx) * j_
                        + 4 * static_cast<std::size_t>(nx) * ny * k_] = depth(i_, j_, k_);
            }
        }
    }
    return g;
}

TestGrdecl sampleGrid(const std::array<int,3>& dims)
{
    auto depth = [](int i_, int j_, int k_) {
        const double x = cellOf(i_) + sideOf(i_);
        const double y = cellOf(j_) + sideOf(j_);
        return 2.0 * (cellOf(k_) + sideOf(k_)) + 0.3 * std::sin(0.9 * x) + 0.2 * std::cos(0.6 * y);
    };
    return makeVerticalPillarGrid(dims, depth);
}


double totalVolume(const Dune::CpGrid& grid)
{
    double vol = 0.0;
    for (const auto& e : Dune::elements(grid.leafGridView())) {
        vol += e.geometry().volume();
    }
    return vol;
}

std::vector<double> sortedLeafVolumes(const Dune::CpGrid& grid)
{
    std::vector<double> v;
    for (const auto& e : Dune::elements(grid.leafGridView())) {
        v.push_back(e.geometry().volume());
    }
    std::ranges::sort(v);
    return v;
}

// Every leaf cell is closed (sum of area-weighted outward normals is zero).
void checkConformalLeaf(const Dune::CpGrid& grid)
{
    for (const auto& e : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> closure(0.0);
        for (const auto& is : Dune::intersections(grid.leafGridView(), e)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            closure += n;
        }
        BOOST_CHECK_SMALL(closure.two_norm(), 1e-9);
    }
}

// Build the reference grid the *static* CARFIN path produces: a coarse grid
// refined in place by addLgrsUpdateLeafView. This is exactly what flow does for
// a CARFIN deck (boxes from the deck), and is our equivalence oracle.
std::shared_ptr<Dune::CpGrid>
staticRefined(const TestGrdecl& g,
              const std::vector<std::array<int,3>>& cellsPerDim,
              const std::vector<std::array<int,3>>& startIJK,
              const std::vector<std::array<int,3>>& endIJK,
              const std::vector<std::string>& names)
{
    auto grid = std::make_shared<Dune::CpGrid>();
    auto raw = g.raw();
    grid->processEclipseFormat(raw, false);
    grid->setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid->setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));
    grid->addLgrsUpdateLeafView(cellsPerDim, startIJK, endIJK, names);
    return grid;
}

// The two leaves must be the same grid: same level structure, same cell counts,
// and the same set of cell geometries (sorted volumes). Same builder + same box
// + same macro input => identical discretisation => identical simulation.
void checkSameLeaf(const Dune::CpGrid& adaptive, const Dune::CpGrid& reference)
{
    BOOST_REQUIRE_EQUAL(adaptive.maxLevel(), reference.maxLevel());
    BOOST_REQUIRE_EQUAL(adaptive.size(0), reference.size(0));        // leaf cells
    BOOST_CHECK_EQUAL(adaptive.size(0, 0), reference.size(0, 0));    // level-0 cells
    for (int lvl = 1; lvl <= reference.maxLevel(); ++lvl) {
        BOOST_CHECK_EQUAL(adaptive.size(lvl, 0), reference.size(lvl, 0));
    }
    BOOST_CHECK_CLOSE(totalVolume(adaptive), totalVolume(reference), 1e-9);

    const auto va = sortedLeafVolumes(adaptive);
    const auto vr = sortedLeafVolumes(reference);
    BOOST_REQUIRE_EQUAL(va.size(), vr.size());
    for (std::size_t i = 0; i < va.size(); ++i) {
        BOOST_CHECK_CLOSE(va[i], vr[i], 1e-9);
    }
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

// The central claim: refining a marked box AFTER constructing a coarse grid
// (AdaptiveCpGrid, no deck CARFIN) yields the same leaf as refining the same box
// AT construction (the static CARFIN path).
BOOST_AUTO_TEST_CASE(adaptiveSingleBoxMatchesStatic)
{
    const std::array<int,3> dims{4, 3, 3};
    auto g = sampleGrid(dims);

    // One CARFIN-equivalent box: parents [1,3)x[1,2)x[1,2) refined 2x2x2.
    auto reference = staticRefined(g, {{2,2,2}}, {{1,1,1}}, {{3,2,2}}, {"LGR1"});

    Opm::AdaptiveCpGrid adaptive(dims, g.coord, g.zcorn, g.actnum);
    BOOST_CHECK(!adaptive.refined());
    BOOST_CHECK_EQUAL(adaptive.grid().maxLevel(), 0);          // coarse before adapt
    adaptive.markBox({1,1,1}, {3,2,2}, {2,2,2});
    BOOST_CHECK_EQUAL(adaptive.markCount(), 2u);   // box [1,3)x[1,2)x[1,2) = 2 cells
    adaptive.adapt();

    BOOST_CHECK(adaptive.refined());
    checkConformalLeaf(adaptive.grid());
    checkSameLeaf(adaptive.grid(), *reference);
}

// Multiple separated boxes (the multi-LGR / model2 pattern) marked and applied
// in one adapt() reproduce the multi-box static refinement.
BOOST_AUTO_TEST_CASE(adaptiveTwoBoxesMatchStatic)
{
    const std::array<int,3> dims{6, 3, 3};
    auto g = sampleGrid(dims);

    auto reference = staticRefined(g,
        {{2,2,2}, {3,3,3}},
        {{0,0,0}, {4,1,1}},
        {{1,1,1}, {5,2,2}},
        {"LGR1", "LGR2"});

    Opm::AdaptiveCpGrid adaptive(dims, g.coord, g.zcorn, g.actnum);
    adaptive.markBox({0,0,0}, {1,1,1}, {2,2,2});
    adaptive.markBox({4,1,1}, {5,2,2}, {3,3,3});
    adaptive.adapt();

    checkConformalLeaf(adaptive.grid());
    checkSameLeaf(adaptive.grid(), *reference);
}

// markCell is just a one-cell box.
BOOST_AUTO_TEST_CASE(markCellMatchesSingleCellBox)
{
    const std::array<int,3> dims{4, 3, 3};
    auto g = sampleGrid(dims);

    auto reference = staticRefined(g, {{3,3,3}}, {{2,1,1}}, {{3,2,2}}, {"LGR1"});

    Opm::AdaptiveCpGrid adaptive(dims, g.coord, g.zcorn, g.actnum);
    adaptive.markCell({2,1,1}, {3,3,3});
    adaptive.adapt();

    checkSameLeaf(adaptive.grid(), *reference);
}

// A purely Cartesian macro grid (uniform pillars) refined by a CARFIN-style box:
// volume conserved and the leaf is the same as the static path. This is the
// shape of the SPE1 CARFIN cases at the grid level.
BOOST_AUTO_TEST_CASE(adaptiveCartesianBox)
{
    const std::array<int,3> dims{10, 10, 3};
    auto flat = [](int, int, int k_) { return 1.0 * (cellOf(k_) + sideOf(k_)); };
    auto g = makeVerticalPillarGrid(dims, flat);

    auto reference = staticRefined(g, {{3,3,3}}, {{8,8,0}}, {{10,10,3}}, {"LGR1"});

    Opm::AdaptiveCpGrid adaptive(dims, g.coord, g.zcorn, g.actnum);
    const double coarseVol = totalVolume(adaptive.grid());
    adaptive.markBox({8,8,0}, {10,10,3}, {3,3,3});
    adaptive.adapt();

    BOOST_CHECK_CLOSE(totalVolume(adaptive.grid()), coarseVol, 1e-9);
    checkConformalLeaf(adaptive.grid());
    checkSameLeaf(adaptive.grid(), *reference);
}

// Re-adapt of an already-refined grid is not implemented in this first cut.
// Cell-by-cell refinement of a contiguous region == one-go static refinement of
// the same region (adjacent same-factor cells merge into one box).
BOOST_AUTO_TEST_CASE(cellByCellMatchesOneGo)
{
    const std::array<int,3> dims{4, 3, 3};
    auto g = sampleGrid(dims);
    // one CARFIN box over the region [1,3)x[1,2)x[1,2) -> cells (1,1,1),(2,1,1).
    auto reference = staticRefined(g, {{2,2,2}}, {{1,1,1}}, {{3,2,2}}, {"LGR1"});

    // (a) mark each cell, one adapt(): the two adjacent marks merge into the box.
    Opm::AdaptiveCpGrid a(dims, g.coord, g.zcorn, g.actnum);
    a.markCell({1,1,1}, {2,2,2});
    a.markCell({2,1,1}, {2,2,2});
    a.adapt();
    checkSameLeaf(a.grid(), *reference);

    // (b) incremental: refine one cell, then re-adapt with the second added.
    Opm::AdaptiveCpGrid b(dims, g.coord, g.zcorn, g.actnum);
    b.markCell({1,1,1}, {2,2,2});
    b.adapt();
    BOOST_CHECK(b.refined());
    b.markCell({2,1,1}, {2,2,2});
    b.adapt();                                  // re-adapt: union -> merged box
    checkSameLeaf(b.grid(), *reference);
}

// Re-adaptation of two separated regions == refining both in one go.
BOOST_AUTO_TEST_CASE(reAdaptRefinesUnion)
{
    const std::array<int,3> dims{6, 3, 3};
    auto g = sampleGrid(dims);
    auto reference = staticRefined(g,
        {{2,2,2}, {2,2,2}}, {{0,0,0}, {4,1,1}}, {{1,1,1}, {5,2,2}}, {"L1", "L2"});

    Opm::AdaptiveCpGrid adaptive(dims, g.coord, g.zcorn, g.actnum);
    adaptive.markBox({0,0,0}, {1,1,1}, {2,2,2});
    adaptive.adapt();
    BOOST_CHECK(adaptive.refined());
    adaptive.markBox({4,1,1}, {5,2,2}, {2,2,2});
    adaptive.adapt();                           // re-adapt: now both regions
    checkSameLeaf(adaptive.grid(), *reference);
}

// Two boxes stacked on top of each other (touching in K) with DIFFERENT vertical
// (K) refinement but the SAME lateral (I,J) refinement build fine: the shared
// face lies in the I-J plane, so only the in-face (I,J) subdivisions must match;
// the perpendicular touch direction (K) is unconstrained. The merge keeps the two
// different-factor slabs as two boxes, and the builder accepts the K-only
// difference. (This is the adaptive analogue of the TLGR_VSTACK_42 deck.)
BOOST_AUTO_TEST_CASE(verticalStackDifferentKMatchesStatic)
{
    const std::array<int,3> dims{2, 2, 4};
    auto g = sampleGrid(dims);

    // Lower slab k in [0,2): x4 in K; upper slab k in [2,4): x2 in K. Same I,J (x2).
    // Box order follows the merge's factor-sorted order ({2,2,2} before {2,2,4}) so
    // the two LGR *levels* line up with the adaptive grid's (the leaf is identical
    // either way; checkSameLeaf compares level-by-level).
    auto reference = staticRefined(g,
        {{2,2,2}, {2,2,4}},
        {{0,0,2}, {0,0,0}},
        {{2,2,4}, {2,2,2}},
        {"VHI", "VLO"});

    Opm::AdaptiveCpGrid adaptive(dims, g.coord, g.zcorn, g.actnum);
    adaptive.markBox({0,0,0}, {2,2,2}, {2,2,4});   // lower slab, fine in K
    adaptive.markBox({0,0,2}, {2,2,4}, {2,2,2});   // upper slab, coarse in K
    adaptive.adapt();

    BOOST_CHECK(adaptive.refined());
    checkConformalLeaf(adaptive.grid());
    checkSameLeaf(adaptive.grid(), *reference);
}

// A2 compatible sub-face mosaic through the adaptive mark API: marking two
// stacked regions with DIFFERENT but compatible in-face (I) factors -- x4 over
// x2 -- builds a conformal leaf via the coarser-side sub-face mosaic. A
// horizontal (K) interface carries no fault, so this is the pure mosaic case.
BOOST_AUTO_TEST_CASE(adaptiveStackedCompatibleInFaceMosaic)
{
    const std::array<int,3> dims{2, 2, 2};
    auto flat = [](int, int, int k_) { return 1.0 * (cellOf(k_) + sideOf(k_)); };
    auto g = makeVerticalPillarGrid(dims, flat);

    // Static reference: the same two stacked boxes refined at construction. The
    // merge orders boxes by factor ({2,2,2} before {4,2,2}), so list HI first.
    auto reference = staticRefined(g,
        {{2,2,2}, {4,2,2}}, {{0,0,1}, {0,0,0}}, {{2,2,2}, {2,2,1}}, {"HI", "LO"});

    Opm::AdaptiveCpGrid adaptive(dims, g.coord, g.zcorn, g.actnum);
    adaptive.markBox({0,0,0}, {2,2,1}, {4,2,2});   // lower slab, fine in I
    adaptive.markBox({0,0,1}, {2,2,2}, {2,2,2});   // upper slab, coarse in I
    adaptive.adapt();

    BOOST_CHECK(adaptive.refined());
    checkConformalLeaf(adaptive.grid());
    checkSameLeaf(adaptive.grid(), *reference);     // == refined-at-construction
}

// No marks => adapt() is a no-op and the grid stays coarse.
BOOST_AUTO_TEST_CASE(adaptWithoutMarksIsNoop)
{
    const std::array<int,3> dims{4, 3, 3};
    auto g = sampleGrid(dims);
    Opm::AdaptiveCpGrid adaptive(dims, g.coord, g.zcorn, g.actnum);
    adaptive.adapt();
    BOOST_CHECK(!adaptive.refined());
    BOOST_CHECK_EQUAL(adaptive.grid().maxLevel(), 0);
}
