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

#define BOOST_TEST_MODULE CornerPointCoarseningTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>

#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

using namespace Opm::Coarsening;

namespace
{

std::size_t cellIdx(const std::array<int,3>& d, int i, int j, int k)
{
    return static_cast<std::size_t>(i) + static_cast<std::size_t>(d[0])*j
        + static_cast<std::size_t>(d[0])*d[1]*k;
}

std::size_t cornerIdx(const std::array<int,3>& d, int i, int j, int k)
{
    return static_cast<std::size_t>(i) + 2*static_cast<std::size_t>(d[0])*j
        + 4*static_cast<std::size_t>(d[0])*d[1]*k;
}

/// Tensor grid with the given cell sizes; layer tops stacked from z0.
Grdecl tensorGrid(const std::vector<double>& dx, const std::vector<double>& dy,
                  const std::vector<double>& dz, double z0 = 1000.0)
{
    Grdecl g;
    g.dims = {int(dx.size()), int(dy.size()), int(dz.size())};

    std::vector<double> xs{0.0}, ys{0.0}, zs{z0};
    for (double d : dx) { xs.push_back(xs.back() + d); }
    for (double d : dy) { ys.push_back(ys.back() + d); }
    for (double d : dz) { zs.push_back(zs.back() + d); }

    for (double y : ys) {
        for (double x : xs) {
            for (double z : {zs.front(), zs.back()}) {
                g.coord.insert(g.coord.end(), {x, y, z});
            }
        }
    }
    // COORD holds top then bottom point per pillar: fix the interleaving above.
    std::vector<double> coord;
    for (std::size_t p = 0; p < xs.size()*ys.size(); ++p) {
        const double x = xs[p % xs.size()], y = ys[p / xs.size()];
        coord.insert(coord.end(), {x, y, zs.front(), x, y, zs.back()});
    }
    g.coord = coord;

    g.zcorn.assign(8*static_cast<std::size_t>(g.dims[0])*g.dims[1]*g.dims[2], 0.0);
    for (int k = 0; k < g.dims[2]; ++k) {
        for (int dk = 0; dk < 2; ++dk) {
            for (int j = 0; j < 2*g.dims[1]; ++j) {
                for (int i = 0; i < 2*g.dims[0]; ++i) {
                    g.zcorn[cornerIdx(g.dims, i, j, 2*k + dk)] = zs[k + dk];
                }
            }
        }
    }
    g.actnum.assign(static_cast<std::size_t>(g.dims[0])*g.dims[1]*g.dims[2], 1);
    return g;
}

Grdecl uniformGrid(int nx, int ny, int nz)
{
    return tensorGrid(std::vector<double>(nx, 100.0),
                      std::vector<double>(ny, 100.0),
                      std::vector<double>(nz, 10.0));
}

CoarsenRequest wholeGrid(const Grdecl& g, int cx, int cy, int cz)
{
    CoarsenRequest r;
    r.startIJK = {0, 0, 0};
    r.endIJK = g.dims;
    r.cellsPerDim = {cx, cy, cz};
    return r;
}

double bulkVolume(const Grdecl& g)
{
    // Tensor grids only: sum of dx*dy*dz taken from the corner data.
    double total = 0.0;
    for (int k = 0; k < g.dims[2]; ++k) {
        for (int j = 0; j < g.dims[1]; ++j) {
            for (int i = 0; i < g.dims[0]; ++i) {
                const double* p0 = &g.coord[6*(static_cast<std::size_t>(j)*(g.dims[0] + 1) + i)];
                const double* p1 = &g.coord[6*(static_cast<std::size_t>(j)*(g.dims[0] + 1) + i + 1)];
                const double* p2 = &g.coord[6*(static_cast<std::size_t>(j + 1)*(g.dims[0] + 1) + i)];
                const double top = g.zcorn[cornerIdx(g.dims, 2*i, 2*j, 2*k)];
                const double bot = g.zcorn[cornerIdx(g.dims, 2*i, 2*j, 2*k + 1)];
                total += (p1[0] - p0[0])*(p2[1] - p0[1])*(bot - top);
            }
        }
    }
    return total;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(NoRequestsIsTheIdentity)
{
    const auto fine = uniformGrid(4, 3, 5);
    const auto result = coarsenCornerPoint(fine, {});

    BOOST_CHECK(result.grid.dims == fine.dims);
    BOOST_CHECK(result.grid.coord == fine.coord);
    BOOST_CHECK(result.grid.zcorn == fine.zcorn);
    BOOST_CHECK(result.grid.actnum == fine.actnum);
    for (std::size_t c = 0; c < result.fineToCoarse.size(); ++c) {
        BOOST_CHECK_EQUAL(result.fineToCoarse[c], int(c));
    }
    BOOST_CHECK(refinementBackToFine(result).empty());
}

BOOST_AUTO_TEST_CASE(UniformCoarseningKeepsVolumeAndTilesTheFineGrid)
{
    const auto fine = uniformGrid(4, 4, 6);
    const auto result = coarsenCornerPoint(fine, {wholeGrid(fine, 2, 2, 3)});

    BOOST_CHECK((result.grid.dims == std::array<int,3>{2, 2, 3}));
    BOOST_CHECK_CLOSE(bulkVolume(result.grid), bulkVolume(fine), 1e-9);

    // Every fine cell belongs to exactly one coarse cell, and the blocks tile.
    std::vector<int> children(result.blocks.size(), 0);
    for (int c : result.fineToCoarse) {
        BOOST_REQUIRE(c >= 0 && c < int(children.size()));
        ++children[c];
    }
    for (std::size_t b = 0; b < result.blocks.size(); ++b) {
        const auto& blk = result.blocks[b];
        const int n = (blk.endIJK[0] - blk.startIJK[0])*(blk.endIJK[1] - blk.startIJK[1])
            * (blk.endIJK[2] - blk.startIJK[2]);
        BOOST_CHECK_EQUAL(children[b], n);
        BOOST_CHECK_EQUAL(n, 2*2*2);
    }

    const auto lgr = refinementBackToFine(result);
    BOOST_REQUIRE_EQUAL(lgr.size(), 1u);
    BOOST_CHECK((lgr[0].cellsPerDim == std::array<int,3>{2, 2, 2}));
    BOOST_CHECK((lgr[0].startIJK == std::array<int,3>{0, 0, 0}));
    BOOST_CHECK((lgr[0].endIJK == result.grid.dims));
}

BOOST_AUTO_TEST_CASE(UnevenSplit)
{
    const auto fine = uniformGrid(5, 1, 4);
    CoarsenRequest r;
    r.startIJK = {0, 0, 0};
    r.endIJK = {5, 1, 4};
    r.cellsPerDim = {2, 1, 1};

    const auto result = coarsenCornerPoint(fine, {r});
    // 5 cells into 2 groups: 3 + 2.
    BOOST_CHECK((result.grid.dims == std::array<int,3>{2, 1, 1}));
    BOOST_CHECK_EQUAL(result.blocks[0].endIJK[0] - result.blocks[0].startIJK[0], 3);
    BOOST_CHECK_EQUAL(result.blocks[1].endIJK[0] - result.blocks[1].startIJK[0], 2);
    BOOST_CHECK_CLOSE(bulkVolume(result.grid), bulkVolume(fine), 1e-9);
}

BOOST_AUTO_TEST_CASE(PartialKBoxCoarsensLayersOnly)
{
    const auto fine = uniformGrid(1, 1, 4);
    CoarsenRequest r;
    r.startIJK = {0, 0, 0};
    r.endIJK = {1, 1, 2};
    r.cellsPerDim = {1, 1, 1};

    const auto result = coarsenCornerPoint(fine, {r});
    BOOST_CHECK((result.grid.dims == std::array<int,3>{1, 1, 3}));
    BOOST_CHECK_CLOSE(bulkVolume(result.grid), bulkVolume(fine), 1e-9);
}

BOOST_AUTO_TEST_CASE(LateralCoarseningOfPartOfTheColumnIsRefused)
{
    // Pillars run through every layer, so a grdecl cannot coarsen laterally in
    // one k-range only; that case needs the CpGridData merge pass.
    const auto fine = uniformGrid(4, 1, 4);
    CoarsenRequest r;
    r.startIJK = {0, 0, 0};
    r.endIJK = {4, 1, 2};
    r.cellsPerDim = {2, 1, 1};
    BOOST_CHECK_THROW(coarsenCornerPoint(fine, {r}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(LateralGroupingMustBeTheSameEverywhere)
{
    const auto fine = uniformGrid(4, 4, 2);
    CoarsenRequest r;              // merges i only in the lower half of j
    r.startIJK = {0, 0, 0};
    r.endIJK = {4, 2, 2};
    r.cellsPerDim = {2, 2, 1};
    BOOST_CHECK_THROW(coarsenCornerPoint(fine, {r}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(FaultInsideABlockIsRefused)
{
    auto fine = uniformGrid(4, 2, 2);
    // Throw at the i-pillar between cells 1 and 2, which a 4->2 merge drops.
    for (int k = 0; k < 2*fine.dims[2]; ++k) {
        for (int j = 0; j < 2*fine.dims[1]; ++j) {
            fine.zcorn[cornerIdx(fine.dims, 2, j, k)] += 3.0;   // right side of cell 0
            fine.zcorn[cornerIdx(fine.dims, 3, j, k)] += 3.0;
        }
    }
    BOOST_CHECK_THROW(coarsenCornerPoint(fine, {wholeGrid(fine, 2, 1, 1)}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(VerticalGapNeedsTheFlagAndIsThenAbsorbed)
{
    auto fine = uniformGrid(1, 1, 2);
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
            fine.zcorn[cornerIdx(fine.dims, i, j, 2)] += 5.0;   // top of layer 2
            fine.zcorn[cornerIdx(fine.dims, i, j, 3)] += 5.0;   // bottom of layer 2
        }
    }
    const auto request = wholeGrid(fine, 1, 1, 1);
    BOOST_CHECK_THROW(coarsenCornerPoint(fine, {request}), std::invalid_argument);

    Options options;
    options.allowVerticalGaps = true;
    const auto result = coarsenCornerPoint(fine, {request}, options);
    BOOST_CHECK_GT(result.report.absorbedGapVolume, 0.0);
    // The coarse cell spans from the top of layer 1 to the base of layer 2.
    const double top = result.grid.zcorn[cornerIdx(result.grid.dims, 0, 0, 0)];
    const double bot = result.grid.zcorn[cornerIdx(result.grid.dims, 0, 0, 1)];
    BOOST_CHECK_CLOSE(bot - top, 25.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(HolesAreFilledForMechanics)
{
    auto fine = uniformGrid(2, 2, 2);
    fine.actnum[cellIdx(fine.dims, 1, 1, 1)] = 0;

    const auto any = coarsenCornerPoint(fine, {wholeGrid(fine, 1, 1, 1)});
    BOOST_CHECK_EQUAL(any.grid.actnum[0], 1);
    BOOST_CHECK_EQUAL(any.report.holesFilled, 1);

    Options strict;
    strict.activity = Activity::AllChildrenActive;
    BOOST_CHECK_EQUAL(coarsenCornerPoint(fine, {wholeGrid(fine, 1, 1, 1)}, strict).grid.actnum[0], 0);

    // A block whose cells are all inactive is still rock under FillHoles.
    std::fill(fine.actnum.begin(), fine.actnum.end(), 0);
    Options fill;
    fill.activity = Activity::FillHoles;
    const auto filled = coarsenCornerPoint(fine, {wholeGrid(fine, 1, 1, 1)}, fill);
    BOOST_CHECK_EQUAL(filled.grid.actnum[0], 1);
    BOOST_CHECK_EQUAL(filled.report.holesFilled, 8);
    BOOST_CHECK_CLOSE(filled.report.holeVolume, 8*100.0*100.0*10.0, 1e-9);

    const auto children = childActnum(fine, filled.blocks[0]);
    BOOST_CHECK_EQUAL(children.size(), 8u);
    BOOST_CHECK_EQUAL(std::accumulate(children.begin(), children.end(), 0), 0);
}

BOOST_AUTO_TEST_CASE(ColumnsWithDifferentLayerGroupingSplitFaces)
{
    const auto fine = uniformGrid(2, 1, 4);
    CoarsenRequest r;              // merge all four layers of the first column
    r.startIJK = {0, 0, 0};
    r.endIJK = {1, 1, 4};
    r.cellsPerDim = {1, 1, 1};

    const auto result = coarsenCornerPoint(fine, {r});
    BOOST_CHECK((result.grid.dims == std::array<int,3>{2, 1, 4}));
    BOOST_CHECK_EQUAL(result.report.maxSubFacesSeen, 4);

    Options capped;
    capped.maxSubFaces = 2;
    BOOST_CHECK_THROW(coarsenCornerPoint(fine, {r}, capped), std::invalid_argument);

    // The short column is padded with collapsed cells at its base.
    const auto& d = result.grid.dims;
    for (int kc = 1; kc < d[2]; ++kc) {
        const double top = result.grid.zcorn[cornerIdx(d, 0, 0, 2*kc)];
        const double bot = result.grid.zcorn[cornerIdx(d, 0, 0, 2*kc + 1)];
        BOOST_CHECK_CLOSE(top, bot, 1e-12);
        BOOST_CHECK_EQUAL(result.grid.actnum[cellIdx(d, 0, 0, kc)], 0);
    }
}

BOOST_AUTO_TEST_CASE(NonNestedGroupingsAreRefused)
{
    const auto fine = uniformGrid(2, 1, 6);
    CoarsenRequest a, b;
    a.startIJK = {0, 0, 0};
    a.endIJK = {1, 1, 6};
    a.cellsPerDim = {1, 1, 2};     // 3 + 3
    b.startIJK = {1, 0, 0};
    b.endIJK = {2, 1, 6};
    b.cellsPerDim = {1, 1, 3};     // 2 + 2 + 2, which does not nest with 3 + 3

    BOOST_CHECK_THROW(coarsenCornerPoint(fine, {a, b}), std::invalid_argument);

    Options ungraded;
    ungraded.requireGradedColumns = false;
    BOOST_CHECK_NO_THROW(coarsenCornerPoint(fine, {a, b}, ungraded));
}

BOOST_AUTO_TEST_CASE(OverlappingRequestsAreRefused)
{
    const auto fine = uniformGrid(4, 4, 2);
    CoarsenRequest a, b;
    a.startIJK = {0, 0, 0};
    a.endIJK = {4, 4, 2};
    a.cellsPerDim = {2, 2, 1};
    b = a;
    BOOST_CHECK_THROW(coarsenCornerPoint(fine, {a, b}), std::invalid_argument);
}
