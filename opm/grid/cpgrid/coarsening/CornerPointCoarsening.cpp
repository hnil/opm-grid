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
#include "config.h"

#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using Opm::Coarsening::Block;
using Opm::Coarsening::Grdecl;

constexpr double Z_TOL = 1e-6;     // m, ZCORN agreement
constexpr double PILLAR_TOL = 1e-6; // m, collinearity of a dropped pillar

std::size_t cellIndex(const std::array<int,3>& dims, int i, int j, int k)
{
    return static_cast<std::size_t>(i)
        + static_cast<std::size_t>(dims[0])*j
        + static_cast<std::size_t>(dims[0])*dims[1]*k;
}

/// ZCORN index of corner (i,j,k) in doubled corner coordinates.
std::size_t cornerIndex(const std::array<int,3>& dims, int i, int j, int k)
{
    return static_cast<std::size_t>(i)
        + 2*static_cast<std::size_t>(dims[0])*j
        + 4*static_cast<std::size_t>(dims[0])*dims[1]*k;
}

double cellZ(const Grdecl& g, int i, int j, int k, int di, int dj, int dk)
{
    return g.zcorn[cornerIndex(g.dims, 2*i + di, 2*j + dj, 2*k + dk)];
}

const double* pillar(const Grdecl& g, int I, int J)
{
    return &g.coord[6*(static_cast<std::size_t>(J)*(g.dims[0] + 1) + I)];
}

[[noreturn]] void fail(const std::string& what)
{
    throw std::invalid_argument("corner-point coarsening: " + what);
}

std::string ijk(int i, int j, int k)
{
    std::ostringstream os;
    os << '(' << i + 1 << ',' << j + 1 << ',' << k + 1 << ')';
    return os.str();
}

/// n cells into m groups, sizes differing by at most one.
std::vector<int> evenSplit(int n, int m)
{
    std::vector<int> sizes(m, n/m);
    for (int r = 0; r < n % m; ++r) {
        ++sizes[r];
    }
    return sizes;
}

/// Corner coordinates of a cell, ordered (di,dj,dk) = 0..7.
std::array<std::array<double,3>,8> cellCorners(const Grdecl& g, int i, int j, int k)
{
    std::array<std::array<double,3>,8> c{};
    for (int dk = 0; dk < 2; ++dk) {
        for (int dj = 0; dj < 2; ++dj) {
            for (int di = 0; di < 2; ++di) {
                const double* p = pillar(g, i + di, j + dj);
                const double z = cellZ(g, i, j, k, di, dj, dk);
                const double zt = p[2], zb = p[5];
                const double t = (std::abs(zb - zt) < 1e-12) ? 0.0 : (z - zt)/(zb - zt);
                c[4*dk + 2*dj + di] = {p[0] + t*(p[3] - p[0]),
                                       p[1] + t*(p[4] - p[1]),
                                       z};
            }
        }
    }
    return c;
}

double tetVolume(const std::array<double,3>& a, const std::array<double,3>& b,
                 const std::array<double,3>& c, const std::array<double,3>& d)
{
    const double u[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
    const double v[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
    const double w[3] = {d[0]-a[0], d[1]-a[1], d[2]-a[2]};
    const double det = u[0]*(v[1]*w[2] - v[2]*w[1])
                     - u[1]*(v[0]*w[2] - v[2]*w[0])
                     + u[2]*(v[0]*w[1] - v[1]*w[0]);
    return std::abs(det)/6.0;
}

/// Hexahedron volume, split into five tetrahedra. Exact for the trilinear
/// cells we build here; good enough for the activity and report numbers.
double hexVolume(const std::array<std::array<double,3>,8>& c)
{
    return tetVolume(c[0], c[1], c[2], c[4])
         + tetVolume(c[1], c[2], c[4], c[7])
         + tetVolume(c[1], c[4], c[5], c[7])
         + tetVolume(c[1], c[2], c[3], c[7])
         + tetVolume(c[2], c[4], c[6], c[7]);
}

/// Distance of the fourth corner from the plane of the first three, relative
/// to the face size.
double nonPlanarity(const std::array<double,3>& a, const std::array<double,3>& b,
                    const std::array<double,3>& c, const std::array<double,3>& d)
{
    const double u[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
    const double v[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
    const double n[3] = {u[1]*v[2] - u[2]*v[1],
                         u[2]*v[0] - u[0]*v[2],
                         u[0]*v[1] - u[1]*v[0]};
    const double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (len < 1e-12) {
        return 0.0;
    }
    const double dist = std::abs(n[0]*(d[0]-a[0]) + n[1]*(d[1]-a[1]) + n[2]*(d[2]-a[2]))/len;
    const double scale = std::sqrt(len);
    return (scale < 1e-12) ? 0.0 : dist/scale;
}

/// Cumulative group start indices, plus the end: sizes {2,1,1} -> {0,2,3,4}.
std::vector<int> groupStarts(const std::vector<int>& sizes)
{
    std::vector<int> starts(sizes.size() + 1, 0);
    std::partial_sum(sizes.begin(), sizes.end(), starts.begin() + 1);
    return starts;
}

} // anonymous namespace

namespace Opm
{
namespace Coarsening
{

namespace
{

/// How the fine cells are merged: tensor-product laterally, per coarse column
/// vertically.
struct Grouping
{
    std::vector<int> iSizes, jSizes;
    /// One entry per coarse column (jc*ncx + ic); each sums to nz.
    std::vector<std::vector<int>> kSizes;

    int ncx() const { return static_cast<int>(iSizes.size()); }
    int ncy() const { return static_cast<int>(jSizes.size()); }
};

/// Per-cell coarsening factors along one axis, taken from the requests.
/// Returns, for each fine cell of that axis, the group sizes seen by the
/// column (i.e. the grouping the requests prescribe there).
std::vector<int> axisGrouping(int n, int axis,
                              const std::vector<CoarsenRequest>& requests,
                              const std::vector<int>& requestOfCell)
{
    std::vector<int> sizes;
    int i = 0;
    while (i < n) {
        const int r = requestOfCell[i];
        if (r < 0) {
            sizes.push_back(1);
            ++i;
            continue;
        }
        const auto& req = requests[r];
        const int len = req.endIJK[axis] - req.startIJK[axis];
        for (int s : evenSplit(len, req.cellsPerDim[axis])) {
            sizes.push_back(s);
        }
        i = req.endIJK[axis];
    }
    return sizes;
}

/// The request covering each cell of one line through the grid, or -1.
std::vector<int> requestsAlong(int n, int axis, const std::array<int,3>& other,
                               const std::vector<CoarsenRequest>& requests)
{
    std::vector<int> of(n, -1);
    for (std::size_t r = 0; r < requests.size(); ++r) {
        const auto& q = requests[r];
        bool inside = true;
        for (int d = 0; d < 3; ++d) {
            if (d == axis) {
                continue;
            }
            inside = inside && other[d] >= q.startIJK[d] && other[d] < q.endIJK[d];
        }
        if (!inside) {
            continue;
        }
        for (int i = q.startIJK[axis]; i < q.endIJK[axis]; ++i) {
            of[i] = static_cast<int>(r);
        }
    }
    return of;
}

void checkRequests(const std::array<int,3>& dims, const std::vector<CoarsenRequest>& requests)
{
    std::vector<int> owner(static_cast<std::size_t>(dims[0])*dims[1]*dims[2], -1);
    for (std::size_t r = 0; r < requests.size(); ++r) {
        const auto& q = requests[r];
        for (int d = 0; d < 3; ++d) {
            if (q.startIJK[d] < 0 || q.endIJK[d] > dims[d] || q.startIJK[d] >= q.endIJK[d]) {
                fail("request " + std::to_string(r + 1) + " is outside the grid");
            }
            const int len = q.endIJK[d] - q.startIJK[d];
            if (q.cellsPerDim[d] < 1 || q.cellsPerDim[d] > len) {
                fail("request " + std::to_string(r + 1) + ": cannot make "
                     + std::to_string(q.cellsPerDim[d]) + " cells from "
                     + std::to_string(len));
            }
        }
        for (int k = q.startIJK[2]; k < q.endIJK[2]; ++k) {
            for (int j = q.startIJK[1]; j < q.endIJK[1]; ++j) {
                for (int i = q.startIJK[0]; i < q.endIJK[0]; ++i) {
                    int& o = owner[cellIndex(dims, i, j, k)];
                    if (o >= 0) {
                        fail("requests " + std::to_string(o + 1) + " and "
                             + std::to_string(r + 1) + " overlap at " + ijk(i, j, k));
                    }
                    o = static_cast<int>(r);
                }
            }
        }
    }
}

/// Lateral grouping must be the same in every column: a dropped pillar line is
/// dropped over the whole grid, which is what keeps the result a grdecl.
Grouping buildGrouping(const std::array<int,3>& dims, const std::vector<CoarsenRequest>& requests)
{
    const auto& d = dims;
    Grouping g;

    for (int axis = 0; axis < 2; ++axis) {
        std::vector<int> reference;
        for (int k = 0; k < d[2]; ++k) {
            for (int j = 0; j < d[axis == 0 ? 1 : 0]; ++j) {
                std::array<int,3> other{};
                other[2] = k;
                other[axis == 0 ? 1 : 0] = j;
                const auto sizes = axisGrouping(d[axis], axis, requests,
                                                requestsAlong(d[axis], axis, other, requests));
                if (reference.empty()) {
                    reference = sizes;
                } else if (sizes != reference) {
                    fail(std::string(axis == 0 ? "i" : "j")
                         + "-grouping differs between columns (first difference near "
                         + ijk(0, j, k) + "); a pillar line must be dropped everywhere "
                         "or nowhere for the result to be a corner-point grid");
                }
            }
        }
        (axis == 0 ? g.iSizes : g.jSizes) = reference;
    }

    const auto iStart = groupStarts(g.iSizes);
    const auto jStart = groupStarts(g.jSizes);
    g.kSizes.resize(static_cast<std::size_t>(g.ncx())*g.ncy());

    for (int jc = 0; jc < g.ncy(); ++jc) {
        for (int ic = 0; ic < g.ncx(); ++ic) {
            std::vector<int> reference;
            for (int j = jStart[jc]; j < jStart[jc + 1]; ++j) {
                for (int i = iStart[ic]; i < iStart[ic + 1]; ++i) {
                    const std::array<int,3> other{i, j, 0};
                    const auto sizes = axisGrouping(d[2], 2, requests,
                                                    requestsAlong(d[2], 2, other, requests));
                    if (reference.empty()) {
                        reference = sizes;
                    } else if (sizes != reference) {
                        fail("layer grouping differs inside the coarse column ("
                             + std::to_string(ic + 1) + "," + std::to_string(jc + 1)
                             + "), at " + ijk(i, j, 0));
                    }
                }
            }
            g.kSizes[static_cast<std::size_t>(jc)*g.ncx() + ic] = reference;
        }
    }
    return g;
}

/// A dropped pillar line may not hide a fault, and it must lie on the line
/// between the pillars that survive, or the coarse cell is not the union of
/// its fine cells.
void checkLateralContinuity(const Grdecl& fine, const Grouping& g)
{
    const auto& d = fine.dims;
    const auto iStart = groupStarts(g.iSizes);
    const auto jStart = groupStarts(g.jSizes);

    for (std::size_t gi = 0; gi < g.iSizes.size(); ++gi) {
        for (int I = iStart[gi] + 1; I < iStart[gi + 1]; ++I) {
            for (int k = 0; k < d[2]; ++k) {
                for (int j = 0; j < d[1]; ++j) {
                    for (int dj = 0; dj < 2; ++dj) {
                        for (int dk = 0; dk < 2; ++dk) {
                            if (std::abs(cellZ(fine, I - 1, j, k, 1, dj, dk)
                                         - cellZ(fine, I, j, k, 0, dj, dk)) > Z_TOL) {
                                fail("fault across the dropped i-pillar line " + std::to_string(I + 1)
                                     + " at " + ijk(I, j, k));
                            }
                        }
                    }
                }
            }
            for (int J = 0; J <= d[1]; ++J) {
                const double* a = pillar(fine, iStart[gi], J);
                const double* b = pillar(fine, I, J);
                const double* c = pillar(fine, iStart[gi + 1], J);
                for (int end = 0; end < 2; ++end) {
                    const int o = 3*end;
                    const double dx = c[o] - a[o], dy = c[o + 1] - a[o + 1];
                    const double len = std::hypot(dx, dy);
                    const double cross = (len < 1e-12) ? std::hypot(b[o] - a[o], b[o + 1] - a[o + 1])
                        : std::abs(dx*(b[o + 1] - a[o + 1]) - dy*(b[o] - a[o]))/len;
                    if (cross > PILLAR_TOL) {
                        fail("dropped i-pillar line " + std::to_string(I + 1) + " at J="
                             + std::to_string(J + 1) + " is not on the line between the "
                             "retained pillars; the coarse cell would not be the union "
                             "of its fine cells");
                    }
                }
            }
        }
    }

    for (std::size_t gj = 0; gj < g.jSizes.size(); ++gj) {
        for (int J = jStart[gj] + 1; J < jStart[gj + 1]; ++J) {
            for (int k = 0; k < d[2]; ++k) {
                for (int i = 0; i < d[0]; ++i) {
                    for (int di = 0; di < 2; ++di) {
                        for (int dk = 0; dk < 2; ++dk) {
                            if (std::abs(cellZ(fine, i, J - 1, k, di, 1, dk)
                                         - cellZ(fine, i, J, k, di, 0, dk)) > Z_TOL) {
                                fail("fault across the dropped j-pillar line " + std::to_string(J + 1)
                                     + " at " + ijk(i, J, k));
                            }
                        }
                    }
                }
            }
            for (int I = 0; I <= d[0]; ++I) {
                const double* a = pillar(fine, I, jStart[gj]);
                const double* b = pillar(fine, I, J);
                const double* c = pillar(fine, I, jStart[gj + 1]);
                for (int end = 0; end < 2; ++end) {
                    const int o = 3*end;
                    const double dx = c[o] - a[o], dy = c[o + 1] - a[o + 1];
                    const double len = std::hypot(dx, dy);
                    const double cross = (len < 1e-12) ? std::hypot(b[o] - a[o], b[o + 1] - a[o + 1])
                        : std::abs(dx*(b[o + 1] - a[o + 1]) - dy*(b[o] - a[o]))/len;
                    if (cross > PILLAR_TOL) {
                        fail("dropped j-pillar line " + std::to_string(J + 1) + " at I="
                             + std::to_string(I + 1) + " is not on the line between the "
                             "retained pillars");
                    }
                }
            }
        }
    }
}

/// Merged layers must touch, unless gaps are explicitly allowed.
void checkVerticalGaps(const Grdecl& fine, const Grouping& g, const Options& options,
                       Report& report)
{
    const auto iStart = groupStarts(g.iSizes);
    const auto jStart = groupStarts(g.jSizes);

    for (int jc = 0; jc < g.ncy(); ++jc) {
        for (int ic = 0; ic < g.ncx(); ++ic) {
            const auto kStart = groupStarts(g.kSizes[static_cast<std::size_t>(jc)*g.ncx() + ic]);
            for (std::size_t gk = 0; gk + 1 < kStart.size(); ++gk) {
                for (int k = kStart[gk]; k + 1 < kStart[gk + 1]; ++k) {
                    for (int j = jStart[jc]; j < jStart[jc + 1]; ++j) {
                        for (int i = iStart[ic]; i < iStart[ic + 1]; ++i) {
                            double gap = 0.0;
                            for (int dj = 0; dj < 2; ++dj) {
                                for (int di = 0; di < 2; ++di) {
                                    gap = std::max(gap,
                                                   std::abs(cellZ(fine, i, j, k, di, dj, 1)
                                                            - cellZ(fine, i, j, k + 1, di, dj, 0)));
                                }
                            }
                            if (gap > Z_TOL) {
                                if (!options.allowVerticalGaps) {
                                    fail("gap of " + std::to_string(gap) + " m inside a merged "
                                         "column at " + ijk(i, j, k) + "; set allowVerticalGaps "
                                         "to absorb it as rock");
                                }
                                // gap times the cell's footprint
                                const auto c = cellCorners(fine, i, j, k);
                                report.absorbedGapVolume += gap*hexVolume(c)
                                    / std::max(1e-12, std::abs(c[4][2] - c[0][2]));
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Neighbouring coarse columns: how many fine faces sit behind one coarse
/// face, and whether the two layer groupings nest.
void checkColumnGrading(const Grouping& g, const Options& options, Report& report)
{
    const auto boundaries = [&](int ic, int jc) {
        return groupStarts(g.kSizes[static_cast<std::size_t>(jc)*g.ncx() + ic]);
    };

    for (int jc = 0; jc < g.ncy(); ++jc) {
        for (int ic = 0; ic < g.ncx(); ++ic) {
            const auto a = boundaries(ic, jc);
            for (int dir = 0; dir < 2; ++dir) {
                const int ni = ic + (dir == 0 ? 1 : 0);
                const int nj = jc + (dir == 1 ? 1 : 0);
                if (ni >= g.ncx() || nj >= g.ncy()) {
                    continue;
                }
                const auto b = boundaries(ni, nj);
                for (std::size_t m = 0; m + 1 < a.size(); ++m) {
                    int sub = 0;
                    for (std::size_t n = 0; n + 1 < b.size(); ++n) {
                        if (b[n] < a[m + 1] && b[n + 1] > a[m]) {
                            ++sub;
                        }
                    }
                    report.maxSubFacesSeen = std::max(report.maxSubFacesSeen, sub);
                    if (sub > options.maxSubFaces) {
                        fail("coarse cell (" + std::to_string(ic + 1) + ","
                             + std::to_string(jc + 1) + "," + std::to_string(m + 1)
                             + ") has " + std::to_string(sub) + " fine faces behind one "
                             "face, over the limit of " + std::to_string(options.maxSubFaces));
                    }
                }
                if (options.requireGradedColumns) {
                    const bool aInB = std::includes(b.begin(), b.end(), a.begin(), a.end());
                    const bool bInA = std::includes(a.begin(), a.end(), b.begin(), b.end());
                    if (!aInB && !bInA) {
                        fail("layer groupings of neighbouring coarse columns ("
                             + std::to_string(ic + 1) + "," + std::to_string(jc + 1) + ") and ("
                             + std::to_string(ni + 1) + "," + std::to_string(nj + 1)
                             + ") do not nest");
                    }
                }
            }
        }
    }
}

/// Coarse dimensions, the blocks and the fine -> coarse map. No geometry.
void layout(const std::array<int,3>& dims, const Grouping& g,
            std::array<int,3>& coarseDims, std::vector<Block>& blocks,
            std::vector<int>& fineToCoarse)
{
    std::size_t maxGroups = 0;
    for (const auto& ks : g.kSizes) {
        maxGroups = std::max(maxGroups, ks.size());
    }
    coarseDims = {g.ncx(), g.ncy(), static_cast<int>(maxGroups)};

    const auto iStart = groupStarts(g.iSizes);
    const auto jStart = groupStarts(g.jSizes);
    blocks.assign(static_cast<std::size_t>(coarseDims[0])*coarseDims[1]*coarseDims[2], Block{});
    fineToCoarse.assign(static_cast<std::size_t>(dims[0])*dims[1]*dims[2], -1);

    for (int jc = 0; jc < coarseDims[1]; ++jc) {
        for (int ic = 0; ic < coarseDims[0]; ++ic) {
            const auto& kSizes = g.kSizes[static_cast<std::size_t>(jc)*coarseDims[0] + ic];
            const auto kStart = groupStarts(kSizes);
            for (int kc = 0; kc < coarseDims[2]; ++kc) {
                const std::size_t coarse = cellIndex(coarseDims, ic, jc, kc);
                auto& block = blocks[coarse];
                block.coarseIndex = static_cast<int>(coarse);
                if (static_cast<std::size_t>(kc) >= kSizes.size()) {
                    // Column with fewer groups: a collapsed cell at its base.
                    block.startIJK = block.endIJK = {iStart[ic], jStart[jc], dims[2] - 1};
                    continue;
                }
                block.startIJK = {iStart[ic], jStart[jc], kStart[kc]};
                block.endIJK = {iStart[ic + 1], jStart[jc + 1], kStart[kc + 1]};
                for (int k = block.startIJK[2]; k < block.endIJK[2]; ++k) {
                    for (int j = block.startIJK[1]; j < block.endIJK[1]; ++j) {
                        for (int i = block.startIJK[0]; i < block.endIJK[0]; ++i) {
                            fineToCoarse[cellIndex(dims, i, j, k)] = static_cast<int>(coarse);
                        }
                    }
                }
            }
        }
    }
}

} // anonymous namespace

BlockLayout blockLayout(const std::array<int,3>& fineDims,
                        const std::vector<CoarsenRequest>& requests)
{
    checkRequests(fineDims, requests);

    BlockLayout out;
    out.blockOfCartesian.assign(static_cast<std::size_t>(fineDims[0])*fineDims[1]*fineDims[2], -1);

    const auto addBlock = [&out, &fineDims](const std::array<int,6>& box) {
        const int id = static_cast<int>(out.boxes.size());
        out.boxes.push_back(box);
        for (int k = box[2]; k <= box[5]; ++k) {
            for (int j = box[1]; j <= box[4]; ++j) {
                for (int i = box[0]; i <= box[3]; ++i) {
                    out.blockOfCartesian[cellIndex(fineDims, i, j, k)] = id;
                }
            }
        }
    };

    for (const auto& q : requests) {
        std::array<std::vector<int>,3> starts;
        for (int d = 0; d < 3; ++d) {
            const auto sizes = evenSplit(q.endIJK[d] - q.startIJK[d], q.cellsPerDim[d]);
            starts[d] = groupStarts(sizes);
            for (auto& v : starts[d]) {
                v += q.startIJK[d];
            }
        }
        for (std::size_t kc = 0; kc + 1 < starts[2].size(); ++kc) {
            for (std::size_t jc = 0; jc + 1 < starts[1].size(); ++jc) {
                for (std::size_t ic = 0; ic + 1 < starts[0].size(); ++ic) {
                    addBlock({starts[0][ic], starts[1][jc], starts[2][kc],
                              starts[0][ic + 1] - 1, starts[1][jc + 1] - 1,
                              starts[2][kc + 1] - 1});
                }
            }
        }
    }

    // Everything the requests left alone is a block of one cell.
    for (int k = 0; k < fineDims[2]; ++k) {
        for (int j = 0; j < fineDims[1]; ++j) {
            for (int i = 0; i < fineDims[0]; ++i) {
                if (out.blockOfCartesian[cellIndex(fineDims, i, j, k)] < 0) {
                    addBlock({i, j, k, i, j, k});
                }
            }
        }
    }
    return out;
}

CartesianMap cartesianMap(const std::array<int,3>& fineDims,
                          const std::vector<CoarsenRequest>& requests)
{
    checkRequests(fineDims, requests);
    const Grouping g = buildGrouping(fineDims, requests);

    CartesianMap out;
    std::vector<Block> blocks;
    layout(fineDims, g, out.coarseDims, blocks, out.fineToCoarse);
    return out;
}

Result coarsenCornerPoint(const Grdecl& fine,
                          const std::vector<CoarsenRequest>& requests,
                          const Options& options)
{
    const auto& d = fine.dims;
    const std::size_t nCells = static_cast<std::size_t>(d[0])*d[1]*d[2];
    if (fine.coord.size() != 6*static_cast<std::size_t>(d[0] + 1)*(d[1] + 1)
        || fine.zcorn.size() != 8*nCells
        || (!fine.actnum.empty() && fine.actnum.size() != nCells)) {
        fail("input arrays do not match the given dimensions");
    }

    checkRequests(d, requests);
    const Grouping g = buildGrouping(d, requests);

    Result result;
    checkLateralContinuity(fine, g);
    checkVerticalGaps(fine, g, options, result.report);
    checkColumnGrading(g, options, result.report);

    const auto iStart = groupStarts(g.iSizes);
    const auto jStart = groupStarts(g.jSizes);

    auto& out = result.grid;
    layout(d, g, out.dims, result.blocks, result.fineToCoarse);
    const std::size_t nCoarse = result.blocks.size();

    out.coord.resize(6*static_cast<std::size_t>(out.dims[0] + 1)*(out.dims[1] + 1));
    for (int Jc = 0; Jc <= out.dims[1]; ++Jc) {
        for (int Ic = 0; Ic <= out.dims[0]; ++Ic) {
            const double* p = pillar(fine, iStart[Ic], jStart[Jc]);
            std::copy(p, p + 6,
                      out.coord.begin() + 6*(static_cast<std::size_t>(Jc)*(out.dims[0] + 1) + Ic));
        }
    }

    out.zcorn.assign(8*nCoarse, 0.0);
    out.actnum.assign(nCoarse, 0);

    for (int kc = 0; kc < out.dims[2]; ++kc) {
        for (int jc = 0; jc < out.dims[1]; ++jc) {
            for (int ic = 0; ic < out.dims[0]; ++ic) {
                const std::size_t coarse = cellIndex(out.dims, ic, jc, kc);
                const auto& block = result.blocks[coarse];
                const bool padding = block.endIJK[2] == block.startIJK[2];
                const int i1 = block.startIJK[0], i2 = block.endIJK[0] - 1;
                const int j1 = block.startIJK[1], j2 = block.endIJK[1] - 1;
                const int k1 = block.startIJK[2];
                const int k2 = padding ? k1 : block.endIJK[2] - 1;

                for (int dk = 0; dk < 2; ++dk) {
                    // A padded cell is collapsed onto the column's base.
                    const int kSrc = (padding || dk == 1) ? k2 : k1;
                    const int dkSrc = (padding || dk == 1) ? 1 : 0;
                    for (int dj = 0; dj < 2; ++dj) {
                        for (int di = 0; di < 2; ++di) {
                            out.zcorn[cornerIndex(out.dims, 2*ic + di, 2*jc + dj, 2*kc + dk)]
                                = cellZ(fine, di ? std::max(i1, i2) : i1,
                                        dj ? std::max(j1, j2) : j1, kSrc, di, dj, dkSrc);
                        }
                    }
                }
                if (padding) {
                    continue;
                }

                bool anyActive = false, allActive = true;
                double inactiveVolume = 0.0;
                int inactiveCells = 0;
                for (int k = k1; k <= k2; ++k) {
                    for (int j = j1; j <= j2; ++j) {
                        for (int i = i1; i <= i2; ++i) {
                            const std::size_t f = cellIndex(d, i, j, k);
                            const bool active = fine.actnum.empty() || fine.actnum[f] != 0;
                            anyActive = anyActive || active;
                            allActive = allActive && active;
                            if (!active) {
                                ++inactiveCells;
                                inactiveVolume += hexVolume(cellCorners(fine, i, j, k));
                            }
                        }
                    }
                }

                const auto corners = cellCorners(out, ic, jc, kc);
                const bool hasVolume = hexVolume(corners) > 0.0;
                switch (options.activity) {
                case Activity::AnyChildActive:   out.actnum[coarse] = anyActive ? 1 : 0; break;
                case Activity::AllChildrenActive: out.actnum[coarse] = allActive ? 1 : 0; break;
                case Activity::FillHoles:        out.actnum[coarse] = hasVolume ? 1 : 0; break;
                }
                if (out.actnum[coarse] != 0 && inactiveCells > 0) {
                    result.report.holesFilled += inactiveCells;
                    result.report.holeVolume += inactiveVolume;
                }

                for (int dk = 0; dk < 2; ++dk) {
                    result.report.maxFaceNonPlanarity
                        = std::max(result.report.maxFaceNonPlanarity,
                                   nonPlanarity(corners[4*dk], corners[4*dk + 1],
                                                corners[4*dk + 3], corners[4*dk + 2]));
                }
            }
        }
    }

    std::ostringstream os;
    os << "coarsened " << nCells << " cells to " << nCoarse
       << " (" << out.dims[0] << 'x' << out.dims[1] << 'x' << out.dims[2] << ')';
    result.report.notes.push_back(os.str());
    if (result.report.holesFilled > 0) {
        result.report.notes.push_back(
            std::to_string(result.report.holesFilled) + " inactive fine cells ("
            + std::to_string(result.report.holeVolume) + " m3) are inside active coarse cells");
    }
    return result;
}

std::vector<Refinement::BlockRefinement> refinementBackToFine(const Result& result)
{
    // Child counts per coarse cell; a padded (collapsed) cell has none.
    const auto& d = result.grid.dims;
    const std::size_t n = static_cast<std::size_t>(d[0])*d[1]*d[2];
    std::vector<std::array<int,3>> counts(n, {0,0,0});
    for (std::size_t c = 0; c < n; ++c) {
        const auto& b = result.blocks[c];
        for (int x = 0; x < 3; ++x) {
            counts[c][x] = b.endIJK[x] - b.startIJK[x];
        }
    }

    std::vector<char> taken(n, 0);
    std::vector<Refinement::BlockRefinement> out;
    const auto at = [&](int i, int j, int k) { return cellIndex(d, i, j, k); };

    for (int k = 0; k < d[2]; ++k) {
        for (int j = 0; j < d[1]; ++j) {
            for (int i = 0; i < d[0]; ++i) {
                const std::size_t c = at(i, j, k);
                const auto cnt = counts[c];
                if (taken[c] || cnt == std::array<int,3>{1,1,1} || cnt[0]*cnt[1]*cnt[2] == 0) {
                    continue;
                }
                int iEnd = i;
                while (iEnd + 1 < d[0] && !taken[at(iEnd + 1, j, k)]
                       && counts[at(iEnd + 1, j, k)] == cnt) {
                    ++iEnd;
                }
                int jEnd = j;
                while (jEnd + 1 < d[1]) {
                    bool ok = true;
                    for (int x = i; x <= iEnd && ok; ++x) {
                        ok = !taken[at(x, jEnd + 1, k)] && counts[at(x, jEnd + 1, k)] == cnt;
                    }
                    if (!ok) {
                        break;
                    }
                    ++jEnd;
                }
                int kEnd = k;
                while (kEnd + 1 < d[2]) {
                    bool ok = true;
                    for (int y = j; y <= jEnd && ok; ++y) {
                        for (int x = i; x <= iEnd && ok; ++x) {
                            ok = !taken[at(x, y, kEnd + 1)] && counts[at(x, y, kEnd + 1)] == cnt;
                        }
                    }
                    if (!ok) {
                        break;
                    }
                    ++kEnd;
                }
                for (int z = k; z <= kEnd; ++z) {
                    for (int y = j; y <= jEnd; ++y) {
                        for (int x = i; x <= iEnd; ++x) {
                            taken[at(x, y, z)] = 1;
                        }
                    }
                }

                Refinement::BlockRefinement req;
                req.name = "COARSEN" + std::to_string(out.size() + 1);
                req.cellsPerDim = cnt;
                req.startIJK = {i, j, k};
                req.endIJK = {iEnd + 1, jEnd + 1, kEnd + 1};
                out.push_back(std::move(req));
            }
        }
    }
    return out;
}

std::vector<int> childActnum(const Grdecl& fine, const Block& block)
{
    std::array<int,3> n{};
    for (int x = 0; x < 3; ++x) {
        n[x] = block.endIJK[x] - block.startIJK[x];
    }
    std::vector<int> out(static_cast<std::size_t>(n[0])*n[1]*n[2], 1);
    if (fine.actnum.empty()) {
        return out;
    }
    for (int k = 0; k < n[2]; ++k) {
        for (int j = 0; j < n[1]; ++j) {
            for (int i = 0; i < n[0]; ++i) {
                out[cellIndex(n, i, j, k)]
                    = fine.actnum[cellIndex(fine.dims, block.startIJK[0] + i,
                                            block.startIJK[1] + j, block.startIJK[2] + k)];
            }
        }
    }
    return out;
}

} // namespace Coarsening
} // namespace Opm
