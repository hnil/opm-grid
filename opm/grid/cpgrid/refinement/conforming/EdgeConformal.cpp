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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/cpgrid/refinement/conforming/EdgeConformal.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/common/utility/SparseTable.hpp>

#include <dune/common/fvector.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>

namespace Opm
{
namespace Refinement
{

namespace
{

using Coord = Dune::FieldVector<double,3>;

// The cell's faces other than f run from a to b over edges each uses alone:
// those nodes hang on f's edge (a,b). Only the order along a->b is checked,
// and the path must be unique.
std::vector<int> hangingPath(int a, int b, int f,
                             const std::multimap<int, std::pair<int,int>>& single,
                             const std::vector<Coord>& P)
{
    Coord d = P[b];
    d -= P[a];
    const double len2 = d.two_norm2();
    const auto param = [&](int c) {
        Coord ac = P[c];
        ac -= P[a];
        return (ac * d) / len2;
    };
    std::vector<int> path;
    std::vector<int> found;
    int complete = 0;
    const auto search = [&](auto&& self, int cur, double tCur) -> void {
        const auto [lo, hi] = single.equal_range(cur);
        for (auto it = lo; it != hi && complete < 2; ++it) {
            const auto [other, face] = it->second;
            if (face == f) {
                continue;
            }
            if (other == b) {
                if (!path.empty()) {
                    ++complete;
                    found = path;
                }
                continue;
            }
            const double t = param(other);
            if (t > tCur && t < 1.0) {
                path.push_back(other);
                self(self, other, t);
                path.pop_back();
            }
        }
    };
    search(search, a, 0.0);
    return (complete == 1) ? found : std::vector<int>{};
}

// Insert nodes that hang on a cell's face edge, found from the cell's other
// faces rather than by distance; repeated until every cell sharing the edge has them.
void closeCellSurfaces(Dune::cpgrid::CpGridData& leaf, const std::vector<Coord>& P)
{
    auto& faceToPoint = GridStateWriter::faceToPoint(leaf);
    const auto& cellToFace = GridStateWriter::cellToFace(leaf);
    const auto& faceToCell = GridStateWriter::faceToCell(leaf);
    const int numCells = cellToFace.size();

    std::vector<std::vector<int>> rows(faceToPoint.size());
    for (int face = 0; face < faceToPoint.size(); ++face) {
        const auto row = faceToPoint[face];
        rows[face].assign(row.begin(), row.end());
    }
    bool changed = false;
    std::vector<int> work;
    const auto revisit = [&](int face) {
        const auto fc = faceToCell[Dune::cpgrid::EntityRep<1>(face, true)];
        for (int q = 0; q < fc.size(); ++q) {
            work.push_back(fc[q].index());
        }
        changed = true;
    };

    // One insertion per visit; the cells on the changed face are visited again.
    const auto insertHanging = [&](int cell) {
        const auto cf = cellToFace[Dune::cpgrid::EntityRep<0>(cell, true)];
        std::map<std::pair<int,int>, std::vector<int>> uses; // undirected edge -> faces
        for (int e = 0; e < cf.size(); ++e) {
            const int face = cf[e].index();
            const auto& row = rows[face];
            for (std::size_t n = 0; n < row.size(); ++n) {
                const int u = row[n];
                const int v = row[(n + 1) % row.size()];
                uses[{std::min(u, v), std::max(u, v)}].push_back(face);
            }
        }
        std::multimap<int, std::pair<int,int>> single; // node -> (other node, face)
        for (const auto& [edge, faces] : uses) {
            if (faces.size() == 1) {
                single.emplace(edge.first, std::pair{edge.second, faces[0]});
                single.emplace(edge.second, std::pair{edge.first, faces[0]});
            }
        }
        for (const auto& [edge, faces] : uses) {
            if (faces.size() != 1) {
                continue;
            }
            const auto [a, b] = edge;
            const auto path = hangingPath(a, b, faces[0], single, P);
            if (path.empty()) {
                continue;
            }
            auto& row = rows[faces[0]];
            const int n = static_cast<int>(row.size());
            for (int k = 0; k < n; ++k) {
                if (row[k] == a && row[(k + 1) % n] == b) {
                    row.insert(row.begin() + k + 1, path.begin(), path.end());
                    break;
                }
                if (row[k] == b && row[(k + 1) % n] == a) {
                    row.insert(row.begin() + k + 1, path.rbegin(), path.rend());
                    break;
                }
            }
            revisit(faces[0]);
            return;
        }
    };

    // A node in only one face of every cell using it is no vertex of a closed
    // cell: a fault split left over from a dropped parent face.
    const auto removeOrphans = [&]() {
        std::map<int, bool> orphan;
        for (int cell = 0; cell < numCells; ++cell) {
            const auto cf = cellToFace[Dune::cpgrid::EntityRep<0>(cell, true)];
            std::map<int, int> faceCount;
            for (int e = 0; e < cf.size(); ++e) {
                for (const int node : rows[cf[e].index()]) {
                    ++faceCount[node];
                }
            }
            for (const auto& [node, count] : faceCount) {
                const auto [it, fresh] = orphan.try_emplace(node, count == 1);
                it->second = it->second && (count == 1);
            }
        }
        const auto isOrphan = [&](int node) { return orphan[node]; };
        for (int face = 0; face < static_cast<int>(rows.size()); ++face) {
            auto& row = rows[face];
            const auto drop = std::count_if(row.begin(), row.end(), isOrphan);
            if (drop > 0 && drop + 3 <= static_cast<long>(row.size())) {
                std::erase_if(row, isOrphan);
                revisit(face);
            }
        }
    };

    for (int cell = 0; cell < numCells; ++cell) {
        work.push_back(cell);
    }
    while (!work.empty()) {
        while (!work.empty()) {
            std::sort(work.begin(), work.end());
            work.erase(std::unique(work.begin(), work.end()), work.end());
            std::vector<int> visit;
            visit.swap(work);
            for (const int cell : visit) {
                insertHanging(cell);
            }
        }
        removeOrphans();
    }

    if (changed) {
        Opm::SparseTable<int> rebuilt;
        for (const auto& row : rows) {
            rebuilt.appendRow(row.begin(), row.end());
        }
        faceToPoint.swap(rebuilt);
    }
}

} // anonymous namespace

void edgeConformalizeLeaf(Dune::cpgrid::CpGridData& leaf)
{
    using Coord = Dune::FieldVector<double,3>;

    const auto& cornerGeom =
        GridStateWriter::geometry(leaf).geomVector(std::integral_constant<int,3>());
    const int numCorners = static_cast<int>(cornerGeom.size());
    if (numCorners == 0) {
        return;
    }

    std::vector<Coord> P(numCorners);
    Coord lo = cornerGeom.get(0).center();
    Coord hi = lo;
    for (int i = 0; i < numCorners; ++i) {
        P[i] = cornerGeom.get(i).center();
        for (int c = 0; c < 3; ++c) {
            lo[c] = std::min(lo[c], P[i][c]);
            hi[c] = std::max(hi[c], P[i][c]);
        }
    }
    double diag = 0.0;
    for (int c = 0; c < 3; ++c) {
        diag += (hi[c] - lo[c]) * (hi[c] - lo[c]);
    }
    diag = std::sqrt(diag);
    if (diag <= 0.0) {
        return;
    }

    // Spatial hash so each edge only tests nearby nodes.
    const double h = diag / 256.0 + 1e-30;
    const double tol = 1e-9 * diag;
    using Bucket = std::tuple<long long, long long, long long>;
    const auto bucketOf = [&](const Coord& p) {
        return Bucket{ std::llround(std::floor((p[0] - lo[0]) / h)),
                       std::llround(std::floor((p[1] - lo[1]) / h)),
                       std::llround(std::floor((p[2] - lo[2]) / h)) };
    };
    std::map<Bucket, std::vector<int>> hash;
    for (int i = 0; i < numCorners; ++i) {
        hash[bucketOf(P[i])].push_back(i);
    }

    auto& faceToPoint = GridStateWriter::faceToPoint(leaf);
    Opm::SparseTable<int> rebuilt;
    std::vector<int> newRow;
    std::vector<std::pair<double,int>> mids;

    for (int face = 0; face < faceToPoint.size(); ++face) {
        auto row = faceToPoint[face];
        const int n = row.size();
        newRow.clear();
        for (int e = 0; e < n; ++e) {
            const int a = row[e];
            const int b = row[(e + 1) % n];
            newRow.push_back(a);

            Coord d = P[b];
            d -= P[a];
            const double len2 = d.two_norm2();
            if (len2 <= tol * tol) {
                continue;
            }

            // On slanted pillars refined nodes sit ~1e-8 of the edge length off
            // the straight coarse edge.
            const double lineTol = std::max(tol, 1e-7 * std::sqrt(len2));
            // Candidate nodes: all buckets overlapping the edge's bounding box.
            const long long x0 = std::llround(std::floor((std::min(P[a][0], P[b][0]) - lo[0] - lineTol) / h));
            const long long x1 = std::llround(std::floor((std::max(P[a][0], P[b][0]) - lo[0] + lineTol) / h));
            const long long y0 = std::llround(std::floor((std::min(P[a][1], P[b][1]) - lo[1] - lineTol) / h));
            const long long y1 = std::llround(std::floor((std::max(P[a][1], P[b][1]) - lo[1] + lineTol) / h));
            const long long z0 = std::llround(std::floor((std::min(P[a][2], P[b][2]) - lo[2] - lineTol) / h));
            const long long z1 = std::llround(std::floor((std::max(P[a][2], P[b][2]) - lo[2] + lineTol) / h));

            mids.clear();
            for (long long bx = x0; bx <= x1; ++bx) {
                for (long long by = y0; by <= y1; ++by) {
                    for (long long bz = z0; bz <= z1; ++bz) {
                        const auto it = hash.find({bx, by, bz});
                        if (it == hash.end()) {
                            continue;
                        }
                        for (const int c : it->second) {
                            if (c == a || c == b) {
                                continue;
                            }
                            Coord ac = P[c];
                            ac -= P[a];
                            const double t = (ac * d) / len2;
                            if (t <= 0.0 || t >= 1.0) {
                                continue; // not strictly between the endpoints
                            }
                            // Distance from the node to the line a-b.
                            Coord closest = d;
                            closest *= t;
                            closest += P[a];
                            Coord diff = P[c];
                            diff -= closest;
                            if (diff.two_norm() < lineTol) {
                                mids.emplace_back(t, c);
                            }
                        }
                    }
                }
            }
            std::sort(mids.begin(), mids.end());
            for (const auto& [t, c] : mids) {
                newRow.push_back(c);
            }
        }
        rebuilt.appendRow(newRow.begin(), newRow.end());
    }

    faceToPoint.swap(rebuilt);
    closeCellSurfaces(leaf, P);
}

} // namespace Refinement
} // namespace Opm
