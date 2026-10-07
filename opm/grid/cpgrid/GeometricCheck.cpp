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
#include <config.h>

#include <opm/grid/cpgrid/GeometricCheck.hpp>

#include <opm/grid/CpGrid.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <sstream>
#include <unordered_map>

namespace Dune::cpgrid
{

namespace
{
using Point = FieldVector<double,3>;

std::uint64_t edgeKey(int a, int b)
{
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32)
        | static_cast<std::uint32_t>(b);
}

// The face's nodes in the order that runs with its stored normal, and its
// vector area (Newell), which depends only on that loop.
struct OrientedFace
{
    std::vector<int> nodes;
    Point area = Point(0.0);
    bool misoriented{false};
};

OrientedFace orientedFace(const CpGrid& grid, int face)
{
    OrientedFace out;
    const int n = grid.numFaceVertices(face);
    out.nodes.resize(n);
    for (int v = 0; v < n; ++v) {
        out.nodes[v] = grid.faceVertex(face, v);
    }
    for (int v = 0; v < n; ++v) {
        const auto& p = grid.vertexPosition(out.nodes[v]);
        const auto& q = grid.vertexPosition(out.nodes[(v + 1) % n]);
        out.area[0] += 0.5*(p[1]*q[2] - p[2]*q[1]);
        out.area[1] += 0.5*(p[2]*q[0] - p[0]*q[2]);
        out.area[2] += 0.5*(p[0]*q[1] - p[1]*q[0]);
    }
    if (out.area*grid.faceNormal(face) < 0.0) {
        out.misoriented = true;
        std::reverse(out.nodes.begin(), out.nodes.end());
        out.area *= -1.0;
    }
    return out;
}

int find(std::vector<int>& parent, int i)
{
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}
} // namespace

bool GeometricCheck::ok() const
{
    return nonPositiveCells == 0 && openCells == 0
        && nonConformingCells == 0 && misorientedFaces == 0
        && unpairedBoundaryEdges == 0 && boundaries.size() <= 1;
}

std::string GeometricCheck::summary(const int maxListed) const
{
    std::ostringstream os;
    os << cells << " cells, " << faces << " faces";
    if (nncFaces > 0) {
        os << " (" << nncFaces << " of them NNC, without geometry)";
    }
    os << ": ";
    if (ok()) {
        os << "geometric (closed cells, edge-conformal, one closed boundary)";
        return os.str();
    }
    const auto item = [&os](int count, const char* what) {
        if (count > 0) {
            os << count << ' ' << what << "; ";
        }
    };
    item(nonPositiveCells, "cells without positive volume");
    item(openCells, "cells whose faces do not close");
    item(nonConformingCells, "cells with a face edge no other face shares (hanging node)");
    if (!unpairedEdges.empty()) {
        os << "e.g. edge (" << unpairedEdges.front().first << ")-(" << unpairedEdges.front().second
           << "); ";
    }
    item(misorientedFaces, "faces with node order against their normal");
    item(unpairedBoundaryEdges, "unpaired boundary edges");
    if (boundaries.size() > 1) {
        os << boundaries.size() - 1 << " closed surfaces besides the outer one (voids or cracks):";
        for (std::size_t s = 1; s < boundaries.size() && static_cast<int>(s) <= maxListed; ++s) {
            const auto& b = boundaries[s];
            os << " [" << b.faces << " faces, area " << b.area << " m2, volume " << b.volume
               << " m3, at " << b.centroid << ']';
        }
    }
    return os.str();
}

GeometricCheck checkGeometric(const CpGrid& grid, const double tolerance)
{
    GeometricCheck check;
    check.cells = grid.numCells();
    check.faces = grid.numFaces();

    std::vector<OrientedFace> oriented(check.faces);
    for (int f = 0; f < check.faces; ++f) {
        if (grid.numFaceVertices(f) == 0) {
            ++check.nncFaces;
            continue;
        }
        oriented[f] = orientedFace(grid, f);
        check.misorientedFaces += oriented[f].misoriented ? 1 : 0;
    }

    // Each cell: positive volume, closed, and every directed edge of its
    // outward faces met once by the reverse edge of another.
    std::unordered_map<std::uint64_t, int> edges;
    for (int c = 0; c < check.cells; ++c) {
        if (!(grid.cellVolume(c) > 0.0)) {
            ++check.nonPositiveCells;
        }
        edges.clear();
        Point closure(0.0);
        double surface = 0.0;
        for (int lf = 0; lf < grid.numCellFaces(c); ++lf) {
            const int f = grid.cellFace(c, lf);
            if (oriented[f].nodes.empty()) {
                continue;
            }
            const bool outward = grid.faceCell(f, 0) == c;
            const auto& nodes = oriented[f].nodes;
            const int n = static_cast<int>(nodes.size());
            for (int v = 0; v < n; ++v) {
                const int a = nodes[v], b = nodes[(v + 1) % n];
                ++edges[outward ? edgeKey(a, b) : edgeKey(b, a)];
            }
            Point area = oriented[f].area;
            area *= outward ? 1.0 : -1.0;
            closure += area;
            surface += area.two_norm();
        }
        if (closure.two_norm() > tolerance*surface) {
            ++check.openCells;
        }
        bool conforming = true;
        for (const auto& [key, count] : edges) {
            const auto reverse = edges.find((key << 32) | (key >> 32));
            if (count != 1 || reverse == edges.end() || reverse->second != 1) {
                conforming = false;
                if (check.unpairedEdges.size() < 10) {
                    check.unpairedEdges.emplace_back(grid.vertexPosition(static_cast<int>(key >> 32)),
                                                     grid.vertexPosition(static_cast<int>(key & 0xffffffffu)));
                }
            }
        }
        check.nonConformingCells += conforming ? 0 : 1;
    }

    // The boundary: outward faces with one cell, joined into surfaces by shared edges.
    std::vector<int> boundaryFaces;
    for (int f = 0; f < check.faces; ++f) {
        if (!oriented[f].nodes.empty() && (grid.faceCell(f, 0) < 0 || grid.faceCell(f, 1) < 0)) {
            boundaryFaces.push_back(f);
        }
    }
    std::vector<int> parent(boundaryFaces.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::unordered_map<std::uint64_t, int> directed, firstOnEdge;
    for (std::size_t i = 0; i < boundaryFaces.size(); ++i) {
        const int f = boundaryFaces[i];
        const bool outward = grid.faceCell(f, 0) >= 0;
        const auto& nodes = oriented[f].nodes;
        const int n = static_cast<int>(nodes.size());
        for (int v = 0; v < n; ++v) {
            const int a = nodes[v], b = nodes[(v + 1) % n];
            ++directed[outward ? edgeKey(a, b) : edgeKey(b, a)];
            const auto [it, isNew] = firstOnEdge.emplace(edgeKey(std::min(a, b), std::max(a, b)),
                                                         static_cast<int>(i));
            if (!isNew) {
                parent[find(parent, static_cast<int>(i))] = find(parent, it->second);
            }
        }
    }
    for (const auto& [key, count] : directed) {
        const auto reverse = directed.find((key << 32) | (key >> 32));
        if (count != 1 || reverse == directed.end() || reverse->second != 1) {
            ++check.unpairedBoundaryEdges;
        }
    }

    Point reference(0.0);
    for (int v = 0; v < grid.numVertices(); ++v) {
        reference += grid.vertexPosition(v);
    }
    if (grid.numVertices() > 0) {
        reference /= grid.numVertices();
    }
    std::unordered_map<int, GeometricCheck::Surface> surfaces;
    for (std::size_t i = 0; i < boundaryFaces.size(); ++i) {
        const int f = boundaryFaces[i];
        Point area = oriented[f].area;
        area *= (grid.faceCell(f, 0) >= 0) ? 1.0 : -1.0;
        auto& s = surfaces[find(parent, static_cast<int>(i))];
        const double a = area.two_norm();
        ++s.faces;
        s.area += a;
        s.volume += (grid.faceCentroid(f) - reference)*area/3.0;
        Point weighted = grid.faceCentroid(f);
        weighted *= a;
        s.centroid += weighted;
    }
    for (auto& [root, s] : surfaces) {
        if (s.area > 0.0) {
            s.centroid /= s.area;
        }
        check.boundaries.push_back(s);
    }
    std::sort(check.boundaries.begin(), check.boundaries.end(),
              [](const auto& a, const auto& b) { return a.volume > b.volume; });
    return check;
}

} // namespace Dune::cpgrid
