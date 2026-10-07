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
#ifndef OPM_CPGRID_GEOMETRIC_CHECK_HPP
#define OPM_CPGRID_GEOMETRIC_CHECK_HPP

#include <dune/common/fvector.hh>

#include <string>
#include <utility>
#include <vector>

namespace Dune { class CpGrid; }

namespace Dune::cpgrid
{

/// Whether a grid is a geometric body, as a mechanics discretisation needs:
/// every cell a closed polyhedron whose faces meet edge to edge, and one
/// closed outer surface, so no voids, cracks or hanging nodes.
struct GeometricCheck
{
    struct Surface
    {
        int faces{0};
        double area{0.0};
        double volume{0.0};          // enclosed, by the divergence theorem
        FieldVector<double,3> centroid = FieldVector<double,3>(0.0);
    };

    int cells{0};
    int faces{0};
    int nncFaces{0};                 // connections without geometry, not part of the body
    int nonPositiveCells{0};
    int openCells{0};                // outward face areas do not sum to zero
    int nonConformingCells{0};       // a face edge no other face of the cell shares
    int misorientedFaces{0};         // node order against the stored normal
    int unpairedBoundaryEdges{0};
    std::vector<Surface> boundaries; // outer surface first, then voids and cracks
    /// The first few unpaired cell edges, as end points (where hanging nodes are).
    std::vector<std::pair<FieldVector<double,3>, FieldVector<double,3>>> unpairedEdges;

    bool ok() const;
    std::string summary(int maxListed = 5) const;
};

/// Checks the current view of the grid.
GeometricCheck checkGeometric(const CpGrid& grid, double tolerance = 1e-9);

} // namespace Dune::cpgrid

#endif // OPM_CPGRID_GEOMETRIC_CHECK_HPP
