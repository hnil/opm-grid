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
#ifndef OPM_ADAPTIVE_CPGRID_HEADER_INCLUDED
#define OPM_ADAPTIVE_CPGRID_HEADER_INCLUDED

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Opm
{

/// First cut of the dynamically / locally refinable corner-point grid -- the
/// `AdaptiveCpGrid` heir of the forest-of-octrees design
/// (docs/DESIGN-parallel-octree.md Sec.10-14, docs/REDISTRIBUTION-requirements.md).
///
/// The whole point: a refinement that the static-LGR (CARFIN) path achieves by
/// refining *at grid construction* is here achieved by building the COARSE grid
/// first and refining the marked region *after* construction -- so a deck needs
/// no initial LGR. The refined leaf is bit-for-bit the same grid the static path
/// produces (same builder, same box), hence the discretisation -- and therefore
/// the simulation -- is the same.
///
/// Representation (design Sec.10):
///   * Layer A -- the persistent macro corner-point description (dims + COORD +
///     ZCORN + ACTNUM), kept so refinement can be (re)computed from the macro
///     geometry rather than communicated.
///   * Layer B -- the derived Dune::CpGrid leaf (level zero immediately; refined
///     levels after adapt()).
///
/// This first cut is the *correctness oracle* path (design Sec.12): each adapt()
/// performs a full conforming rebuild through the existing refinement builder
/// (Opm::Refinement::ConformingBlockBuilder), the same one the static CARFIN
/// path drives. The fast in-place local adapt (append + free-list, Sec.12), the
/// factor-2 dynamic re-adapt levels (Sec.13) and the parallel root-tree
/// migration (Sec.5/14) are deliberately deferred; this class establishes the
/// mark -> adapt -> equivalent-leaf seam and proves the static/adaptive
/// equivalence.
///
/// Restrictions (inherited from the oracle builder, throw): serial; "GLOBAL"
/// parent only (no nested). Marked regions of the *same* factor merge and may be
/// adjacent. Two regions of *different* factor become touching boxes; the builder
/// accepts them when their subdivisions match in the two in-face directions (the
/// touch/perpendicular direction is always free, e.g. vertically-stacked boxes
/// may differ in K), OR when they differ by a *compatible* multiple with one box
/// uniformly the finer side on the shared face -- then the coarser side is tiled
/// by a sub-face mosaic (LGR_GAPS A2). Still rejected: incompatible (non-multiple)
/// in-face factors, edge/corner contact, and mixed nesting.
class AdaptiveCpGrid
{
public:
    /// Build a coarse adaptive grid from its macro corner-point description.
    /// The level-zero Dune::CpGrid is constructed immediately; refinement is
    /// applied by adapt(). @p actnum may be empty (all cells active).
    AdaptiveCpGrid(const std::array<int,3>& dims,
                   std::vector<double> coord,
                   std::vector<double> zcorn,
                   std::vector<int> actnum = {});

    /// Mark the half-open Cartesian box [startIJK, endIJK) of level-zero cells
    /// for refinement by @p cellsPerDim sub-cells per parent cell (the
    /// CARFIN-exact first anisotropic split, design Sec.13). @p name defaults to
    /// "ADAPT<n>". This is exactly one CARFIN box, expressed programmatically.
    void markBox(const std::array<int,3>& startIJK,
                 const std::array<int,3>& endIJK,
                 const std::array<int,3>& cellsPerDim,
                 const std::string& name = {});

    /// Mark a single level-zero cell (i,j,k) for refinement by @p cellsPerDim.
    void markCell(const std::array<int,3>& ijk,
                  const std::array<int,3>& cellsPerDim);

    /// Apply the accumulated marks: refine every marked level-zero cell and
    /// rebuild the leaf so it is consistent (parent/child maps, ids, index sets).
    /// Adjacent marked cells with the same refinement factor are merged into
    /// maximal boxes first, so a contiguous region marked cell-by-cell refines
    /// like one CARFIN box (the no-touching-boxes builder rule is satisfied).
    ///
    /// Re-adaptable: marks persist, so calling markCell()/markBox() again then
    /// adapt() refines the union of all marks. This first cut re-adapts via the
    /// full-rebuild oracle (rebuild coarse from Layer A, then refine all marks);
    /// the fast in-place mutation (design Sec.12) is the next step. Throws if a
    /// mark is invalid or the (merged) refinement is unsupported (e.g. two
    /// different-factor regions that end up touching -- gap A2); on throw the grid
    /// is left as the coarse grid.
    void adapt();

    //! \brief Whether adapt() has produced refined levels.
    bool refined() const { return grid_->maxLevel() > 0; }

    //! \brief Number of marked (not-yet-merged) level-zero cells.
    std::size_t markCount() const { return marks_.size(); }

    //! \brief The underlying Dune::CpGrid (coarse before adapt(), refined after).
    Dune::CpGrid&       grid()       { return *grid_; }
    const Dune::CpGrid& grid() const { return *grid_; }

    const std::array<int,3>& dims() const { return dims_; }

private:
    // (Re)build the coarse level-zero grid from Layer A into grid_.
    void buildCoarse_();
    // Greedily merge same-factor marked cells into maximal axis-aligned boxes.
    std::vector<Refinement::BlockRefinement> mergeMarksIntoBoxes_() const;

    // Layer A: persistent macro corner-point description.
    std::array<int,3>   dims_;
    std::vector<double> coord_;
    std::vector<double> zcorn_;
    std::vector<int>    actnum_;

    // Layer B: the derived leaf grid (held by pointer so buildCoarse_ can reset
    // it cleanly for a re-adapt).
    std::unique_ptr<Dune::CpGrid> grid_;

    // Accumulated marks: level-zero cell (i,j,k) -> refinement factor per dir.
    std::map<std::array<int,3>, std::array<int,3>> marks_;
};

} // namespace Opm

#endif // OPM_ADAPTIVE_CPGRID_HEADER_INCLUDED
