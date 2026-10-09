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
#ifndef OPM_GRID_CORNER_POINT_COARSENING_HEADER_INCLUDED
#define OPM_GRID_CORNER_POINT_COARSENING_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace Opm
{
namespace Coarsening
{

/// Corner-point (grdecl) description, same layout and conventions as
/// Refinement::RefinedBlockGrdecl.
struct Grdecl
{
    std::array<int,3> dims{};
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;
};

/// One COARSEN record: the box becomes cellsPerDim coarse cells. IJK are
/// 0-based and the box is half-open, as in Refinement::BlockRefinement.
struct CoarsenRequest
{
    std::array<int,3> startIJK{};
    std::array<int,3> endIJK{};
    std::array<int,3> cellsPerDim{1,1,1};
};

/// Which coarse cells are active.
enum class Activity
{
    AnyChildActive,    ///< flow: active where at least one fine cell is
    AllChildrenActive, ///< strict: active only where every fine cell is
    FillHoles,         ///< mechanics: active wherever the block has volume
};

struct Options
{
    Activity activity{Activity::AnyChildActive};
    /// Merge layers with a gap between them (the gap becomes part of the
    /// coarse cell). Rock-filling for mechanics; a volume change for flow.
    bool allowVerticalGaps{false};
    /// Cap on the fine faces behind one coarse face; too many nodes on a cell
    /// destabilises VEM.
    int maxSubFaces{16};
    /// Neighbouring columns' layer groupings must nest.
    bool requireGradedColumns{true};
};

/// The index side of a coarsening: no geometry, so it can be computed wherever
/// only the dimensions are at hand (a load balancer, for instance).
struct CartesianMap
{
    std::array<int,3> coarseDims{};
    /// Fine Cartesian index -> coarse Cartesian index.
    std::vector<int> fineToCoarse;
};

/// The coarse dimensions and the fine -> coarse Cartesian map the requests
/// imply. Validates the requests as coarsenCornerPoint() does, apart from the
/// checks that need geometry.
CartesianMap cartesianMap(const std::array<int,3>& fineDims,
                          const std::vector<CoarsenRequest>& requests);

/// The blocks the requests define, without asking whether they can be written
/// as a corner-point description: what the topological merge needs.
struct BlockLayout
{
    /// Block of each fine Cartesian cell.
    std::vector<int> blockOfCartesian;
    /// Each block's box as {i1, j1, k1, i2, j2, k2}, inclusive.
    std::vector<std::array<int,6>> boxes;
};

/// Split the requested boxes into blocks. Every cell belongs to exactly one:
/// cells no request covers are blocks of their own.
BlockLayout blockLayout(const std::array<int,3>& fineDims,
                        const std::vector<CoarsenRequest>& requests);

/// One coarse cell and the fine box it covers (half-open, fine indices).
struct Block
{
    std::array<int,3> startIJK{};
    std::array<int,3> endIJK{};
    int coarseIndex{};
};

struct Report
{
    int holesFilled{};
    double holeVolume{};
    double absorbedGapVolume{};
    int maxSubFacesSeen{1};
    double maxFaceNonPlanarity{};
    std::vector<std::string> notes;
};

struct Result
{
    Grdecl grid;
    /// One per coarse cell, in coarse Cartesian order.
    std::vector<Block> blocks;
    /// Fine Cartesian index -> coarse Cartesian index.
    std::vector<int> fineToCoarse;
    Report report;
};

/// Merge boxes of a corner-point grid into single cells.
///
/// The result is a corner-point grid whose cells are unions of the input's:
/// coarse pillars are input pillars, and a coarse cell's corners are the
/// corners of the input cells at the block's corners. The fine grid is then a
/// (graded) refinement of it — see refinementBackToFine().
///
/// Throws std::invalid_argument when the requests cannot produce a
/// corner-point grid: overlapping boxes, a pillar line dropped in one region
/// and kept in another, coarsening across a fault, a gap inside a merged
/// column when allowVerticalGaps is false, or a violated face-count or
/// grading limit. The message names the offending indices.
Result coarsenCornerPoint(const Grdecl& fine,
                          const std::vector<CoarsenRequest>& requests,
                          const Options& options = {});

/// The refinement that takes the coarse grid back to the fine one, as
/// CARFIN-style requests: one per maximal box of coarse cells with the same
/// child counts. Together they cover every coarsened cell, so the fine grid
/// can be viewed as an LGR of the coarse grid.
std::vector<Refinement::BlockRefinement>
refinementBackToFine(const Result& result);

/// The fine grid's own ACTNUM for one block, in refined Cartesian order.
/// Needed because a coarse cell may be active while some children are not.
std::vector<int> childActnum(const Grdecl& fine, const Block& block);

} // namespace Coarsening
} // namespace Opm

#endif // OPM_GRID_CORNER_POINT_COARSENING_HEADER_INCLUDED
