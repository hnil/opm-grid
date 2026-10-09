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

#include <opm/grid/cpgrid/AdaptiveCpGrid.hpp>

#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{

AdaptiveCpGrid::AdaptiveCpGrid(const std::array<int,3>& dims,
                               std::vector<double> coord,
                               std::vector<double> zcorn,
                               std::vector<int> actnum)
    : dims_(dims)
    , coord_(std::move(coord))
    , zcorn_(std::move(zcorn))
    , actnum_(std::move(actnum))
{
    buildCoarse_();
}

void AdaptiveCpGrid::buildCoarse_()
{
    // Reset to a fresh coarse (level-zero) leaf built from the macro input.
    grid_ = std::make_unique<Dune::CpGrid>();
    grdecl g{};
    g.dims[0] = dims_[0];
    g.dims[1] = dims_[1];
    g.dims[2] = dims_[2];
    g.coord  = coord_.data();
    g.zcorn  = zcorn_.data();
    g.actnum = actnum_.empty() ? nullptr : actnum_.data();
    grid_->processEclipseFormat(g, false);
    grid_->setLgrBackend(Refinement::Backend::Conforming);
    grid_->setRefinementBuilder(std::make_shared<Refinement::ConformingBlockBuilder>(
        dims_, coord_, zcorn_, actnum_));
}

void AdaptiveCpGrid::markBox(const std::array<int,3>& startIJK,
                             const std::array<int,3>& endIJK,
                             const std::array<int,3>& cellsPerDim,
                             const std::string& /*name*/)
{
    for (int k = startIJK[2]; k < endIJK[2]; ++k)
        for (int j = startIJK[1]; j < endIJK[1]; ++j)
            for (int i = startIJK[0]; i < endIJK[0]; ++i)
                marks_[{i, j, k}] = cellsPerDim;
}

void AdaptiveCpGrid::markCell(const std::array<int,3>& ijk,
                              const std::array<int,3>& cellsPerDim)
{
    marks_[ijk] = cellsPerDim;
}

std::vector<Refinement::BlockRefinement>
AdaptiveCpGrid::mergeMarksIntoBoxes_() const
{
    using Cell = std::array<int,3>;
    using CellSet = std::set<Cell>;

    // Group marked cells by their refinement factor; only same-factor cells can
    // share a box.
    std::map<std::array<int,3>, CellSet> byFactor;
    for (const auto& [cell, fac] : marks_) {
        byFactor[fac].insert(cell);
    }

    auto rowPresent = [](const CellSet& s, int i0, int i1, int j, int k) {
        for (int i = i0; i <= i1; ++i)
            if (s.find({i, j, k}) == s.end()) return false;
        return true;
    };
    auto slabPresent = [&](const CellSet& s, int i0, int i1, int j0, int j1, int k) {
        for (int j = j0; j <= j1; ++j)
            if (!rowPresent(s, i0, i1, j, k)) return false;
        return true;
    };

    std::vector<Refinement::BlockRefinement> boxes;
    int idx = 0;
    for (const auto& [fac, cellsConst] : byFactor) {
        CellSet cells = cellsConst;
        while (!cells.empty()) {
            const Cell seed = *cells.begin();
            const int i0 = seed[0], j0 = seed[1], k0 = seed[2];
            // Greedily grow a maximal box: extend in i, then j, then k while the
            // whole new row/slab is marked.
            int i1 = i0;
            while (cells.find({i1 + 1, j0, k0}) != cells.end()) ++i1;
            int j1 = j0;
            while (rowPresent(cells, i0, i1, j1 + 1, k0)) ++j1;
            int k1 = k0;
            while (slabPresent(cells, i0, i1, j0, j1, k1 + 1)) ++k1;

            for (int k = k0; k <= k1; ++k)
                for (int j = j0; j <= j1; ++j)
                    for (int i = i0; i <= i1; ++i)
                        cells.erase({i, j, k});

            Refinement::BlockRefinement r;
            r.name = "ADAPT" + std::to_string(++idx);
            r.parentGridName = "GLOBAL";
            r.cellsPerDim = fac;
            r.startIJK = {i0, j0, k0};
            r.endIJK   = {i1 + 1, j1 + 1, k1 + 1};
            boxes.push_back(std::move(r));
        }
    }
    return boxes;
}

void AdaptiveCpGrid::adapt()
{
    if (marks_.empty()) {
        return;
    }
    const auto boxes = mergeMarksIntoBoxes_();

    // Full-rebuild oracle: refine the union of all (merged) marks from the coarse
    // grid. On a re-adapt the grid is already refined, so reset it to coarse
    // first (rebuild from Layer A); on the first adapt the ctor's coarse grid is
    // reused (no rebuild). Either way marks persist, so a later markCell()/
    // markBox() + adapt() refines the union. The fast in-place mutation (design
    // Sec.12) replaces this rebuild and is the next step.
    if (grid_->maxLevel() > 0) {
        buildCoarse_();
    }

    std::vector<std::array<int,3>> cellsPerDim, startIJK, endIJK;
    std::vector<std::string>       names;
    cellsPerDim.reserve(boxes.size());
    startIJK.reserve(boxes.size());
    endIJK.reserve(boxes.size());
    names.reserve(boxes.size());
    for (const auto& b : boxes) {
        cellsPerDim.push_back(b.cellsPerDim);
        startIJK.push_back(b.startIJK);
        endIJK.push_back(b.endIJK);
        names.push_back(b.name);
    }

    grid_->addLgrsUpdateLeafView(cellsPerDim, startIJK, endIJK, names);
}

} // namespace Opm
