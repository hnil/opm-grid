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
#ifndef OPM_GRID_RETAINED_CORNER_POINT_INPUT_HEADER_INCLUDED
#define OPM_GRID_RETAINED_CORNER_POINT_INPUT_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RetainedCornerPointInput.hpp>

#include <memory>

namespace Opm
{

/// COORD/ZCORN after MINPV and PINCH processing: what a second grid on the
/// same geometry (a mechanics coarsening, say) has to be built from.
using RetainedCornerPointInput = Refinement::RetainedCornerPointInput;

/// Ask grid processing to keep that description. Off by default; processing
/// happens on rank 0, so what is kept here is the whole, undistributed grid.
class RetainCornerPointInput
{
public:
    static void enable(bool on = true) { wanted() = on; }
    static bool enabled() { return wanted(); }

    static void store(std::shared_ptr<const RetainedCornerPointInput> input)
    { held() = std::move(input); }

    /// The description of the last grid processed, or nullptr.
    static const RetainedCornerPointInput* get() { return held().get(); }

private:
    static bool& wanted()
    {
        static bool on = false;
        return on;
    }

    static std::shared_ptr<const RetainedCornerPointInput>& held()
    {
        static std::shared_ptr<const RetainedCornerPointInput> input;
        return input;
    }
};

} // namespace Opm

#endif // OPM_GRID_RETAINED_CORNER_POINT_INPUT_HEADER_INCLUDED
