/*
  Copyright 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE LgrOnFaultedGridTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <opm/grid/CpGridLGR.hpp>
#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <array>
#include <cmath>
#include <string>
#include <vector>

struct Fixture {
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

// Faulted 4x1x1 grid with COORD/ZCORN.
// Cells 1-2 have top/bottom at 2000/2010, cells 3-4 at 1995/2005.
// This creates a 5m fault between cells 2 and 3.
// The LGR 'LGR1' covers cells I=2 to I=3 (spanning the fault) with 2x1x2 refinement.
const std::string faulted_deck_string = R"(
RUNSPEC
DIMENS
 4 1 1 /
TITLE
  FAULTED EXAMPLE WITH LOCAL GRID REFINEMENT
START
  11 FEB 2026 /
OIL
WATER
METRIC
GRID
COORD
  0.0   0.0 0.0    0.0   0.0 0.0
100.0   0.0 0.0  100.0   0.0 0.0
200.0   0.0 0.0  200.0   0.0 0.0
300.0   0.0 0.0  300.0   0.0 0.0
400.0   0.0 0.0  400.0   0.0 0.0
--
  0.0 100.0 0.0    0.0 100.0 0.0
100.0 100.0 0.0  100.0 100.0 0.0
200.0 100.0 0.0  200.0 100.0 0.0
300.0 100.0 0.0  300.0 100.0 0.0
400.0 100.0 0.0  400.0 100.0 0.0
/
ZCORN
 2000.0  2000.0   2000.0  2000.0    1995.0  1995.0    1995.0  1995.0
 2000.0  2000.0   2000.0  2000.0    1995.0  1995.0    1995.0  1995.0
 2010.0  2010.0   2010.0  2010.0    2005.0  2005.0    2005.0  2005.0
 2010.0  2010.0   2010.0  2010.0    2005.0  2005.0    2005.0  2005.0
/
CARFIN
  'LGR1'  2 3  1 1  1 1  2 1 2 /
ENDFIN
EQUALS
  PORO    0.3 /
  PERMX 200.0 /
  PERMY 200.0 /
  PERMZ  25.0 /
/
PROPS
SOLUTION
SCHEDULE
)";

BOOST_AUTO_TEST_CASE(faulted_grid_diagnose)
{
    // First, diagnose the faulted grid structure before LGR
    Opm::Parser parser;
    const auto deck = parser.parseString(faulted_deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();

    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);

    const auto& data = *grid.currentData().back();
    BOOST_TEST_MESSAGE("Grid: " << data.size(0) << " cells, " << data.numFaces() << " faces, " << data.size(3) << " vertices");

    for (int cell = 0; cell < data.size(0); ++cell) {
        auto cellFaces = data.cellToFace(cell);
        std::string facesStr;
        for (const auto& face : cellFaces) {
            auto tag = data.faceTag(face.index());
            std::string tagStr;
            switch(tag) {
                case 0: tagStr = "I"; break;
                case 1: tagStr = "J"; break;
                case 2: tagStr = "K"; break;
                case 3: tagStr = "NNC"; break;
                default: tagStr = "?"; break;
            }
            facesStr += "[" + std::to_string(face.index()) + " " + tagStr + (face.orientation() ? "+" : "-") + "] ";
        }
        BOOST_TEST_MESSAGE("Cell " << cell << " has " << cellFaces.size() << " faces: " << facesStr);
    }

    // Check that cells 1 and 2 (which will be refined) exist and characterize them
    BOOST_CHECK_EQUAL(data.size(0), 4);
}

BOOST_AUTO_TEST_CASE(faulted_grid_lgr_basic)
{
    // Create the faulted grid with LGR spanning across the fault
    Dune::CpGridLGR grid;
    Opm::createGridAndAddLgrs(grid,
                              faulted_deck_string,
                              {{2, 1, 2}}, // 2x1x2 refinement per parent cell
                              {{1, 0, 0}}, // startIJK (I=2 -> index 1, J=1 -> index 0, K=1 -> index 0)
                              {{3, 1, 1}}, // endIJK (I=3, J=1, K=1)
                              {"LGR1"});

    // Basic checks that the grid was created successfully
    const auto& data = grid.currentData();

    // Level 0 should have 4 cells
    BOOST_CHECK_EQUAL(data[0]->size(0), 4);

    // LGR level (level 1) should have 8 cells (2 parent cells * 2*1*2 children)
    const int lgr1_level = grid.getLgrNameToLevel().at("LGR1");
    BOOST_CHECK_EQUAL(data[lgr1_level]->size(0), 8);

    // Leaf grid: 2 unrefined + 8 refined = 10 cells
    const auto& leafData = grid.currentData().back();
    BOOST_CHECK_EQUAL(leafData->size(0), 10);

    // Run standard grid checks
    Opm::checkGridWithLgrs(grid, {{2, 1, 2}}, {"LGR1"});
}

BOOST_AUTO_TEST_CASE(faulted_grid_lgr_face_geometry)
{
    // Create the faulted grid with LGR spanning across the fault
    Dune::CpGridLGR grid;
    Opm::createGridAndAddLgrs(grid,
                              faulted_deck_string,
                              {{2, 1, 2}}, // 2x1x2 refinement per parent cell
                              {{1, 0, 0}},
                              {{3, 1, 1}},
                              {"LGR1"});

    // Check face properties on the leaf grid
    const auto& leafView = grid.leafGridView();
    for (const auto& element : Dune::elements(leafView)) {
        for (const auto& intersection : Dune::intersections(leafView, element)) {
            // Every face must have a positive area
            const double faceArea = intersection.geometry().volume();
            BOOST_CHECK_GT(faceArea, 0.0);

            // Face normals should be unit normals (length ~1)
            const auto& normal = intersection.centerUnitOuterNormal();
            const double normalLength = normal.two_norm();
            BOOST_CHECK_CLOSE(normalLength, 1.0, 1e-6);

            // If the intersection has a neighbor, verify that inside and outside cells
            // are different
            if (intersection.neighbor()) {
                BOOST_CHECK(intersection.inside().index() != intersection.outside().index());
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(faulted_grid_lgr_cell_volumes)
{
    // Create the faulted grid with LGR spanning across the fault
    Dune::CpGridLGR grid;
    Opm::createGridAndAddLgrs(grid,
                              faulted_deck_string,
                              {{2, 1, 2}}, // 2x1x2 refinement per parent cell
                              {{1, 0, 0}},
                              {{3, 1, 1}},
                              {"LGR1"});

    // Sum of child cell volumes should equal parent cell volume
    // Parent cell 1 (I=2): 100*100*10 = 100000
    // Parent cell 2 (I=3): 100*100*10 = 100000
    double parentVol = 100.0 * 100.0 * 10.0;
    double childVolSum_cell1 = 0.0;
    double childVolSum_cell2 = 0.0;

    const int lgr1_level = grid.getLgrNameToLevel().at("LGR1");
    const auto& levelView = grid.levelGridView(lgr1_level);

    int childCount = 0;
    for (const auto& element : Dune::elements(levelView)) {
        double vol = element.geometry().volume();
        BOOST_CHECK_GT(vol, 0.0);

        // Each parent cell has 2*1*2=4 children
        if (childCount < 4) {
            childVolSum_cell1 += vol;
        } else {
            childVolSum_cell2 += vol;
        }
        ++childCount;
    }

    BOOST_CHECK_EQUAL(childCount, 8);
    BOOST_CHECK_CLOSE(childVolSum_cell1, parentVol, 1e-6);
    BOOST_CHECK_CLOSE(childVolSum_cell2, parentVol, 1e-6);
}

BOOST_AUTO_TEST_CASE(faulted_grid_lgr_intersection_across_fault)
{
    // Create the faulted grid with LGR spanning across the fault
    Dune::CpGridLGR grid;
    Opm::createGridAndAddLgrs(grid,
                              faulted_deck_string,
                              {{2, 1, 2}}, // 2x1x2 refinement per parent cell
                              {{1, 0, 0}},
                              {{3, 1, 1}},
                              {"LGR1"});

    const auto& leafView = grid.leafGridView();

    // Count the number of interior faces (faces with two neighbors) on the leaf grid
    int interiorFaceCount = 0;
    int boundaryFaceCount = 0;
    for (const auto& element : Dune::elements(leafView)) {
        for (const auto& intersection : Dune::intersections(leafView, element)) {
            if (intersection.neighbor()) {
                ++interiorFaceCount;
            } else {
                ++boundaryFaceCount;
            }
        }
    }

    // Each interior face is counted twice (once from each side)
    // Expected interior faces (unique):
    // - 2 unrefined cells: each has connections
    // - Within each refined parent cell: internal faces
    // - Between the two refined parent cells: faces across the fault
    //
    // For 2x1x2 refinement of cell I=2:
    //   I-faces interior: 1*1*2 = 2
    //   K-faces interior: 2*1*1 = 2
    //   Total interior faces in cell I=2: 4
    // Same for cell I=3: 4
    // Between cell I=2 and cell I=3: 1*2 = 2 faces across the fault
    // Between unrefined cell I=1 and LGR (I=2): 1*2 = 2 faces
    // Between LGR (I=3) and unrefined cell I=4: 1*2 = 2 faces
    // Between unrefined cell I=1 and I=2 (not applicable - I=2 is refined)
    // Total unique interior faces:
    //   Within LGR cell I=2: 4
    //   Within LGR cell I=3: 4
    //   Between LGR cells I=2 and I=3 (fault): 2
    //   Between unrefined I=1 and LGR boundary: 2
    //   Between LGR boundary and unrefined I=4: 2
    // Total unique: 14
    // Each counted twice from each side: 28
    // But we also need to account for the fact that each face is seen from each element,
    // so interiorFaceCount = 2 * number_of_unique_interior_faces
    // Boundary faces:
    //   Cell I=1: I_min, J_min, J_max, K_min, K_max = 5 - 1 (shared with LGR) = up to 5 boundary
    //   Cell I=4: I_max, J_min, J_max, K_min, K_max = 5 - 1 (shared with LGR)
    //   LGR cells have boundary faces on J_min, J_max, K_min, K_max boundaries + external I faces

    // The most important check: there should be faces connecting cells across the fault
    // Count connections between cells on different sides of the fault
    // The fault is between the refined cells of parent I=2 and parent I=3
    bool foundCrossFaultConnection = false;
    for (const auto& element : Dune::elements(leafView)) {
        for (const auto& intersection : Dune::intersections(leafView, element)) {
            if (intersection.neighbor()) {
                const auto& insideCenter = element.geometry().center();
                const auto& outsideCenter = intersection.outside().geometry().center();

                // Check if the connection crosses the fault (x around 200)
                // Inside cell is on one side (x < 200) and outside cell is on the other (x > 200)
                // or vice versa
                if ((insideCenter[0] < 200.0 && outsideCenter[0] > 200.0) ||
                    (insideCenter[0] > 200.0 && outsideCenter[0] < 200.0)) {
                    foundCrossFaultConnection = true;
                    // The face connecting cells across the fault should have reasonable geometry
                    const double faceArea = intersection.geometry().volume();
                    BOOST_CHECK_GT(faceArea, 0.0);
                }
            }
        }
    }
    BOOST_CHECK_MESSAGE(foundCrossFaultConnection,
                        "Expected to find face connections across the fault between refined cells");
}

BOOST_AUTO_TEST_CASE(faulted_grid_lgr_face_has_4_vertices)
{
    // Create the faulted grid with LGR spanning across the fault
    Dune::CpGridLGR grid;
    Opm::createGridAndAddLgrs(grid,
                              faulted_deck_string,
                              {{2, 1, 2}}, // 2x1x2 refinement per parent cell
                              {{1, 0, 0}},
                              {{3, 1, 1}},
                              {"LGR1"});

    // Check that all faces have 4 vertices and at most 2 neighboring cells
    Opm::checkFaceHas4VerticesAndMax2NeighboringCells(grid, grid.currentData());
}
