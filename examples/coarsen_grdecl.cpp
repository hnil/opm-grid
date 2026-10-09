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

/*
  Coarsen a corner-point grid description (SPECGRID/COORD/ZCORN/ACTNUM) with
  COARSEN-style records and write the result as a grdecl include file.

    coarsen_grdecl in.grdecl out.grdecl \
        --coarsen "1 4 1 4 1 6 2 2 3" [--coarsen ...] \
        [--fill-holes] [--allow-vertical-gaps] [--max-sub-faces N] [--lgr]

  Record items are 1-based and inclusive, as in the deck: I1 I2 J1 J2 K1 K2
  NX NY NZ.
*/

#include "config.h"

#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

/// Read one keyword's numbers, skipping comments; stops at the terminating '/'.
template <class T>
std::vector<T> readKeyword(std::istream& is)
{
    std::vector<T> values;
    std::string token;
    while (is >> token) {
        if (token.rfind("--", 0) == 0) {
            std::getline(is, token);
            continue;
        }
        if (token == "/") {
            break;
        }
        const auto star = token.find('*');
        if (star != std::string::npos) {          // ECL repeat count, n*value
            const int n = std::stoi(token.substr(0, star));
            const T v = static_cast<T>(std::stod(token.substr(star + 1)));
            values.insert(values.end(), n, v);
            continue;
        }
        if (token == "F" || token == "T") {       // SPECGRID's trailing flag
            continue;
        }
        values.push_back(static_cast<T>(std::stod(token)));
    }
    return values;
}

Opm::Coarsening::Grdecl readGrdecl(const std::string& path)
{
    std::ifstream is(path);
    if (!is) {
        throw std::runtime_error("cannot open " + path);
    }

    Opm::Coarsening::Grdecl g;
    std::string token;
    while (is >> token) {
        if (token.rfind("--", 0) == 0) {
            std::getline(is, token);
        } else if (token == "SPECGRID" || token == "DIMENS") {
            const auto v = readKeyword<int>(is);
            g.dims = {v.at(0), v.at(1), v.at(2)};
        } else if (token == "COORD") {
            g.coord = readKeyword<double>(is);
        } else if (token == "ZCORN") {
            g.zcorn = readKeyword<double>(is);
        } else if (token == "ACTNUM") {
            g.actnum = readKeyword<int>(is);
        }
    }
    if (g.coord.empty() || g.zcorn.empty()) {
        throw std::runtime_error(path + ": no COORD/ZCORN found");
    }
    return g;
}

void writeVector(std::ostream& os, const char* keyword, const std::vector<double>& v,
                 int perLine)
{
    os << keyword << '\n' << std::setprecision(10);
    for (std::size_t i = 0; i < v.size(); ++i) {
        os << ' ' << v[i];
        if ((i + 1) % perLine == 0) {
            os << '\n';
        }
    }
    os << "\n/\n\n";
}

Opm::Coarsening::CoarsenRequest parseRecord(const std::string& text)
{
    std::istringstream is(text);
    std::vector<int> v;
    for (int x = 0; is >> x;) {
        v.push_back(x);
    }
    if (v.size() != 9) {
        throw std::runtime_error("a COARSEN record needs 9 numbers: " + text);
    }
    Opm::Coarsening::CoarsenRequest r;
    for (int d = 0; d < 3; ++d) {
        r.startIJK[d] = v[2*d] - 1;
        r.endIJK[d] = v[2*d + 1];
        r.cellsPerDim[d] = v[6 + d];
    }
    return r;
}

} // anonymous namespace

int main(int argc, char** argv)
try {
    std::vector<std::string> args(argv + 1, argv + argc);
    std::vector<Opm::Coarsening::CoarsenRequest> requests;
    Opm::Coarsening::Options options;
    std::string in, out;
    bool lgr = false;

    for (std::size_t a = 0; a < args.size(); ++a) {
        if (args[a] == "--coarsen") {
            requests.push_back(parseRecord(args.at(++a)));
        } else if (args[a] == "--fill-holes") {
            options.activity = Opm::Coarsening::Activity::FillHoles;
        } else if (args[a] == "--allow-vertical-gaps") {
            options.allowVerticalGaps = true;
        } else if (args[a] == "--max-sub-faces") {
            options.maxSubFaces = std::stoi(args.at(++a));
        } else if (args[a] == "--lgr") {
            lgr = true;
        } else if (in.empty()) {
            in = args[a];
        } else {
            out = args[a];
        }
    }
    if (in.empty() || out.empty()) {
        std::cerr << "usage: coarsen_grdecl in.grdecl out.grdecl --coarsen \"I1 I2 J1 J2 K1 K2 "
                     "NX NY NZ\" [...]\n";
        return EXIT_FAILURE;
    }

    const auto fine = readGrdecl(in);
    const auto result = coarsenCornerPoint(fine, requests, options);

    std::ofstream os(out);
    os << "-- coarsened from " << in << " by coarsen_grdecl\n";
    for (const auto& note : result.report.notes) {
        os << "-- " << note << '\n';
    }
    os << "SPECGRID\n " << result.grid.dims[0] << ' ' << result.grid.dims[1] << ' '
       << result.grid.dims[2] << " 1 F /\n\n";
    writeVector(os, "COORD", result.grid.coord, 6);
    writeVector(os, "ZCORN", result.grid.zcorn, 8);
    os << "ACTNUM\n";
    for (std::size_t i = 0; i < result.grid.actnum.size(); ++i) {
        os << ' ' << result.grid.actnum[i];
        if ((i + 1) % 20 == 0) {
            os << '\n';
        }
    }
    os << "\n/\n";

    for (const auto& note : result.report.notes) {
        std::cout << note << '\n';
    }
    std::cout << "max sub-faces per face: " << result.report.maxSubFacesSeen
              << ", max face non-planarity: " << result.report.maxFaceNonPlanarity << '\n';

    if (lgr) {
        std::cout << "refinement back to the fine grid:\n";
        for (const auto& r : refinementBackToFine(result)) {
            std::cout << "  " << r.name << " CARFIN " << r.startIJK[0] + 1 << '-' << r.endIJK[0]
                      << ' ' << r.startIJK[1] + 1 << '-' << r.endIJK[1]
                      << ' ' << r.startIJK[2] + 1 << '-' << r.endIJK[2]
                      << "  cells/dim " << r.cellsPerDim[0] << ' ' << r.cellsPerDim[1] << ' '
                      << r.cellsPerDim[2] << '\n';
        }
    }
    return EXIT_SUCCESS;
} catch (const std::exception& e) {
    std::cerr << "coarsen_grdecl: " << e.what() << '\n';
    return EXIT_FAILURE;
}
