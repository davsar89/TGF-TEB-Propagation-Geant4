////////////////////////////////////////////////////////////////////////////////

// /* GEANT4 code for propagation of gamma-rays, electron and positrons in
// Earth's atmosphere */
//
// //
// // ********************************************************************
// // * License and Disclaimer                                           *
// // *                                                                  *
// // * The  Geant4 software  is  copyright of the Copyright Holders  of *
// // * the Geant4 Collaboration.  It is provided  under  the terms  and *
// // * conditions of the Geant4 Software License,  included in the file *
// // * LICENSE and available at  http://cern.ch/geant4/license .  These *
// // * include a list of copyright holders.                             *
// // *                                                                  *
// // * Neither the authors of this software system, nor their employing *
// // * institutes,nor the agencies providing financial support for this *
// // * work  make  any representation or  warranty, express or implied, *
// // * regarding  this  software system or assume any liability for its *
// // * use.  Please see the license in the file  LICENSE  and URL above *
// // * for the full disclaimer and the limitation of liability.         *
// // *                                                                  *
// // * This  code  implementation is the result of  the  scientific and *
// // * technical work of the GEANT4 collaboration.                      *
// // * By using,  copying,  modifying or  distributing the software (or *
// // * any work based  on the software)  you  agree  to acknowledge its *
// // * use  in  resulting  scientific  publications,  and indicate your *
// // * acceptance of all terms of the Geant4 Software license.          *
// // ********************************************************************
////////////////////////////////////////////////////////////////////////////////

#include <Analysis.hh>
#include <Settings.hh>
#include <fstream>
#include <iomanip>

#include "G4SteppingManager.hh"
#include "G4Track.hh"

// singleton pattern initialization
Analysis *Analysis::instance = nullptr;

// constructor
Analysis::Analysis() {
    const G4double ALT_MAX_RECORDED = *std::max_element(
            settings->record_altitudes.begin(), settings->record_altitudes.end());
    //
    const G4String output_filename_second_part =
            std::to_string(settings->RANDOM_SEED) + "_" +
            std::to_string(int(ALT_MAX_RECORDED)) + "_" +
            std::to_string(int(settings->SOURCE_ALT)) + "_" +
            std::to_string(int(settings->OPENING_ANGLE)) + "_" +
            settings->BEAMING_TYPE + "_" +
            std::to_string(int(settings->SOURCE_SIGMA_TIME)) + ".out";
    //
    asciiFileName2 = "./output_ascii/detParticles_" + output_filename_second_part;
    std::ofstream asciiFile00(asciiFileName2,
                              std::ios::trunc); // to clean the output file
    asciiFile00.close();
    //
    output_lines.clear();
    //

    if ((settings->BEAMING_TYPE == "Uniform") ||
        (settings->BEAMING_TYPE == "uniform")) {
        number_beaming = 0;
    } else if ((settings->BEAMING_TYPE == "Gaussian") ||
               (settings->BEAMING_TYPE == "gaussian") ||
               (settings->BEAMING_TYPE == "normal") ||
               (settings->BEAMING_TYPE == "Normal")) {
        number_beaming = 1;
    }

    //
    if (settings->RECORD_ELEC_POSI_ONLY) { // avoid having a too big buffer size
        // if record is only for leptons
        if (output_buffer_size > 100) {
            output_buffer_size = 100;
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Analysis::get_NB_RECORDED() const { return NB_RECORDED; }

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Analysis::~Analysis() = default;

Analysis *Analysis::getInstance() {
    if (instance == nullptr) {
        instance = new Analysis;
    }

    return instance;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Analysis::save_in_output_buffer(
        const G4int PDG_NB, const G4double &time, const G4double &energy,
        const G4double &dist_rad, const G4int ID, const G4double &ecef_x,
        const G4double &ecef_y, const G4double &ecef_z, const G4double &mom_x,
        const G4double &mom_y, const G4double &mom_z, const G4double& lat, const G4double& lon, const G4double& alt) {

    //
    double alt2 = alt / 1000.0; // m to km

    bool record_or_not = false;

    if (settings->RECORD_ONLY_IN_WINDOW) {
        const bool is_inside_record_window = (lat > settings->RECORD_WIN.min_lat && lat < settings->RECORD_WIN.max_lat)
                                             && (lon > settings->RECORD_WIN.min_lon && lon < settings->RECORD_WIN.max_lon);

        if (is_inside_record_window) record_or_not = true;
    } else {
        record_or_not = true;
    }

    if (record_or_not) {

        // ASCII OUTPUT
        if (settings->OUTPUT_TO_ASCII_FILE) {
            std::stringstream buffer;
            buffer << std::scientific
                   << std::setprecision(7); // scientific notation with
            // 5 significant digits
            //   asciiFile << name;
            //   asciiFile << ' ';
            buffer << settings->RANDOM_SEED;
            buffer << ' ';
            buffer << settings->SOURCE_ALT;
            buffer << ' ';
            buffer << settings->OPENING_ANGLE;
            buffer << ' ';
            buffer << settings->TILT_ANGLE;
            buffer << ' ';
            buffer << settings->NB_EVENT; // 5
            buffer << ' ';
            buffer << ID;
            buffer << ' ';
            buffer << PDG_NB;
            buffer << ' ';
            buffer << time;
            buffer << ' ';
            buffer << energy;
            buffer << ' ';
            buffer << alt2; // 10
            buffer << ' ';
            buffer << lat;
            buffer << ' ';
            buffer << lon;
            buffer << ' ';
            buffer << dist_rad;
            buffer << ' ';
            buffer << ecef_x;
            buffer << ' ';
            buffer << ecef_y; // 15
            buffer << ' ';
            buffer << ecef_z;
            buffer << ' ';
            buffer << mom_x;
            buffer << ' ';
            buffer << mom_y;
            buffer << ' ';
            buffer << mom_z;
            buffer << ' ';
            buffer << number_beaming; // 20 // number_beaming == 0 for uniform and 1 for // gaussian
            buffer << ' ';
            buffer << settings->SOURCE_LAT;
            buffer << ' ';
            buffer << settings->SOURCE_LONG;
            buffer << ' ';
            buffer << G4endl;
            //
            NB_RECORDED++;
            //

            output_lines.push_back(buffer.str());
            //
            write_in_output_file();
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Analysis::write_in_output_file() {

    if (output_lines.size() <= output_buffer_size) {
        return;
    }

    std::ofstream asciiFile2;
    asciiFile2.open(asciiFileName2, std::ios::app);

    if (asciiFile2.is_open()) {
        for (G4String &line : output_lines) {
            asciiFile2 << line;
        }

        asciiFile2.close();
        output_lines.clear();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Analysis::write_in_output_file_endOfRun() {
    if (output_lines.empty()) {
        return;
    }

    std::ofstream asciiFile1;
    asciiFile1.open(asciiFileName2, std::ios::app);

    if (asciiFile1.is_open()) {
        for (G4String &line : output_lines) {
            asciiFile1 << line;
        }

        asciiFile1.close();
        output_lines.clear();
    }
}
