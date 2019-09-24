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
#pragma once

#include "G4ios.hh"
#include "globals.hh"

#include "EarthMagField_WMM.hh"

#include "Settings.hh"
#include "geodetic_converter.hh"

// following a singleton pattern

class G4Track;

class Analysis {
public:
    ~Analysis();

    void save_in_output_buffer(const G4int PDG_NB, const G4double &time,
                               const G4double &energy, const G4double &dist_rad,
                               const G4int ID, const G4double &ecef_x,
                               const G4double &ecef_y, const G4double &ecef_z,
                               const G4double &mom_x, const G4double &mom_y,
                               const G4double &mom_z, const G4double& lat, const G4double& lon, const G4double& alt);

    static Analysis *getInstance();

    G4int get_NB_RECORDED() const;

    Analysis(Analysis const &) = delete; // Don't Implement copy constructor
    void
    operator=(Analysis const &) = delete; // don't implement equality operator

    void write_in_output_file();

    void write_in_output_file_endOfRun();

private:
    Analysis();

    static Analysis *instance;

    Settings *settings = Settings::getInstance();

    std::vector<G4String> output_lines;

    G4String asciiFileName2;

    uint output_buffer_size = 3;
    //

    G4int NB_RECORDED = 0;

    G4int number_beaming = 0;

};
