////////////////////////////////////////////////////////////////////////////////

// /* GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's atmosphere */
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

#include <Settings.hh>
#include <vector>

#include <Analysis.hh>
#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4VSensitiveDetector.hh"
#include <fstream>
#include "Settings.hh"
#include "Analysis.hh" // singleton
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

class G4Step;

struct record_coords {
    G4ThreeVector position;
    G4double time;
};

class SteppingAction : public G4UserSteppingAction {
public:

    SteppingAction();


    ~SteppingAction() override;

    void UserSteppingAction(const G4Step *aStep) override;

private:

    Settings *settings = Settings::getInstance();
    Analysis *analysis = Analysis::getInstance();

    G4StepPoint *thePrePoint = nullptr;
    G4StepPoint *thePostPoint = nullptr;

    const G4int PDG_positron = -11;
    const G4int PDG_electron = 11;
    const G4int PDG_photon = 22;

    const G4double earth_radius = 6371.137 * km;
    const G4double photon_max_altitude = 800.0 * km;

    G4double Get_dist_rad(const G4double &lat, const G4double &lon,
                          const G4double &alt_record);
};
