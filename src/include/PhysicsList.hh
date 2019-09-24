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

#include "G4SystemOfUnits.hh"
#include "G4VUserPhysicsList.hh"
#include "Settings.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option4_dr.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4RadioactiveDecay.hh"

class G4VPhysicsConstructor;

class TGF_PhysicsList : public G4VUserPhysicsList {
public:

    TGF_PhysicsList();

    ~TGF_PhysicsList() override;

    void SetCuts() override;

    void ConstructParticle() override;

    void ConstructProcess() override;

private:

    Settings *settings = Settings::getInstance();

    G4double cutForGamma = 0;
    G4double cutForElectron = 0;
    G4double cutForPositron = 0;

    G4VPhysicsConstructor *emPhysicsList = nullptr;
    G4VPhysicsConstructor *radDecPhysicsList = nullptr;

    void Add_StepMax_for_record_regions();

    void AddStepMax_GLOBAL(G4double stepMaxVal_elec, G4double stepMaxVal_gamma);
};
