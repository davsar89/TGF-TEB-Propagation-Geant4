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

#include "G4FieldManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"

#if defined(__linux__) && defined(__GNUG__)

#include "EarthMagField_IGRF.hh"

#endif

#include "EarthMagField_WMM.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManagerStore.hh"
#include "EarthMagField_GeographicLib.hh"

class myG4FieldManager : public G4FieldManager {
public:

    myG4FieldManager();

    ~myG4FieldManager() override;

    void ConfigureForTrack(const G4Track *) override;

    inline double Get_Larmor_Radius() { return larmor_radius_m; };

    inline void Set_Larmor_Radius(const double &lr) { larmor_radius_m = lr; };

    double CalculateLarmorRadius(const G4StepPoint *preStepPoint) const;

private:

    bool mag_field_model_is_igrf = false;
    bool mag_field_model_is_geolib = false;

    const G4int PDG_positron = -11;
    const G4int PDG_electron = 11;

    const double e_c = 1.60217662e-19;
    const double sol = 299792458.0;
    const double m_e = 9.10938356e-31;
    const double ere = 510.9989500; // in keV
    const double factor = 5.685630e-12; // m_e/e_c

    double mag_tmp;

    int previous_ID = -87;

    G4double larmor_radius_m = 1.0e10;

#if defined(__linux__) && defined(__GNUG__)
    EarthMagField_IGRF *myEarthMagField_IGRF = nullptr;
#endif

    EarthMagField_WMM *myEarthMagField_WMM = nullptr;
    EarthMagField_GeographicLib *myEarthMagField_GeographicLib = nullptr;

};