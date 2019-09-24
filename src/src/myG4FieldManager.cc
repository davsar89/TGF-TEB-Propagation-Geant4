//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// -------------------------------------------------------------------

#include "myG4FieldManager.hh"

myG4FieldManager::myG4FieldManager() {

    if (Settings::MAGNETIC_FIELD_MODEL == "IGRF") {
#if defined(__linux__) && defined(__GNUG__)
        mag_field_model_is_igrf = true;
        mag_field_model_is_geolib = false;
        myEarthMagField_IGRF = new EarthMagField_IGRF();
#endif
    } else if (Settings::MAGNETIC_FIELD_MODEL == "WMM") {
        mag_field_model_is_igrf = false;
        mag_field_model_is_geolib = false;
        myEarthMagField_WMM = new EarthMagField_WMM();
    } else if (Settings::MAGNETIC_FIELD_MODEL == "GEOLIB") {
        mag_field_model_is_geolib = true;
        mag_field_model_is_igrf = false;
        myEarthMagField_GeographicLib = new EarthMagField_GeographicLib();
    } else {
        G4cout << "Magnetic model must be IGRF or WMM. Aborting." << G4endl;
        std::abort();
    }

}

myG4FieldManager::~myG4FieldManager() {
    if (Settings::MAGNETIC_FIELD_MODEL == "IGRF") {
#if defined(__linux__) && defined(__GNUG__)
        delete myEarthMagField_IGRF;
#endif
    } else if (Settings::MAGNETIC_FIELD_MODEL == "WMM") {
        delete myEarthMagField_WMM;
    } else if (Settings::MAGNETIC_FIELD_MODEL == "GEOLIB") {
        delete myEarthMagField_GeographicLib;
    }
}

void myG4FieldManager::ConfigureForTrack(const G4Track *aTrack) {
// if a new track is being processes, put the maximum step for charged particles to a small value compared to larmor radius

    const G4int PDG = aTrack->GetParticleDefinition()->GetPDGEncoding();

    //#ifndef NDEBUG // debug mode only
//    if (settings.CHECK_ELECTRON_TRACKING_IN_MAG_FIELD) {
    // for leptons : check step size compared to Larmor radius
    if (PDG == PDG_electron || PDG == PDG_positron) {

        if (aTrack->GetKineticEnergy() >= Settings::MIN_ENERGY_OUTPUT) {

            larmor_radius_m = CalculateLarmorRadius(aTrack->GetStep()->GetPreStepPoint());
            // affects only charged particles
            G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(
                    larmor_radius_m * m / 10.0);
//                static_cast<G4FieldManager *>(this)->SetDeltaOneStep(larmor_radius_m / 10.0);
//                static_cast<G4FieldManager *>(this)->SetDeltaIntersection(larmor_radius_m / 50.0);
//                fDefault_Delta_One_Step_Value(0.01),    // mm
//                fDefault_Delta_Intersection_Val(0.001), // mm

        }
    }

}

double myG4FieldManager::CalculateLarmorRadius(const G4StepPoint *preStepPoint) const {

    const G4ThreeVector &position = preStepPoint->GetPosition();
    const G4ThreeVector &momentum_dir = preStepPoint->GetMomentumDirection();

    G4ThreeVector local_mag_field{0, 0, 0};

    if (mag_field_model_is_igrf) {
#if defined(__linux__) && defined(__GNUG__)
        local_mag_field = myEarthMagField_IGRF->GetFieldComponents(position / meter);
#endif
    } else if (mag_field_model_is_geolib) {
        local_mag_field = myEarthMagField_GeographicLib->GetFieldComponents(position / meter);
    } else {
        local_mag_field = myEarthMagField_WMM->GetFieldComponents(position / meter);
    }

    const double Kenergy = preStepPoint->GetKineticEnergy() / keV;
    const double gamma = Kenergy / ere + 1.0;
    const double velocity_mag = sol * sqrt(1.0 - powf(gamma, -2.0));
    const G4ThreeVector velocity_vector = momentum_dir * velocity_mag;
    const G4ThreeVector mag_field_dir = local_mag_field / local_mag_field.mag();
    const G4ThreeVector v_par = mag_field_dir * mag_field_dir.dot(velocity_vector);
    const G4ThreeVector v_perp = velocity_vector - v_par;

    return factor * gamma * v_perp.mag() / local_mag_field.mag();
}
