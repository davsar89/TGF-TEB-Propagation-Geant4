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

#include <RunAction.hh>
#include <SteppingAction.hh>
#include "G4Track.hh"
#include "G4RunManager.hh"

SteppingAction::SteppingAction() {
    //    asciiFileName = "./extra_outputs/created_electrons";
    //    std::ofstream asciiFile00(asciiFileName, std::ios::trunc); // to clean the output file
    //    asciiFile00.close();
    //    asciiFileName_phot = "./extra_outputs/created_photons";
    //    std::ofstream asciiFile11(asciiFileName_phot, std::ios::trunc); // to clean the output file
    //    asciiFile11.close();
}

SteppingAction::~SteppingAction() = default;

void SteppingAction::UserSteppingAction(const G4Step *aStep) {

    const G4int PDG = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

#ifndef NDEBUG // debug mode only
    if (settings->USE_STEP_MAX_GLOBAL) {
        if ((PDG == PDG_electron) || (PDG == PDG_positron)) {
            if ((aStep->GetStepLength() > settings->GLOBAL_STEP_MAX_VAL_LEPTONS)) {
                G4cout << "Current step length : " << aStep->GetStepLength() / meter << " meter" << G4endl;
                G4cout << "Error in SteppingAction : step is lager than the maximum allowed." << G4endl;
                std::abort();
            }
        }
    }
#endif // end debug mode only
    //
    G4Track *track = aStep->GetTrack();
    G4double global_time = track->GetGlobalTime();

    if (global_time > settings->TIME_LIMIT) {
        track->SetTrackStatus(fStopAndKill);
        return;
    }

    if (track->GetKineticEnergy() < settings->MIN_ENERGY_OUTPUT) { // killing particles with energy lower than threshold
        if (PDG != PDG_positron) { // don't kill positrons to make sure they do annihilation before disappearing
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    const G4double rough_altitude = aStep->GetPreStepPoint()->GetPosition().mag() - earth_radius;

    if (PDG == PDG_photon) {

        if (rough_altitude > photon_max_altitude) {
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    if (rough_altitude > 300.0 * km) {

        const G4int ID = track->GetTrackID();

        thePrePoint = aStep->GetPreStepPoint();
        double ecef_x = thePrePoint->GetPosition().x() / m;
        double ecef_y = thePrePoint->GetPosition().y() / m;
        double ecef_z = thePrePoint->GetPosition().z() / m;
        double pre_lat, pre_lon, pre_alt;
        geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x, ecef_y, ecef_z, pre_lat, pre_lon, pre_alt);

        thePostPoint = aStep->GetPostStepPoint();
        ecef_x = thePostPoint->GetPosition().x() / m;
        ecef_y = thePostPoint->GetPosition().y() / m;
        ecef_z = thePostPoint->GetPosition().z() / m;
        double post_lat, post_lon, post_alt;
        geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x, ecef_y, ecef_z, post_lat, post_lon, post_alt);

        if ((pre_alt/1000.0 <= settings->record_altitudes[0] && post_alt/1000.0 > settings->record_altitudes[0]) ||
            (pre_alt/1000.0 > settings->record_altitudes[0] && post_alt/1000.0 <= settings->record_altitudes[0])) {

            double avg_alt = (pre_alt + post_alt) / 2.0;
            double avg_lon = (pre_lon + post_lon) / 2.0;
            double avg_lat = (pre_lat + post_lat) / 2.0;

            double dist_rad = Get_dist_rad(avg_lat, avg_lon,
                                           avg_alt * meter);

            const double energy = thePrePoint->GetKineticEnergy() / keV;
            const double time1 =
                    (thePrePoint->GetGlobalTime()) / microsecond;
            const double time2 =
                    (thePostPoint->GetGlobalTime()) / microsecond;
            const double time = (time1 + time2) / 2.0;

            const double mom_x = thePrePoint->GetMomentumDirection().getX();
            const double mom_y = thePrePoint->GetMomentumDirection().getY();
            const double mom_z = thePrePoint->GetMomentumDirection().getZ();

            analysis->save_in_output_buffer(
                    PDG, time, energy, dist_rad / km, ID, ecef_x / 1000.0,
                    ecef_y / 1000.0, ecef_z / 1000.0, // from m to km
                    mom_x, mom_y, mom_z, avg_lat,avg_lon,avg_alt);

        }
    }
}

///////////////////////////////////////////////////////////////////////////

G4double SteppingAction::Get_dist_rad(const G4double &lat, const G4double &lon,
                                      const G4double &alt_record) {
    // getting the radial distance (along curve parrallel to Earth surface)
    // Equirectangular approximation -> leads to "effective speed" of the output
    // data that can be slightly greater than the speed of light
    //    G4double RR = (settings->EarthRadius() + alt_record) ;
    //    G4double delta_lambda = (settings->SourceLong() - lon) * degree;
    //    G4double sum_phi = (lat + settings->SourceLat()) * degree;
    //
    //    G4double xx = delta_lambda * cos((sum_phi) / 2.);
    //    G4double yyy = (settings->SourceLat() - lat) * degree;
    //    G4double dist_rad = sqrt(pow(xx, 2) + pow(yyy, 2)) * RR;
    // haversine formula (better)
    const G4double RR = (settings->earthRadius + alt_record);
    const G4double phi1 = (lat) * degree;
    const G4double phi2 = (settings->SOURCE_LAT) * degree;
    const G4double delta_phi = (phi2 - phi1);
    const G4double sin_delta_phi_over_2 = sin(delta_phi / 2.);
    const G4double delta_lambda = (settings->SOURCE_LONG - lon) * degree;
    const G4double sin_delta_lambda_over_2 = sin(delta_lambda / 2.);
    const G4double aa =
            sin_delta_phi_over_2 * sin_delta_phi_over_2 +
            cos(phi1) * cos(phi2) * sin_delta_lambda_over_2 * sin_delta_lambda_over_2;
    const G4double cc = 2.0 * atan2(sqrt(aa), sqrt(1.0 - aa));
    // G4double cc = 2. * std::asin(std::min(1., sqrt(aa)));
    const G4double dist_rad = RR * cc;
    return dist_rad;
}