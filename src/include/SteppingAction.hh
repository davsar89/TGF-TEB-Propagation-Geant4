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

#include "Analysis.hh"
#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4VSensitiveDetector.hh"
#include <fstream>
#include "Settings.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "EarthMagField_WMM.hh"

#if defined(__linux__) && defined(__GNUG__)

#include "EarthMagField_IGRF.hh"

#endif

#include <RunAction.hh>
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "myG4FieldManager.hh"
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Constants.hpp>

class G4Step;

struct record_coords {
    G4ThreeVector position;
    G4double time;
};

class SteppingAction : public G4UserSteppingAction {
public:

    explicit SteppingAction();

    ~SteppingAction() override;

    void UserSteppingAction(const G4Step *aStep) override;

private:

    GeographicLib::Geocentric *earth = nullptr;

    myG4FieldManager *myG4FieldMan = nullptr;

    const double e_c = 1.60217662e-19;
    const double sol = 299792458.0;
    const double m_e = 9.10938356e-31;
    const double ere = 510.9989500; // in keV
    const double factor = 5.685630e-12; // m_e/e_c

    double mag_tmp;
    double velocity_mag;
    double Kenergy;
    double gamma;

    G4ThreeVector momentum_dir, local_mag_field, velocity_vector, mag_field_dir;
    G4ThreeVector v_par, v_perp, field_components;

    int previous_ID = -87;

    G4ThreeVector previous_pos{0, 0, 0};

    Analysis& analysis = Analysis::getInstance();

    double sol_SI = 299792458.0;

    struct min_max_altitudes {
        double min_val;
        double max_val;
    };

    min_max_altitudes min_max_alt{0.0, 0.0};

    void record_particles(const G4Step *aStep, const int PDG);

    double project_to_record_alt(
            const G4int type_number, const double &kinE,
            double &ecef_x, double &ecef_y, double &ecef_z,
            const double &lat_in, const double &lon_in, const double &alt_in);

    double ALT_RECORD_in_meter = Settings::record_altitude * 1000.0;

    const G4int PDG_positron = -11;
    const G4int PDG_electron = 11;
    const G4int PDG_photon = 22;

    const G4double earth_radius = 6371.137 * km;
    const G4double photon_max_altitude = 800.0 * km;

    G4double Get_dist_rad(const G4double &lat, const G4double &lon,
                          const G4double &alt_record);

    const static int nb = 512;
    std::vector<double> grid_latitude;
    std::vector<double> grid_longitude;

    min_max_altitudes find_min_max_rough_altitude();

    double find_rough_altitude_for_real_altitude(const double &lat, const double &lon, const double &target_alt);

    /////////////////////

    template<typename T>
    std::vector<double> linspace(T start_in, T end_in, int num_in) {

        std::vector<double> linspaced;

        double start = static_cast<double>(start_in);
        double end = static_cast<double>(end_in);
        double num = static_cast<double>(num_in);

        if (num == 0) { return linspaced; }
        if (num == 1) {
            linspaced.push_back(start);
            return linspaced;
        }

        double delta = (end - start) / (num - 1);

        for (int i = 0; i < num - 1; ++i) {
            linspaced.push_back(start + delta * i);
        }
        linspaced.push_back(end); // I want to ensure that start and end
        // are exactly the same as the input
        return linspaced;
    }
};
