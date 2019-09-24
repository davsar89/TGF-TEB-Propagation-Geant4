//////////////////////////////////////////////////////////////////////////////////

//// /* GEANT4 code for propagation of gamma-rays, electron and positrons in
/// Earth's atmosphere */
////
//// //
//// // ********************************************************************
//// // * License and Disclaimer                                           *
//// // *                                                                  *
//// // * The  Geant4 software  is  copyright of the Copyright Holders  of *
//// // * the Geant4 Collaboration.  It is provided  under  the terms  and *
//// // * conditions of the Geant4 Software License,  included in the file *
//// // * LICENSE and available at  http://cern.ch/geant4/license .  These *
//// // * include a list of copyright holders.                             *
//// // *                                                                  *
//// // * Neither the authors of this software system, nor their employing *
//// // * institutes,nor the agencies providing financial support for this *
//// // * work  make  any representation or  warranty, express or implied, *
//// // * regarding  this  software system or assume any liability for its *
//// // * use.  Please see the license in the file  LICENSE  and URL above *
//// // * for the full disclaimer and the limitation of liability.         *
//// // *                                                                  *
//// // * This  code  implementation is the result of  the  scientific and *
//// // * technical work of the GEANT4 collaboration.                      *
//// // * By using,  copying,  modifying or  distributing the software (or *
//// // * any work based  on the software)  you  agree  to acknowledge its *
//// // * use  in  resulting  scientific  publications,  and indicate your *
//// // * acceptance of all terms of the Geant4 Software license.          *
//// // ********************************************************************
//////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <vector>

// singleton pattern with parameters set as public to avoid the overhead of
// using getters and setters

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class Settings
{
private:
    Settings()
    {
        record_altitudes.clear();
    } // Private so that it can not be called

    // copy constructor is private
    // assignment operator is private
    static Settings *instance;

public:
    Settings(Settings const &) = delete; // Don't Implement copy constructor
    void
    operator=(Settings const &) = delete; // don't implement equality operator

    static Settings *getInstance();

public:
    //////////////// Parameters are listed below ////////////////
    // date
    const int year = 2018;
    const double day = 24;
    const double month = 3;

    // to record particles only inside a specific latitude and longitude window
    const bool RECORD_ONLY_IN_WINDOW = true;
    struct record_window
    {
        double min_lat; // degrees
        double max_lat;
        double min_lon;
        double max_lon;
    };

    const record_window RECORD_WIN{-3.5, 3.5, 52.5, 58.5}; // degrees

    // Earth radius
    const double earthRadius = 6378.137 * km;

    const G4bool MAG_FIELD_ON = true;

    const G4bool OUTPUT_ALT_LAYERS_TO_FILE =
        false; // output list of altitude and densities of layer to file (for
    // debug)

    // parameters : initialization values, can be modified in the main code
    G4long RANDOM_SEED = 12345; // dummy value that will be replaced
    G4int NB_EVENT = 0;

    const double CACHED_LENGTH =
        10.0 * meter; // for magnetic field solver, in meters

    // Source parameters, geodetic coordinates ( = geographic = GPS)

    // ASIM TEB 180916
    // -6.998426914215088, 55.92176818847656 -> magnetic field line tracing
    // -7.0492, 55.9129 −> VAISALA
    const double SOURCE_LAT = -7.0492;  // degree
    const double SOURCE_LONG = 55.9129; // degree

    double SOURCE_ALT = 15.;    // km
    double OPENING_ANGLE = 30.; // degree
    G4String BEAMING_TYPE = "Gaussian";
    double TILT_ANGLE = 0.0;

    const double TIME_LIMIT = 2.0 * second;
    const double MIN_ENERGY_OUTPUT = 10.0 * keV;
    G4double SOURCE_SIGMA_TIME = 0.; // microsecond
    std::vector<double>
        record_altitudes; // ! : geodetic altitudes (remark: when building the
    // geometry, geocentric altitudes are used)

    // force max step only for layers where particles are recorded
    const G4bool USE_STEP_MAX_for_record = false;
    const G4bool USE_STEP_MAX_GLOBAL = true; // force max step everywhere
    const double GLOBAL_STEP_MAX_VAL_LEPTONS = 200.0 * meter;
    const double GLOBAL_STEP_MAX_VAL_PHOTONS = 400.0 * meter;
    const double STEP_MAX_RECORD_AREA = 15.0 * meter;

    const G4bool RECORD_ELEC_POSI_ONLY = false; // record only electron and positrons
    const G4bool RECORD_PHOT_ONLY = false;      // record only photons

    const G4bool OUTPUT_TO_ASCII_FILE = true;

    const G4String MAGNETIC_FIELD_MODEL = "WMM"; // "WMM" or "IGRF" ; "IGRF" only possible on LINUX. Windows only has WWM

    const double dr_over_R = 0.001; // stepping parameter for ionization, default is 0.2, that may be too high
};

// #pragma once

// #include "G4RunManager.hh"
// #include "G4UnitsTable.hh"
// #include "G4SystemOfUnits.hh"
// #include "G4PhysicalConstants.hh"

// #include <vector>

////
///....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// #include <string>

// class Settings
// {
//    private:

//        Settings() {} // Private so that it can not be called

//        Settings(Settings const &) {}

//        // copy constructor is private
//        // assignment operator is private
//        static Settings *instance;

//        // Earth radius
//        const G4double earthRadius = 6378.137 * km;

//        // parameters : initialization values, can be modified in the main
//        code G4int Rand_seed = 78;

//        G4int NB_EVENT = 0;

//        G4String CACHED_LENGTH = "10";   // for magnetic field solver

//        // Source parameters, geodetic coordinates ( = geographic = GPS)
//        G4double SOURCE_LAT  = -13.;     // degree
//        G4double SOURCE_LONG = 32.;      // degree
//        G4double SOURCE_ALT  = 15;       // km
//        G4double ALT_MAX_RECORDED     = 500.;     // initial value

//        G4double OPENING_ANGLE = 45.;    // degree
//        G4String BEAMING_TYPE  = "Uniform";

//        G4double SOURCE_SIGMA_TIME = 0.; // microsecond

//        // output altitudes
//        G4int nb_altitude_record = 0;
//        std::vector < G4double > record_altitudes; // ! : geodetic altitudes
//        (remark: when building the geometry, geocentric altitudes are used)

//        void set_AltMax_recorded();

//        G4bool MAG_FIELD_ON_ = false;

//        G4bool USE_STEP_MAX_ = false;

//    public:

//        static Settings *get_Instance()
//        {
//            if (instance == 0)   // Only allow one instance of class to be
//            generated (lazy initialization)
//                {
//                    instance = new
//                    Settings();
//                }

//            return instance;
//        }

//        G4double AltMax_recorded() const
//        {
//            return ALT_MAX_RECORDED;
//        }

//        G4double AltMax_recorded_times_km() const
//        {
//            return ALT_MAX_RECORDED * km;
//        }

//        const G4String &CachedLength() const
//        {
//            return CACHED_LENGTH;
//        }

//        void set_CachedLength(const G4String &cachedLength)
//        {
//            CACHED_LENGTH = cachedLength;
//        }

//        G4int NbAltitudeRecord() const
//        {
//            return nb_altitude_record;
//        }

//        G4int NbEvent() const
//        {
//            return NB_EVENT;
//        }

//        void set_NbEvent(G4int nbEvent)
//        {
//            NB_EVENT = nbEvent;
//        }

//        G4double OpeningAngle() const
//        {
//            return OPENING_ANGLE;
//        }

//        void set_OpeningAngle(G4double openingAngle)
//        {
//            OPENING_ANGLE = openingAngle;
//        }

//        const std::vector < G4double > &Output_Altitudes() const
//        {
//            return record_altitudes;
//        }

//        void set_OutputAltitudes(const std::vector < G4double >
//        &outputAltitudes)
//        {
//            record_altitudes = outputAltitudes;
//        }

//        G4int RandSeed() const
//        {
//            return Rand_seed;
//        }

//        void set_RandSeed(G4int randSeed)
//        {
//            Rand_seed = randSeed;
//        }

//        G4double SourceAlt() const
//        {
//            return SOURCE_ALT;
//        }

//        void set_SourceAlt(G4double sourceAlt)
//        {
//            SOURCE_ALT = sourceAlt;
//        }

//        G4double SourceLat() const
//        {
//            return SOURCE_LAT;
//        }

//        void set_SourceLat(G4double sourceLat)
//        {
//            SOURCE_LAT = sourceLat;
//        }

//        G4double SourceLong() const
//        {
//            return SOURCE_LONG;
//        }

//        void set_SourceLong(G4double sourceLong)
//        {
//            SOURCE_LONG = sourceLong;
//        }

//        G4double SourceSigmaTime() const
//        {
//            return SOURCE_SIGMA_TIME;
//        }

//        void set_SourceSigmaTime(G4double sourceSigmaTime)
//        {
//            SOURCE_SIGMA_TIME = sourceSigmaTime;
//        }

//        void add_OutputAltitude(G4double altitude)
//        {
//            record_altitudes.push_back(altitude);
//            nb_altitude_record++;
//            set_AltMax_recorded();
//        }

//        std::vector < G4double > Output_Altitudes()
//        {
//            return record_altitudes;
//        }

//        G4double EarthRadius() const
//        {
//            return earthRadius;
//        }

//        void diplay_settings();

//        G4String BeamingType() const;
//        void     set_BeamingType(const G4String &value);
//        G4bool MAG_FIELD_ON() const;
//        void set_MAG_FIELD_ON(const G4bool MAG_FIELD_BOOL);
//        G4bool USE_STEP_MAX() const;
// };
