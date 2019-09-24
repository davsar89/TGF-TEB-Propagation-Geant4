#pragma once

#include "G4Track.hh"
#include "G4VSensitiveDetector.hh"
#include <fstream>
#include "Settings.hh"
#include "Analysis.hh" // singleton
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "geodetic_converter.hh"

extern "C" {
#include <C/coordinates_conversions.h>
}

using namespace std;

class G4Step;

class G4HCofThisEvent;

class G4TouchableHistory;

class SensitiveDet : public G4VSensitiveDetector {
public:

    SensitiveDet(const G4String &name, const G4int ID, const G4double &ALT_RECORD_km);

    ~SensitiveDet() override;

    void Initialize(G4HCofThisEvent *HCE) override;

    G4bool ProcessHits(G4Step *aStep,
                       G4TouchableHistory *ROhist) override;

    void EndOfEvent(G4HCofThisEvent *HCE) override;

    G4int getID_SD() const;

private:

    Settings *settings = Settings::getInstance();

    int evtctr = 0;
    G4StepPoint *thePrePoint = nullptr;

    Analysis *analysis = Analysis::getInstance();

    void given_altitude_particle_record(const G4Step *step);

    const G4int PDG_phot = 22;
    const G4int PDG_elec = 11;
    const G4int PDG_posi = -11;

    G4int ID_SD = 0;
    G4double ALT_RECORD_in_km = 0;
    G4double ALT_RECORD_in_meter = 0;

    std::vector<int> recorded_ID_inside_event_and_SD;

    bool IDpart_not_recorded_yet(G4int ID);

    G4double Get_dist_rad(const G4double &lat,
                          const G4double &lon,
                          const G4double &alt_record);

    G4int nb_event = 0;

    G4double project_to_record_alt(const G4int type_number,
                                   const G4double &kinE,
                                   G4double &ecef_x,
                                   G4double &ecef_y,
                                   G4double &ecef_z,
                                   const G4double &alt_in,
                                   const G4double &lat_in,
                                   const G4double &lon_in);

    G4double c_light = 299792458.0; // m/s
};
