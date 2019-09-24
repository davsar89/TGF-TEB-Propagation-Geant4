#include "SD.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"

// ONE DIFFERENT SENSITIVE DETECTOR IS SET UP BY RECORD ALTITUDE

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDet::SensitiveDet(const G4String &name, const G4int ID,
                           const G4double &ALT_RECORD_km)
        : G4VSensitiveDetector(name) {
    this->ID_SD = ID;
    this->ALT_RECORD_in_km = ALT_RECORD_km;
    this->ALT_RECORD_in_meter = ALT_RECORD_km * 1000.0;
    // initialization
    recorded_ID_inside_event_and_SD.clear();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDet::~SensitiveDet() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDet::Initialize(
        G4HCofThisEvent *) // executed at begin of each event
{
    // if it is a new event, clear the recorded ID vector
//    recorded_ID_inside_event_and_SD.clear();
//    nb_event++;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SensitiveDet::ProcessHits(G4Step *aStep,
                                 G4TouchableHistory * /*ROhist*/) {
//    given_altitude_particle_record(aStep);
    return true;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDet::EndOfEvent(G4HCofThisEvent * /*HCE*/) {
    // RK : see EndOfEventAction method in EventAction.cc
    recorded_ID_inside_event_and_SD.clear(); // redundant, but jsut to be safe
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDet::given_altitude_particle_record(const G4Step *step) {
//    G4Track *track = step->GetTrack();
//    const G4int type_number = track->GetParticleDefinition()->GetPDGEncoding();
//
//    if (settings->RECORD_ELEC_POSI_ONLY) {
//        if (type_number == 22) {
//            return;
//        } // WARNING: no elec posi is recorded
//    }
//
//    if (settings->RECORD_PHOT_ONLY) {
//        if (type_number != 22) {
//            return;
//        } // WARNING: no photon is recorded
//    }
//
//    //    const G4int PDG_nb =
//    //    aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
//
//    if (!((type_number == 22) || (type_number == 11) || (type_number == -11))) {
//        return;
//    }
//
//    thePrePoint = step->GetPreStepPoint();
//
//    if (thePrePoint->GetStepStatus() ==
//        fGeomBoundary) { // if the particle has just entered the volume
//        //            const G4double geodetic_alt_record = ALT_RECORD_in_km * km;
//        const G4int ID = track->GetTrackID();
//        //            if (IDpart_not_recorded_yet(ID))
//        //                {
//        // getting geodetic coordinates
//        G4double ecef_x = thePrePoint->GetPosition().x() / m;
//        G4double ecef_y = thePrePoint->GetPosition().y() / m;
//        G4double ecef_z = thePrePoint->GetPosition().z() / m;
//        G4double pre_lat, pre_lon, pre_alt;
//        geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x, ecef_y, ecef_z, pre_lat,
//                                                    pre_lon, pre_alt);
//        G4double delta_time =
//                project_to_record_alt(type_number, track->GetKineticEnergy(), ecef_x,
//                                      ecef_y, ecef_z, pre_alt, pre_lat, pre_lon);
//        G4double dist_rad = -1.0;
//        dist_rad =
//                Get_dist_rad(pre_lat, pre_lon,
//                             pre_alt * meter); // pre_alt is converted to geant4 unit
//        // PDG encoding is 22 for gamma, 11 for electrons, -11 for positrons
//        const G4double energy = thePrePoint->GetKineticEnergy() / keV;
//        const G4double time =
//                (thePrePoint->GetGlobalTime() + delta_time) / microsecond;
//        //    if ((type_number != gamma_number) || (energy <= 0.))
//        //        {
//        //            return;
//        //        }
//        //            G4String creator_process = "initial";
//        //            if (track->GetCreatorProcess())
//        //                {
//        //                    creator_process =
//        //                    track->GetCreatorProcess()->GetProcessName();
//        //                }
//        //                    if (creator_process == "annihil" && energy > 511.0)
//        //                    return;
//        const double mom_x = thePrePoint->GetMomentumDirection().getX();
//        const double mom_y = thePrePoint->GetMomentumDirection().getY();
//        const double mom_z = thePrePoint->GetMomentumDirection().getZ();
//        recorded_ID_inside_event_and_SD.push_back(ID);
//        analysis->save_in_output_buffer(
//                type_number, time, energy, dist_rad / km, ID, ecef_x / 1000.0,
//                ecef_y / 1000.0, ecef_z / 1000.0, // from m to km
//                mom_x, mom_y, mom_z, thePrePoint->GetVelocity());
//    }

    //        }
}

//////////////////

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SensitiveDet::getID_SD() const { return ID_SD; }

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SensitiveDet::IDpart_not_recorded_yet(G4int ID) {
    if (recorded_ID_inside_event_and_SD.empty()) {
        return true;
    }

    return !(std::find(recorded_ID_inside_event_and_SD.begin(),
                       recorded_ID_inside_event_and_SD.end(),
                       ID) != recorded_ID_inside_event_and_SD.end());
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SensitiveDet::Get_dist_rad(const G4double &lat, const G4double &lon,
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

// // linear interpolation to get the coordinates at the output altitude
// (alt_tmp)

// added delta in time
G4double SensitiveDet::project_to_record_alt(
        const G4int type_number, const G4double &kinE, G4double &ecef_x,
        G4double &ecef_y, G4double &ecef_z, const G4double &alt_in,
        const G4double &lat_in, const G4double &lon_in)

// input and output are in meters
{
    G4double sin_lon, cos_lon, sin_lat, cos_lat;
    sincos(lon_in * degree, &sin_lon, &cos_lon);
    sincos(lat_in * degree, &sin_lat, &cos_lat);
    const G4double local_vertical_x = cos_lon * cos_lat;
    const G4double local_vertical_y = sin_lon * cos_lat;
    const G4double local_vertical_z = sin_lat;
    G4double delta_alt = alt_in - ALT_RECORD_in_meter;
    ecef_x += -local_vertical_x * delta_alt;
    ecef_y += -local_vertical_y * delta_alt;
    ecef_z += -local_vertical_z * delta_alt;
    G4double delta_time = 0;

    if (type_number == 22) {
        delta_time = -delta_alt / c_light * second;
    } else {
        G4double gamma = kinE / (510.9989461 * keV) + 1.0;
        G4double beta = std::sqrt(1.0 - std::pow(gamma, -2.0));
        delta_time = -delta_alt / (beta * c_light) * second;
    }

#ifndef NDEBUG // debug mode

    if (std::isnan(alt_in)) {
        G4cout << G4endl
               << "Error: alt_in is Nan in SensitiveDet::project_to_record_alt. "
                  "Aborting"
               << G4endl;
        std::abort();
    }

#endif // ifndef NDEBUG
    return delta_time;
}
