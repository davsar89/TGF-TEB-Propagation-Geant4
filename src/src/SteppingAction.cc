#include <utility>

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

#include <SteppingAction.hh>

// If you want to identify the first step in a volume, pick fGeomBoudary status in PreStepPoint.
// If you want to identify a step getting out of a volume, pick fGeomBoundary status in PostStepPoint

SteppingAction::SteppingAction()
{

    earth = new GeographicLib::Geocentric(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());

    //    asciiFileName = "./extra_outputs/created_electrons";
    //    std::ofstream asciiFile00(asciiFileName, std::ios::trunc); // to clean the output file
    //    asciiFile00.close();
    //    asciiFileName_phot = "./extra_outputs/created_photons";
    //    std::ofstream asciiFile11(asciiFileName_phot, std::ios::trunc); // to clean the output file
    //    asciiFile11.close();

    min_max_alt = find_min_max_rough_altitude();
}

///

SteppingAction::~SteppingAction()
{
    delete analysis;
    delete earth;
}

///

void SteppingAction::UserSteppingAction(const G4Step *aStep)
{

    G4Track *track = aStep->GetTrack();
    const G4StepPoint *preStepPoint = aStep->GetPreStepPoint();

    const G4int PDG = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

    //	if (!((PDG == PDG_electron) || (PDG == PDG_positron) || (PDG == PDG_photon))) return;

#ifndef NDEBUG // debug mode only
    if (Settings::USE_STEP_MAX_for_record)
    {

        if (preStepPoint)
        {

            const G4TouchableHandle &preStepTouch = preStepPoint->GetTouchableHandle();
            const G4VPhysicalVolume *volume = preStepTouch->GetVolume();

            const G4LogicalVolume *lVolume = volume->GetLogicalVolume();
            const G4Region *region = lVolume->GetRegion();
            const G4String reg_name = region->GetName();

            if (reg_name == "RECORD_REGION")
            {
                //                const double rough_altitude = (aStep->GetPreStepPoint()->GetPosition().mag() - earth_radius) / km;
                //                if (PDG==22) G4cout << aStep->GetStepLength()/m << G4endl;
                if ((PDG == PDG_electron) || (PDG == PDG_positron) || (PDG == PDG_photon))
                {
                    if ((aStep->GetStepLength() > Settings::STEP_MAX_RECORD_AREA))
                    {
                        G4cout << "Current step length : " << aStep->GetStepLength() / meter << " meter" << G4endl;
                        G4cout << "Error in SteppingAction : step is lager than the maximum allowed." << G4endl;
                        std::abort();
                    }
                }
            }
        }
    }
#endif // end debug mode only
    //

    G4double global_time = track->GetGlobalTime();

    if (!Settings::RECORD_PHOT_ONLY)
    {
        if (global_time > Settings::TIME_LIMIT)
        {
            track->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    if (track->GetKineticEnergy() < Settings::MIN_ENERGY_OUTPUT)
    { // killing particles with energy lower than threshold
        if (PDG != PDG_positron)
        { // don't kill positrons to make sure they do annihilation before disappearing
            track->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    // #ifndef NDEBUG // debug mode only
    //     if (Settings::CHECK_ELECTRON_TRACKING_IN_MAG_FIELD) {
    //  for leptons : check step size compared to Larmor radius
    if (!Settings::RECORD_PHOT_ONLY)
    {
        if (PDG == PDG_electron || PDG == PDG_positron)
        {

            if (track->GetKineticEnergy() >= Settings::MIN_ENERGY_OUTPUT)
            {

                const G4double rough_altitude_km = (aStep->GetPreStepPoint()->GetPosition().mag() - earth_radius) / km;

                if (rough_altitude_km > 30.)
                {

                    if ((previous_pos - track->GetPosition()).mag() > 200.0 * meter)
                    {
                        double larmor_radius_m = 0;

                        // the global field manager is actually a myG4FieldManager type (set un FieldSetup.cc)
                        // so it can be casted to this type, because we want to use the method to calculate the Larmor radius
                        // and to set (update) the larmor_radius member variable
                        myG4FieldMan = dynamic_cast<myG4FieldManager *>(G4TransportationManager::GetTransportationManager()->GetFieldManager());

                        larmor_radius_m = myG4FieldMan->CalculateLarmorRadius(preStepPoint);

#ifndef NDEBUG // debug mode only
                        const double step_length_m = aStep->GetStepLength() / m;

                        if (step_length_m > larmor_radius_m / 8.0)
                        {
                            G4cout << rough_altitude_km << G4endl;
                            G4cout << step_length_m << G4endl;
                            G4cout << larmor_radius_m / 8.0 << G4endl;
                            G4cout << "ERROR: step size is lager than 1/8 of Larmor radius, the integration has too much error. Aborting." << G4endl;
                            std::abort();
                        }
#endif // end debug mode only

                        // affects only charged particles
                        G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(
                            larmor_radius_m * m / 10.0);

                        myG4FieldMan->Set_Larmor_Radius(larmor_radius_m);
                    }
                }
            }
        }
    }
    //    }
    // #endif // end debug mode only
    // kill photons that reached a too high altitude
    // kill also electrons and positrons if the setting RECORD_PHOT_ONLY is set.
    if (PDG == PDG_photon || ((PDG == PDG_electron || PDG == PDG_positron) && Settings::RECORD_PHOT_ONLY))
    {

        const double rough_altitude = aStep->GetPreStepPoint()->GetPosition().mag() - earth_radius;

        if (rough_altitude > photon_max_altitude)
        {
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    if (Settings::RECORD_ELEC_POSI_ONLY)
    {
        if (PDG == PDG_electron || PDG == PDG_positron)
        {
            record_particles(aStep, PDG);
        }
    }
    else if (Settings::RECORD_PHOT_ONLY)
    {
        if (PDG == PDG_photon)
        {
            record_particles(aStep, PDG);
        }
    }
    else
    {
        record_particles(aStep, PDG);
    }
}

///////////////////////////////////////////////////////////////////////////

void SteppingAction::record_particles(const G4Step *aStep, const int PDG)
{

    const double rough_altitude = aStep->GetPreStepPoint()->GetPosition().mag() - earth_radius;

    if (rough_altitude > min_max_alt.min_val * km * 0.9 && rough_altitude < min_max_alt.max_val * km * 1.1)
    {

        const G4Track *track = aStep->GetTrack();

        const G4int ID = track->GetTrackID();

        G4StepPoint *thePrePoint = aStep->GetPreStepPoint();
        double ecef_x = thePrePoint->GetPosition().x() / m;
        double ecef_y = thePrePoint->GetPosition().y() / m;
        double ecef_z = thePrePoint->GetPosition().z() / m;
        double pre_lat, pre_lon, pre_alt;
        //        geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x, ecef_y, ecef_z, pre_lat, pre_lon, pre_alt);
        earth->Reverse(ecef_x, ecef_y, ecef_z,
                       pre_lat, pre_lon, pre_alt);

        G4StepPoint *thePostPoint = aStep->GetPostStepPoint();
        double ecef_x2 = thePostPoint->GetPosition().x() / m;
        double ecef_y2 = thePostPoint->GetPosition().y() / m;
        double ecef_z2 = thePostPoint->GetPosition().z() / m;
        double post_lat, post_lon, post_alt;

        // skip if it is a point that is fixed (should not happen, but just in case)
        if (ecef_x == ecef_x2 && ecef_y == ecef_y2 && ecef_z == ecef_z2)
        {
            return;
        }

        //        geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x2, ecef_y2, ecef_z2, post_lat, post_lon, post_alt);
        earth->Reverse(ecef_x2, ecef_y2, ecef_z2,
                       post_lat, post_lon, post_alt);

        if ((pre_alt / 1000.0 <= Settings::record_altitude && post_alt / 1000.0 > Settings::record_altitude) ||
            (pre_alt / 1000.0 > Settings::record_altitude && post_alt / 1000.0 <= Settings::record_altitude))
        {

            const double energy = thePrePoint->GetKineticEnergy();

            const double delta_time = project_to_record_alt(
                PDG,                        // input
                energy,                     // input
                ecef_x, ecef_y, ecef_z,     // modified (input and output)
                pre_lat, pre_lon, pre_alt); // input

            double pre_lat2, pre_lon2, pre_alt2;
            //            geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x, ecef_y, ecef_z, // inputs
            //                                                        pre_lat2, pre_lon2, pre_alt2); // outputs
            earth->Reverse(ecef_x, ecef_y, ecef_z,
                           pre_lat2, pre_lon2, pre_alt2);

            double dist_rad = Get_dist_rad(pre_lat2, pre_lon2,
                                           pre_alt2 * meter);

            const double time =
                (thePrePoint->GetGlobalTime() + delta_time * second) / microsecond;

            const double mom_x = thePrePoint->GetMomentumDirection().getX();
            const double mom_y = thePrePoint->GetMomentumDirection().getY();
            const double mom_z = thePrePoint->GetMomentumDirection().getZ();

            G4int event_nb = 0;
            const G4Event *evt = G4RunManager::GetRunManager()->GetCurrentEvent();
            if (evt)
                event_nb = evt->GetEventID();

            analysis->save_in_output_buffer(
                PDG, time, energy / keV, dist_rad / km, ID, ecef_x / 1000.0,
                ecef_y / 1000.0, ecef_z / 1000.0, // from m to km
                mom_x, mom_y, mom_z, pre_lat2, pre_lon2, pre_alt2, event_nb);

            // kill after record if photon, to make 100% sure there is no double record
            if (PDG == PDG_photon)
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);

            if (Settings::RECORD_PHOT_ONLY)
            {
                if (PDG == PDG_electron || PDG == PDG_positron)
                    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////

G4double SteppingAction::Get_dist_rad(const G4double &lat, const G4double &lon,
                                      const G4double &alt_record)
{
    // getting the radial distance (along curve parrallel to Earth surface)
    // Equirectangular approximation -> leads to "effective speed" of the output
    // data that can be slightly greater than the speed of light
    //    G4double RR = (Settings::EarthRadius() + alt_record) ;
    //    G4double delta_lambda = (Settings::SourceLong() - lon) * degree;
    //    G4double sum_phi = (lat + Settings::SourceLat()) * degree;
    //
    //    G4double xx = delta_lambda * cos((sum_phi) / 2.);
    //    G4double yyy = (Settings::SourceLat() - lat) * degree;
    //    G4double dist_rad = sqrt(pow(xx, 2) + pow(yyy, 2)) * RR;
    // haversine formula (better)
    const G4double RR = (Settings::earthRadius + alt_record);
    const G4double phi1 = (lat)*degree;
    const G4double phi2 = (Settings::SOURCE_LAT)*degree;
    const G4double delta_phi = (phi2 - phi1);
    const G4double sin_delta_phi_over_2 = sin(delta_phi / 2.);
    const G4double delta_lambda = (Settings::SOURCE_LONG - lon) * degree;
    const G4double sin_delta_lambda_over_2 = sin(delta_lambda / 2.);
    const G4double aa =
        sin_delta_phi_over_2 * sin_delta_phi_over_2 +
        cos(phi1) * cos(phi2) * sin_delta_lambda_over_2 * sin_delta_lambda_over_2;
    const G4double cc = 2.0 * atan2(sqrt(aa), sqrt(1.0 - aa));
    // G4double cc = 2. * std::asin(std::min(1., sqrt(aa)));
    const G4double dist_rad = RR * cc;
    return dist_rad;
}

///////////////////////////////////////////////////////////////////////////

SteppingAction::min_max_altitudes SteppingAction::find_min_max_rough_altitude()
{

    grid_latitude = linspace(-90.0, 90.0, nb);
    grid_longitude = linspace(-180.0, 180.0, nb);

    const double target_alt = Settings::record_altitude;

    double rough_alt;

    double min_rough_alt = 1000.0;
    double max_rough_alt = 0.0;

    for (int i_lat = 0; i_lat < nb; ++i_lat)
    {
        for (int i_lon = 0; i_lon < nb; ++i_lon)
        {

            double lat = grid_latitude[i_lat];
            double lon = grid_longitude[i_lon];

            rough_alt = find_rough_altitude_for_real_altitude(
                lat,
                lon,
                target_alt);

            if (rough_alt < min_rough_alt)
            {
                min_rough_alt = rough_alt;
            }
            if (rough_alt > max_rough_alt)
            {
                max_rough_alt = rough_alt;
            }
        }
    }

    min_max_altitudes min_max_alt_out;
    min_max_alt_out.max_val = max_rough_alt;
    min_max_alt_out.min_val = min_rough_alt;
    return min_max_alt_out;
    //    G4cout << min_rough_alt << G4endl;
}

///////////////////////////////////////////////////////////////////////////

double SteppingAction::find_rough_altitude_for_real_altitude(const double &lat, const double &lon, const double &target_alt)
{

    double ecef_x, ecef_y, ecef_z;
    G4ThreeVector ECEF_position{0, 0, 0};

    //    geod_conv::GeodeticConverter::geodetic2ecef(lat, lon, target_alt * 1000.0, ecef_x, ecef_y, ecef_z);
    earth->Forward(lat, lon, target_alt * 1000.0, ecef_x, ecef_y, ecef_z);
    ECEF_position.setX(ecef_x);
    ECEF_position.setY(ecef_y);
    ECEF_position.setZ(ecef_z); // in meters

    const double rough_altitude_km = (ECEF_position.mag() - earth_radius / m) / 1000.0; // in km

    return rough_altitude_km;
}

///////////////////////////////////////////////////////////////////////////

double SteppingAction::project_to_record_alt(
    const G4int type_number, const double &kinE,
    double &ecef_x, double &ecef_y, double &ecef_z,
    const double &lat_in, const double &lon_in, const double &alt_in)
// input and output lengths are in meters
// input energy is in keV
// output time is in seconds
{
    double beta;
    double gammaa;
    double sin_lon, cos_lon, sin_lat, cos_lat;

#if defined(__linux__) && defined(__GNUG__) // if linux and GCC
    sincos(lat_in * degree, &sin_lat, &cos_lat);
    sincos(lon_in * degree, &sin_lon, &cos_lon);
#else
    sin_lat = sin(lat_in * degree);
    cos_lat = cos(lat_in * degree);
    sin_lon = sin(lon_in * degree);
    cos_lon = cos(lon_in * degree);
#endif

    const double local_vertical_x = cos_lon * cos_lat;
    const double local_vertical_y = sin_lon * cos_lat;
    const double local_vertical_z = sin_lat;
    double delta_alt = alt_in - ALT_RECORD_in_meter; // in meters
    ecef_x += -local_vertical_x * delta_alt;
    ecef_y += -local_vertical_y * delta_alt;
    ecef_z += -local_vertical_z * delta_alt;
    double delta_time = 0;

    if (type_number == 22)
    {
        delta_time = -delta_alt / sol_SI;
        beta = 1.0;
        gammaa = 1.e6;
    }
    else
    {
        gammaa = kinE / (510.9989461 * keV) + 1.0;
        beta = std::sqrt(1.0 - std::pow(gammaa, -2.0));
        delta_time = -delta_alt / (beta * sol_SI);
    }

#ifndef NDEBUG // debug mode

    if (std::isnan(alt_in))
    {
        G4cout << G4endl
               << "Error: alt_in is Nan in SensitiveDet::project_to_record_alt. "
                  "Aborting"
               << G4endl;
        std::abort();
    }

#endif // ifndef NDEBUG
    return delta_time;
}
