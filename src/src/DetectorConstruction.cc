////////////////////////////////////////////////////////////////////////////////

// /* GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's environment */
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

#include <DetectorConstruction.hh>
//#include <fortran.hh>

using namespace std;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// COPY PASTED COMMENT FROM MSIS FORTRAN CODE
//C     INPUT VARIABLES:
//C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
//C              (Year ignored in current model)
//C        SEC - UT(SEC)
//C        ALT - ALTITUDE(KM)
//C        GLAT - GEODETIC LATITUDE(DEG)
//C        GLONG - GEODETIC LONGITUDE(DEG)
//C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
//C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
//C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
//C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
//C           - ARRAY CONTAINING:
//C             (1) DAILY AP
//C             (2) 3 HR AP INDEX FOR CURRENT TIME
//C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
//C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
//C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
//C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
//C                    TO CURRENT TIME
//C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
//C                    TO CURRENT TIME
//C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
//C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
//C                 MASS 17 IS Anomalous O ONLY.)
//C
//C     NOTES ON INPUT VARIABLES:
//C        UT, Local Time, and Longitude are used independently in the
//C        model and are not of equal importance for every situation.
//C        For the most physically realistic calculation these three
//C        variables should be consistent (STL=SEC/3600+GLONG/15).
//C        The Equation of Time departures from the above formula
//C        for apparent local time can be included if available but
//C        are of minor importance.
//c
//C        F107 and F107A values used to generate the model correspond
//C        to the 10.7 cm radio flux at the actual distance of the Earth
//C        from the Sun rather than the radio flux at 1 AU. The following
//C        site provides both classes of values:
//C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
//C
//C        F107, F107A, and AP effects are neither large nor well
//C        established below 80 km and these parameters should be set to
//C        150., 150., and 4. respectively.
//C
//C     OUTPUT VARIABLES:
//C        D(1) - HE NUMBER DENSITY(CM-3)
//C        D(2) - O NUMBER DENSITY(CM-3)
//C        D(3) - N2 NUMBER DENSITY(CM-3)
//C        D(4) - O2 NUMBER DENSITY(CM-3)
//C        D(5) - AR NUMBER DENSITY(CM-3)
//C        D(6) - TOTAL MASS DENSITY(GM/CM3)
//C        D(7) - H NUMBER DENSITY(CM-3)
//C        D(8) - N NUMBER DENSITY(CM-3)
//C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
//C        T(1) - EXOSPHERIC TEMPERATURE
//C        T(2) - TEMPERATURE AT ALT

// IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T

// extrernal fortran subroutine to get MSIS atmospheric densities
//extern "C" {
//void gtd7_(INTEGER &IYD, // YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
//           REAL &SEC, // UT(SEC)
//           REAL &ALT, // ALTITUDE(KM)
//           REAL &GLAT, // GEODETIC LATITUDE(DEG)
//           REAL &GLONG, // GEODETIC LONGITUDE(DEG)
//           REAL &STL, // LOCAL APPARENT SOLAR TIME
//           REAL &F107A, // 81 day AVERAGE OF F10.7 FLUX (centered on day DDD
//           REAL &F107, // DAILY F10.7 FLUX FOR PREVIOUS DAY
//           REAL &AP,  // MAGNETIC INDEX(DAILY)
//           INTEGER &MASS, // MASS NUMBER
//           REAL *D, REAL *T); // OUTPUT VARIABLES temperatures
//}
// END COPY PASTED COMMENT FROM MSIS FORTRAN CODE

typedef unsigned int uint;

TGFDetectorConstruction::TGFDetectorConstruction() {
    logicalWorld = nullptr;
    physicalWorld = nullptr;

    globalfieldMgr = nullptr;

    if (settings->OUTPUT_ALT_LAYERS_TO_FILE) {
        asciiFile.open("alt_dens_debug.txt", std::ios::trunc);

        asciiFile << "altitude (km) // density (g/cm2)" << G4endl;
    }

}

TGFDetectorConstruction::~TGFDetectorConstruction() = default;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume *TGFDetectorConstruction::Construct() {
    //    G4FieldManager *null_field = nullptr;

    // cleaning geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    Airs.clear();
    Airs = GENERATE_AIR_MATERIALS(sea_level_density, km_150_density, number_of_AIR_materials);

    vaccum = man->FindOrBuildMaterial("G4_Galactic");

    // construct field managers (for magnetic fields)
    if (settings->MAG_FIELD_ON) {
        Construct_MagField_Managers();
    }

    // generating the world layers of constant density
    calculate_radii_list();

    G4NistManager *man = G4NistManager::Instance();
    vac = man->FindOrBuildMaterial("G4_Galactic");

    // World solid
    G4Sphere *solidWorld;
    solidWorld = new G4Sphere("world_S", base_radius, (base_radius + world_max_altitude),
                              0 * degree, 360 * degree, 0 * degree, 180 * degree);
// World logical

    logicalWorld = new G4LogicalVolume(solidWorld,  // solid
                                       vac,         // material
                                       "world_L");

    if (settings->MAG_FIELD_ON) {
        logicalWorld->SetFieldManager(globalfieldMgr, true);
    }

    // Physical volume
    physicalWorld = new G4PVPlacement(nullptr, G4ThreeVector(), "world_P", // name (2nd constructor)
                                      logicalWorld,             // logical volume
                                      nullptr,                      // mother volume
                                      false,              // no boolean operation
                                      0);                          // copy number

    G4VisAttributes *VisAttWorld = new G4VisAttributes(G4Colour(204 / 255., 255 / 255., 255 / 255.));
    logicalWorld->SetVisAttributes(VisAttWorld);

    // setting default (world) region info
    G4Region *defaultRegion = (*(G4RegionStore::GetInstance()))[0]; // the default (world) region is index 0 in the region store
    auto *defaultRInfo = new RegionInformation();
    defaultRInfo->Set_World();
    defaultRegion->SetUserInformation(defaultRInfo);

    // Make Invisible
    //  logicalWorld -> SetVisAttributes(G4VisAttributes::Invisible);

    G4double innerRad = 0;
    G4double outerRad = 0;

    G4double DELTA_PHI = 10. * degree;

    G4double DELTA_THETA = 10. * degree;

    int nb_total_angles = int((360.0 * degree) / DELTA_PHI) * int(180.0 * degree / DELTA_THETA);

    // atmosphere construction

    for (int i_phi = 0; i_phi < int((360.0 * degree) / DELTA_PHI); ++i_phi) {
        for (int i_theta = 0; i_theta < int(180.0 * degree / DELTA_THETA); ++i_theta) {
            // loop on "altitudes"
            for (int jj = 0; jj < radius_list.size() - 1; jj++) {
//    for (int i_phi = 1; i_phi < 2; ++i_phi) {
//        for (int i_theta = 6; i_theta < 7; ++i_theta) {
//            for (unsigned int jj = 200; jj < 201; jj++) {

                // here "radius" is the radius of spheres in the ECEF space that should overlap the geodetic altitudes
                innerRad = base_radius + radius_list[jj];
                outerRad = base_radius + radius_list[jj + 1];

                double center_theta = DELTA_THETA * i_theta + DELTA_THETA / 2.0;
                double center_phi = i_phi * DELTA_PHI + DELTA_PHI / 2.0;
                double avg_radius = (innerRad + outerRad) / 2.;

                double ecef_x = avg_radius * std::sin(center_theta) * std::cos(center_phi);
                double ecef_y = avg_radius * std::sin(center_theta) * std::sin(center_phi);
                double ecef_z = avg_radius * std::cos(center_theta);

                double ecef_x_in_m = ecef_x / m;
                double ecef_y_in_m = ecef_y / m;
                double ecef_z_in_m = ecef_z / m;

                double geod_lat = 0, geod_lon = 0, geod_alt = 0;
                geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x_in_m, ecef_y_in_m, ecef_z_in_m, geod_lat,
                                                            geod_lon, geod_alt);

                double geod_alt_m = geod_alt;

                // test box to see if conversion of G4Sphere radius/theta/phi to Geant4 X,Y,Z coordinates (that is used as ECEF) is OK
//                G4Box *Abox = new G4Box('test_box',100000*meter,100000*meter,100000*meter);
//                G4LogicalVolume* G4box_logical = new G4LogicalVolume(Abox, Airs[jj], "box_test_" + std::to_string(jj), nullptr,
//                                                    nullptr, nullptr);
//                G4VPhysicalVolume* positioned_box = new G4PVPlacement(nullptr, G4ThreeVector(ecef_x,ecef_y,ecef_z), "box_test_PV" + std::to_string(jj),
//                                                                      G4box_logical, physicalWorld, false, 0,
//                false);
                double geod_alt_km = geod_alt / 1000.0;

                // first altitude layer
                // check if we covered altitudes down to 0
                if (jj == 0) {
                    done_for_this_theta_phi = false;
                    double local_min_alt = geod_alt_km;
                    if (local_min_alt > 0.0) {
                        G4cout << "ERROR: first (lower) altitude for the theta/phi segment is above 0; so the atmosphere is incomplete. ABORTING.";
                        std::abort();
                    }
                }
                // last altitude layer
                // check if we covered all altitudes up to the highest record altitude
                // settings->record_altitudes must be in ascending order
//                if (jj == radius_list.size() - 2) {
//                    double local_alt = geod_alt_km;
//                    if (local_alt < settings->record_altitudes[settings->record_altitudes.size() - 1]) {
//                        G4cout << "ERROR: last (higher) altitude for the theta/phi segment is not reaching the recording altitude. ABORTING.";
//                        std::abort();
//                    }
//                }

                // skip some layers in the middle
                if (geod_alt_km > 150.0) {
                    continue;
                }

                atmosLayers_S.push_back(
                        new G4Sphere("atmosLayer_S_" + std::to_string(jj), innerRad, outerRad, i_phi * DELTA_PHI, DELTA_PHI,
                                     DELTA_THETA * i_theta, DELTA_THETA));

                if (geod_alt_km > 0.0 && geod_alt_km < 150.0) {

                    // this function takes input altitude in meters
                    int idx_mat = find_atmosphere_part_material(geod_lat, geod_lon, geod_alt_m);

                    atmosLayers_LV.push_back(
                            new G4LogicalVolume(atmosLayers_S.back(), Airs[idx_mat], "atmosphere_LV_" + std::to_string(jj), nullptr,
                                                nullptr, nullptr));
                }
                    // vaccum if altitude > ALT_MAX or altitude < 0
                else {
                    atmosLayers_LV.push_back(
                            new G4LogicalVolume(atmosLayers_S.back(), vac, "atmosphere_LV_" + std::to_string(jj), nullptr,
                                                nullptr, nullptr));
                }

                // assigning magnetic field
                if (settings->MAG_FIELD_ON) {
                    atmosLayers_LV.back()->SetFieldManager(globalfieldMgr, true);

                    // null magnetic field if altitude < 45 km
                    if (geod_alt_km < 35.0) {
                        atmosLayers_LV.back()->SetFieldManager(Null_FieldManager, true);
                    }
                }

                G4String name_PV = "atmosphere_PV_" + std::to_string(jj);
                atmosLayers_PV.push_back(
                        new G4PVPlacement(nullptr, G4ThreeVector(), name_PV, atmosLayers_LV.back(), physicalWorld, false, 0, false));

            }
        }
    }

//    if (id_SD == 0 || id_SD > nb_total_angles) {
//        G4cout << "ERROR with the number of sensitive detectors : there is none of them, or too much of them. ABORTING." << G4endl;
//        std::abort();
//    }

    G4cout << G4endl << "Geometry built successfully." << G4endl << G4endl;

    return physicalWorld;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool TGFDetectorConstruction::CHECK_ASCENDING(const std::vector<double> &alt_list) {

    if (alt_list.size() == 1) return true;

    for (int i = 0; i < alt_list.size() - 1; ++i) {

        if (alt_list[i] > alt_list[i + 1]) {
            return false;
        }

    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void TGFDetectorConstruction::calculate_radii_list()
// fills the vector altitudes
{

    if (!CHECK_ASCENDING(settings->record_altitudes)) {
        G4cout << "ERROR : settings->record_altitudes is not sorted from low to high. Aborting." << G4endl;
        std::abort();
    }

    const G4double alt_max_construction = 180.0 * km;

    // defining the altitude vector
    for (G4int jj = 0; jj < nb_altitudes; jj++) {
        radius_list.push_back(
                exp(log(alt_min) + (log(alt_max_construction) - log(alt_min)) * double(jj) / double(nb_altitudes - 1)));
    }

    if (hasDuplicates(radius_list)) { // ERROR if there is any duplicates
        G4cout << "ERROR : There are duplicates values in the altitude list. Aborting." << G4endl;
        std::abort();
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool TGFDetectorConstruction::hasDuplicates(const std::vector<G4double> &arr) {
    for (uint i = 0; i < arr.size(); ++i) {
        for (uint j = i + 1; j < arr.size(); ++j) {
            if (arr[i] == arr[j]) {
                return true;
            }
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
G4bool TGFDetectorConstruction::not_contains(G4double x, const std::vector<G4double> &v) {
    return !(std::find(v.begin(), v.end(), x) != v.end());
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Caltulating materials of the atmopsheric layers, based on the MSIS C++ model integrated to this code
// ref : https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html
// simplified version where it is standard G4_AIR material with varying density
// It is also possible to include all MSIS-given elements with the corresponding densities; but in practice it give very similar results and runs slower
int
TGFDetectorConstruction::find_atmosphere_part_material(double lat, double lon, double alt) {

    const double altitude_in_km = alt / 1000.0; // input alt is in meters

    struct nrlmsise_output output{};
    struct nrlmsise_input input{};
    struct nrlmsise_flags flags{};
    struct ap_array aph{};

    /* input values */
    for (int i = 0; i < 7; i++) {
        aph.a[i] = 100;
    }

    flags.switches[0] = 0;
    for (int i = 1; i < 24; i++) {
        flags.switches[i] = 1;
    }

    input.doy = int(datetools::day_of_year(settings->year, settings->month, settings->day));
    input.year = settings->year; /* without effect */
    input.sec = 29000;
    input.alt = altitude_in_km;
    input.g_lat = lat;
    input.g_long = lon;
    input.lst = 16;
    input.f107A = 150;
    input.f107 = 150;
    input.ap = 4;
    input.ap_a = &aph;

    gtd7(&input, &flags, &output);

    // G4cout << altitude << G4endl;

//            gtd7_(input_iyd, input_sec, input_alt, input_g_lat, input_g_long, input_lst, input_f107A, input_f107, input_ap, input_mass, output_D, output_T); // MSIS, fortran function call

    if (std::isnan(output.d[5]) || std::isinf(output.d[5])) {
        G4cout << "ERROR : density from gtd7_ is NaN of inf. Aborting" << G4endl;
        std::abort();
    }

    // getting density and converting it to the GEANT4 system of unit
    G4double density_air = output.d[5] * g / cm3;

    // find the index of the closest density air material

    std::vector<double> diffs(number_of_AIR_materials, 0);

    for (int jj = 0; jj < density_grid.size(); ++jj) {
        diffs[jj] = std::abs((density_air - density_grid[jj]) / density_grid[jj]);
    }

    int minElementIndex = std::min_element(diffs.begin(), diffs.end()) - diffs.begin();

//      Airs.push_back(man->BuildMaterialWithNewDensity("Air_" + std::to_string(idx_air_materials), "G4_AIR", density_air));

    return minElementIndex;

}

///////////////////////////////////////////////////////////////////////

std::vector<G4Material *> TGFDetectorConstruction::GENERATE_AIR_MATERIALS(double min_density, double max_density, int number) {


    const double sea_level_density2 = min_density * 1.3;
    const double km_150_density2 = max_density * 0.7;

//    density_grid = linspace(km_150_density2, sea_level_density2, number);
    density_grid = pyLogspace(std::log10(km_150_density2), std::log10(sea_level_density2), number);

    std::vector<G4Material *> Airs_out;
    Airs_out.reserve(number);

    for (int jj = 0; jj < density_grid.size(); ++jj) {

        double density_air = density_grid[jj];

        G4int ncomponents = 2;
        G4double fractionmass;
        Airs_out.push_back(new G4Material("Air_" + std::to_string(jj), density_air, ncomponents = 2));
        Airs_out.back()->AddElement(elN, fractionmass = 0.755 + 0.0128 / 2.0);
        Airs_out.back()->AddElement(elO, fractionmass = 0.2322 + 0.0128 / 2.0);

    }

    return Airs_out;

}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void TGFDetectorConstruction::Construct_MagField_Managers() {
    /////////// Magnetic field

    // myEarthMagField = new EarthMagField_alt;

    G4double distanceConst = settings->CACHED_LENGTH;

    globalfieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

    if (distanceConst == 0.0) {

#ifdef _WIN32
        G4cout << "Windows system detected : magnetic field model must be WMM." << G4endl;

        pMagFldEquation = new G4Mag_UsualEqRhs(new EarthMagField_WMM);
        fStepper = new G4ClassicalRK4(pMagFldEquation);
        fChordFinder = new G4ChordFinder(new EarthMagField_WMM, fMinStep, fStepper);
        globalfieldMgr->SetChordFinder(fChordFinder);
        globalfieldMgr->SetDetectorField(new EarthMagField_WMM);
#elif  __unix__
        if (settings->MAGNETIC_FIELD_MODEL == "IGRF") {
            G4cout << "Linux system detected, and magnetic field model set to IGRF." << G4endl;
            pMagFldEquation = new G4Mag_UsualEqRhs(new EarthMagField_IGRF);
            fStepper = new G4ClassicalRK4(pMagFldEquation);
            fChordFinder = new G4ChordFinder(new EarthMagField_IGRF, fMinStep, fStepper);
            globalfieldMgr->SetChordFinder(fChordFinder);
            globalfieldMgr->SetDetectorField(new EarthMagField_IGRF);
        } else if (settings->MAGNETIC_FIELD_MODEL == "WMM") {
            G4cout << "Linux system detected, and magnetic field model set to WMM." << G4endl;
            pMagFldEquation = new G4Mag_UsualEqRhs(new EarthMagField_WMM);
            fStepper = new G4ClassicalRK4(pMagFldEquation);
            fChordFinder = new G4ChordFinder(new EarthMagField_WMM, fMinStep, fStepper);
            globalfieldMgr->SetChordFinder(fChordFinder);
            globalfieldMgr->SetDetectorField(new EarthMagField_WMM);
        } else {
            G4cout << "ERROR : settings->MAGNETIC_FIELD_MODEL is not a string of value 'IGRF' or 'WMM' . Aborting."
                   << G4endl;
            std::abort();
        }
#endif

    } else {
#ifdef _WIN32
        myCachedEarthMagField = new G4CachedMagneticField(new EarthMagField_WMM, distanceConst);
#elif  __unix__
        if (settings->MAGNETIC_FIELD_MODEL == "IGRF") {
            myCachedEarthMagField = new G4CachedMagneticField(new EarthMagField_IGRF, distanceConst);
        } else if (settings->MAGNETIC_FIELD_MODEL == "WMM") {
            myCachedEarthMagField = new G4CachedMagneticField(new EarthMagField_WMM, distanceConst);
        } else {
            G4cout << "ERROR : settings->MAGNETIC_FIELD_MODEL is not a string of value 'IGRF' or 'WMM' . Aborting."
                   << G4endl;
            std::abort();
        }
#endif
        pMagFldEquation = new G4Mag_UsualEqRhs(myCachedEarthMagField);
        fStepper = new G4ClassicalRK4(pMagFldEquation);
        fChordFinder = new G4ChordFinder(myCachedEarthMagField, fMinStep, fStepper);
        globalfieldMgr->SetChordFinder(fChordFinder);
        globalfieldMgr->SetDetectorField(myCachedEarthMagField);
    }

    //    G4cout << "CACHED_LENGTH = " << settings->CachedLength() << G4endl;
    //    !!! RQ : avoid G4NystromRK4 , bug with use of G4CachedMagneticField
    //           avoid G4ConstRK4 and G4ImplicitEuler : bug with G4CachedMagneticField
    //               avoid G4DormandPrinceRK78
    //    G4DormandPrince745 is the best, from first tests
    //    CACHED_LENGTH at 1000 m gives bad results
    //   this was tested with G4 10.02 may be fixed with further updates.
    //
    globalfieldMgr->SetMinimumEpsilonStep(minEps);
    globalfieldMgr->SetMaximumEpsilonStep(maxEps);
    globalfieldMgr->SetDeltaOneStep(70. * cm);
    globalfieldMgr->SetDeltaIntersection(50. * cm);
    globalfieldMgr->GetChordFinder()->SetDeltaChord(50. * cm);

    //// NUll field manager
    ///
    Null_FieldManager = new G4FieldManager();
    magField_null = new G4UniformMagField(G4ThreeVector(0., 0., 0.));
    Null_FieldManager->SetDetectorField(magField_null);
    Null_FieldManager->CreateChordFinder(magField_null);
    Null_FieldManager->SetMinimumEpsilonStep(minEps);
    Null_FieldManager->SetMaximumEpsilonStep(maxEps);
    Null_FieldManager->SetDeltaOneStep(70. * cm);
    Null_FieldManager->SetDeltaIntersection(50. * cm);
    Null_FieldManager->GetChordFinder()->SetDeltaChord(50. * cm);

    // set maximum acceptable step everywhere, only for charged particles (i.e. affected by EM-field)
//    G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(
//            settings->STEP_MAX_DetConst);

}
