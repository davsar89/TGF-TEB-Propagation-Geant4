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

#include "PrimaryGeneratorAction.hh"
#include "EarthMagField_WMM.hh"

using namespace std;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField_WMM::EarthMagField_WMM() {

    /* Memory allocation */
    earth = new GeographicLib::Geocentric(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());

    strncpy(VersionDate, VERSIONDATE_LARGE + 39, 11);
    VersionDate[11] = '\0';

    char filename[] = "./mag_data/WMM.COF";

    if (!MAG_robustReadMagModels(filename, &MagneticModels, epochs)) {
        printf("\n WMM.COF not found.  Aborting... \n ");
        std::abort();
    }

    if (nMax < MagneticModels[0]->nMax)
        nMax = MagneticModels[0]->nMax;
    NumTerms = ((nMax + 1) * (nMax + 2) / 2);
    TimedMagneticModel = MAG_AllocateModelMemory(NumTerms); /* For storing the time modified WMM Model parameters */

    MAG_SetDefaults(&Ellip, &Geoid); /* Set default values and constants */

    Geoid.Geoid_Initialized = 1;
    /* Set EGM96 Geoid parameters END */

    UserDate.Year = Settings::dt_year;
    UserDate.Day = int(Settings::dt_day);
    UserDate.Month = int(Settings::dt_month);
    UserDate.DecimalYear = double(UserDate.Year) + double(UserDate.Month) / 12.0 + double(UserDate.Day) / 365.24;

    MAG_TimelyModifyMagneticModel(UserDate, MagneticModels[0],
                                  TimedMagneticModel); /* Time adjust the coefficients, Equation 19, WMM Technical report */

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField_WMM::~EarthMagField_WMM() {
    delete earth;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EarthMagField_WMM::GetFieldValue(const double Point[3], double *Bfield) const {
    //  geodetic_converter::GeodeticConverter g_geodetic_converter;

    xx = Point[0] / m;
    yy = Point[1] / m;
    zz = Point[2] / m;

    rough_alt = sqrt(xx * xx + yy * yy + zz * zz) - earth_radius_m;

    if (rough_alt < 30000.0) {
        Bfield[0] = 0;
        Bfield[1] = 0;
        Bfield[2] = 0;
        return;
    }

    const G4ThreeVector mag_field = GetFieldComponents(G4ThreeVector{xx, yy, zz});

    Bfield[0] = mag_field.getX() * tesla;
    Bfield[1] = mag_field.getY() * tesla;
    Bfield[2] = mag_field.getZ() * tesla;

} // EarthMagField_WMM::GetFieldValue

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// input ECEF coordinates in meters
// output field in TESLA
G4ThreeVector EarthMagField_WMM::GetFieldComponents(const G4ThreeVector &position_ECEF) const {

    // conversion ECEF to geodetic

    // input x y z in meters, output altitude in km, output lat lon in degrees
//    geod_conv::GeodeticConverter::ecef2Geodetic(position_ECEF.x(), position_ECEF.y(), position_ECEF.z(),
//                                                lat, lon, alt);

    earth->Reverse(position_ECEF.getX(), position_ECEF.getY(), position_ECEF.getZ(),
                   lat, lon, alt);

    elong = lon;                 // degrees
    colat = (180.0 / 2.0 - lat); // degrees
    alt_km = alt / 1000.0;

    //

    CoordGeodetic.lambda = lon;
    CoordGeodetic.phi = lat;
    CoordGeodetic.HeightAboveEllipsoid = alt_km;

    //    if (MAG_GetUserInput(MagneticModels[0], &Geoid, &CoordGeodetic, &UserDate) == 1) /*Get User Input */
    //    {

    MAG_GeodeticToSpherical(Ellip, CoordGeodetic,
                            &CoordSpherical); /*Convert from geodetic to Spherical Equations: 17-18, WMM Technical report*/

    MAG_Geomag(Ellip, CoordSpherical, CoordGeodetic, TimedMagneticModel,
               &GeoMagneticElements); /* Computes the geoMagnetic field elements and their time change*/
    //    }

    Bx = GeoMagneticElements.X; // X Y Z output in nanoTesla
    By = GeoMagneticElements.Y;
    Bz = GeoMagneticElements.Z;

    // NED (North East Down) to ECEF convertion

    mag_field_ecef = myUtils::ned2ecef(Bx, By, Bz, lat, lon);

    G4ThreeVector field_out{mag_field_ecef.u * 1e-9,
                            mag_field_ecef.v * 1e-9,
                            mag_field_ecef.w * 1e-9};// output in TESLA

#ifndef NDEBUG // debug mode only
    if (std::isnan(field_out.x()) || std::isinf(field_out.x())
        || std::isnan(field_out.y()) || std::isinf(field_out.y())
        || std::isnan(field_out.z()) || std::isinf(field_out.z())) {
        std::abort();
    }
#endif // end debug mode only

    return field_out;
}