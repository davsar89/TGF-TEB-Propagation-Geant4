#include <G4ThreeVector.hh>
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

#if defined(__linux__) && defined(__GNUG__) // if linux and GCC

#include "PrimaryGeneratorAction.hh"
#include "EarthMagField_IGRF.hh"

using namespace std;

// interface with igrf12 Fortran code from : http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
extern "C" {
void igrf12syn_(int *isv, double *date, int *itype, double *alt, double *colat, double *elong, double *x, double *y,
                double *z, double *f, int *ier);
}

void cart2sphr(double x, double y, double z, double *r, double *theta, double *phi) {
    double hypotxy = hypot(x, y);
    *r = hypot(hypotxy, z);
    *theta = atan2(z, hypotxy) * 180.0 / CLHEP::pi;
    *phi = atan2(y, x) * 180.0 / CLHEP::pi;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField_IGRF::EarthMagField_IGRF() {

    date = double(Settings::dt_year) + double(Settings::dt_month) / 12.0 + double(Settings::dt_day) / 365.25;
    xx = 0;
    yy = 0;
    zz = 0;
    earth = new GeographicLib::Geocentric(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());

//#ifndef NDEBUG
//    run_sanity_check();
//#endif

//    G4int thread_ID = G4Threading::G4GetThreadId();
//
//    G4cout << "Thread ID" << thread_ID << G4endl;

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField_IGRF::~EarthMagField_IGRF() {
    delete earth;
};


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EarthMagField_IGRF::GetFieldValue(const double Point[3], double *Bfield) const {
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

//    G4ThreeVector EMF = GetFieldComponents(G4ThreeVector{xx, yy, zz}); // output in Tesla, input in meters
     G4ThreeVector EMF = GetFieldComponents_cached(G4ThreeVector{xx, yy, zz}); // output in Tesla, input in meters

    Bfield[0] = EMF.x() * tesla;
    Bfield[1] = EMF.y() * tesla;
    Bfield[2] = EMF.z() * tesla;

} // EarthMagField_WMM::GetFieldValue

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector EarthMagField_IGRF::GetFieldComponents(const G4ThreeVector &position_ECEF) const {
    // INPUT in meters

    // conversion ECEF to geodetic

    // input x y z in meters, output altitude in km, output lat lon in degrees
//    geod_conv::GeodeticConverter::ecef2Geodetic(position_ECEF.getX(), position_ECEF.getY(), position_ECEF.getZ(),
//                                                lat, lon, alt);

/// METHOD 1 and 2 should give same values. METHOD 2 should be a bit faster.

///// METHOD 1 : GEODETIC LATITUDE AND LONGITUDE
    earth->Reverse(position_ECEF.getX(), position_ECEF.getY(), position_ECEF.getZ(),
                   lat, lon, alt);
    elong = lon; // degrees
    colat = (90.0 - lat); // degrees
    alt_km = alt / 1000.0;
    igrf12syn_(&isv, &date, &itype, &alt_km, &colat, &elong, &Bx, &By, &Bz, &f,
               &ier); // gives NED (North East Down) components, in nT
    mag_field_ecef = myUtils::ned2ecef(Bx, By, Bz, lat, lon);

    /// METHOD 2 : GEOCENTRIC LATITUDE AND LONGITUDE
//    cart2sphr(position_ECEF.getX(), position_ECEF.getY(), position_ECEF.getZ(), &RR, &lat, &lon);
//    elong = lon; // degrees
//    colat = (90.0 - lat); // degrees
//    RR_km = RR / 1000.0;
//    itype = 2;
//    igrf12syn_(&isv, &date, &itype, &RR_km, &colat, &elong, &Bx, &By, &Bz, &f,
//               &ier); // gives NED (North East Down) components, in nT
//    mag_field_ecef = myUtils::ned2ecef(Bx, By, Bz, lat, lon);


    G4ThreeVector field_out{mag_field_ecef.u * 1e-9,
                            mag_field_ecef.v * 1e-9,
                            mag_field_ecef.w * 1e-9};


#ifndef NDEBUG // debug mode only
    if (std::isnan(field_out.x()) || std::isinf(field_out.x())
        || std::isnan(field_out.y()) || std::isinf(field_out.y())
        || std::isnan(field_out.z()) || std::isinf(field_out.z())) {
        std::abort();
    }
#endif // end debug mode only

    return field_out;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector EarthMagField_IGRF::GetFieldComponents_cached(const G4ThreeVector &position_ECEF) const {

    double distance_squared = (position_ECEF - position_cached_meters).mag2();

    if (distance_squared < cache_distance_squared) {
        return cached_FieldComponents;
    } else {
        position_cached_meters = position_ECEF;
        cached_FieldComponents = GetFieldComponents(position_ECEF);
        return cached_FieldComponents;
    }

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EarthMagField_IGRF::run_sanity_check() {
    /// Check with pre-calculated values of the magnetic field

    /// read the file
    std::string filename = "./mag_data/sanity_check/values_to_check.txt";

    std::ifstream file(filename);

    double XX, YY, ZZ, U_ecef, V_ecef, W_ecef;
    std::vector<double> Xs, Ys, Zs, U_ecefs, V_ecefs, W_ecefs;

    int nb = 0;

    if (file) {
        while (!file.eof()) {
            file >> XX >> YY >> ZZ >> U_ecef >> V_ecef >> W_ecef;

            Xs.push_back(XX);
            Ys.push_back(YY);
            Zs.push_back(ZZ);
            U_ecefs.push_back(U_ecef);
            V_ecefs.push_back(V_ecef);
            W_ecefs.push_back(W_ecef);

            nb++;
        }
    } else {
        G4cout << "ERROR : impossible to read the input file." << G4endl;
    }

    double max_diff = 0;
    double min_value = 100;
    double max_thres = 16; // maximum acceptable discrepency between reference data and the model used here (in nT)
    // using relative values does not work because they can be up to 80 % for low field values.

    double wt0 = myUtils::get_wall_time();
    for (int kk = 0; kk < 400; ++kk) { // number of loops can be increased for benchmarking

        for (int ii = 0; ii < nb; ++ii) {

            G4ThreeVector field_out = GetFieldComponents(G4ThreeVector{Xs[ii], Ys[ii], Zs[ii]});

            double diff1 = std::abs(U_ecefs[ii] - field_out.x()) * 1e9; // to nT to make it easier
            double diff2 = std::abs(V_ecefs[ii] - field_out.y()) * 1e9;
            double diff3 = std::abs(W_ecefs[ii] - field_out.z()) * 1e9;

            if (diff1 > max_diff) max_diff = diff1;
            if (diff2 > max_diff) max_diff = diff2;
            if (diff3 > max_diff) max_diff = diff3;
            if (std::abs(U_ecefs[ii]) < min_value) min_value = std::abs(U_ecefs[ii]) * 1e9;
            if (std::abs(V_ecefs[ii]) < min_value) min_value = std::abs(V_ecefs[ii]) * 1e9;
            if (std::abs(W_ecefs[ii]) < min_value) min_value = std::abs(W_ecefs[ii]) * 1e9;
//            G4cout << ii << "; maximum difference found: " << max_diff << " nT; minimum value found:  " << min_value << " nT" << G4endl;

            if (diff1 > max_thres || diff2 > max_thres || diff3 > max_thres) {
                G4cout << diff1 << " " << diff2 << " " << diff3 << G4endl;
                std::abort();
            }
        }
    }
    double wt1 = myUtils::get_wall_time();
    G4cout << "Magnetic Field sanity check OK " << G4endl;
    G4cout << " maximum difference found: " << max_diff << " nT; minimum value found:  " << min_value << " nT" << G4endl;
    G4cout << "Wall time taken: " << (wt1 - wt0) / 1.e6 << " s" << G4endl;
//    std::abort();
}

#endif


///