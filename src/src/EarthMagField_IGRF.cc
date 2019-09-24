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

#ifndef _WIN32 // not usable on Windows

#include "PrimaryGeneratorAction.hh"
#include "EarthMagField_IGRF.hh"

using namespace std;

// interface with igrf12 Fortran code from : http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
extern "C" {
	void igrf12syn_(int *isv, double *date, int *itype, double *alt, double *colat, double *elong, double *x, double *y,
		double *z, double *f, int *ier);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField_IGRF::EarthMagField_IGRF() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField_IGRF::~EarthMagField_IGRF() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EarthMagField_IGRF::GetFieldValue(const double Point[3], double *Bfield) const {
	//  geodetic_converter::GeodeticConverter g_geodetic_converter;

	xx = Point[0] / m;
	yy = Point[1] / m;
	zz = Point[2] / m;

	// conversion ECEF to geodetic

	// input x y z in meters, output altitude in km, output lat lon in degrees
	geod_conv::GeodeticConverter::ecef2Geodetic(xx, yy, zz, lat, lon, alt);

	if (alt < 45000.0) {
		Bfield[0] = 0;
		Bfield[1] = 0;
		Bfield[2] = 0;
		return;
	}

	elong = lon; // degrees
	colat = (180.0 / 2.0 - lat); // degrees
	alt_km = alt / 1000.0;

	//
	// // Get altitude from latitude
	// G4double altt = p*coslat + (z + e2*Nphi*sinlat)*sinlat - Nphi;

	//     itype = 1 if geodetic (spheroid)
	//     itype = 2 if geocentric (sphere)
	//     alt   = height in km above sea level if itype = 1
	//           = distance from centre of Earth in km if itype = 2 (>3485 km)                                              // degrees

	// IGRF magnetic field
	igrf12syn_(&isv, &date, &itype, &alt_km, &colat, &elong, &Bx, &By, &Bz, &f,
		&ier); // gives NED (North East Down) components, in nT

//     G4cout << Bx << " " << By << " " << Bz << G4endl;

// NED (North East Down) to ECEF convertion

// NOT USING EIGEN TO GET MORE PERFORMANCE

// Local north direction      //Local east direction     // Local vertical (down)
//        MAT_enu_ecef << -coslon*sinlat,            -sinlon,                  -coslon*coslat,
//                        -sinlon*sinlat,             coslon,                  -sinlon*coslat,
//                  coslat,                      0.   ,                     -sinlat;

//        Bvec=MAT_enu_ecef*Bvec;

	lat = lat * degree;   // degree to radians
	lon = lon * degree;   // degree to radians

	sincos(lon, &sinlon, &coslon); // important for performance
	sincos(lat, &sinlat, &coslat);

	Bfield_ecef_x = -coslon * sinlat * Bx - sinlon * By - coslon * coslat * Bz;
	Bfield_ecef_y = -sinlon * sinlat * Bx + coslon * By - sinlon * coslat * Bz;
	Bfield_ecef_z = coslat * Bx - sinlat * Bz;

	Bfield[0] = Bfield_ecef_x * nano_tesla_to_G4;
	Bfield[1] = Bfield_ecef_y * nano_tesla_to_G4;
	Bfield[2] = Bfield_ecef_z * nano_tesla_to_G4;

} // EarthMagField_WMM::GetFieldValue


#endif