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

#include <vector>

#include "globals.hh"
#include "G4MagneticField.hh"
#include <vector>
#include <geodetic_converter.hh>
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include <string>

extern "C"
{
#include "GeomagnetismHeader.h"
}

#include "Settings.hh"

class EarthMagField_WMM : public G4MagneticField {
public:
    EarthMagField_WMM();

    ~EarthMagField_WMM() override;

    void GetFieldValue(const double Point[3],
                       double *Bfield) const override;

private:

    Settings *settings = Settings::getInstance();

    mutable double altitude_MagField_off_in_meters = 35000.0;

    mutable int isv = 0;
    mutable int itype = 1;
    mutable int ier = 0;
    mutable double date = double(settings->year) + double(settings->month) / 12.0 + double(settings->day) / 365.25;
    mutable double alt = 0;
    mutable double Bx = 0, By = 0, Bz = 0, f = 0, lat = 0, lon = 0;
    mutable double xx, yy, zz;
    mutable double coslon = 0, sinlon = 0, sinlat = 0, coslat = 0;

    mutable double elong = 0;  // degrees
    mutable double colat = 0;  // degrees
    mutable double alt_km = 0; // km

    mutable double Bfield_ecef_x = 0;
    mutable double Bfield_ecef_y = 0;
    mutable double Bfield_ecef_z = 0;

    const G4double nano_tesla_to_G4 = tesla * 1.e-9;

    mutable MAGtype_MagneticModel *MagneticModels[1], *TimedMagneticModel;
    mutable MAGtype_Ellipsoid Ellip;
    mutable MAGtype_CoordSpherical CoordSpherical;
    mutable MAGtype_CoordGeodetic CoordGeodetic;
    mutable MAGtype_Date UserDate;
    mutable MAGtype_GeoMagneticElements GeoMagneticElements, Errors;
    mutable MAGtype_Geoid Geoid;

    mutable char ans[20], b;

    mutable char VersionDate[12];
    mutable int NumTerms, Flag = 1, nMax = 0;
    mutable int epochs = 1;

};
