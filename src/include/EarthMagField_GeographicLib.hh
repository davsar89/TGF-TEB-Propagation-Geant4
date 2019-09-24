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

#ifndef _WIN32 // not usable on Windows

#include <vector>

#include "globals.hh"
#include "G4MagneticField.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <Settings.hh>
#include "G4TransportationManager.hh"
#include <GeographicLib/MagneticModel.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Math.hpp>
#include "myUtils.hh"

class EarthMagField_GeographicLib : public G4MagneticField {
public:

    EarthMagField_GeographicLib();

    ~EarthMagField_GeographicLib() override;

    void GetFieldValue(const double Point[3],
                       double *Bfield) const override;

    G4ThreeVector GetFieldComponents(const G4ThreeVector &position_ECEF) const;

    G4ThreeVector GetFieldComponents_cached(const G4ThreeVector &position_ECEF) const;

private:
    void run_sanity_check();

    mutable myUtils::vector3D mag_field_ecef{0, 0, 0};

    mutable int isv = 0;
    mutable int itype = 1;
    mutable int ier = 0;
    mutable double date;
    mutable double alt = 0;
    mutable double Bx = 0, By = 0, Bz = 0, f = 0, lat = 0, lon = 0;
    mutable double Bx2 = 0, By2 = 0, Bz2 = 0;
    mutable double xx, yy, zz;
    mutable double coslon = 0, sinlon = 0, sinlat = 0, coslat = 0;

    mutable double elong = 0; // degrees
    mutable double colat = 0; // degrees
    mutable double alt_m = 0; // km

    mutable double Bfield_ecef_x = 0;
    mutable double Bfield_ecef_y = 0;
    mutable double Bfield_ecef_z = 0;

    mutable double rough_alt = 0;
    const double earth_radius_m = 6371000.137;

    const double cache_distance = 25.0; // in meters
    const double cached_distance2 = cache_distance * cache_distance;
    mutable G4ThreeVector center_pos_cached{0, 0, 0};
    mutable G4ThreeVector field_cached{0, 0, 0};

    GeographicLib::MagneticModel *mag = nullptr;
    mutable GeographicLib::MagneticModel::field_comp field_compo;

    mutable GeographicLib::Geocentric *earth = nullptr;

};


#endif