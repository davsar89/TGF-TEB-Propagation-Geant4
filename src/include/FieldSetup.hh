//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file field/field03/include/F03FieldSetup.hh
/// \brief Definition of the F03FieldSetup class
//
//
// $Id$
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#pragma once

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4DormandPrince745.hh"
#include "EarthMagField_WMM.hh"
#include "EarthMagField_GeographicLib.hh"

#if defined(__linux__) && defined(__GNUG__)

#include "EarthMagField_IGRF.hh"

#endif

#include "G4CachedMagneticField.hh"


class G4FieldManager;

class G4ChordFinder;

class G4Mag_UsualEqRhs;

class G4MagIntegratorStepper;

class F03FieldMessenger;

///  A class for setting up the Magnetic Field
///
///  It also creates the necessary classes to control accuracy of propagation.
///  In this example
///    - There is a global field for most of the setup;
///    - A local field overides it for some volume(s) and it assumed to be
///      uniform.

class FieldSetup {
public:
    explicit FieldSetup();           //  A zero field
    virtual ~FieldSetup();

    void SetStepper();

    void SetMinStep(G4double sss) { fMinStep = sss; }

    void UpdateField();

    G4FieldManager *GetLocalFieldManager() { return fLocalFieldManager; }

private:

protected:

    // Find the global Field Manager

    G4FieldManager *GetGlobalFieldManager();

    G4FieldManager *fFieldManager;
    G4FieldManager *fLocalFieldManager;
    G4ChordFinder *fChordFinder;
    G4ChordFinder *fLocalChordFinder;
    G4Mag_UsualEqRhs *fEquation;
    G4Mag_UsualEqRhs *fLocalEquation;
    G4MagneticField *fMagneticField;
    G4MagneticField *fLocalMagneticField;

    G4MagIntegratorStepper *fStepper;
    G4MagIntegratorStepper *fLocalStepper;

    G4double fMinStep;

};
