#include <utility>

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
/// \file field/field03/src/F03FieldSetup.cc
/// \brief Implementation of the F03FieldSetup class
//
//
// $Id$
//
//
//   Field Setup class implementation.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "FieldSetup.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "myG4FieldManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FieldSetup::FieldSetup()
        : fFieldManager(0),
          fLocalFieldManager(0),
          fChordFinder(0),
          fLocalChordFinder(0),
          fEquation(0),
          fLocalEquation(0),
          fMagneticField(0),
          fLocalMagneticField(0),
          fStepper(0),
          fLocalStepper(0) {

    if (Settings::MAGNETIC_FIELD_MODEL == "IGRF") {
#if defined(__linux__) && defined(__GNUG__) // if linux and GCC
        fMagneticField = new EarthMagField_IGRF();
#endif
    } else if (Settings::MAGNETIC_FIELD_MODEL == "WMM") {
        fMagneticField = new EarthMagField_WMM();
    } else if (Settings::MAGNETIC_FIELD_MODEL == "GEOLIB") {
        fMagneticField = new EarthMagField_GeographicLib();
    } else {
        G4cout << "Magnetic model must be IGRF or WMM or GEOLIB. Aborting." << G4endl;
        std::abort();
    }

    fEquation = new G4Mag_UsualEqRhs(fMagneticField);

    fMinStep = 0.1 * mm; // minimal step of 1 mm is default

    G4TransportationManager::GetTransportationManager()
            ->SetFieldManager(new myG4FieldManager());

    fFieldManager = GetGlobalFieldManager();

    UpdateField();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FieldSetup::~FieldSetup() {
    delete fMagneticField;
    delete fChordFinder;
    delete fStepper;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FieldSetup::UpdateField() {
    // It must be possible to call 'again' - e.g. to choose an alternative stepper
    //   has been chosen, or in case other changes have been made.

    // 1. First clean up previous state.
    delete fChordFinder;
    fChordFinder = nullptr;
    delete fLocalChordFinder;
    fLocalChordFinder = nullptr;


    // 2. Create the steppers ( Note: this also deletes the previous ones. )
    SetStepper();

    double minEps = 1.0e-6; //   Minimum & value for smallest steps
    double maxEps = 1.0e-5; //   Maximum & value for largest steps

    // 3. Create the chord finder(s)
    fChordFinder = new G4ChordFinder(fMagneticField, fMinStep, fStepper);

    fFieldManager->SetChordFinder(fChordFinder);

    fFieldManager->SetMinimumEpsilonStep(minEps);
    fFieldManager->SetMaximumEpsilonStep(maxEps);
    fFieldManager->SetDeltaOneStep(100.0 * cm);
    fFieldManager->SetDeltaIntersection(10.0 * cm);

    // 4. Ensure that the field is updated (in Field manager & equation)
    fFieldManager->SetDetectorField(fMagneticField);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FieldSetup::SetStepper() {
    delete fStepper;
    fStepper = nullptr;

    fStepper = new G4DormandPrince745(fEquation);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FieldManager *FieldSetup::GetGlobalFieldManager() {
    return G4TransportationManager::GetTransportationManager()
            ->GetFieldManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
