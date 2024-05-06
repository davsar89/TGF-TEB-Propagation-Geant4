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
#include "PhysicsList.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TGF_PhysicsList::TGF_PhysicsList() : G4VUserPhysicsList() {

//    emPhysicsList = new G4EmStandardPhysics_option4_dr(0);
    emPhysicsList = new G4EmStandardPhysics_option1_dr(0);
    emPhysicsList->SetVerboseLevel(0);
//    this->DumpCutValuesTable();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TGF_PhysicsList::~TGF_PhysicsList() {
    delete emPhysicsList;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::ConstructParticle() {
    emPhysicsList->ConstructParticle();
    //   G4GenericIon::GenericIon();
    //   G4NeutrinoE::NeutrinoEDefinition();
    //   G4AntiNeutrinoE::AntiNeutrinoEDefinition();
    //   G4Alpha::AlphaDefinition();
    //
    //   G4Geantino::GeantinoDefinition();
    //   G4ChargedGeantino::ChargedGeantinoDefinition();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::ConstructProcess() {
    AddTransportation();
    emPhysicsList->ConstructProcess();
    //   G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
    //   radioactiveDecay->SetICM(true);                //Internal Conversion
    //   radioactiveDecay->SetARM(true);               //Atomic Rearangement
    // radioactiveDecay->SetFBeta(true);
    // radioactiveDecay->SetBRBias(true);
    //   G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    //
    //   ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

    /// STEP LIMITATION for photons in record area

    if (Settings::USE_STEP_MAX_for_record) {
        auto myParticleIterator = GetParticleIterator();
        myParticleIterator->reset();

        auto *stepLimiter = new G4StepLimiter();

        while ((*myParticleIterator)()) {
            G4ParticleDefinition *particle = myParticleIterator->value();
            G4ProcessManager *pmanager = particle->GetProcessManager();

            if (!particle->IsShortLived()) {
                if ((particle->GetPDGEncoding() == 22) || (particle->GetPDGEncoding() == 11) || (particle->GetPDGEncoding() == -11)) {
                    // All particles should have a step limiter
                    // to make sure that the steps do not get too long.
                    pmanager->AddDiscreteProcess(stepLimiter);
                }
            }
        }
    }

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::SetCuts() {
    defaultCutValue = 10. * cm;
    //
    cutForGamma = defaultCutValue;
    cutForElectron = defaultCutValue;
    cutForPositron = defaultCutValue;
    //
    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");

//    G4double lowlimit = Settings::MIN_ENERGY_OUTPUT;
    G4ProductionCutsTable *aPCTable = G4ProductionCutsTable::GetProductionCutsTable();
    aPCTable->SetEnergyRange(10.0 * keV, 100 * CLHEP::GeV);

}