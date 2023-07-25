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
// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <Analysis.hh>
#include <EventAction.hh>
#include "G4Event.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction() {

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *evt) {

    G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(0.1 * mm);

    // initializations
    evtNb++;
//#ifndef NDEBUG // debug mode
//
//    if (evtNb % printModulo == 0) G4cout << "\n---> Begin Of Event: " << evtNb << G4endl;
//
//#endif // ifndef NDEBUG

    G4int thread_ID = G4Threading::G4GetThreadId();

    //if (thread_ID == 0) {
        if (evtNb == 1) G4cout << "\n---> Begin Of Event: " << evtNb << G4endl;
        if (evtNb % printModulo == 0) G4cout << "\n---> Begin Of Event: " << evtNb << G4endl;
    //}
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *) {
    //     xgreAnalysis* analysis = xgreAnalysis::getInstance();
    //     if (fTotalEnergyDeposit>5.*keV) analysis->analyseDeposit(fTotalEnergyDeposit);
    //     if (fTotalEnergyDepositPlas>5.*keV) analysis->analyseDepositPlas(fTotalEnergyDepositPlas);
    //     if (fOpticalCount>4) analysis->analyseOptical(fOpticalCount);
    //   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //   analysisManager->FillH1(1, fTotalEnergyDeposit/keV);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
