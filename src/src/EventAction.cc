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

EventAction::EventAction() : G4UserEventAction()
{
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *evt)
{

    G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(0.1 * mm);

    // initializations
    evtNb++;
    // #ifndef NDEBUG // debug mode
    //
    //     if (evtNb % printModulo == 0) G4cout << "\n---> Begin Of Event: " << evtNb << G4endl;
    //
    // #endif // ifndef NDEBUG

    G4int thread_ID = G4Threading::G4GetThreadId();

    // if (thread_ID == 0) {
    if (evtNb == 1)
    {
        std::string dt_str = myUtils::getCurrentDateTimeStr();
        G4cout << "-- " << dt_str << ":---> Begin Of Event: " << evtNb << G4endl;
    }
    if (evtNb % printModulo == 0)
    {
        const std::string dt_str = myUtils::getCurrentDateTimeStr();
        const std::string remaining_time_str = estimateRemainingTime(evtNb, Settings::NB_EVENT_TOTAL, start_time);
        G4cout << "-- " << dt_str << ":---> Begin Of Event: " << evtNb << G4endl;
        G4cout << "                   "
               << "Remaining duration: " << remaining_time_str << " (for " << Settings::NB_EVENT_TOTAL << " events)" << G4endl;
    }
    //}
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *)
{
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Function to estimate remaining time based on completed iterations
std::string EventAction::estimateRemainingTime(const int evtNb, int const NB_EVENT_TOTAL, const std::chrono::steady_clock::time_point &start_time)
{
    if (evtNb == 0)
    {
        return "99:99:99"; // Return very long time if no events have been processed
    }

    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    std::chrono::steady_clock::duration time_elapsed = now - start_time;
    std::chrono::steady_clock::duration average_time_per_event = time_elapsed / evtNb;
    int events_left = NB_EVENT_TOTAL - evtNb;
    std::chrono::steady_clock::duration remaining_time = average_time_per_event * events_left;

    // Convert remaining_time to hours, minutes, and seconds
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(remaining_time).count();
    int hours = seconds / 3600;
    int minutes = (seconds % 3600) / 60;
    seconds = seconds % 60;

    // Format time into a string "HH:MM:SS"
    std::ostringstream stream;
    stream << std::setw(2) << std::setfill('0') << hours << ":"
           << std::setw(2) << std::setfill('0') << minutes << ":"
           << std::setw(2) << std::setfill('0') << seconds;

    return stream.str();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......