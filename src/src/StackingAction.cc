//
// ********************************************************************

//

#include "StackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4StackManager.hh"

extern const double time_step_length;
extern double current_time_step;

StackingAction::StackingAction()
{
}

StackingAction::~StackingAction()
{  }


G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track *aTrack)
{
    //    if (aTrack->GetTrackStatus() == fSuspend)
    //        {
    //            return fWaiting;
    //        }

    return fUrgent;

}

void StackingAction::NewStage()
// called when particles in 'waiting' stack are upgraded to particles in 'urgent' stack

{

}

void StackingAction::PrepareNewEvent()
{

}



