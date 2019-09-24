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

#pragma once

#include "globals.hh"
#include "G4VUserRegionInformation.hh"

class RegionInformation : public G4VUserRegionInformation {
public:

    RegionInformation();


    ~RegionInformation() override;

    void Print() const override;

    inline void Set_World(G4bool v = true) {
        fis_world = v;
    }

    inline G4bool is_World() const {
        return fis_world;
    }

    inline void set_considered_atmosphere(G4bool v = true) {
        fis_considered_atmosphere = v;
    }

    inline G4bool is_considered_atmosphere() const {
        return fis_considered_atmosphere;
    }

private:

    G4bool fis_world;
    G4bool fis_considered_atmosphere;
};
