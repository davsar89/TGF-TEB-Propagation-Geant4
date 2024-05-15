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

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "Settings.hh"
#include <vector>

#include <GeographicLib/Geocentric.hpp>

struct DataPoint
{
    double x;
    double y;
};

class G4ParticleGun;

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    explicit PrimaryGeneratorAction();

    ~PrimaryGeneratorAction();

    // methods
    virtual void GeneratePrimaries(G4Event *);

private:
    // data members
    G4ParticleGun *fParticleGun; // pointer a to G4 service class
    G4int nofParticles = 1;

    GeographicLib::Geocentric *earth = nullptr;

    std::vector<DataPoint> Celestin_160MV_data;
    std::vector<DataPoint> Celestin_60MV_data;

    double max_60MV; // max y values for rejection sampling
    double max_160MV;

    double Sample_one_RREA_gammaray_energy(const double &MinEner, const double &MaxEner, const double &cut_ener);

    double Sample_one_BowersFormula_gammaray_energy(const double &MinEner, const double &MaxEner);

    double BrokenPL(double &p1, double &p2, double &ec, double &x);

    std::vector<DataPoint> loadCSV(const std::string &filename);
    void normalize(std::vector<DataPoint> &data);
    double logLogInterpolate(const std::vector<DataPoint> &data, const double x);
    double rejectionSampling(const std::vector<DataPoint> &data, const double maxY);
    double computeMaxY(const std::vector<DataPoint> &data);
};
