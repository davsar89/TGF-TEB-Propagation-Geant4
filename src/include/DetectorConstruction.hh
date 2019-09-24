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

#define GEANT_VERSION 4103

#include <Settings.hh>

#include "EarthMagField_WMM.hh"
#include "EarthMagField_IGRF.hh"

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Sphere.hh"
#include <vector>
#include "G4ClassicalRK4.hh"
#include "G4VPVParameterisation.hh"
#include "G4Element.hh"
#include "G4HelixImplicitEuler.hh"

extern "C" {
#include <nrlmsise-00.h>
}

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>
#include <fstream>
#include <string>
#include "G4UserLimits.hh"
#include "G4SDManager.hh"

#include "RegionInformation.hh"
#include "G4RegionStore.hh"

#include <PrimaryGeneratorAction.hh>
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4String.hh"
#include "G4KM_DummyField.hh"
#include "G4UniformMagField.hh"
#include "G4CachedMagneticField.hh"
#include "G4PVParameterised.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "EarthMagField_WMM.hh"
#include <string>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"

#if GEANT_VERSION == 4104
#include "G4IntegrationDriver.hh"
#include "G4VIntegrationDriver.hh"
#endif

#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4DormandPrince745.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include <G4ExtendedMaterial.hh>
#include "FieldSetup.hh"
#include "G4Cache.hh"
#include "G4AutoDelete.hh"

class G4UserLimits;

using namespace std;

class TGFDetectorConstruction : public G4VUserDetectorConstruction {
public:

    TGFDetectorConstruction();

    ~TGFDetectorConstruction() override;

    G4VPhysicalVolume *Construct() override;

    void ConstructSDandField() override;

private:

    GeographicLib::Geocentric *earth = nullptr;

    G4Region *RECORD_REGION = new G4Region("RECORD_REGION");

    G4Cache<FieldSetup *> fEmFieldSetup;

    ///

    double base_radius;

    ///

    double find_radius_for_record_altitude(const double &center_theta, const double &center_phi, const double &rec_alt);

    bool CHECK_ASCENDING(const std::vector<double> &alt_list);

    std::vector<G4Material *> GENERATE_AIR_MATERIALS(double min_density, double max_density, int number);

    void calculate_radii_list();

    int find_atmosphere_part_material(double lat, double lon, double alt);

    std::vector<G4Material *> Airs;

    //    void ConstructAtmosMats2();
    //    void ConstructAtmosMats3();
    //    void ReadInputAtmosFile();

    //    double interp1(vector < double >,
    //                   vector < double >,
    //                   double);
    //    int    findNearestNeighbourIndex(double,
    //                                     vector < double >);

    G4LogicalVolume *logicalWorld;
    G4VPhysicalVolume *physicalWorld;
    G4Material *vac = nullptr;

    std::vector<G4Sphere *> atmosLayers_S;
    std::vector<G4LogicalVolume *> atmosLayers_LV;
    std::vector<G4VPhysicalVolume *> atmosLayers_PV;

    std::vector<G4Sphere *> det_layers_S;
    std::vector<G4LogicalVolume *> det_layers_LV;
    std::vector<G4VPhysicalVolume *> det_layers_PV;

    std::vector<double> radius_list; // (geodetic)altitudes intervals of the layers
    G4int nb_altitudes = 256;
    double radius_min = 0.1 * km; //
    double min_alt_to_build = 1.0 * km; // just inititialization value
    const int nb_theta = 45;
    const int nb_phi = 6;

    // because generating more than thousands of air materials uses a lot of memory
    const int number_of_AIR_materials = 256;

    //
    double world_max_altitude = 15000. * km;

    double saved_radius_test = 0; // initialization

    double fMinStep = 0.010 * mm;
    double minEps = 1.0e-7; //   Minimum & value for smallest steps
    double maxEps = 1.0e-6; //   Maximum & value for largest steps

    G4bool not_contains(double value, const std::vector<double> &vec);

    std::ofstream asciiFile;

    bool hasDuplicates(const std::vector<double> &arr);

    // Vaccum
    G4NistManager *man = G4NistManager::Instance();
    G4Material *vaccum = nullptr;

    G4Element *elN = new G4Element("Nitrogen", "N", 7., 14.01 * g / mole);
    G4Element *elO = new G4Element("Oxygen", "O", 8., 16.00 * g / mole);

    const double sea_level_density = 1.304E-03 * g / cm3;
    const double km_150_density = 2.190E-12 * g / cm3;
    std::vector<double> density_grid;

    template<typename T>
    std::vector<double> linspace(T start_in, T end_in, int num_in) {

        std::vector<double> linspaced;

        double start = static_cast<double>(start_in);
        double end = static_cast<double>(end_in);
        double num = static_cast<double>(num_in);

        if (num == 0) { return linspaced; }
        if (num == 1) {
            linspaced.push_back(start);
            return linspaced;
        }

        double delta = (end - start) / (num - 1);

        for (int i = 0; i < num - 1; ++i) {
            linspaced.push_back(start + delta * i);
        }
        linspaced.push_back(end); // I want to ensure that start and end
        // are exactly the same as the input
        return linspaced;
    }

    template<typename T>
    class Logspace {
    private:
        T curValue, base;

    public:
        Logspace(T first, T base) : curValue(first), base(base) {}

        T operator()() {
            T retval = curValue;
            curValue *= base;
            return retval;
        }
    };

    std::vector<double> pyLogspace(double start, double stop, int num = 50, double base = 10) {
        double realStart = pow(base, start);
        double realBase = pow(base, (stop - start) / num);

        std::vector<double> retval;
        retval.reserve(num);
        std::generate_n(std::back_inserter(retval), num + 1, Logspace<double>(realStart, realBase));
        return retval;
    }

};

namespace datetools {
    namespace details {
        constexpr unsigned int days_to_month[2][12] =
                {
                        // non-leap year
                        {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
                        // leap year
                        {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335},
                };
    }

    constexpr bool is_leap(int const year) noexcept {
        return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    }

    constexpr unsigned int day_of_year(int const year, unsigned int const month, unsigned int const day) {
        return details::days_to_month[is_leap(year)][month - 1] + day;
    }
}