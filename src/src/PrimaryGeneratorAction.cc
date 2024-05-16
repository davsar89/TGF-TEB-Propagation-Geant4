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
#include <PrimaryGeneratorAction.hh>
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

using namespace std;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr)
{

    fParticleGun = new G4ParticleGun(nofParticles);

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);

    earth = new GeographicLib::Geocentric(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());

    Celestin_60MV_data = loadCSV("Celestin_60MV.csv");
    Celestin_160MV_data = loadCSV("Celestin_160MV.csv");

    normalize(Celestin_60MV_data);
    normalize(Celestin_160MV_data);

    max_60MV = computeMaxY(Celestin_60MV_data);
    max_160MV = computeMaxY(Celestin_160MV_data);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    delete earth;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    // this function is called at the begining of event
    // for 1/E spectrum: G4double energy=MinEner*(pow(((MaxEner)/(MinEner)),G4UniformRand()))

    ////////////// ENERGY /////////////////
    double MaxEner = 60. * CLHEP::MeV;  // Max energy
    double MinEner = 300. * CLHEP::keV; // Min energy
    double energy = 0;

    if (Settings::SPECTRUM_MODEL == 0)
    {
        double cut_ener = 7.3 * MeV; // exponential cut-off factor
        energy = Sample_one_RREA_gammaray_energy(MinEner, MaxEner, cut_ener);
    }
    else if (Settings::SPECTRUM_MODEL == 1)
    {
        energy = Sample_one_BowersFormula_gammaray_energy(MinEner, MaxEner);
    }
    else if (Settings::SPECTRUM_MODEL == 2)
    {
        energy = rejectionSampling(Celestin_60MV_data, max_60MV) * CLHEP::electronvolt;
    }
    else if (Settings::SPECTRUM_MODEL == 3)
    {
        energy = rejectionSampling(Celestin_160MV_data, max_160MV) * CLHEP::electronvolt;
    }
    else
    {
        G4cout << "Invalid SPECTRUM_MODEL value (should be 0,1,2 or 3)" << G4endl;
        std::abort();
    }

    ////////////// POSITION / ANGLE /////////////////
    // ! : sampling theta uniformly between 0 and SOURCE_OPENING_ANGLE*degree does not sample uniformly over the area
    //  G4double theta = G4UniformRand()*shared_var::SOURCE_OPENING_ANGLE*degree;
    //     G4double theta = CLHEP::RandGauss::shoot(0., shared_var::SOURCE_OPENING_ANGLE * degree);
    //    G4double theta = Settings::OpeningAngle() * degree * G4UniformRand();
    //    photons.alpha=acos(1.-rand*(1-cos(isot/180.*!pi)))/!pi*180.
    //    G4double theta;
    G4double X_try = 0, Y_try = 0;
    G4double R_max;
    G4double R_try;
    G4double sigma_sample_R;

    if ((Settings::BEAMING_TYPE == "Uniform") || (Settings::BEAMING_TYPE == "uniform"))
    {
        R_max = std::tan(Settings::SOURCE_OPENING_ANGLE * degree);
        R_try = R_max + 10.; // just for initialization

        while (R_try > R_max)
        {
            X_try = 2. * (G4UniformRand() - 0.5) * R_max;
            Y_try = 2. * (G4UniformRand() - 0.5) * R_max;
            R_try = sqrt(X_try * X_try + Y_try * Y_try);
        }

        //         theta = std::acos(1. - G4UniformRand()*(1. - std::cos(Settings::OpeningAngle() * degree))); // uniform over spherical(i.e a part of
        // sphere) area
    }
    else if ((Settings::BEAMING_TYPE == "Gaussian") || (Settings::BEAMING_TYPE == "gaussian") ||
             (Settings::BEAMING_TYPE == "normal") || (Settings::BEAMING_TYPE == "Normal"))
    {
        R_max = 10000.; // -> maximum angle is atan(10000) = 89.9943 degrees
        sigma_sample_R = std::tan(Settings::SOURCE_OPENING_ANGLE * degree);
        R_try = R_max + 10.; // just for initialization

        while (R_try > R_max)
        {
            X_try = CLHEP::RandGauss::shoot(0., sigma_sample_R); // gaussian position sample
            Y_try = CLHEP::RandGauss::shoot(0., sigma_sample_R); // gaussian position sample
            R_try = sqrt(X_try * X_try + Y_try * Y_try);
            //  G4cout << X_try << G4endl;
        }
    }
    else
    {
        G4cout << "ERROR : Beaming type is not Gaussian or Uniform. Aborting." << G4endl;
        std::abort();
    }

    G4double lat = Settings::SOURCE_LAT;
    G4double lon = Settings::SOURCE_LONG;
    G4double alt = Settings::SOURCE_ALT * 1000.0; // km to m
    G4double ecef_x, ecef_y, ecef_z;

    earth->Forward(lat, lon, alt, ecef_x, ecef_y, ecef_z);

    //    geod_conv::GeodeticConverter::geodetic2ecef(lat, lon, alt, ecef_x, ecef_y, ecef_z);

    G4ThreeVector position;
    position.setX(ecef_x * m);
    position.setY(ecef_y * m);
    position.setZ(ecef_z * m);

    G4ThreeVector localVertical, localVertical_perp1, localVertical_perp2;

    localVertical = position / position.mag();
    // computing two vector perpendicular to the local vertical and perpendicular to each other
    G4double ux = localVertical[0];
    G4double uy = localVertical[1];
    G4double uz = localVertical[2];
    G4double norme = sqrt(pow(ux, 2) + pow(uy, 2));
    G4double vx = -uy / norme;
    G4double vy = ux / norme;
    G4double vz = 0.;
    G4double wx = -uz * vy;
    G4double wy = uz * vx;
    G4double wz = ux * vy - uy * vx;
    localVertical_perp1 = G4ThreeVector(vx, vy, vz);
    localVertical_perp2 = G4ThreeVector(wx, wy, wz);

    // adding tilt angle and recomputing the two perpendicular vectors
    if (Settings::TILT_ANGLE != 0.0)
    {
        G4ThreeVector tilt_shift =
            localVertical_perp1 * std::sin(Settings::TILT_ANGLE * degree); // we could also use localVertical_perp2
        localVertical = localVertical + tilt_shift;
        localVertical = position / position.mag();
        ux = localVertical[0];
        uy = localVertical[1];
        uz = localVertical[2];
        norme = sqrt(pow(ux, 2) + pow(uy, 2));
        vx = -uy / norme;
        vy = ux / norme;
        vz = 0.;
        wx = -uz * vy;
        wy = uz * vx;
        wz = ux * vy - uy * vx;
        localVertical_perp1 = G4ThreeVector(vx, vy, vz);
        localVertical_perp2 = G4ThreeVector(wx, wy, wz);
    }

    G4ThreeVector position_virtual_particle;
    position_virtual_particle = position + localVertical + localVertical_perp1 * X_try + localVertical_perp2 * Y_try;
    G4ThreeVector momentumDirection = (position_virtual_particle - position);
    momentumDirection = momentumDirection / momentumDirection.mag();
    //     theta = std::acos(momentumDirection.dot(localVertical))/degree;
    //     G4cout << theta << G4endl;
    ////////////// TIME /////////////////
    G4double time;

    if (Settings::SOURCE_SIGMA_TIME == 0.)
    {
        time = 0.;
    }
    else
    {
        time = CLHEP::RandGauss::shoot(0., Settings::SOURCE_SIGMA_TIME) * microsecond;
    }

    ////////////// assignments /////////////////
    fParticleGun->SetParticleTime(time);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->SetParticleMomentumDirection(momentumDirection);
    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::BrokenPL(double &p1, double &p2, double &ec, double &x)
{
    double broken_PL = 0;

    if (x < ec)
    {
        broken_PL = std::pow(x, p1);
    }
    else if (x >= ec)
    {
        broken_PL = std::pow(ec, p1 - p2) * std::pow(x, p2);
    }

    return broken_PL;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::Sample_one_RREA_gammaray_energy(const double &MinEner, const double &MaxEner, const double &cut_ener)
{
    // random samples the energy of one RREA gamma ray
    // "classical" fully developped RREA gamma spectrum is approximately = 1/E * exp (-E / 7.3MeV)
    // (rejection method)

    auto RREA_spec = [&cut_ener](const double x)
    {
        return 1.0 / x * exp(-x / cut_ener);
    };

    const double pMax = RREA_spec(MinEner) * 1.01;
    const double pMin = RREA_spec(MaxEner);
    double pOfeRand = 0.0;
    double pRand = 1.0;
    double eRand = 0.0;

    while (pOfeRand < pRand)
    { // rejection
        pRand = pMin + (pMax - pMin) * G4UniformRand();
        eRand = MinEner + (MaxEner - MinEner) * G4UniformRand();
        pOfeRand = RREA_spec(eRand);
    }

    return eRand;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::Sample_one_BowersFormula_gammaray_energy(const double &MinEner, const double &MaxEner)
{
    // Sample Bowers et al. 2018 formula for reverse positron beam TGF (much harder than classical TGF)
    // (rejection method)

    const double E0 = 100 * CLHEP::MeV;
    const double power = 0.85;

    auto Bowers_spec = [&E0, &power](const double E)
    {
        return 1.0 / E * exp(-E0 / (std::pow(E0, power) - std::pow(E, power)));
    };

    const double pMax = Bowers_spec(MinEner) * 1.01;
    const double pMin = Bowers_spec(MaxEner);
    double pOfeRand = 0.0;
    double pRand = 1.0;
    double eRand = 0.0;

    while (pOfeRand < pRand)
    { // rejection
        pRand = pMin + (pMax - pMin) * G4UniformRand();
        eRand = MinEner + (MaxEner - MinEner) * G4UniformRand();
        pOfeRand = Bowers_spec(eRand);
    }

    return eRand;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<DataPoint> PrimaryGeneratorAction::loadCSV(const std::string &filename)
{
    std::vector<DataPoint> data;
    data.reserve(60);

    std::ifstream file(filename);
    std::string line;

    // Skip header line
    std::getline(file, line);

    double value1, value2;
    char comma; // To read and discard the comma

    // Read data
    while (std::getline(file, line))
    {
        while (file >> value1 >> comma >> value2)
        {
            if (comma != ',')
            {
                std::cerr << "Error: Malformed CSV file." << std::endl;
                std::abort();
            }
            const DataPoint point = {value1, value2 * 1.0e10};
            data.push_back(point);
        }
    }

    return data;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::normalize(std::vector<DataPoint> &data)
{
    double sum = 0.0;
    for (const auto &point : data)
    {
        sum += point.y;
    }
    for (auto &point : data)
    {
        point.y /= sum;
    }
}

double PrimaryGeneratorAction::logLogInterpolate(const std::vector<DataPoint> &data, const double x)
{
    // Use lower_bound to find the first element in data where p.x >= x
    auto it = std::lower_bound(data.begin(), data.end(), x, [](const DataPoint &p, double val)
                               { return p.x < val; });

    // Check if the iterator is at the beginning or end of the vector
    if (it == data.end() || it == data.begin())
    {
        throw std::out_of_range("Interpolation x value out of range");
    }

    // Get the previous iterator (it - 1) and current iterator (it)
    auto it1 = it - 1;
    auto it2 = it;

    // Calculate the logarithms of the x and y values for interpolation
    double x1 = std::log(it1->x);
    double y1 = std::log(it1->y);
    double x2 = std::log(it2->x);
    double y2 = std::log(it2->y);

    // Perform the linear interpolation in log-log space
    double y = y1 + (y2 - y1) * (std::log(x) - x1) / (x2 - x1);

    // Return the exponentiated result to convert back from log space
    return std::exp(y);
}

double PrimaryGeneratorAction::computeMaxY(const std::vector<DataPoint> &data)
{
    double maxY = 0.0;
    for (const auto &point : data)
    {
        if (point.y > maxY)
        {
            maxY = point.y;
        }
    }
    return maxY;
}

double PrimaryGeneratorAction::rejectionSampling(const std::vector<DataPoint> &data, const double maxY)
{
    while (true)
    {
        // Generate a random x value within the range of the data set
        double x = data.front().x + (data.back().x - data.front().x) * G4UniformRand();

        // Generate a random y value between 0 and maxY
        double y = maxY * G4UniformRand();

        // Interpolate the y value for the generated x using log-log interpolation
        double interpolatedY = logLogInterpolate(data, x);

        // Accept the sample if the generated y is less than the interpolated y
        if (y < interpolatedY)
        {
            return x;
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
