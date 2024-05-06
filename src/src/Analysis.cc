////////////////////////////////////////////////////////////////////////////////

// /* GEANT4 code for propagation of gamma-rays, electron and positrons in
// Earth's atmosphere */
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

#include <Analysis.hh>

// constructor
Analysis::Analysis()
{

    ///

    const double ALT_MAX_RECORDED = Settings::record_altitude;

    ///

    const G4String output_filename_second_part =
        std::to_string(Settings::RAND_SEED) + "_" +
        std::to_string(int(ALT_MAX_RECORDED)) + "_" +
        std::to_string(int(Settings::SOURCE_ALT)) + "_" +
        std::to_string(int(Settings::SOURCE_OPENING_ANGLE)) + "_" +
        Settings::BEAMING_TYPE + "_" +
        std::to_string(int(Settings::SOURCE_SIGMA_TIME)) + ".out";
    //

    //

    if ((Settings::BEAMING_TYPE == "Uniform") ||
        (Settings::BEAMING_TYPE == "uniform"))
    {
        number_beaming = 0;
    }
    else if ((Settings::BEAMING_TYPE == "Gaussian") ||
             (Settings::BEAMING_TYPE == "gaussian") ||
             (Settings::BEAMING_TYPE == "normal") ||
             (Settings::BEAMING_TYPE == "Normal"))
    {
        number_beaming = 1;
    }

    ///

    const std::string base_outputfolder = "./output_ascii/";

    createDirectory(base_outputfolder);

    std::string output_folder_name = base_outputfolder + std::to_string(int(ALT_MAX_RECORDED)) + "_" +
                              std::to_string(int(Settings::SOURCE_ALT)) + "_" +
                              std::to_string(int(Settings::SOURCE_OPENING_ANGLE)) + "_" +
                              Settings::BEAMING_TYPE + "_" +
                              std::to_string(int(Settings::SOURCE_SIGMA_TIME))+"/";

    createDirectory(output_folder_name);

    asciiFileName = output_folder_name + "/detParticles_" + output_filename_second_part;

    std::ofstream asciiFile(asciiFileName,
                              std::ios::trunc); // to clean the output file

    asciiFile.close();

    output_lines.clear();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Analysis::get_NB_RECORDED() const { return NB_RECORDED; }

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Analysis::save_in_output_buffer(
    const G4int PDG_NB, const G4double &time, const G4double &energy,
    const G4double &dist_rad, const G4int ID, const G4double &ecef_x,
    const G4double &ecef_y, const G4double &ecef_z, const G4double &mom_x,
    const G4double &mom_y, const G4double &mom_z, const G4double &lat, const G4double &lon, const G4double &alt, const int event_nb)
{

    //
    double alt2 = alt / 1000.0; // m to km

    bool record_or_not = false;

    if (Settings::RECORD_PHOT_ONLY && (PDG_NB == 11 || PDG_NB == -11))
    {
        return;
    }

    if (Settings::RECORD_ELEC_POSI_ONLY && PDG_NB == 22)
    {
        return;
    }

    if (Settings::RECORD_ONLY_IN_WINDOW)
    {
        const bool is_inside_record_window = (lat > Settings::RECORD_WIN.min_lat && lat < Settings::RECORD_WIN.max_lat) && (lon > Settings::RECORD_WIN.min_lon && lon < Settings::RECORD_WIN.max_lon);

        if (is_inside_record_window)
            record_or_not = true;
    }
    else
    {
        record_or_not = true;
    }

    if (record_or_not)
    {

        // ASCII OUTPUT
        if (Settings::OUTPUT_TO_ASCII_FILE)
        {
            std::stringstream buffer;
            buffer << std::scientific
                   << std::setprecision(5); // scientific notation with
            // 5 significant digits
            //   asciiFile << name;
            //   asciiFile << ' ';
            buffer << Settings::RAND_SEED;
            buffer << ' ';
            buffer << Settings::SOURCE_ALT;
            buffer << ' ';
            buffer << Settings::SOURCE_OPENING_ANGLE;
            buffer << ' ';
            buffer << Settings::TILT_ANGLE;
            buffer << ' ';
            buffer << event_nb; // 5
            buffer << ' ';
            buffer << ID;
            buffer << ' ';
            buffer << PDG_NB;
            buffer << ' ';
            buffer << time;
            buffer << ' ';
            buffer << energy;
            buffer << ' ';
            buffer << alt2; // 10
            buffer << ' ';
            buffer << lat;
            buffer << ' ';
            buffer << lon;
            buffer << ' ';
            buffer << dist_rad;
            buffer << ' ';
            buffer << ecef_x;
            buffer << ' ';
            buffer << ecef_y; // 15
            buffer << ' ';
            buffer << ecef_z;
            buffer << ' ';
            buffer << mom_x;
            buffer << ' ';
            buffer << mom_y;
            buffer << ' ';
            buffer << mom_z;
            buffer << ' ';
            buffer << number_beaming; // 20 // number_beaming == 0 for uniform and 1 for // gaussian
            buffer << ' ';
            buffer << Settings::SOURCE_LAT;
            buffer << ' ';
            buffer << Settings::SOURCE_LONG;
            buffer << ' ';
            buffer << G4endl;
            //
            NB_RECORDED++;
            //

            output_lines.push_back(buffer.str());
            //
            write_in_output_file();
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Analysis::write_in_output_file()
{

    if (output_lines.size() <= output_buffer_size)
    {
        return;
    }

    std::ofstream asciiFile;
    asciiFile.open(asciiFileName, std::ios::app);

    if (asciiFile.is_open())
    {
        for (const G4String &line : output_lines)
        {
            asciiFile << line;
        }

        asciiFile.close();
        output_lines.clear();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Analysis::write_in_output_file_endOfRun()
{
    if (output_lines.empty())
    {
        return;
    }

    std::ofstream asciiFile;
    asciiFile.open(asciiFileName, std::ios::app);

    if (asciiFile.is_open())
    {
        for (G4String &line : output_lines)
        {
            asciiFile << line;
        }

        asciiFile.close();
        output_lines.clear();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Function to create a directory if it does not already exist
bool Analysis::createDirectory(const std::string &dir)
{
    std::cout << "Creating directory:" << dir << std::endl;

    struct stat info;
    if (stat(dir.c_str(), &info) != 0)
    {
        // Directory does not exist, attempt to create it
        if (mkdir(dir.c_str(), 0755) == 0)
        {
            return true; // Directory created successfully
        }
        else
        {
            std::cerr << "Failed to create directory. Error: " << strerror(errno) << std::endl;
            return false; // Failed to create directory
        }
    }
    else if (info.st_mode & S_IFDIR)
    {
        std::cout << "Directory already exists." << std::endl;
        return false; // Directory already exists
    }
    else
    {
        std::cerr << "Path exists but is not a directory." << std::endl;
        return false; // Path exists but is not a directory
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....