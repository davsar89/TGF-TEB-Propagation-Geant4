#include "Settings.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Settings *Settings::instance = nullptr;

Settings *Settings::getInstance()
{
    if (instance == nullptr)
    {
        instance = new Settings;
    }

    return instance;
}

////// Global static pointer used to ensure a single instance of the class.
///// (singleton pattern)
// Settings *Settings::instance = 0;

// void
// Settings::diplay_settings()
// {
//    G4cout << G4endl;
//    G4cout << "*************************************************************" << G4endl;
//    G4cout << " SIMULATION SETTINGS : " << G4endl;
//    G4cout << "    Random Number Seed : " << Rand_seed << G4endl;
//    G4cout << "    Cached length for Integrator : " << CACHED_LENGTH << " m" << G4endl;
//    G4cout << "    Initial altitude : " << SOURCE_ALT << " km" << G4endl;
//    G4cout << "    Beaming : " << G4endl;
//    G4cout << "        Type : " << BEAMING_TYPE << G4endl;
//    G4cout << "        Angle : " << OPENING_ANGLE << " degrees; corresponds to" << G4endl;
//    G4cout << "          - max angle if Type = uniform;" << G4endl;
//    G4cout << "          - sigma if Type = gaussian (normal)" << G4endl;

//    G4cout << "    Source timimg (gaussian sigma) : " << SOURCE_SIGMA_TIME << " microsecond" << G4endl;

//    G4cout << "    Output Altitudes : " << G4endl;
//    G4cout << "        ";

//    for (unsigned int ii = 0; ii < record_altitudes.size(); ++ii)
//        {
//            G4cout << record_altitudes[ii] << " km, ";
//        }

//    G4cout << G4endl;

//    G4cout << "*************************************************************" << G4endl;
//    G4cout << G4endl;
// }

// G4String Settings::BeamingType() const
// {
//    return BEAMING_TYPE;
// }

// void Settings::set_BeamingType(const G4String &value)
// {
//    BEAMING_TYPE = value;
// }

// void Settings::set_AltMax_recorded()
// {
//    G4double max_alt = 0.;

//    for (int ii = 0; ii < nb_altitude_record; ++ii)
//        {
//            if (record_altitudes[ii] > max_alt)
//                {
//                    max_alt = record_altitudes[ii];
//                }
//        }

//    ALT_MAX_RECORDED = max_alt;
// }

// G4bool Settings::USE_STEP_MAX() const
// {
//    return USE_STEP_MAX_;
// }

// void Settings::set_MAG_FIELD_ON(const G4bool MAG_FIELD_BOOL)
// {
//    MAG_FIELD_ON_ = MAG_FIELD_BOOL;
// }

// G4bool Settings::MAG_FIELD_ON() const
// {
//    return MAG_FIELD_ON_;
// }
