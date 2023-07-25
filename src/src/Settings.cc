#include "Settings.hh"

namespace Settings {

    long RAND_SEED = 777 ;

    double SOURCE_LAT = -7.0492;  // degree
    double SOURCE_LONG = 55.9129; // degree
    double SOURCE_ALT = 15.;    // km
    double record_altitude = 408.0;

    double SOURCE_OPENING_ANGLE = 30.; // degree
    G4String BEAMING_TYPE = "Gaussian";

    double TILT_ANGLE = 0.0;

    double SOURCE_SIGMA_TIME = 0.; // microsecond

    /// Below variables should stay CONSTANT

    G4bool RECORD_PHOT_ONLY = false;      // record only photons

    // force max step only for layers where particles are recorded
    G4bool USE_STEP_MAX_for_record = true;
    double STEP_MAX_RECORD_AREA = 100.0 * meter;

    bool CHECK_ELECTRON_TRACKING_IN_MAG_FIELD = true; // for leptons : check step size compared to Larmor radius

    G4bool OUTPUT_TO_ASCII_FILE = true;

    G4String MAGNETIC_FIELD_MODEL = "IGRF"; // "WMM" or "IGRF" or "GEOLIB", should be WMM on MAC
    G4String MAGNETIC_FIELD_MODEL_FOR_GEOLIB = "IGRF"; // "WMM" or "IGRF" or "EMM"
    // WARNING : EMM is about 1000 times slower.

    double dr_over_R = 0.05; // stepping parameter for ionization, default is 0.2, that may be too high
    double TIME_LIMIT = 0.2 * second;
    double MIN_ENERGY_OUTPUT = 10.0 * keV;

    // date
    int dt_year = 2019;
    double dt_month = 3;
    double dt_day = 24;
    double dt_hour = 0;
    double dt_minute = 31;
    double dt_second = 53;
    double dt_microsecond = 135444;

    // Earth radius
    double earthRadius = 6378.137 * km;

    G4bool OUTPUT_ALT_LAYERS_TO_FILE =
            false; // output list of altitude and densities of layer to file (for debug)

    // to record particles only inside a specific latitude and longitude window

    bool RECORD_ONLY_IN_WINDOW = false;

    double iss_lat = 0.1575;
    double iss_lon = 55.3015;

    record_window RECORD_WIN{iss_lat - 5.0, iss_lat + 5.0, iss_lon - 5.0, iss_lon + 5.0}; // degrees
};
