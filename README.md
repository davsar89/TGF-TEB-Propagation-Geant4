TGF-TEB-Propagation-Geant4
=====
Model for Terrestrial Gamma-ray Flashes (TGF) and associated electrons and positrons (TEB) propagation in Earth atmosphere and environment (incl. magnetic field). Based on [GEANT4](https://geant4.web.cern.ch/).
=====

contact : <david (dot) sarria (at) uib (dot) no>

## Generalities
- Propagation of photons, electrons and positron in Earth's environment (atmosphere, ionosphere, magnetosphere), in the context of Terrestrial Gamma-ray Flashes (TGF) and associated electrons and positrons beams.
- Uses mostly Geant4 features. See [documentation](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/index.html "Geant4 documentation").
- Integrates the [NRL-MSISE-00 model](https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/) for the atmosphere and the [IGRF-12](http://wdc.kugi.kyoto-u.ac.jp/igrf/index.html) or WMM model for the magnetic field.
- The coordinate system is set so that the Geant4 X,Y,Z coordinates correspond to the ECEF X,Y,Z coordinates (earth-centered, earth-fixed).
- This code is probably not perfect. Feel free to suggest improvements.

## Compilation, installation
- The source code of `TGF-TEB-Propagation-Geant4` is located in `src/` and the build should be done in the folder `build/`.
- Requires [Geant4](https://geant4.web.cern.ch/) compiled and [installed](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/index.html) properly. Recommended to use geant4-10-07-patch-03 (19-November-2021). Minor changes in the source code may be required for other versions.
### Linux
- Geant4 needs to be compiled, installed and set-up (environement) properly. Easy installation scripts of Geant4 for Linux are provided [here](https://github.com/DavidSarria89/GEANT4-easy-install-scripts)
- Once Geant4 is installed and set-up properly, open a terminal in the `build/` folder and run `cmake ../` (to generate the `makefile` using CMake) and then `make` to compile. It will make the executable file `build/tgf_propa`. 

## Simulation Settings:
Most of settings can be adjusted in `src/Settings.cc`. In particular:
- `Settings::RECORD_PHOT_ONLY` = boolean to record only photons (not recording electrons and positrons, and also turning OFF the magnetic field).
- `Settings::RECORD_ELEC_POSI_ONLY` = boolean to record only electrons and positrons (not recording photons).
- `Settings::record_altitude` = record altitude (in km), default is 408 km.
- `Settings::SOURCE_ALT` = TGF source altitude in km, default is 15 km
- `Settings::SOURCE_LAT` = TGF source latitude in degrees
- `Settings::SOURCE_LONG` = TGF source longitude in degrees
- `Settings::SOURCE_OPENING_ANGLE` = half-cone TGF opening angle in degrees if "Uniform" is selected for `Settings::BEAMING_TYPE`. If "Gaussian" is selected for `Settings::BEAMING_TYPE`, it is the sigma of the gaussian distribution.
- `Settings::TILT_ANGLE` = TGF tilt angle in degrees. Default is 0 degrees.
- `Settings::BEAMING_TYPE` = TGF beaming type, that is a string that values "Uniform" or "Gaussian" for isotropic (within half cone angle) or gaussian distribution.
- `Settings::SOURCE_SIGMA_TIME` = TGF sigma time. Assumes the TGF has an intrinsic duration, that has Gaussian (=normal) distribution. The parameter is the sigma of this distribution, in microseconds. Default is 0.
### Other settings:
- Two modes are possible: `visualization` and `run`. `visualization` will show the 3D geometry (simplified Earth) and particle track. `run` will not show any 3D visualization, to run the code as quickly as possible. By default, the mode is set to `visu` if no input argument for the executable is specified and `run` otherwise. This can be changed by editing the `G4String` variable `Mode` in the main function located in the source file `src/tgf_propa.cc`, that can be set to `"visu"` or `"run"`. Default mode is `visu`.
- Primary Generator is a point source, with adjustable altitude and geometry. See `src/src/PrimaryGeneratorAction.hh` and `src/src/PrimaryGeneratorAction.cc`
- Atmosphere density is not constant with altitude, it evolves ~exponentially. However, Geant4 can only handle volumes with constant density, therefore the atmosphere is simulated by 256 exponentially spaced layers, each with constant density, covering altitude from 1 km to 150 km (negligible above). The Earth is also separed into 45 "latitude" and 6 "longitude" parts (which are actually spherical coordinates `theta` and `phi`). This can be changed in the source code, with the `src/src/DetectorConstruction.hh` and `src/src/DetectorConstruction.cc` files.
- If required, Magnetic Field can be turned ON with the Setting : `Settings::MAG_FIELD_ON` set to `true`. Magnetic field is always turned OFF below 30 km altitude (where it is negligible), for performance.
- `Settings::MAGNETIC_FIELD_MODEL` = magnetic field model to use. Should be set to `"IGRF"`.

## The program executable can accept input parameters in this order: 
- `number_st` = number of particles (TGF photons) to shoot, the program is stopped after they are shot.
- `Settings::SOURCE_ALT` = TGF source altitude in km
- `Settings::SOURCE_LAT` = TGF source latitude in degrees
- `Settings::SOURCE_LONG` = TGF source longitude in degrees
- `Settings::SOURCE_SIGMA_TIME` = TGF sigma time. Assumes the TGF has an intrinsic duration, that has Gaussian (=normal) distribution. The parameter is the sigma of this distribution, in microseconds
- `Settings::SOURCE_OPENING_ANGLE` = half-cone TGF opening angle in degrees. If "Gaussian" is selected for `Settings::BEAMING_TYPE`, it is the sigma of the gaussian distribution.
- `Settings::TILT_ANGLE` = TGF tilt angle in degrees
- `Settings::BEAMING_TYPE` = TGF beaming type, that is a string that values `"Uniform"` or `"Gaussian"` for isotropic or gaussian distribution
- `Settings::record_altitude` = record altitude (in km) of the TGF (and secondary electron and positrons), default is 408 km.

## Program output:
- Recorded particles are outputed as a list (one by one, line by line) in files located in the `build/output_ascii/`. See `build/README_output.txt` to find which quantity is in which column. The same electron/positron can be recorded several time since it may cross several times the limit altitude (it is the same particle if it has the same random seed and the same ID).

## Additional information:
- By default, the code uses the `G4EmStandardPhysics_option1` physics list. This can be changed inside the source file `src/src/PhysicsList.cc`.
