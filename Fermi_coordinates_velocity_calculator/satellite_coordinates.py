import sys

###### check python version, required >= 3.6
if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("This script requires Python version >= 3.6")
    sys.exit(1)
######
import calendar
import re
import datetime
import numpy as np
import skyfield.sgp4lib as sgp4lib
from pathlib import Path
from skyfield import api
import pyproj
import math

######

class satellite_coordinates:
    ##

    def __init__(self, name):
        """
        :param name: satellite name, must be 'ISS' or 'Fermi'
        """
        self.name = name

        self.ISS_data_file = Path("./dataFiles/ISS_orbital_info.txt")
        self.Fermi_data_file = Path("./dataFiles/Fermi_GLAST_orbital_info.txt")
        self.AGILE_data_file = Path("./dataFiles/AGILE_orbital_info.txt")

        if self.name == "ISS":
            self.data_file_path = self.ISS_data_file
        elif self.name == "Fermi":
            self.data_file_path = self.Fermi_data_file
        elif self.name == "AGILE":
            self.data_file_path = self.AGILE_data_file
        else:
            raise Exception("name input argument should be 'Fermi' or 'ISS' or 'AGILE' ")

        self.datetimes, self.TLE_line_1, self.TLE_line_2 = self.read_satellite_TLE_data(self.data_file_path)

        self.time_value=[]
        for dt in self.datetimes:
            self.time_value.append(calendar.timegm(dt.utctimetuple()))
            # print(calendar.timegm(dt.utctimetuple()))
        self.time_value = np.array(self.time_value)


        self.ts = api.load.timescale()

        self.ecef = pyproj.Proj(proj='geocent', ellps='GRS80')
        self.lla = pyproj.Proj(proj='latlong', ellps='GRS80')

    ##

    def nearestDate(self, base):
        base_timestamp = calendar.timegm(base.utctimetuple())
        differences = np.abs(base_timestamp - self.time_value)
        # print(differences)
        return np.argmin(differences)

    ##

    def read_satellite_TLE_data(self,datafile_path):
        """
        convert textfile (txt) containing list of TLE of a given satellite to python lists of times and TLE
        :return: datetimes : list of datetime read in the file
                 dates_year : list of corresponding years, redundant information with datetimes
                 dates_doy : list of corresponding days of year, redundant information with datetimes
                 TLE_line_1 : first line of TLE at given 'datetimes' date
                 TLE_line_2 : second line of TLE at given 'datetimes' date
        """


        TLE_line_1 = []
        TLE_line_2 = []
        datetimes = []

        ## raw text data with TLE text lines
        with open(self.data_file_path, "r") as f:
            self.datafile = f.read()

        splited_txt = self.datafile.split('\n')
        splited_txt.pop(-1)

        # print(splited_txt[-1])

        for ii in range(0, int(len(splited_txt)), 2):

            TLE_line_1.append(splited_txt[ii])
            TLE_line_2.append(splited_txt[ii + 1])

            epoch_str = TLE_line_1[-1][18:32]

            year_tle = float(epoch_str[0:2])

            if year_tle > 70:
                 year_tle = year_tle + 1900
            else:
                 year_tle = year_tle + 2000

            integ = np.floor(float(epoch_str))
            day_fract = float(epoch_str) - integ
            DOY = float(epoch_str[2:5])

            date = datetime.datetime(year=int(year_tle), month=1,day=1) + datetime.timedelta(days=DOY-1) + datetime.timedelta(days=day_fract)

            datetimes.append(date)
            # print(epoch_str)
            # print(date)

        return datetimes, TLE_line_1, TLE_line_2

    ##

    def get_satellite_coordinates(self, input_datetime):
        """
    Calculates the longitude (deg), latitude (deg), altitude (km) and velocity vector (normalized) of a satellite at a given time
        :param input_datetime: python datetime structure identifying the time where the coordinates of the satellite is wanted
        :return: longitude (deg), latitude (deg), altitude (km), velocity vector (normalized)
        """

        # input_datetime = datetime.datetime(year=2018, month=9, day=16, hour=13,minute=15,second=40)

        ## finding which TLE is the closest to the time we want

        closest_idx = self.nearestDate(input_datetime)

        line1 = self.TLE_line_1[closest_idx]
        line2 = self.TLE_line_2[closest_idx]

        # print(closest_idx)
        # print(line1)
        # print(line2)

        satellite = sgp4lib.EarthSatellite(line1, line2)

        # print(satellite)

        tttt = self.ts.utc(input_datetime.year, input_datetime.month, input_datetime.day, input_datetime.hour,
                      input_datetime.minute, input_datetime.second + input_datetime.microsecond/1.0e6)

        # tttt = ts.utc(input_datetime.year, input_datetime.month, input_datetime.day, input_datetime.hour,
                      # input_datetime.minute, input_datetime.second)

        # print(tttt)

        position, velocity, error = satellite.ITRF_position_velocity_error(tttt)
        au_to_Km = 149597870.700

        v_vec_itrf = velocity / math.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2) # unit vector

        position = np.asarray(position) * au_to_Km * 1000.0    # to meters

        lon, lat, alt = pyproj.transform(self.ecef, self.lla, position[0], position[1],
                                         position[2], radians=False)
        alt = alt / 1000.0

        return lon, lat, alt, v_vec_itrf
