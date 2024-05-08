import datetime
from satellite_coordinates import satellite_coordinates
import os
import pyproj
import math
import numpy as np

#Fermi TEB: dt = '2009-12-14T11:53:27.83'

# test for ISS:
# 2019-Mar-24 00:31:53.135444
# should be close to 0.157454678525265    55.3015430661528

def gps_to_ecef_pyproj(lat, lon, alt):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)

    return x, y, z

if __name__ == "__main__":

    #sat_coordinates = satellite_coordinates(name='ISS')
    sat_coordinates = satellite_coordinates(name='Fermi')

    #input_datetime = datetime.datetime(year=2019, month=3, day=24, hour=0, minute=31, second=53, microsecond=135444)
    input_datetime = datetime.datetime(year=2009, month=12, day=14, hour=11, minute=53, second=27, microsecond=830000)
    
    lon, lat, alt,v_vec = sat_coordinates.get_satellite_coordinates(input_datetime)

    X,Y,Z = gps_to_ecef_pyproj(lat, lon, alt)

    # Z: normalized inverted position vector
    LVLH_baseZ = np.array([-1.0*X,-1.0*Y,-1.0*Z])
    LVLH_baseZ = LVLH_baseZ/np.linalg.norm(LVLH_baseZ)

    v_vec_val = v_vec[0]**2 + v_vec[1]**2 + v_vec[2]**2

    v_vec[0] = v_vec[0]/ v_vec_val
    v_vec[1] = v_vec[1]/ v_vec_val
    v_vec[2] = v_vec[2]/ v_vec_val

    # Y: Z cross velocity
    LVLH_baseY = np.cross(LVLH_baseZ, v_vec)
    LVLH_baseY = LVLH_baseY/np.linalg.norm(LVLH_baseY)
    # X: Y cross Z
    LVLH_baseX = np.cross(LVLH_baseY, LVLH_baseZ)
    LVLH_baseX = LVLH_baseX/np.linalg.norm(LVLH_baseX)

    print('Satellite name : ', sat_coordinates.name)
    print('Requested time : ', input_datetime)
    print(f'Longitude (deg), Latitude (deg), Altitude (km) : {round(lon,4)}, {round(lat,4)}, {round(alt,4)}')
    print('Velocity vector (ECEF, unit) : ', round(v_vec[0],8), round(v_vec[1],8), round(v_vec[2],8))
    print('Position vector (ECEF, unit) : ', round(LVLH_baseZ[0],8), round(LVLH_baseZ[1],8), round(LVLH_baseZ[2],8))
    print('     ')
    print('LVLH base vector X (ECEF) : ', round(LVLH_baseX[0],8), round(LVLH_baseX[1],8), round(LVLH_baseX[2],8))
    print('LVLH base vector Y (ECEF) : ', round(LVLH_baseY[0],8), round(LVLH_baseY[1],8), round(LVLH_baseY[2],8))
    print('LVLH base vector Z (ECEF) : ', round(LVLH_baseZ[0],8), round(LVLH_baseZ[1],8), round(LVLH_baseZ[2],8))
    print('\n')

