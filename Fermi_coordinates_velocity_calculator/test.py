import datetime
from satellite_coordinates import satellite_coordinates

if __name__ == "__main__":
    fermi_coordinates = satellite_coordinates(name='Fermi')

    #TGF-190324-00:31:52-33518
    input_datetime = datetime.datetime(year=2019, month=3, day=24, hour=0, minute=31, second=52, microsecond=335180)

    print('\n')

    lon, lat, alt, v_vec = fermi_coordinates.get_satellite_coordinates(input_datetime)
    print('Satellite name : ', fermi_coordinates.name)
    print('Requested time : ', input_datetime)
    print('Longitude (deg), Latitude (deg), Altitude (km) : ', round(lon,3), round(lat,3), round(alt,3))
    print('\n')
