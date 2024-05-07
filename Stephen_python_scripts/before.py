from gbm.data import PosHist
from gbm.time import Met

import math
 
def euler_from_quaternion(quaternion):
        """
        https://automaticaddison.com/how-to-convert-a-quaternion-into-euler-angles-in-python/
        Convert a quaternion into euler angles (roll, pitch, yaw)
        roll is rotation around x in radians (counterclockwise)
        pitch is rotation around y in radians (counterclockwise)
        yaw is rotation around z in radians (counterclockwise)
        """

        w, x, y, z = quaternion

        t0 = +2.0 * (w * x + y * z)
        t1 = +1.0 - 2.0 * (x * x + y * y)
        roll_x = math.atan2(t0, t1)

        t2 = +2.0 * (w * y - z * x)
        t2 = +1.0 if t2 > +1.0 else t2
        t2 = -1.0 if t2 < -1.0 else t2
        pitch_y = math.asin(t2)

        t3 = +2.0 * (w * z + x * y)
        t4 = +1.0 - 2.0 * (y * y + z * z)
        yaw_z = math.atan2(t3, t4)

        return roll_x, pitch_y, yaw_z # in radians

file = './glg_poshist_all_091214_v00.fit' #'./glg_poshist_all_191214_v01.fit'
dt = '2009-12-14T11:53:27.83'

met = Met.from_iso(dt)

poshist = PosHist.open(file)

lat = poshist.get_latitude(met.met)
Elon = poshist.get_longitude(met.met)
alt = poshist.get_altitude(met.met)

print('')
print('lat:', lat)
print('lon:', Elon)
print('alt:', alt*1e-3)
print('')

q = poshist.get_quaternions(met.met)
print('quaternion:', q)
print('')

roll, pitch, yaw = euler_from_quaternion(q)

print('[in degrees]')
print('roll (x):', math.degrees(roll))
print('pitch (y):', math.degrees(pitch))
print('yaw (z):', math.degrees(yaw))






