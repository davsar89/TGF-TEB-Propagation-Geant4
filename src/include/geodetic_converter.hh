#pragma once

#include <cmath>

// A bunch of static functions to convert coordinates

// x y z are in meters, lat lon are in degrees, alt in m
namespace geod_conv {
    static double kSemimajorAxis = 6378137;
    static double kSemiminorAxis = 6356752.3142;
    static double kFirstEccentricitySquared = 6.69437999014 * 0.001;
    static double kSecondEccentricitySquared = 6.73949674228 * 0.001;
//static double kFlattening                = 1 / 298.257223563;

    class GeodeticConverter {
    public:

        GeodeticConverter() {
        }

        ~GeodeticConverter() {
        }

        // Default copy constructor and assignment operator are OK.

        static void
        geodetic2ecef(const double &latitude, const double &longitude, const double &altitude, double &x, double &y,
                      double &z)

        // x y z are in meters, lat lon are in degrees, alt in m
        {
            // Convert geodetic coordinates to ECEF.
            // http://code.google.com/p/pysatel/source/browse/trunk/coord.py?r=22
            double lat_rad = deg2Rad(latitude);
            double lon_rad = deg2Rad(longitude);
            double coslat = 0, sinlat = 0, sinlon = 0, coslon = 0;
            sincos(lat_rad, &sinlat, &coslat);
            sincos(lon_rad, &sinlon, &coslon);
            double xi = sqrt(1 - kFirstEccentricitySquared * sinlat * sinlat);
            x = (kSemimajorAxis / xi + altitude) * coslat * coslon;
            y = (kSemimajorAxis / xi + altitude) * coslat * sinlon;
            z = (kSemimajorAxis / xi * (1 - kFirstEccentricitySquared) + altitude) * sinlat;
        }

        static void
        ecef2Geodetic(const double &x, const double &y, const double &z, double &latitude, double &longitude,
                      double &altitude) {
            // Convert ECEF coordinates to geodetic coordinates.
            // J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates
            // to geodetic coordinates," IEEE Transactions on Aerospace and
            // Electronic Systems, vol. 30, pp. 957-961, 1994.
            double r = sqrt(x * x + y * y);
            double Esq = kSemimajorAxis * kSemimajorAxis - kSemiminorAxis * kSemiminorAxis;
            double F = 54 * kSemiminorAxis * kSemiminorAxis * z * z;
            double G = r * r + (1 - kFirstEccentricitySquared) * z * z - kFirstEccentricitySquared * Esq;
            double C = (kFirstEccentricitySquared * kFirstEccentricitySquared * F * r * r) / pow(G, 3);
            double S = cbrt(1 + C + sqrt(C * C + 2 * C));
            double P = F / (3 * pow((S + 1 / S + 1), 2) * G * G);
            double Q = sqrt(1 + 2 * kFirstEccentricitySquared * kFirstEccentricitySquared * P);
            double r_0 = -(P * kFirstEccentricitySquared * r) / (1 + Q) + sqrt(
                    0.5 * kSemimajorAxis * kSemimajorAxis * (1 + 1.0 / Q) -
                    P * (1 - kFirstEccentricitySquared) * z * z / (Q * (1 + Q)) - 0.5 * P * r * r);
            double U = sqrt(pow((r - kFirstEccentricitySquared * r_0), 2) + z * z);
            double V = sqrt(pow((r - kFirstEccentricitySquared * r_0), 2) + (1 - kFirstEccentricitySquared) * z * z);
            double Z_0 = kSemiminorAxis * kSemiminorAxis * z / (kSemimajorAxis * V);
            altitude = U * (1 - kSemiminorAxis * kSemiminorAxis / (kSemimajorAxis * V));
            latitude = rad2Deg(atan((z + kSecondEccentricitySquared * Z_0) / r));
            longitude = rad2Deg(atan2(y, x));
        }

    private:

        static inline double rad2Deg(const double &radians) {
            return (radians / M_PI) * 180.0;
        }

        static inline double deg2Rad(const double &degrees) {
            return (degrees / 180.0) * M_PI;
        }
    };

// class GeodeticConverter
} // namespace geodetic_conv
