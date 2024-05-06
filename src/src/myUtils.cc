

#include "myUtils.hh"

/////////////////////////////////////////////////////////

namespace myUtils
{
    ///
    double get_wall_time()
    {
        std::chrono::high_resolution_clock m_clock;
        double time = std::chrono::duration_cast<std::chrono::microseconds>(
                          m_clock.now().time_since_epoch())
                          .count();
        return time;
    }
    ///

    ///
    long reversDigits(long num)
    {
        long rev_num = 0;
        while (num > 0)
        {
            rev_num = rev_num * 10 + num % 10;
            num = num / 10;
        }
        return rev_num;
    }

    bool addOvf(long *result, long a, long b)
    {
        // adding two long integers but checking if the sum does not overflow
        if (a > std::numeric_limits<long>::max() - b)
            return false;
        else
        {
            *result = a + b;
            return true;
        }
    }

    /////////////////////////////////////////////////////////
    // returns ISS altitude (km) at a given time
    double get_ISS_altitude_at_time(const int year,
                                    const int month,
                                    const int day,
                                    const int hour, // UTC
                                    const int minute,
                                    const int second,
                                    const int microsecond)
    {
        const double pii = 3.14159265359;
        double dateTime_diff;
        //
        sgp4::DateTime dt_TGF(year, month, day, hour, minute, second, microsecond);

        int64_t time_tgf_in_ticks = dt_TGF.Ticks();

        // READ ISS TLE (Two Line Element) data info and find the TLE with closest epoch (time)
        std::ifstream inFile("./TLE_data/ISS_orbital_info.txt");

        std::string a_line;
        std::vector<std::string> lines1;
        std::vector<std::string> lines2;
        std::vector<std::string> Epoch_Year; // last two digits of year
        std::vector<std::string> Epoch_day;  // day of the year and fractional portion of the day
        std::vector<sgp4::DateTime> Epochs_dt;

        int64_t smallest_dt = 6418429131356360; // in ticks
        int idx_of_smallest_dt0 = 0;
        int idx = 0;

        bool flag = true; // flag to alternate file reading to store in lines1 and lines2 vectors

        if (inFile)
        {
            while (getline(inFile, a_line, '\n'))
            {
                if (flag)
                {
                    lines1.push_back(a_line);
                    Epoch_Year.push_back(a_line.substr(18, 2));
                    Epoch_day.push_back(a_line.substr(20, 12));
                    flag = false;

                    unsigned int extract_year = std::stoi(Epoch_Year.back());
                    double extract_day = std::stod(Epoch_day.back());

                    if (extract_year < 57)
                        extract_year += 2000;
                    else
                        extract_year += 1900;

                    int64_t time_tle_in_ticks = sgp4::DateTime(extract_year, extract_day).Ticks();

                    dateTime_diff = std::abs(time_tle_in_ticks - time_tgf_in_ticks);

                    if (dateTime_diff < smallest_dt)
                    {
                        smallest_dt = dateTime_diff;
                        idx_of_smallest_dt0 = idx;
                    }
                    idx++;
                }
                else
                {
                    lines2.push_back(a_line);
                    flag = true;
                }
            }
        }
        else
        {
            G4cout << "ERROR : impossible to read the input file." << G4endl;
        }

        ///

        Tle tle = Tle("ISS                ",
                      lines1[idx_of_smallest_dt0],
                      lines2[idx_of_smallest_dt0]);
        SGP4 sgp4(tle);

        //        std::cout << lines1[idx_of_smallest_dt0] << std::endl; // print for debug
        //        std::cout << lines2[idx_of_smallest_dt0] << std::endl; // print for debug
        //    std::cout << tle << std::endl; // print for debug

        /*
 * calculate satellite position
 */
        Eci eci = sgp4.FindPosition(dt_TGF);
        /*
 * convert satellite position to geodetic coordinates
 */
        CoordGeodetic geo = eci.ToGeodetic(); // warning: output latitude and longitude are in radians

        geo.latitude = geo.latitude * 180.0 / pii;
        geo.longitude = geo.longitude * 180.0 / pii;

        std::cout << "ISS position : " << geo.altitude << " " << geo.latitude << " " << geo.longitude << std::endl; // print for debug

        return geo.altitude;
    }

    /////////////////////////////////////////////////////////

    vector3D enu2ecefvFormula(const double uEast, const double vNorth, const double wUp, const double lat0, const double lon0)
    {

        double cosPhi;
        double sinPhi;
        double cosLambda;
        double sinLambda;

        GeographicLib::Math::sincosd(lat0, sinPhi, cosPhi); // input in degrees
        GeographicLib::Math::sincosd(lon0, sinLambda, cosLambda);

        const double t = cosPhi * wUp - sinPhi * vNorth;

        return vector3D{cosLambda * t - sinLambda * uEast, sinLambda * t + cosLambda * uEast, sinPhi * wUp + cosPhi * vNorth};
    }

    vector3D ned2ecef(const double uNorth, const double vEast, const double wDown, const double lat0, const double lon0)
    {
        return enu2ecefvFormula(vEast, uNorth, -wDown, lat0, lon0);
    }

    ///

    std::string convertToString(char a[37], int size)
    {
        int i;
        std::string s = "";
        for (i = 0; i < size; i++)
        {
            s = s + a[i];
        }
        return s;
    }

    std::string change_letters_to_numbers(std::string my_str)
    {
        std::replace(my_str.begin(), my_str.end(), 'a', '1');
        std::replace(my_str.begin(), my_str.end(), 'b', '2');
        std::replace(my_str.begin(), my_str.end(), 'c', '3');
        std::replace(my_str.begin(), my_str.end(), 'd', '4');
        std::replace(my_str.begin(), my_str.end(), 'e', '5');
        std::replace(my_str.begin(), my_str.end(), 'f', '6');
        std::replace(my_str.begin(), my_str.end(), 'g', '7');
        std::replace(my_str.begin(), my_str.end(), 'h', '8');
        std::replace(my_str.begin(), my_str.end(), 'i', '9');
        std::replace(my_str.begin(), my_str.end(), 'j', '1');
        std::replace(my_str.begin(), my_str.end(), 'j', '2');
        std::replace(my_str.begin(), my_str.end(), 'l', '3');
        std::replace(my_str.begin(), my_str.end(), 'm', '4');
        std::replace(my_str.begin(), my_str.end(), 'n', '5');
        std::replace(my_str.begin(), my_str.end(), 'o', '6');
        std::replace(my_str.begin(), my_str.end(), 'p', '7');
        std::replace(my_str.begin(), my_str.end(), 'k', '8');
        std::replace(my_str.begin(), my_str.end(), 'r', '9');
        std::replace(my_str.begin(), my_str.end(), 's', '1');
        std::replace(my_str.begin(), my_str.end(), 't', '2');
        std::replace(my_str.begin(), my_str.end(), 'u', '3');
        std::replace(my_str.begin(), my_str.end(), 'v', '4');
        std::replace(my_str.begin(), my_str.end(), 'w', '5');
        std::replace(my_str.begin(), my_str.end(), 'x', '6');
        std::replace(my_str.begin(), my_str.end(), 'y', '7');
        std::replace(my_str.begin(), my_str.end(), 'z', '8');
        std::replace(my_str.begin(), my_str.end(), '-', '0');
        return my_str;
    }

    long generate_a_unique_ID()
    {

        uuid_t uu;
        uuid_generate(uu);
        char uuid[37];
        uuid_unparse(uu, uuid);
        //        std::cout << uuid << std::endl;
        int mys = sizeof(uuid) / sizeof(char);
        std::string my_str = std::string(convertToString(uuid, mys));

        my_str = change_letters_to_numbers(my_str);

        my_str = my_str.substr(1, my_str.size() - 27);
        long output = std::stol(my_str);
        //        int dummy = 1 + 1;
        std::cout << output << std::endl;
        return output;
    }
    
    long generateUniqueRandomLong() // alternative to generate_a_unique_ID
    {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = now.time_since_epoch();
        auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
        auto pid = getpid();
        long base = static_cast<long>((nanoseconds % 1000000000L) ^ (pid % 100000L));
        std::random_device rd;
        std::mt19937 eng(rd());
        std::uniform_int_distribution<long> distr;
        return base ^ distr(eng);
    }
} // namespace myUtils
