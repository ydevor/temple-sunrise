# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:07:59 2020
@author: Jonathan Devor

Compare with HORIZONS (https://ssd.jpl.nasa.gov/horizons.cgi)
"""

import math
from skyfield import api, almanac, timelib   
# Skyfield Documentation: https://rhodesmill.org/skyfield/api.html
# To install the latest development version:  pip install https://github.com/skyfielders/python-skyfield/archive/master.zip


# Uses Bennett's empirical formula for calculating the refraction angle from the apparent altitude
# Reference https://en.wikipedia.org/wiki/Atmospheric_refraction
# Bennett, G.G. (1982), The Calculation of Astronomical Refraction in Marine Navigation, Journal of Navigation, 35 (2), 255–259
def Bennett_refraction(apparent_altitude_degrees):
    ang_degrees = apparent_altitude_degrees + (7.31 / (apparent_altitude_degrees + 4.4))
    return (1.0 / (60 * math.tan(ang_degrees * math.pi / 180)))


# date_tuple is a tuple of int or real numbers (year, month, day)
def print_date (date_tuple):
    ts = api.load.timescale()
    date_strings = ts.utc(*date_tuple).utc_jpl().split(' ')
    return (date_strings[0] + ' ' + date_strings[1])


# Returns a time string, adding +2 hours for the local Israel time zone
# time_tuple is a tuple of int or real numbers (hour, minute, seconds)
# Hour and minute must be int, seconds can be int or real
# Note: hours are assume to be befor 22:00 
def print_time_plus2 (time_tuple):
    # hours : minutes : seconds
    return ('%d:%02d:%08.5f' % (time_tuple[0] + 2, time_tuple[1], time_tuple[2]))


# Returns a date string
# The input date structure is assumed to contain 3 int
def print_proleptic_Gregorian_date(t):
    # year , month , day of month
    return ("%d,%02d,%02d" % (t.utc.year, t.utc.month, t.utc.day))


# ----------------

# Local ephemeride files
# Info: https://ssd.jpl.nasa.gov/?planet_eph_export
# Download from:  ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp
    
#path_ephem = r"C:\Users\jonathan.d\Desktop\jdevor\astro\de422.bsp"    # goes back to  -3000
#path_ephem = r"C:\Users\jonathan.d\Desktop\jdevor\astro\de431t.bsp"   # goes back to -13000
#path_ephem = r"C:\Users\jdevor\Desktop\Sky calc\de422.bsp"            # goes back to  -3000
path_ephem = r"C:\Users\jdevor\Desktop\Sky calc\de431t.bsp"            # goes back to -13000
ephem = api.load(path_ephem)

site_local = api.Topos('31.7781 N', '35.236 E')  # Jerusalem Temple coordinates
site_abs = ephem['earth'] + site_local

alt_of_first_light = 3.258    # apparent altitude of first light in degrees ; Mt. of Olives measurment
angular_radius_sun = 0.2666   # range: 0.2711 - 0.2621 degrees


goal_alt = alt_of_first_light - Bennett_refraction(alt_of_first_light) - angular_radius_sun
print ('Goal altitude:', goal_alt, '\n') # The astronomical altitude of the center of the Sun's disk, without atmospheric refraction

ts = api.load.timescale()
t_start = ts.utc(-514, 3, 1)  # scan range -- proleptic Gregorian dates (e.g. -959 = BC 960)
t_finsh = ts.utc(-514, 3, 20)

almanac_func = almanac.risings_and_settings(ephem, ephem['sun'], site_local, horizon_degrees=goal_alt)
t, y = almanac.find_discrete(t_start, t_finsh, almanac_func)

print ('The moment of the local sunrise on each day:\n')
print ('Gregorian     Julian Day        Julian Date      Time (UT1 +2)  Altitude   Azimuth')
print ('----------   -------------    ----------------   -------------   -------   -------')
for ti, yi in zip(t,y):
    if yi == True:     # sunrise
        position = site_abs.at(ti).observe(ephem['sun']).apparent()
        alt, az, distance = position.altaz()
        calendar_gregorian = ti.ut1_calendar()
        # Converts to Julian date when earlier than 2299161 = 15 October 1582 ; start of the Gregorian calendar
        calendar_corrected = timelib.compute_calendar_date(ti.ut1, 2299161)
        print(print_proleptic_Gregorian_date(ti), '  %.5f' % ti.ut1, '  ', print_date(calendar_corrected), ' ', print_time_plus2(calendar_gregorian[3:]), '  %.5f' % alt.degrees, '  %.5f' % az.degrees)
