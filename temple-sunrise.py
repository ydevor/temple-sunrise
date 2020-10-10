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


# for altitude set: altaz = 0
# for azimuth set:  altaz = 1
# Returns the value, in degrees, for the given time offset
def get_azimuth_t_offset(site_abs, sun_abs, ti, altaz, offset_days):
    ti_offset = ts.tt_jd(ti.tt + offset_days)
    position = site_abs.at(ti_offset).observe(sun_abs).apparent()    
    return(position.altaz()[altaz].degrees)
    

# for altitude set: altaz = 0
# for azimuth set:  altaz = 1
# Returns the time object closes to the goal value
# Assumes that the value changes monotonically within the search range
# max_err_days = the maximum search error.
def find_goal_altaz(site_abs, sun_abs, ti, search_range_days, altaz, goal_val, max_err_days):
    delta_t_lo = -search_range_days
    delta_t_hi = +search_range_days
    val_lo = get_azimuth_t_offset(site_abs, sun_abs, ti, altaz, delta_t_lo) # used only for the sanity check and swap
    val_hi = get_azimuth_t_offset(site_abs, sun_abs, ti, altaz, delta_t_hi)
    if val_lo > val_hi:
        delta_t_lo, delta_t_hi = delta_t_hi, delta_t_lo
        val_lo, val_hi = val_hi, val_lo
    if (goal_val < val_lo) or (goal_val > val_hi):
        print ("ERROR: goal value is out of range", goal_val, val_lo, val_hi)
        raise
    while abs(delta_t_hi - delta_t_lo) > max_err_days:
        delta_t_mid = (delta_t_lo + delta_t_hi) / 2
        val_mid = get_azimuth_t_offset(site_abs, sun_abs, ti, altaz, delta_t_mid)
        if val_mid < goal_val:
            delta_t_lo = delta_t_mid
        else:
            delta_t_hi = delta_t_mid
    return (ts.tt_jd(ti.tt + delta_t_mid))


# Given the distance to the center of the Sun in kilometers
# Returns the Sun's disk angular radius (in degrees)
# Note: the solar radius is 696,340 km
def calc_sun_disk_radius(distance):
    return (math.asin(696340 / distance.km) * 180 / math.pi)


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

# Ephemeride files
# Info: https://ssd.jpl.nasa.gov/?planet_eph_export
# Download from:  ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp
    
#path_ephem = r"C:\Users\jonathan.d\Desktop\jdevor\astro\de422.bsp"    # goes back to  -3000
#path_ephem = r"C:\Users\jonathan.d\Desktop\jdevor\astro\de431t.bsp"   # goes back to -13000
#path_ephem = r"C:\Users\jdevor\Desktop\Sky calc\de422.bsp"            # goes back to  -3000
path_ephem = r"C:\Users\jdevor\Desktop\Sky calc\de431t.bsp"            # goes back to -13000
ephem = api.load(path_ephem)

site_local = api.Topos('31.7781 N', '35.236 E')  # Jerusalem Temple coordinates
site_abs = ephem['earth'] + site_local  # the absolute location of the observer
sun_abs = ephem['sun']

alt_of_first_light = 3.258    # apparent altitude of first light in degrees ; Mt. of Olives measurment
approx_sun_radius = 0.2666    # angular radius range: 0.2711 - 0.2621 degrees

goal_az = 99.7                # degrees ; for azimuth search
max_err_days = 0.0000001      # = 0.0086 sec ; max search error
search_range_days = 0.1       # = 2.4 hours ; search at either side of the approx sunrise


goal_alt_approx = alt_of_first_light - Bennett_refraction(alt_of_first_light) - approx_sun_radius
print ('Goal altitude:', goal_alt_approx, '\n') # The astronomical altitude of the center of the Sun's disk, without atmospheric refraction

ts = api.load.timescale()
t_start = ts.utc(-954, 10, 9)  # scan dates: proleptic Gregorian dates (e.g. -959 = BC 960)
t_finsh = ts.utc(-954, 10, 20)

almanac_func = almanac.risings_and_settings(ephem, sun_abs, site_local, horizon_degrees = goal_alt_approx)
t, y = almanac.find_discrete(t_start, t_finsh, almanac_func)

print ("For each day in range, we check:")
print (' Approximate local sunrise, with fixed Sun radius')
print (' More accurate local sunrise, with variable Sun radius')
print (' When the Sun is near an azimuth of', goal_az)
print ()

print ('Gregorian     Julian Day        Julian Date      Time (UT1 +2)  Altitude   Azimuth')
print ('----------   -------------    ----------------   -------------   -------   -------')
for ti, yi in zip(t,y):
    if yi == True:     # sunrise
        position = site_abs.at(ti).observe(sun_abs).apparent()
        alt, az, distance = position.altaz()
        calendar_gregorian = ti.ut1_calendar()
        # Converts to Julian date when earlier than 2299161 = 15 October 1582 ; start of the Gregorian calendar
        calendar_corrected = timelib.compute_calendar_date(ti.ut1, 2299161)
        print(print_proleptic_Gregorian_date(ti), '  %.5f' % ti.ut1, '  ', print_date(calendar_corrected), ' ', print_time_plus2(calendar_gregorian[3:]), '  %.5f' % alt.degrees, '  %.5f' % az.degrees)
        #---
        sun_disk_radius = calc_sun_disk_radius(distance)
        goal_alt_accurate = goal_alt_approx + (approx_sun_radius - sun_disk_radius)
        ti_goal = find_goal_altaz(site_abs, sun_abs, ti, search_range_days, 0, goal_alt_accurate, max_err_days)
        position = site_abs.at(ti_goal).observe(sun_abs).apparent()
        alt, az, distance = position.altaz()
        calendar_gregorian = ti_goal.ut1_calendar()
        print(print_proleptic_Gregorian_date(ti_goal), '  %.5f' % ti_goal.ut1, '  ', print_date(calendar_corrected), ' ', print_time_plus2(calendar_gregorian[3:]), '  %.5f' % alt.degrees, '  %.5f' % az.degrees)
        #---
        ti_goal = find_goal_altaz(site_abs, sun_abs, ti, search_range_days, 1, goal_az, max_err_days)
        position = site_abs.at(ti_goal).observe(sun_abs).apparent()
        alt, az, distance = position.altaz()
        calendar_gregorian = ti_goal.ut1_calendar()
        print(print_proleptic_Gregorian_date(ti_goal), '  %.5f' % ti_goal.ut1, '  ', print_date(calendar_corrected), ' ', print_time_plus2(calendar_gregorian[3:]), '  %.5f' % alt.degrees, '  %.5f' % az.degrees)
        print()
