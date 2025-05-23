from astropy.time import Time, TimeDelta
def oldtime2newtime(old_time_mjd, old_time_utc_s):

    """
    This is just a function that might turn out to be useful to check dates and makes sure they correspond to the correct epoch
    Takes a time in mjdn + seconds in a day (old fdets format) as input
    Spits out the corresponding datetime in utc
    """
    epoch = Time(old_time_mjd, format='mjd', scale='utc')
    seconds = TimeDelta(old_time_utc_s, format='sec')
    time = epoch + seconds
    print(time)
    print(time.isot)
    print(time.datetime)
    return time.datetime

# Example
oldtime2newtime(5.6654000000000000e+04, 7.0920500000000000e+04 )