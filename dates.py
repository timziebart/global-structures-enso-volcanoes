

from data_loader import NCEP_NCAR

import datetime as dt
import argparse




default_begin_date = dt.date(NCEP_NCAR.begin_year, 1, 1)
default_end_date = dt.date(NCEP_NCAR.end_year, 12, 31)

def parse_date(s):
    try:
        date = dt.datetime.strptime(s, "%Y-%m-%d").date()
        if not (default_begin_date <= date <= default_end_date):
            raise ValueError
        return date
    except ValueError:
        msg = "Not a valid date '{0}', should be 'yyyy-mm-dd' and between {} and {}".format(s, default_begin_date, default_end_date)
        # TODO: check why the message is not shown correctly
        raise argparse.ArgumentTypeError(msg)

def isleap(year):
    """Return True for leap years, False for non-leap years."""
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

def getNumFeb29s(dat1, dat2, ignore_dat2 = False, full_return = False): # include only dat2 ifboth if they are Feb 29
    # ignore if dat1 is Feb 29, but do not ignore dat2
    # if ignore_dat2 == True, then ignore both
    assert dat1 != dat2
    pre = +1
    if dat1 > dat2:
        pre = -1
        _pre, years = getNumFeb29s(dat2, dat1, ignore_dat2 = True, full_return = True)
        assert _pre == +1
    else:
        years = list(range(dat1.year, dat2.year + 1))
        years = list(filter(isleap, years))
        if dat1.year in years:
            if dat1 >= dt.date(dat1.year, 2, 29):
                years.remove(dat1.year)
        if dat2.year in years:
            if dat2 <= dt.date(dat2.year, 2, 29):
                years.remove(dat2.year)
    nyears = len(years)
    if not ignore_dat2:
        if dat2.month == 2 and dat2.day == 29:
            assert not dat2.year in years
            nyears += 1
            years.append(dat2.year)

    if full_return:
        return pre, years
    #else
    return pre * nyears

def sumdate(date, days):
    if not days:
        return date
    date2 = date + dt.timedelta(days=days)
    num_add = getNumFeb29s(date, date2)
    while num_add:
        ##         print("-> skipping", num_add, end=" ")
        _date2 = date2 + dt.timedelta(days=num_add)
        num_add = getNumFeb29s(date2, _date2)
        date2 = _date2

    return date2

def get_date_pairs_list(begin_date, final_date, *, time_step, time_between):
    date_pairs = []
    start_date = begin_date
    end_date = sumdate(start_date, time_between)
    while end_date <= final_date:
        date_pairs.append((start_date, end_date))
        start_date = sumdate(start_date, time_step)
        end_date = sumdate(start_date, time_between)
    return date_pairs