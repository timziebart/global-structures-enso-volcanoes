
import data_loader as dl

import datetime as dt
import numpy as np

ClimEventsYear = {
        "El Nino": { # always the beginning year
            "weak": [1951, 1952, 1953, 1958, 1968, 1969, 1976, 1977, 1979, 1977, 1979, 1994, 2004, 2006],
            "moderate": [1963, 1986, 1987, 1991, 2002, 2009],
            "strong": [1957, 1965, 1972],
            "very strong": [1982, 1997, 2015],
            },
        "La Nina": { # always the beginning year
            "weak": [1950, 1954, 1964, 1967, 1971, 1974, 1983, 1984, 1995, 2000, 2011],
            "moderate": [1955, 1970, 1998, 1999, 2007, 2010],
            "strong": [1973, 1975, 1988],
            },
        }

YearClimEvents = { year : (event, strength) for event in ClimEventsYear for strength in ClimEventsYear[event] for year in ClimEventsYear[event][strength] }

moderate_plus_strengths = ["moderate", "strong", "very strong"]
strong_plus_strengths = ["strong", "very strong"]

event_years_mark_en_cp = [1953, 1958, 1963, 1968, 1969, 1977, 1979, 1986, 1987, 1991, 1994, 2002, 2004, 2006, 2009]
event_years_mark_en_ep = [1957, 1965, 1972, 1976, 1982, 1997]
event_years_mark_ln_cp = [1954, 1955, 1967, 1971, 1974, 1975, 1984, 1995, 2000, 2001, 2011]
event_years_mark_ln_ep = [1964, 1970, 1973, 1988, 1998, 2007, 2010]

event_years = {}

event_years["all-years"] = list(range(dl.NCEP_NCAR.begin_year, dl.NCEP_NCAR.end_year + 1))

event_years["oni"] = sorted(YearClimEvents)

event_years["oni-moderate+"] = sorted(
    [year for year in YearClimEvents if YearClimEvents[year][1] in moderate_plus_strengths])

event_years["oni-strong+"] = sorted(
    [year for year in YearClimEvents if YearClimEvents[year][1] in strong_plus_strengths])

event_years["mark-EN-cp"] = event_years_mark_en_cp
event_years["mark-EN-ep"] = event_years_mark_en_ep
event_years["marc-LN-cp"] = event_years_mark_ln_cp
event_years["marc-LN-ep"] = event_years_mark_ln_ep


def simple_composite_dates(alldates, events=None):
    # assume allyears is sorted

    if events is None:
        events = list(event_years)

    event_dates = { ev:[] for ev in events }
    event_dates["other"] = []

    year_begin, year_end = alldates[0].year, alldates[-1].year
    for y in range(year_begin, year_end+1):
        christmas_index = np.searchsorted(alldates, dt.date(y, 12, 24))
        if christmas_index == len(alldates):
            break # stopping because end of array
        found_event = False
        for ev in events:
            if y in event_years[ev]:
                if ev != "all-years":
                    found_event = True
                event_dates[ev].append(alldates[christmas_index])
        if not found_event:
            event_dates["other"].append(alldates[christmas_index])

    return event_dates



