
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt


ONI_DATAFILE= "detrend.nino34.ascii.txt"

def _load_oni():
    with open(ONI_DATAFILE, "r") as f:
        data = f.readlines()
    data = data[1:] # first line is the header
    oni_dict = {}
##     print data[0].split()
    for dat in data:
        year, month, _, _, oni = dat.split()
        year = int(year)
        month = int(month)
        oni = float(oni)
        if not year in oni_dict:
            oni_dict[year] = {}
        oni_dict[year][month] = oni
    return oni_dict


def getdate(year, month):
    return  dt.date(year, month, 15)

    
def get_oni(dates = None, ret_type = "dict"):
    onidict = _load_oni()

    if dates is None:
        if ret_type == "dict":
            return onidict
        elif ret_type == "lists":
            ds = [dt.date(year, month, 1) for year in onidict for month in onidict[year] ]
            onis = [onidict[d.year][d.month] for d in ds]
            return ds, onis
        else:
            raise TypeError("unknown ret_type: %s"%repr(ret_type))

    # else
    valid_dates = [dat for dat in dates if dat.year in onidict and dat.month in onidict[dat.year]] 
    return valid_dates, [onidict[dat.year][dat.month] for dat in valid_dates]

def plot_oni(ax = "new", dates = None, with_labels = False, style = {}, style_threshold = {}, with_threshold = False):
    newfig = ( ax == "new")
    if ax == "new":
        ax = plt.figure(figsize = (16, 4)).add_subplot(111)
    elif ax == "current" or ax is None:
        ax = plt.gca()
    fig = ax.get_figure()


    if with_labels:
        ax.set_xlabel("date")
        ax.set_ylabel("oni")
    
    assert isinstance(ax, plt.Axes)

    if dates is None:
        dates, onidata = get_oni(ret_type = "lists")
    else:
        min_date, max_date = min(dates), max(dates)
        dates, onidata = get_oni(dates = dates)

    ax.plot(dates, onidata, **style)
    if with_threshold:
        mindat, maxdat = min(dates), max(dates)
        for thr in np.linspace(-2.0, 2.0, 9, True):
            if thr == 0.0: continue # jump the middle line
            ax.plot([mindat, maxdat], [thr, thr], **style_threshold)

    
    

    if newfig:
        fig.tight_layout()
    return ax




## ONI_DATAFILE= "oni.data"

if __name__ == "__main__":
##     avOpts = ["full-redo", "save"]
    avOpts = ["save"]
    avKeyOpts = {}
    opts = uf.OptDict(avOpts = avOpts, avKeyOpts  = avKeyOpts)
##     if opts.isset("full-redo"):
##         full_redo_oni()
    oniax = plot_oni(with_labels = True, style = {"color" : "green"}, with_threshold = True, style_threshold = {"linestyle":":", "color":"black", "alpha":0.4})
    cn.plotClimEvents(oniax)
    if opts.isset("save"):
        filename = "oni.svg"
        uf.savePlot(filename)
    plt.show()




