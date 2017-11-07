#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK


import data_post_processor as dpp
import events
import fullrun as fr
import graph_analysis as ga
import icosahedral_grid as ico
import locations as locs


import argparse, argcomplete
import datetime as dt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import dates as mdates
from matplotlib import rc
import os

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


volcanos = {
    "agung": dict(
        name="Agung",
        shortname="Agung",
        location=locs.Point(**{"lat": -8.343, "lon": 115.507}),
        eruption=dt.date(1963, 2, 18),
        shift=locs.DegreeDelta(lat=20, lon=-5),
        delay=180,
    ),
    "st-helens": dict(
        name="St. Helens",
        shortname="helens",
        location=locs.Point(**{"lat": 46.191, "lon": -122.196}),
        eruption=dt.date(1980, 5, 18),
        shift=None,
        delay=90,
    ),
    "el-chichon": dict(
        name="El Chichon",
        shortname="Chichon",
        location=locs.Point(**{"lat": 17.359, "lon": -93.231}),
        eruption=dt.date(1982, 3, 15),  # roughly, wikipedia was not precise
        shift=locs.DegreeDelta(lat=15 - 17.359, lon=-120 + 93.231),
        # shift=locs.DegreeDelta(lat=15, lon=8),
        delay=90,
    ),
    "pinatubo": dict(
        name="Pinatubo",
        shortname="Pinatubo",
        location=locs.Point(**{"lat": 15.143, "lon": 120.350}),
        eruption=dt.date(1991, 6, 15),
        shift=locs.DegreeDelta(lat=-15, lon=10),
        delay=540,
    )
}

volcano_radius = 0.12  # unit: length on a unit-sphere
volcano_degrees = 60 * volcano_radius  # for plotting only
volcano_color = "blue"
volcano_color_shift = "red"

# for volc in volcanos.values():
#     print(f"{volc['name']:<12}: {volc['location']!s:<40} eruption {volc['eruption']!s}   shift {volc['shift']!s:<40} {volc['delay']:>3} days")




big2loc = locs.rectangle_from_infsup(dict(
    lat_inf=-30.0, lat_sup=10.0,
    lon_inf=180, lon_sup=-60,
))

def plot_regions_of_interest():
    fig, ax, m = data.create_base_map()

    patch34_style = dict(
        edgecolor="#de2d26",
        fill=False,
        linewidth=5
    )
    patch3_style = dict(
        # edgecolor="brown",
        fill=True,
        color="forestgreen",
        alpha=0.8,
        linewidth=0
    )
    patch4_style = dict(
        # edgecolor="red",
        fill=True,
        color="darkblue",
        alpha=0.8,
        linewidth=0
    )
    patchBig_style = dict(
        edgecolor="green",
        fill=False,
        linewidth=8,
        linestyle="dashed"
    )
    d = locs.DegreeDelta(lat=2)

    fontsize = 24

    data.draw_map_rectangle(ga.AreaCoordinates["nino-3-region"]["location"], m=m, style=patch3_style)
    data.text_on_map("3", m=m, point=ga.AreaCoordinates["nino-3-region"]["location"].lower_center - d, va="top", color=patch3_style["color"], fontsize=fontsize)

    data.draw_map_rectangle(ga.AreaCoordinates["nino-4-region"]["location"], m=m, style=patch4_style)
    data.text_on_map("4", m=m, point=ga.AreaCoordinates["nino-4-region"]["location"].lower_center - d, va="top", color=patch4_style["color"], fontsize=fontsize)

    data.draw_map_rectangle(ga.AreaCoordinates["nino-3-4-region"]["location"], m=m, style=patch34_style)
    data.text_on_map("3.4", m=m, point=ga.AreaCoordinates["nino-3-4-region"]["location"].lower_center - d, va="top", color=patch34_style["edgecolor"], fontsize=fontsize)

    # data.draw_map_rectangle(ga.AreaCoordinates["nino-big-region"]["location"], m=m, style=patchBig_style)
    # data.text_on_map("EN-big", m=m, point=ga.AreaCoordinates["nino-big-region"]["location"].upper_center + d)

    data.draw_map_rectangle(big2loc, m=m, style=patchBig_style)
    data.text_on_map("ENSO-big", m=m, point=big2loc.lower_left - d - d, va="top", color=patchBig_style["edgecolor"], fontsize=fontsize)


    for volc in volcanos:
        if volc == "st-helens":
            continue
        mask = locs.Circle(center=volcanos[volc]["location"], radius=volcano_radius).get_mask(data.grid_obj)
        data.draw_mask_on_map(mask, m=m, color=volcano_color)
        tag_point = locs.Point(**dict(volcanos[volc]["location"]))
        if volc == "el-chichon":
            tag_point.lon -= volcano_degrees + 1
            data.text_on_map(volcanos[volc]["name"], m=m, point=tag_point, color=volcano_color, ha="right", fontsize=fontsize)
        elif volc == "agung":
            tag_point.lon -= volcano_degrees - 2
            tag_point.lat -= volcano_degrees - 4
            data.text_on_map(volcanos[volc]["name"], m=m, point=tag_point, color=volcano_color, ha="right", va="top", fontsize=fontsize)
        else:
            tag_point.lon += volcano_degrees + 1
            data.text_on_map(volcanos[volc]["name"], m=m, point=tag_point, color=volcano_color, fontsize=fontsize)

    if args.save:
        fname = "regions-of-interest.pdf"
        print(f"saving {fname} ... ", end="", flush=True)
        fig.savefig(fname)
        print("done")


TIMESERIES_GLOBAL_META_DATA = {
    "global-transitivity" : ("a", r"$\mathcal{T}$"),
    "modularity-walktrap" : ("b", r"$\mathcal{Q}$"),
    "global-avg-link-length-field" : ("c", r"$\left<\left<d\right>\right>$")
}
def plot_global_timeseries():
    for field_type in data.field_dict:
        name="global-" + field_type
        data.create_timeseries(
            name,
            field_name=field_type,
            location=None,
        )

    for name in TIMESERIES_GLOBAL_META_DATA:
        if name not in data.timeseries:
            print()
            print(f"{name} NOT FOUND!")
            print()
            continue
        fig, ax = data.plot_timeseries(
            name,
            meta_data=TIMESERIES_GLOBAL_META_DATA[name],
            add_label=True
        )
        if args.save:
            fname = name+".pdf"
            print(f"saving {fname} ... ", end="", flush=True)
            fig.savefig(fname)
            print("done")


    for field_type in data.field_dict:
        name="global-" + field_type
        data.delete_timeseries(name)

TIMESERIES_LOCAL_ENSO_META_DATA = {
    "nino-3-4-region-degree-field" : ("a", r"$k_{\textrm{Ni{\~{n}}o3.4}}$"),
    "nino-3-region-degree-field" : ("b", r"$k_{\textrm{Ni{\~{n}}o3}}$"),
    "nino-4-region-degree-field" : ("c", r"$k_{\textrm{Ni{\~{n}}o4}}$"),
    "nino-3-4-region-avg-link-length-field" : ("d", r"$d_{\textrm{Ni{\~{n}}o3.4}}$"),
    "nino-3-region-avg-link-length-field" : ("e", r"$d_{\textrm{Ni{\~{n}}o3}}$"),
    "nino-4-region-avg-link-length-field" : ("f", r"$d_{\textrm{Ni{\~{n}}o4}}$"),
    # "nino-3-4-region-degree-field" : ("a", "ENSO 3.4\nAvg. Degree"),
    # "nino-3-region-degree-field" : ("b", "ENSO 3\nAvg. Degree"),
    # "nino-4-region-degree-field" : ("c", "ENSO 4\nAvg. Degree"),
    # "nino-3-4-region-avg-link-length-field" : ("d", "ENSO 3.4\nAvg. Link Length"),
    # "nino-3-region-avg-link-length-field" : ("e", "ENSO 3\nAvg. Link Length"),
    # "nino-4-region-avg-link-length-field" : ("f", "ENSO 4\nAvg. Link Length"),
}
def plot_local_enso_timeseries():
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    for field_type in data.field_dict:
        for region in ["nino-3-4-region", "nino-3-region", "nino-4-region"]:
            name = region + "-" + field_type
            data.create_timeseries(
                name,
                field_name=field_type,
                location=ga.AreaCoordinates[region]["location"],
            )

    for name in TIMESERIES_LOCAL_ENSO_META_DATA:
        fig, ax = data.plot_timeseries(
            name,
            meta_data=TIMESERIES_LOCAL_ENSO_META_DATA[name],
            add_label=True
        )
        if args.save:
            fname = name+".pdf"
            print(f"saving {fname} ... ", end="", flush=True)
            fig.savefig(fname)
            print("done")
        # break # only do once for testing



    for field_type in data.field_dict:
        for region in ["nino-3-4-region", "nino-3-region", "nino-4-region"]:
            name = region + "-" + field_type
            data.delete_timeseries(name)


TIMESERIES_VOLCANOES_META_DATA = {
    "pinatubo-degree-field" : ("a", "Pinatubo\nAvg. Degree"),
    "el-chichon-degree-field" : ("b", "El Chichon\nAvg. Degree"),
    "agung-degree-field" : ("c", "Mt. Agnung\nAvg. Degree"),
}
def plot_volcano_timeseries():
    for field in data.field_dict:
        for volc, loc in volcanos.items():
            name = volc + "-" + field
            shift_name = name + "-shift"
            data.create_timeseries(
                name,
                field_name=field,
                location=locs.Circle(center=volcanos[volc]["location"], radius=volcano_radius)
            )
            shift = volcanos[volc].get("shift", None)
            if shift is not None:
                data.create_timeseries(
                    shift_name,
                    field_name=field,
                    location=locs.Circle(center=volcanos[volc]["location"]+shift, radius=volcano_radius)
                )

    YLIMS = {
        "pinatubo" : (0, 500),
        "agung" : (0, 250)
    }

    field = "degree-field"
    next_label = "a"
    for volcanoname in ["pinatubo", "agung", "el-chichon"]:

        name = volcanoname + "-" + field
        shift_name = name + "-shift"
        date_eruption = volcanos[volcanoname]["eruption"]
        date_signature = date_eruption + dt.timedelta(days=volcanos[volcanoname]["delay"])
        shift = volcanos[volcanoname]["shift"]

        print(f"plotting {volcanoname}")

        mask = locs.Circle(
            center=volcanos[volcanoname]["location"],
            radius=volcano_radius
        ).get_mask(data.grid_obj)
        mask_shift = locs.Circle(
            center=volcanos[volcanoname]["location"] + shift,
            radius=volcano_radius
        ).get_mask(data.grid_obj)

        fig = plt.figure(volcanoname, figsize=(14, 4))
        # ax = fig.add_axes((0.10, 0.26, 0.895, 0.62))
        ax       = fig.add_axes((0.06, 0.62, 0.50, 0.31))
        ax_shift = fig.add_axes((0.06, 0.13, 0.50, 0.31))
        # ax_map = fig.add_axes((0.60, 0.125, 0.395, 0.75))
        if volcanoname == "pinatubo":
            ax_cb = fig.add_axes((0.6, 0.909, 0.39, 0.09))
            fig.text(0.586, 0.94, r"$k_i$", fontsize=22, ha="center", va="center")
            ax_map = fig.add_axes((0.60, 0.1, 0.395, 0.75))
        else:
            ax_map = fig.add_axes((0.60, 0.125, 0.395, 0.75))
        label       = "(%s)"%next_label
        label_shift = "(%s)"%chr(ord(next_label)+1)
        label_map   = "(%s)"%chr(ord(next_label)+2)
        next_label  = chr(ord(next_label)+3)
        fig.text(0, 0.95, label)
        fig.text(0, 0.45, label_shift)
        if volcanoname == "pinatubo":
            fig.text(0.58, 0.75, label_map)
        else:
            fig.text(0.58, 0.795, label_map)

        marker_eruption_style = dict(color = "purple", lw=2, zorder=0)
        marker_signature_style = dict(color="green", lw=2, zorder=1)

        yeardiff = 10
        begin_date = dt.date(day=date_eruption.day, month=date_eruption.month, year=date_eruption.year-yeardiff)
        end_date = dt.date(day=date_eruption.day, month=date_eruption.month, year=date_eruption.year+yeardiff)

        if volcanoname in YLIMS:
            ax.set_ylim(YLIMS[volcanoname])
        data.plot_timeseries(name, ax=ax, meta_data=("", ""))
        ylim = ax.get_ylim()
        ax.plot([mdates.date2num(date_eruption)]*2, ylim, **marker_eruption_style)
        ax.plot([mdates.date2num(date_signature)]*2, ylim, **marker_signature_style)
        ax.set_ylim(ylim)
        # title = volcanos[volcanoname]["name"] + "\n Avg. Degree"
        title = r"$k_{\textrm{" + volcanos[volcanoname]['shortname'] + r"}}$"
        fig.text(0.005, 0.75, f"{title}", ha="left", va="center", rotation=90, multialignment="center")
        ax.set_xlim((begin_date, end_date))

        if volcanoname in YLIMS:
            ax_shift.set_ylim(YLIMS[volcanoname])
        data.plot_timeseries(shift_name, ax=ax_shift, meta_data=("", ""))
        ylim = ax_shift.get_ylim()
        ax_shift.plot([mdates.date2num(date_eruption)]*2, ylim, **marker_eruption_style)
        ax_shift.plot([mdates.date2num(date_signature)]*2, ylim, **marker_signature_style)
        ax_shift.set_ylim(ylim)
        # line = list(ax.get_lines())[0]
        # line.set_zorder(3)
        # print(line, line.get_zorder())
        # title = volcanos[volcanoname]["name"] + " Shift\n Avg. Degree"
        title = r"$k^\prime_{\textrm{" + volcanos[volcanoname]['shortname'] + r"}}$"
        fig.text(0.005, 0.25, f"{title}", ha="left", va="center", rotation=90, multialignment="center")
        ax_shift.set_xlim((begin_date, end_date))

        _, _, m = data.create_base_map(ax=ax_map)
        _, _, _, mappable = data.plot_field(date_signature, field, m=m, set_title=False)
        data.draw_mask_on_map(mask, m=m, ms=2, color=volcano_color)
        data.draw_mask_on_map(mask_shift, m=m, color=volcano_color_shift, ms=2)

        if volcanoname == "pinatubo":
            cbar = plt.colorbar(mappable=mappable, cax=ax_cb, orientation="horizontal")
            cbar.formatter.set_powerlimits((0, 0))
            # off_txt = cbar.ax.xaxis.get_offset_text()
            # print(off_txt)
            # off_txt.set_position((1, 0))
            cbar.update_ticks()

        patchBig_style = dict(
            fill=True,
            color="white",
            linewidth=0
        )
        data.draw_map_rectangle(big2loc, m=m, style=patchBig_style)

        if args.save:
            # fname = name+".pdf"
            fname = name+".jpg"
            print(f"saving {fname} ... ", end="", flush=True)
            fig.savefig(fname, dpi=200)
            print("done")
        break # to do it only once for testing



    for field in data.field_dict:
        for volc, loc in volcanos.items():
            data.delete_timeseries(volc + "-" + field)


def plot_composites():

    use_events = [
        "mark-EN-ep",
        "mark-EN-cp",
        "marc-LN-ep",
        "marc-LN-cp",
        # "other",
    ]
    other_event = "other"
    use_events_with_other = use_events + [other_event]

    use_fields = ["degree-field", "avg-link-length-field"]

    count = 0
    labela = "a"

    figwidth = 6
    figheight = 3
    cb_figheight = 1.2

    figsize = (figwidth, figheight)
    axsize = (0.06, 0, 0.935, 1)

    cb_figsize = (figwidth, cb_figheight)
    cb_axsize = (0.06, 0.4, 0.88, 0.26)

    event_dates = events.simple_composite_dates(data.dates, events=use_events)

    plotted_cb_already = {field:False for field in use_fields}

    for ev in use_events_with_other:
        dates = event_dates[ev]
        # for field in data.field_dict:
        for field in use_fields:
            composite_name = f"{ev}-{field}"
            print(f"creating composite {composite_name} ... ", flush=True, end = "")
            data.create_composite(composite_name, field=field, dates=dates)
            print("plotting ... ", flush=True, end = "")

            fig = plt.figure(composite_name, figsize=figsize)
            ax = fig.add_axes(axsize)

            # ax = fig.add_subplot(event_count, field_count, count)
            # count += 1

            _, _, m = data.create_base_map(ax=ax)
            _, _, _, mappable = data.plot_composite(composite_name, m=m, add_title=False)

            label = chr(ord(labela) + count)
            count += 1
            fig.text(0, 0.93, f"({label})")
            del label

            print("deleting the composite again ... ", flush=True, end = "")

            data.delete_composite(composite_name)

            print("done")

            if not plotted_cb_already[field]:
                plotted_cb_already[field] = True
                print(f"plotting colorbar for {field} ... ", flush=True, end = "")
                cb_fig = plt.figure(f"colorbar-{field}", figsize=cb_figsize)
                # ax = fig.add_subplot(111)
                cb_ax = cb_fig.add_axes(cb_axsize)
                cbar = plt.colorbar(mappable=mappable, cax=cb_ax, orientation="horizontal")
                cbar.formatter.set_powerlimits((0, 0))
                # off_txt = cbar.ax.xaxis.get_offset_text()
                # print(off_txt)
                # off_txt.set_position((1, 0))
                cbar.update_ticks()
                if field == "degree-field":
                    cb_label = r"Degree $k_i$"
                elif field == "avg-link-length-field":
                    cb_label = r"Average Link Distance $d_i$"
                cb_fig.text(0.5, 0.7, cb_label, ha="center", fontsize=22)
                if args.save:
                    fname = f"colorbar-{field}.jpg"
                    print(f"saving {fname} ... ", end="", flush=True)
                    cb_fig.savefig(fname, dpi=200)
                print("done")

            if args.save:
                fname = f"{composite_name}.jpg"
                print(f"saving {fname} ... ", end="", flush=True)
                fig.savefig(fname, dpi=200)
                print("done")
        #     break # for testing just run it once
        # break # for testing just run it once

    # for comp in sorted(data.composites):
    #     data.plot_composite(comp)


    # fig.add_subplot()

def cmp_modulariy():
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # data.timeseries.plot(ax=ax)

    title = "modularity-comparison"

    modularities = [
        'modularity-fast-greedy',
        'modularity-infomap',
        'modularity-label-propagation',
        'modularity-leading-eigenvector',
    ]

    l = len("modularity-")
    modularities_labels = [ mod[l:].replace("-", " ") for mod in modularities]

    walktrap = 'modularity-walktrap'
    walktraplabel = "walktrap"

    fig = plt.figure(title, figsize=(14, 5))
    ax = fig.add_axes((0.06, 0.16, 0.935, 0.81))
    for modularity, label in zip(modularities, modularities_labels):
        print(f"plotting {modularity}")
        alpha = 0.5
        if "walktrap" in modularity:
            print("FOUND WALKTRAP")
            alpha = 1
        # data.timeseries[modularity].plot(ax=ax, legend=True, alpha = alpha)
        fig, ax = data.plot_timeseries(
            modularity,
            ax=ax,
            show_ENSO=False,
            alpha=0.5,
            color=None,
            add_title=False,
            legend=True,
            dropna=True,
            label=label
        )

    data.plot_timeseries(
        walktrap,
        ax=ax,
        show_ENSO=False,
        add_title=False,
        dropna=True,
        label=walktraplabel,
        xlabel="time"
    )
    ax.legend(loc="lower right", fontsize="small")

    ax.set_ylabel("modularity")
    ax.set_ylim(0.5, 0.86)

    if args.save:
        fname = f"{title}.pdf"
        print(f"saving {fname} ... ", end="", flush=True)
        fig.savefig(fname)
        print("done")

def plot_enso_colorbar():
    title = "enso-colorbar"
    # fig = plt.figure(title, figsize=(16, 1), tight_layout=False)
    # ax = fig.add_axes((0.07, 0.4, 0.86, 0.59))
    # fig.text(0.05, 0.68, "EN", fontsize=28, ha="right", va="center")
    # fig.text(0.95, 0.68, "LN", fontsize=28, ha="left", va="center")

    fig = plt.figure(title, figsize=(16, 1.5), tight_layout=False)

    # ax = fig.add_axes((0.10, 0.59, 0.83, 0.4))
    # fig.text(0.08, 0.68, "EN", fontsize=28, ha="right", va="center")
    # fig.text(0.95, 0.68, "LN", fontsize=28, ha="left", va="center")
    # anomaly_height = 0.1
    # anomaly_xdelta = 0.83 / 9
    # anomaly_xoffset = 0.1 + 0.83/18
    # fig.text(0, anomaly_height, "ONI", fontsize=24)

    ax_width = 0.86
    ax_xoffset = 0.07

    ax = fig.add_axes((ax_xoffset, 0.59, ax_width, 0.4))
    fig.text(0.05, 0.68, "EN", fontsize=28, ha="right", va="center")
    fig.text(0.95, 0.68, "LN", fontsize=28, ha="left", va="center")
    anomaly_height = 0.1
    anomaly_classes_num = 9
    anomaly_xdelta = ax_width / anomaly_classes_num
    anomaly_xoffset = ax_xoffset + ax_width / (2*anomaly_classes_num)
    fig.text(0.05, anomaly_height, "ONI", fontsize=24, ha="right")
    fig.text(anomaly_xoffset + 0*anomaly_xdelta, anomaly_height, r"$>2.0^\circ$", fontsize=24, ha="center")
    fig.text(anomaly_xoffset + 1*anomaly_xdelta, anomaly_height, r"$>1.5^\circ$", fontsize=24, ha="center")
    fig.text(anomaly_xoffset + 2*anomaly_xdelta, anomaly_height, r"$>1.0^\circ$", fontsize=24, ha="center")
    fig.text(anomaly_xoffset + 3*anomaly_xdelta, anomaly_height, r"$>0.5^\circ$", fontsize=24, ha="center")

    fig.text(anomaly_xoffset + 8*anomaly_xdelta, anomaly_height, r"$<2.0^\circ$", fontsize=24, ha="center")
    fig.text(anomaly_xoffset + 7*anomaly_xdelta, anomaly_height, r"$<1.5^\circ$", fontsize=24, ha="center")
    fig.text(anomaly_xoffset + 6*anomaly_xdelta, anomaly_height, r"$<1.0^\circ$", fontsize=24, ha="center")
    fig.text(anomaly_xoffset + 5*anomaly_xdelta, anomaly_height, r"$<0.5^\circ$", fontsize=24, ha="center")



    _en_cmap = mpl.colors.ListedColormap(["#ee4444", "#2266aa"])
    data_picks = np.linspace(-4.5, 4.5, 19, endpoint=True)
    # print(data_picks)
    color_picks = _en_cmap(data_picks)
    # print()
    # print(color_picks)
    # alphas = np.repeat(np.linspace(0, 1, 5), 2)
    alphas = np.repeat([0., 0.2, 0.4, 0.6, 0.8], 2)
    color_picks[8:-1, -1] = alphas
    color_picks[-1, -1] = alphas[-1]
    color_picks[0:10, -1] = alphas[::-1]
    # print()
    # print(color_picks)
    en_cmap = mpl.colors.ListedColormap(color_picks)
    # assert False
    bounds = data_picks
    norm = mpl.colors.BoundaryNorm(bounds, en_cmap.N)
    en_cb = mpl.colorbar.ColorbarBase(ax,
                                      cmap=en_cmap,
                                      norm=norm, boundaries=bounds,
                                      orientation="horizontal")
    # en_cb.ax.majorticks_on()
    # en_cb.minorticks
    en_strengths = ["weak", "moderate", "strong", "very strong"]
    ticklabels = en_strengths[::-1] + [""] + en_strengths
    # ticklabels.insert(1, "")
    # print(ticklabels)
    # ticklabels = np.array([""]*19)
    # ticklabels[:8:2] = en_strengths
    # ticklabels[-2:10:-2] = en_strengths
    en_cb.ax.set_xticklabels(ticklabels, fontsize=24)

    if args.save:
        fname = f"{title}.pdf"
        print(f"saving {fname} ... ", end="", flush=True)
        fig.savefig(fname)
        print("done")


if __name__ == "__main__":

    pix_modes = [
        "ENSO-global",
        "ENSO-local",
        "volcanoes",
        "grid",
        "regions-of-interest",
        "composites",
        "cmp-modularity",
        "enso-colorbar",
    ]

    parser = argparse.ArgumentParser(description="creating the pix for the paper")

    parser.add_argument(
        "input_file", metavar="input-file"
    )

    parser.add_argument(
        "grid", default=fr.RunGrids.icosahedral.name, choices=fr.grid_choices,
        help="choose which grid has been used for the computation"
    )

    parser.add_argument(
        "modes", metavar="mode", nargs="+", choices=pix_modes,
        help="which pix are to be shown, choose from: " + ", ".join(pix_modes)
    )

    parser.add_argument(
        "--no-show", action="store_false", dest="show",
        help="do not show the plots"
    )
    parser.add_argument(
        "--save", action="store_true",
        help="save the pictures"
    )

    argcomplete.autocomplete(parser)

    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        parser.error(f"{args.input_file} is not a file")

    filename = args.input_file

    grid_obj = fr.RunGrids[args.grid].value() # choose the necessary grid from RunGrids (as given in the command line arguments

    print(f"loading input file '{filename}' ... ", flush=True, end="")
    data = dpp.DataPostProcessor(filename, grid_obj=grid_obj)
    print("done")
    if {"teleconnectivity-field", "degree-field"}.issubset(data.field_dict):
        data.field_dict["avg-link-length-field"] = data.field_dict["teleconnectivity-field"] / data.field_dict["degree-field"]
    if "elnino-deg" in data.timeseries:
        data.delete_timeseries("elnino-deg") # already computed during the run, but only for checking
    if "elnino-tele" in data.timeseries:
        data.delete_timeseries("elnino-tele") # already computed during the run, but only for checking

    print("loaded time series:", sorted(data.timeseries))
    print("loaded fields:", sorted(data.field_dict))

    if "grid" in args.modes:
        _, _, m = data.create_base_map()
        mask = locs.WholeWorld().get_mask(data.grid_obj)
        data.draw_mask_on_map(mask, m=m)
        if args.show:
            plt.show()
        plt.close("all")

    if "regions-of-interest" in args.modes:
        plot_regions_of_interest()
        if args.show:
            plt.show()
        plt.close("all")

    if "ENSO-global" in args.modes:
        assert grid_obj.__class__ is ico.IcosahedralGrid
        plot_global_timeseries()
        if args.show:
            plt.show()
        plt.close("all")

    if "ENSO-local" in args.modes:
        assert grid_obj.__class__ is ico.IcosahedralGrid
        plot_local_enso_timeseries()
        if args.show:
            plt.show()
        plt.close("all")

    if "volcanoes" in args.modes:
        assert grid_obj.__class__ is ico.IcosahedralGrid_PartRemoved
        plot_volcano_timeseries()
        if args.show:
            plt.show()
        plt.close("all")

    if "composites" in args.modes:
        assert grid_obj.__class__ is ico.IcosahedralGrid
        plot_composites()
        if args.show:
            plt.show()
        plt.close("all")

    if "cmp-modularity" in args.modes:
        cmp_modulariy()
        if args.show:
            plt.show()
        plt.close("all")

    if "enso-colorbar" in args.modes:
        plot_enso_colorbar()
        if args.show:
            plt.show()
        plt.close("all")





























