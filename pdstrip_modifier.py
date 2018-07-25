# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 12:11:00 2018

@author: Admaren02
"""

import numpy as np


def makeafullgeom(inputfilepath: str, outputfilepath: str):
    linecounter = 0
    ws = " " * 5
    with open(outputfilepath, "w") as _fout:
        with open(inputfilepath, "r") as _fin:
            while True:

                if linecounter == 0:
                    firstline = _fin.readline()
                    ns, *_ = firstline.split()
                    _fout.write("{0}   F  {1}\n".format(ns, _[-1]))
                else:
                    infoline = _fin.readline()

                    loc, np, gap = infoline.split()

                    npts = int(np)
                    _fout.write(
                        "{0:>11}{1:>16}{2:>3}\n".format(loc, (2 * npts) - 1, gap)
                    )

                    ctr = 0
                    xarr = []
                    yarr = []
                    while ctr < npts:
                        _line = _fin.readline().split()
                        for v in _line:
                            xarr.append(v)
                            ctr += 1
                    ctr = 0
                    while ctr < npts:
                        _line = _fin.readline().split()
                        for v in _line:
                            yarr.append(v)
                            ctr += 1

                    xarrnew = [float(v) * -1 for v in xarr[1:][::-1]] + xarr
                    yarrnew = yarr[1:][::-1] + yarr
                    _fout.write(ws.join([str(v) for v in xarrnew]) + "\n")
                    _fout.write(ws.join(yarrnew) + "\n")

                linecounter += 1
                if linecounter == int(ns) + 1:
                    break


def section_area(x=[3, 5, 12, 9, 5], y=[4, 11, 8, 5, 6]):
    a1 = 0

    for i in range(1, len(x) + 1):

        if i == len(x):

            a1 += x[i - 1] * y[0] - (x[0] * y[i - 1])
        else:

            a1 += x[i - 1] * y[i] - (x[i] * y[i - 1])

    return a1 / 2


def section_centroid(x=[0, 1, 2, -1, 0, -1, -2, -1], y=[2, 1, 0, 2, -2, -1, 0, 1]):

    a1 = 0
    cx = 0
    cy = 0

    for i in range(1, len(x) + 1):

        if i == len(x):

            a = (x[i - 1] * y[0]) - (x[0] * y[i - 1])
            a1 += a
            cx += (x[i - 1] + x[0]) * a
            cy += (y[i - 1] + y[0]) * a
        else:

            
            a = (x[i - 1] * y[i]) - (x[i] * y[i - 1])
            a1 += a
            cx += (x[i - 1] + x[i]) * a
            cy += (y[i - 1] + y[i]) * a

    area = a1 / 2.
    return cx / (6 * area), cy / (6 * area)


def output_array(inputfilepath: str):
    linecounter = 0
    sectionss = {}
    with open(inputfilepath, "r") as _fin:
        while True:
            # print(linecounter)
            if linecounter == 0:
                firstline = _fin.readline()
                ns, *_ = firstline.split()

            else:
                infoline = _fin.readline()

                loc, np, gap = infoline.split()

                npts = int(np)

                ctr = 0
                xarr = []
                yarr = []
                while ctr < npts:
                    _line = _fin.readline().split()
                    for v in _line:
                        xarr.append(float(v))
                        ctr += 1
                ctr = 0
                while ctr < npts:
                    _line = _fin.readline().split()
                    for v in _line:
                        yarr.append(float(v))
                        ctr += 1

                sectionss[float(loc)] = [xarr, yarr]
            linecounter += 1
            if linecounter == int(ns) + 1:
                break
    return sectionss


def pdstriphullproperties(geomfilepath: str, midshipindex=None) -> dict:
    """Generate information about geometry file
    
    Assumptions
    #########
    maximum section area is assumed as midship  other wise midship index to be
    explicitly provided
    
    """

    res = output_array(outputfilepath)

    hull_res = {}

    hull_res["sectioninfo"] = []

    for loc in sorted(res):
        sec_res = {}
        sec0x, sec0y = output_array(outputfilepath)[loc]

        bmin = min(sec0x)
        bmax = max(sec0x)

        dmax = max(sec0y)
        dmin = min(sec0y)

        sec_res["x"] = float(loc)
        sec_res["depth"] = dmax - dmin
        sec_res["breadth"] = bmax - bmin
        sec_res["area"] = section_area(sec0x, sec0y)
        sec_res["centroid"] = section_centroid(sec0x, sec0y)
        hull_res["sectioninfo"].append(sec_res)

    arealist = [sectiondict["area"] for sectiondict in hull_res["sectioninfo"]]
    xlist = [sectiondict["x"] for sectiondict in hull_res["sectioninfo"]]
    cxmom = [
        sectiondict["x"] * sectiondict["area"]
        for sectiondict in hull_res["sectioninfo"]
    ]
    cvmom = [
        sectiondict["centroid"][1] * sectiondict["area"]
        for sectiondict in hull_res["sectioninfo"]
    ]

    hull_res["displacement"] = np.trapz(arealist, xlist)
    hull_res["centroid_x"] = np.trapz(cxmom, xlist) / hull_res["displacement"]
    hull_res["centroid_y"] = np.trapz(cvmom) / hull_res["displacement"]

    hull_res["length"] = max(xlist) - min(xlist)
    hull_res["breadth"] = max(hull_res["sectioninfo"], key=lambda x: x["breadth"])[
        "breadth"
    ]

    hull_res["depth"] = max(hull_res["sectioninfo"], key=lambda x: x["depth"])["depth"]
    hull_res["cb"] = hull_res["displacement"] / (
        hull_res["length"] * hull_res["breadth"] * hull_res["depth"]
    )
    if not midshipindex:
        hull_res["cm"] = max(arealist) / (hull_res["breadth"] * hull_res["depth"])
    else:
        hull_res["cm"] = arealist[midshipindex] / (
            hull_res["breadth"] * hull_res["depth"]
        )

    return hull_res


if __name__ == "__main__":

    inputfilepath = r"C:\Users\Admaren02\Desktop\joel\pdstrip_example\geomet_bck.out"
    outputfilepath = r"C:\Users\Admaren02\Desktop\joel\pdstrip_example\jjj.out"
    resdic = pdstriphullproperties(outputfilepath)
    print(resdic)