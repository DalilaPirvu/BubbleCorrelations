import os, sys
import numpy as np
import math
import scipy as scp
import scipy.optimize as sco
import scipy.signal as scs
import scipy.interpolate as intp

from itertools import groupby, cycle
import random
from functools import partial
import scipy.ndimage
from scipy.ndimage import gaussian_filter, gaussian_filter1d

import statistics as stat
from scipy.stats import moment
from collections import OrderedDict

def find_peak_positions(slice, threshold):
    """ Finds x coordinate of peaks in masked field with mask applied at threshold. """
    peak_coord = scs.find_peaks(slice, height = threshold)[0]
    # and since we have periodic boundaries:
    if slice[0]>=threshold or slice[-1]>=threshold:
        aa, bb = slice[-2:]
        cc, dd = slice[:2]
        if bb > cc and bb > aa:
            peak_coord = np.concatenate((peak_coord,[len(slice)-1]))
        elif bb < cc and cc > dd:
            peak_coord = np.concatenate((peak_coord,[0]))
    return peak_coord

def identify_bubble_sites2(simulation, tMinCutoff, tMaxCutoff, tFutureCheck, lightc, crit_thresh, bubble_thr, collision_thr):
    nT, nN = np.shape(simulation)
    vals, times = [], []

    for tt in range(0, tMaxCutoff, tFutureCheck):
        # find new peaks in field and check if they are good bubble candidates
        fldslice = simulation[tt]
        peaks_new = find_peak_positions(fldslice, bubble_thr)
        peaks_new = set(peaks_new) - set(vals)

        for pk in peaks_new:
            # step 2:
            # discard element if future lightcone is not in tv
            overshoot = False

            tmax = min(tt+tFutureCheck+1, nT)
            for tparse in range(tt, tmax):
                sizeCone = int(np.abs(tparse+1-tt)*lightc)
                coordsCone = np.arange(pk - sizeCone, pk + sizeCone + 1)
                conelist = simulation[tparse, coordsCone%nN]
                if np.any(conelist < bubble_thr):# or np.any(conelist > collision_thr):
                    overshoot = True
                    break

            if not overshoot:
                # Step 3:
                # check bubble hasn't been detected already
                bubbleCoords = find_bubbles_at_t(fldslice, bubble_thr, nN)
                for bubble in bubbleCoords:
                    if pk in bubble and not common_member(vals, bubble):
                        vals.append(pk)
                        times.append(tt)
                        break
    return [(ti, va) for ti, va in zip(times, vals) if ti>tMinCutoff]

def identify_bubble_sites(simulation, tMinCutoff, tMaxCutoff, tFutureCheck, lightc, crit_thresh, bubble_thr, collision_thr):
    nT, nN = np.shape(simulation)
    vals, times = [], []

    for tt in range(tMaxCutoff):
        # find new peaks in field and check if they are good bubble candidates
        fldslice = simulation[tt]
        peaks_new = find_peak_positions(fldslice, bubble_thr)
        peaks_new = set(peaks_new) - set(vals)

        for pk in peaks_new:
            # step 1:
            # ensure new peaks are not in the light cone of the previously identified bubble nucleation sites
            target_bounds = []
            for tval, xval in zip(times, vals):
                sizeCone = int(np.abs(tt+1-tval)*lightc)
                boundCoords = np.arange(xval - sizeCone, xval + sizeCone + 1)
                target_bounds += (boundCoords%nN).tolist()
            
            if pk not in target_bounds:
                # step 2:
                # discard element if future lightcone is not in tv
                overshoot = False

                tmax = min(tt+tFutureCheck+1, nT)
                for tparse in range(tt, tmax):
                    sizeCone = int(np.abs(tparse+1-tt)*lightc)
                    coordsCone = np.arange(pk - sizeCone, pk + sizeCone + 1)
                    conelist = simulation[tparse, coordsCone%nN]
                    if np.any(conelist < bubble_thr) or np.any(conelist > collision_thr):
                        overshoot = True
                        break

                if not overshoot:
                    # Step 3:
                    # check bubble hasn't been detected already
                    bubbleCoords = find_bubbles_at_t(fldslice, bubble_thr, nN)
                    for bubble in bubbleCoords:
                        if pk in bubble and not common_member(vals, bubble):
                            vals.append(pk)
                            times.append(tt)
                            break
    return [(ti, va) for ti, va in zip(times, vals) if ti>tMinCutoff]

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    if (a_set & b_set):
        return True
    else:
        return False

def extract_data(nL, path_sim):
    data = np.genfromtxt(path_sim)
    nNnT, nC = np.shape(data)
    return np.asarray([np.reshape(data[:,cc], (nNnT//nL, nL)) for cc in range(nC)])

def get_realisation(nL, sim, path_sim, outcome, phieq):
    data = extract_data(nL, path_sim)

    if outcome == 1:
        data[0] = 2.*phieq - data[0]
        data[1] = - data[1]
    return np.asarray(data), sim

def bubble_counts_at_equal_time(bubble, thresh):
    return np.count_nonzero(bubble > thresh, axis=1)

def reflect_against_equil(bubble, phi_init):
    return abs(bubble-phi_init) + phi_init

def find_t_max_width(bubble, light_cone, phi_init, crit_thresh, crit_rad, t0):
    T, N = np.shape(bubble)
    bubble = bubble[t0:]
    refl_bubble = reflect_against_equil(bubble, phi_init)

    bubble_counts = bubble_counts_at_equal_time(refl_bubble, crit_thresh)
    bubble_diffs = bubble_counts[1:] - bubble_counts[:-1]

    tmax = np.argwhere(bubble_diffs[::-1] >= light_cone*2).flatten()
    out = next((ii for ii, jj in zip(tmax[:-1], tmax[1:]) if jj-ii == 1), 1)
    return T-out

def find_bubbles_at_t(field_slice, crit_thresh, nN):
    vals = np.arange(1, nN+1)*(field_slice>crit_thresh)

    first_zero = next((i for i, x in enumerate(vals) if x == 0.), 0)
    vals = np.roll(vals, -first_zero)

    bubbleCoords = [np.asarray(list(g))-1 for k, g in groupby(vals, lambda x: x != 0) if k]
    return bubbleCoords

def round_to_n(x, n):
    return x if x == 0 else round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

