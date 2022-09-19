"""
This file contains the many helper functions used in IVL_search
"""
import numpy as np
import efel


def evaluator(t_vec, v_vec):
    """
    Take a time vector and voltage vector and then does score calculation as in metric test run
    both vectors has to be numpy arrays. The two vectors should have the same length.
    The time vector is in ms, and the voltage vector in mV.
    """
    trace = {'V': v_vec, 'T': t_vec, 'stim_start': [0], 'stim_end': [t_vec[-1]]}
    ef_list = ["AP_begin_indices", "AP_end_indices", "ISI_CV", "peak_voltage", "peak_indices"]

    efel.api.setThreshold(-20)
    efel.api.setDerivativeThreshold(1)

    features = efel.getFeatureValues([trace], ef_list)[0]
    spike_starts = features['AP_begin_indices']
    spike_ends = features['AP_end_indices']
    isicv = features["ISI_CV"][0]
    peakvm = features['peak_voltage']
    f = len(features['peak_indices']) / t_vec[-1] * 1000

    if len(peakvm) == len(v_vec[spike_starts]):
        ap_amps = peakvm - v_vec[spike_starts]
    else:
        pass

    v_spikeless = spike_ridder(v_vec, spike_starts, spike_ends)

    v_mean = np.mean(v_spikeless)
    v_variance = np.sum((v_spikeless - v_mean) ** 2) / len(v_spikeless)
    peakvm_avg = np.mean(ap_amps)

    print(int(v_mean > -70.588))
    print(int(v_variance > 2.2))
    print(int(isicv > 0.8))
    print(int(3 < f < 25))
    print(int(peakvm_avg < 40))

    score = (int(v_mean > -70.588)) + (int(v_variance > 2.2)) + (int(isicv > 0.8)) + (int(3 < f < 25)) - 5 * (
        int(peakvm_avg < 40))

    return score


def evaluator2(t_vec, v_vec):
    """
    Compare to evaluator, used threshold for ap detection instead of derivative. Derivative has shown to be confusing
    sub-threshold fluctuation as spiking event.
    """
    trace = {'V': v_vec, 'T': t_vec, 'stim_start': [0], 'stim_end': [t_vec[-1]]}
    ef_list = ["ISI_CV", "peak_voltage", "peak_indices"]

    efel.api.setThreshold(-20)
    efel.api.setDerivativeThreshold(1)

    features = efel.getFeatureValues([trace], ef_list)[0]
    if features["ISI_CV"] is None:
        isicv = 0
    else:
        isicv = features["ISI_CV"][0]
    peakvm = features['peak_voltage']
    spike_rate = len(features['peak_indices'])/t_vec[-1] * 1000

    v_spikeless = v_vec[v_vec < -20]

    v_mean = np.mean(v_spikeless)
    v_variance = np.sum((v_spikeless - v_mean) ** 2) / len(v_spikeless)

    if peakvm is None:
        peakvm_avg = 0
    else:
        peakvm_avg = np.mean(peakvm - v_mean)

    # print(int(v_mean > -70.588))
    # print(int(v_variance > 2.2))
    # print(int(isicv > 0.8))
    # print(int(3 < spike_rate < 25))
    # print(int(peakvm_avg < 40))

    score = (int(v_mean > -70.588)) + (int(v_variance > 2.2)) + (int(isicv > 0.8)) + (int(3 < spike_rate < 25)) - 5 * (
        int(peakvm_avg < 40))

    return score


def evaluator3(t_vec, v_vec):
    """
    Compared to evaluator2, added cutting around the -20 mV threshold to get less of the spike. The average above
    -20 mv part is around 0.73 ms, and AP is around 2 ms, so we opt to cut 0.8 ms from either side of the spike.
    """
    trace = {'V': v_vec, 'T': t_vec, 'stim_start': [0], 'stim_end': [t_vec[-1]]}
    ef_list = ["ISI_CV", "peak_voltage", "peak_indices"]

    efel.api.setThreshold(-30)
    efel.api.setDerivativeThreshold(1)

    features = efel.getFeatureValues([trace], ef_list)[0]
    if features["ISI_CV"] is None:
        isicv = 0
    else:
        isicv = features["ISI_CV"][0]
    peakvm = features['peak_voltage']
    spike_rate = len(features['peak_indices']) / t_vec[-1] * 1000

    v_spikeless = spike_ridder3(v_vec, features['peak_indices'], 70)

    v_mean = np.mean(v_spikeless)
    v_variance = np.var(v_spikeless)

    if peakvm is None:
        peakvm_avg = 0
    else:
        peakvm_avg = np.mean(peakvm - v_mean)

    # print(int(v_mean > -70.588))
    # print(int(v_variance > 2.2))
    # print(int(isicv > 0.8))
    # print(int(3 < spike_rate < 25))
    # print(int(peakvm_avg < 40))

    score = (int(v_mean > -70.588)) + (int(v_variance > 2.2)) + (int(isicv > 0.8)) + (int(3 < spike_rate < 25)) - 5 * (
        int(peakvm_avg < 40))

    return score


def evaluator4(t_vec, v_vec):
    """
    Evaluator3 but with smaller range for rate, between 11 and 16 hz
    """
    trace = {'V': v_vec, 'T': t_vec, 'stim_start': [0], 'stim_end': [t_vec[-1]]}
    ef_list = ["ISI_CV", "peak_voltage", "peak_indices"]

    efel.api.setThreshold(-20)
    efel.api.setDerivativeThreshold(1)

    features = efel.getFeatureValues([trace], ef_list)[0]
    if features["ISI_CV"] is None:
        isicv = 0
    else:
        isicv = features["ISI_CV"][0]
    peakvm = features['peak_voltage']
    spike_rate = len(features['peak_indices']) / t_vec[-1] * 1000

    v_spikeless = spike_ridder3(v_vec, features['peak_indices'], 70)

    v_mean = np.mean(v_spikeless)
    v_variance = np.sum((v_spikeless - v_mean) ** 2) / len(v_spikeless)

    if peakvm is None:
        peakvm_avg = 0
    else:
        peakvm_avg = np.mean(peakvm - v_mean)

    # print(int(v_mean > -70.588))
    # print(int(v_variance > 2.2))
    # print(int(isicv > 0.8))
    # print(int(3 < spike_rate < 25))
    # print(int(peakvm_avg < 40))

    score = (int(v_mean > -70.588)) + (int(v_variance > 2.2)) + (int(isicv > 0.8)) + (int(11 < spike_rate < 16)) - 5 * (
        int(peakvm_avg < 40))

    return score


def spike_ridder(v_vec, start_times, end_times):
    for i in range(len(start_times)):
        v_vec = np.hstack((v_vec[:start_times[i]], v_vec[end_times[i] + 1:]))
        off_set = end_times[i] - start_times[i] + 1
        start_times -= off_set
        end_times -= off_set
    return v_vec


def spike_ridder2(v_vec, cut_away):
    """Using -20 mV as threshold, remove everything above it, and 0.8 ms around the point that cross the threshold"""
    v_vec = v_vec + 0   # creates a copy

    v_vec[v_vec > -20] = 0
    i = 1
    while i < (len(v_vec) - 1):
        if v_vec[i] != 0 and v_vec[i + 1] == 0:
            if i < cut_away:
                v_vec[:i + 1] = 0
            else:
                v_vec[i - cut_away:i + 1] = 0
            i += 1
        elif v_vec[i] == 0 and v_vec[i + 1] != 0:
            if i + cut_away >= len(v_vec):
                v_vec[i:] = 0
            else:
                v_vec[i + 1:i + cut_away] = 0
            i += 9

    return v_vec


def spike_ridder3(v_vec, spike_indices, cut_away):
    """
    From each spike index, remove values in cut_away distance on the left and right.
    """

    length = len(v_vec)
    v_vec = v_vec + 0   # creates a copy

    for spike_time in spike_indices:
        if spike_time < cut_away:
            v_vec[:spike_time + cut_away+1] = 20
        elif spike_time + cut_away >= length:
            v_vec[spike_time - (cut_away+1):] = 20
        else:
            v_vec[spike_time - cut_away: spike_time + cut_away + 1] = 20

    return v_vec[v_vec != 20]

