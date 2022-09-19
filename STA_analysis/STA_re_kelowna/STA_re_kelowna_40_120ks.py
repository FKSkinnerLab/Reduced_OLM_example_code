import numpy as np
from neuron import h, gui
import record_1comp as r1
import efel
import matplotlib.pyplot as plt
import multiprocessing as mp
from currents_visualization import *
import itertools

# Instantiate Model
h.load_file("init_1comp.hoc")
h.cvode_active(0)
h.dt = 0.1
h.steps_per_ms = 10
DNQX_Hold = 0.004
TTX_Hold = -0.0290514

tau_e = 12.7
tau_i = 12.05

def one_state_sta_pert(all_paras, lst, runtime, pretime, fi, set_num, seed_seg, numpts, shared_dic):
    """run seed in seed sag. Update sta counts in sta_count_lst"""
    OU_in = h.Gfluct2(h.soma(0.5))
    OU_in.E_i = -87.1
    paras = all_paras[lst[set_num]]
    OU_in.g_e0 = paras[2]
    OU_in.g_i0 = paras[3]
    OU_in.std_e = paras[4]
    OU_in.std_i = paras[5]
    OU_in.tau_e = tau_e
    OU_in.tau_i = tau_i

    h.ic_hold.delay = 0
    h.ic_hold.dur = runtime
    h.tstop = runtime

    # creating synapse
    stim = h.NetStim()  # average number of spikes. convert to ms, then s
    stim.noise = 0  # deterministic
    stim.start = 0

    syn = h.Exp2Syn(h.soma(0.5))
    syn.tau1 = 1.6  # ms rise time
    syn.tau2 = 12.0  # ms decay time
    syn.e = -87.1  # reversal potential

    netcon = h.NetCon(stim, syn)  # threshold is irrelevant with event-based source
    netcon.weight[0] = 0.0054

    stim.interval = 1000 / fi  # interval in ms, but fi are in Hz
    stim.number = fi * (runtime / 1000)

    ef_list = ["peak_indices"]
    efel.api.setThreshold(-30)
    efel.api.setDerivativeThreshold(1)
    t_vec = np.array((range(int(runtime / h.dt) + 1))) * h.dt
    vecs = r1.set_up_full_recording()
    pre_points = int(pretime / h.dt)
    count = 0
    final = []
    i = 0
    while i < seed_seg.shape[0] and count < numpts:
        OU_in.new_seed(seed_seg[i])
        h.run()
        v_vec = np.array(vecs[0])
        all_vecs = np.array([np.array(arr) for arr in vecs])
        trace = {'V': v_vec, 'T': t_vec, 'stim_start': [0], 'stim_end': [t_vec[-1]]}

        spike_idx = efel.getFeatureValues([trace], ef_list)[0]['peak_indices']
        for j in range(len(spike_idx)):
            curr_idx = spike_idx[j]
            slice_start = curr_idx - pre_points
            if j == 0:
                pre_surpass = slice_start < 0
            else:
                pre_surpass = spike_idx[j - 1] > slice_start
            if not pre_surpass:
                interspike_iv = all_vecs[:, slice_start:curr_idx]
                final.append(interspike_iv)
                count += 1
        i += 1

    final = np.array(final)
    shared_dic[(set_num, fi)] = final


def multi_seed_run(all_paras, lst, runtime, pretime, fr_array, final_stfr, seeds, savepath, numpts):
    shared_dic = mp.Manager().dict()
    proc_list = []
    for i in range(len(fr_array)):
        for j in range(10):
            fr = fr_array[i]
            seed_seg = seeds[np.nonzero(final_stfr[:, j] == fr)[0]]
            p = mp.Process(target=one_state_sta_pert, args=(all_paras, lst, runtime, pretime, fr,
                                                            j, seed_seg, numpts, shared_dic))
            p.start()
            proc_list.append(p)
    for p in proc_list:
        p.join()
    all_sta = []
    for set_num in range(10):
        local_sta = []
        for fr in fr_array:
            sta = shared_dic[(set_num, fr)]     # shape of 50, 7, 2000
            local_sta.append(sta)
        local_sta = np.array(local_sta)
        all_sta.append(local_sta)
    all_sta = np.array(all_sta)
    np.save(savepath, all_sta)


def paras_sorter(para_path):
    para = np.load(para_path)
    values = []
    dtype = [('minr', float), ('maxr', float), ('ge', float), ('gi', float), ('stde', float), ('stdi', float)]
    for i in range(para.shape[0]):
        values.append((para[i][0], para[i][1], para[i][2], para[i][3], para[i][4], para[i][5]))
    para_sorted = np.array(values, dtype=dtype)
    para_sorted = np.sort(para_sorted, order=['minr', 'maxr'])
    return para_sorted


if __name__ == "__main__":
    all_paras = paras_sorter("data/ps5_sl3_stats.npy")
    final_stfr = np.load("data/stfr_inh_40_120ks.npy")
    savepath = "data/sta_50pts_40_120ks.npy"
    multi_seed_run(all_paras, lst=[10, 40, 70, 120, 170, 220, 260, 300, 330, 350], runtime=10000, pretime=200,
                   fr_array=np.array([0.5, 1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                      21.0, 22.0, 23.0, 24.0, 25, 30]),
                   final_stfr=final_stfr, 
                   seeds=np.array(range(40000, 120000))+2021,
                   savepath=savepath, numpts=50)
