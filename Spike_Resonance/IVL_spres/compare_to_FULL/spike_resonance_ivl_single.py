import numpy as np
from neuron import h, gui
import record_1comp as r1
import efel
import scipy.signal as signal
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import efel
# import glob
# import IPython, os
import multiprocessing as mp
from currents_visualization import *

### Instantiate Model ###
h.load_file("init_1comp.hoc")
h.cvode_active(0)
h.dt = 0.1
h.steps_per_ms = 10
DNQX_Hold = 0.004
TTX_Hold = -0.0290514
# cell = h.cell

"""
Most of the testing for spike resonance was done while setting up the inhibitory perturbation, so all the related tests
are in spike_resonance_inh.py.
In addition, baseline spiking rates are calculated in spike_resonance_inh.py, so here we use it right away.
"""


def base_rates_finder(runtime, filepath, ivl_paras=None):
    seeds = np.array(range(50)) + 2021
    if ivl_paras is None:
        ivl_paras = np.sort(np.load("data/IVL_select/ps1_agm_sl3_stats.npy"), axis=0)[53]
    base_rates = np.zeros(len(seeds))
    apc = h.APCount(h.soma(0.5))
    apc.thresh = -30
    h.tstop = runtime

    OU_in = h.Gfluct2(h.soma(0.5))
    OU_in.E_i = -87.1
    for k in range(50):
        line = ivl_paras
        OU_in.g_e0 = line[0]
        OU_in.g_i0 = line[1]
        OU_in.std_e = line[2]
        OU_in.std_i = line[3]
        OU_in.tau_e = 12.7
        OU_in.tau_i = 12.05
        OU_in.new_seed(seeds[k])

        spike_times = h.Vector()
        apc.record(spike_times)

        h.run()

        spike_times = np.array(spike_times)
        base_rates[k] = len(spike_times) / (runtime / 1000)  # covert runtime from ms to s
    np.save(filepath, base_rates)


def resonance_exp_vivo(runtime, fi, ivl_paras, savepaths, seeds):
    """This function finds baseline ratios and spike rate at each fi and ampb"""

    h.ic_hold.delay = 0
    h.ic_hold.dur = runtime
    h.tstop = runtime

    # creating synapse
    stim = h.NetStim()  # average number of spikes. convert to ms, then s
    stim.noise = 0  # deterministic
    stim.start = 0

    syn = h.Exp2Syn(h.soma(0.5))
    syn.tau1 = 2.4  # ms rise time
    syn.tau2 = 12.7  # ms decay time
    syn.e = 0  # reversal potential

    netcon = h.NetCon(stim, syn)  # threshold is irrelevant with event-based source
    netcon.weight[0] = 0.006

    apc = h.APCount(h.soma(0.5))
    apc.thresh = -30
    spike_timing = h.Vector()
    apc.record(spike_timing)

    # each row a diff fb, each col a diff fi
    results = np.zeros((50, len(fi)))
    firing_rates = np.zeros((50, len(fi)))

    # set up ivl
    OU_in = h.Gfluct2(h.soma(0.5))
    OU_in.E_i = -87.1
    # seeds = np.array(range(50)) + 2021
    for k in range(50):
        line = ivl_paras
        OU_in.g_e0 = line[0]
        OU_in.g_i0 = line[1]
        OU_in.std_e = line[2]
        OU_in.std_i = line[3]
        OU_in.tau_e = 12.7
        OU_in.tau_i = 12.05
        OU_in.new_seed(seeds[k])
        for i in range(len(fi)):
            if fi[i] == 0:
                stim.number = 0
                stim.interval = 1
            else:
                stim.interval = 1000 / fi[i]   # interval in ms, but fi are in Hz
                stim.number = fi[i] * (runtime/1000)

            h.run()

            spike_times = np.array(spike_timing) / h.dt  # in 0.1 ms
            spike_times = spike_times.astype(int)
            spike_train = np.zeros(int(runtime/h.dt))

            spike_train[spike_times] = 1
            spike_train = spike_train[10001:]   # take away first 1 second

            psd = signal.welch(spike_train, fs=1/(h.dt/1000), scaling='density', nperseg=20000)

            if fi[i] == 0:
                base_psd = psd
                results[k][i] = 1
            else:
                results[k][i] = psd[1][np.argwhere(psd[0] == fi[i])] / base_psd[1][np.argwhere(psd[0] == fi[i])]

            firing_rates[k][i] = np.sum(spike_train) / (runtime/1000 - 1) # since we took away data from 1st second

    # np.save('data/spres_fs4_sl1/ext_baseline_ratio_vivo_single', results)
    # np.save('data/spres_fs4_sl1/ext_firing_rate_vivo_single', firing_rates)

    np.save(savepaths[0], results)
    np.save(savepaths[1], firing_rates)


def resonance_exp_vivo_inh(runtime, fi, ivl_paras, savepaths, seeds):
    """This function finds baseline ratios and spike rate at each fi and ampb"""

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

    apc = h.APCount(h.soma(0.5))
    apc.thresh = -30
    spike_timing = h.Vector()
    apc.record(spike_timing)

    # each row a diff seed, each col a diff fi
    results = np.zeros((50, len(fi)))
    firing_rates = np.zeros((50, len(fi)))

    # set up ivl
    OU_in = h.Gfluct2(h.soma(0.5))
    OU_in.E_i = -87.1
    # seeds = np.array(range(50)) + 2021
    for k in range(50):
        line = ivl_paras
        OU_in.g_e0 = line[0]
        OU_in.g_i0 = line[1]
        OU_in.std_e = line[2]
        OU_in.std_i = line[3]
        OU_in.tau_e = 12.7
        OU_in.tau_i = 12.05
        OU_in.new_seed(seeds[k])

        for i in range(len(fi)):
            if fi[i] == 0:
                stim.number = 0
                stim.interval = 1
            else:
                stim.interval = 1000 / fi[i]   # interval in ms, but fi are in Hz
                stim.number = fi[i] * (runtime/1000)

            h.run()

            spike_times = np.array(spike_timing) / h.dt  # in 0.1 ms
            spike_times = spike_times.astype(int)
            spike_train = np.zeros(int(runtime/h.dt))

            spike_train[spike_times] = 1
            spike_train = spike_train[10001:]   # take away first 1 second

            psd = signal.welch(spike_train, fs=1/(h.dt/1000), scaling='density', nperseg=20000)

            if fi[i] == 0:
                base_psd = psd
                results[k][i] = 1
            else:
                results[k][i] = psd[1][np.argwhere(psd[0] == fi[i])] / base_psd[1][np.argwhere(psd[0] == fi[i])]

            firing_rates[k][i] = np.sum(spike_train) / (runtime/1000 - 1) # since we took away data from 1st second

    # np.save('data/spres_fs4_sl1/inh_baseline_ratio_vivo_single', results)
    # np.save('data/spres_fs4_sl1/inh_firing_rate_vivo_single', firing_rates)

    np.save(savepaths[0], results)
    np.save(savepaths[1], firing_rates)


def baseline_ratio_plotter(fi, fb_path, dot_height, bsr_path, savepath):
    # bsr_path is the file that holds the baseline ratios
    bl_ratio = np.load(bsr_path)     # each row a diff fb, each col a diff fi
    baserates = np.load(fb_path)
    minrate = min(baserates)
    maxrate = max(baserates)
    fi_str = np.array(fi, dtype=str)

    fig = plt.figure()
    ax = fig.add_subplot(1,15,(1,14))
    ax1 = fig.add_subplot(1,15,15)
    colors = plt.cm.plasma((baserates - minrate) / (maxrate - minrate))

    for y in range(len(bl_ratio)):
        ax.plot(fi_str, bl_ratio[y], color=colors[y])

    # find fi with maximum baseline ratio excluding the first col, where no modulation
    res_freqs = fi[np.argmax(bl_ratio[:], axis=1)]
    res_freqs = res_freqs.astype(str)

    cmap = plt.cm.plasma
    norm = matplotlib.colors.Normalize(vmin=np.amin(baserates), vmax=np.amax(baserates))

    ax.scatter(res_freqs, np.zeros(len(res_freqs))+dot_height, c=colors, cmap=cmap, norm=norm)

    # ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax.set_xlim(0, 16)
    ax.set_xticks(fi_str)
    ax.set_xticklabels(fi_str, rotation=45)
    ax.set_xlabel(r'$f_i$')
    ax.set_ylabel('Baseline Ratio')

    cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label(r'$f_B$')
    plt.tight_layout()
    plt.savefig(savepath)


def resonant_freq_histogram_plotter(fi, bsr_path, savepath):
    # bsr_path is the file that holds the baseline ratios
    bl_ratio = np.load(bsr_path)  # each row a diff fb, each col a diff fi

    fig = plt.figure()
    ax = fig.add_subplot(15,1,(1,14))

    res_freqs = fi[np.argmax(bl_ratio[:], axis=1)]
    # res_freqs = res_freqs.astype(str)

    # the way plt.hist works: last bin is [fi[-2], fi[-1]], while other bins are [fi[x], fi[x+1])
    # we want each bin to include only 1 fi value, so we replicate the last fi value
    bins = list(fi)
    bins.append(bins[-1]+1)
    plt.hist(res_freqs, bins=bins, rwidth=0.1, align='left')
    ax.set_xlim(0, 16)
    ax.set_xticks(bins)
    ax.set_xticklabels(bins, rotation=45)
    ax.set_xlabel(r'$f_i$')
    plt.tight_layout()
    plt.savefig(savepath)


def fb_v_fr_plot(fi, fb_path, bsr_path, savepath):
    # bsr_path is the file that holds the baseline ratios
    bl_ratio = np.load(bsr_path)  # each row a diff fb, each col a diff fi
    baserates = np.load(fb_path)

    fig = plt.figure()
    ax = fig.add_subplot(15, 1, (1, 14))
    res_freqs = fi[np.argmax(bl_ratio, axis=1)]
    ax.scatter(baserates, res_freqs, color='b')
    ax.plot([0, 36], [0, 36], linestyle='--')
    ax.set_xlabel(r"$f_B$")
    ax.set_ylabel(r"$f_r$")
    plt.savefig(savepath)


def combined_fb_v_fr_plot(fi, fb_path_lst, bsr_path_lst, savepath):
    fig = plt.figure()
    ax = fig.add_subplot(15, 1, (1, 14))
    ax.plot([0, 36], [0, 36], linestyle='--')
    for i in range(len(fb_path_lst)):
        bsr_path = bsr_path_lst[i]
        fb_path = fb_path_lst[i]
        bl_ratio = np.load(bsr_path)
        baserates = np.load(fb_path)
        res_freqs = fi[np.argmax(bl_ratio, axis=1)]
        ax.scatter(baserates, res_freqs, color='b')
    ax.set_xlabel(r"$f_B$")
    ax.set_ylabel(r"$f_r$")
    plt.savefig(savepath)


def exstat_v_fr_plot(seed_to_fr, exstat, fb_vivo, exstat_name, savepath):
    # bsr_path is the file that holds the baseline ratios
    seed_to_fr = np.reshape(seed_to_fr, (-1, 2))[:, 1]
    exstat = np.reshape(exstat, (-1,))
    cmap = plt.cm.viridis
    norm = matplotlib.colors.Normalize(vmin=np.min(fb_vivo), vmax=np.max(fb_vivo))
    fig = plt.figure()
    ax = fig.add_subplot(1, 15, (1, 14))
    ax1 = fig.add_subplot(1, 15, 15)
    ax.scatter(exstat, seed_to_fr, c=fb_vivo, norm=norm, cmap=cmap)
    ax.set_xlabel(f"{exstat_name}")
    ax.set_ylabel(r"$f_r$")
    cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label(r'$f_B$')
    plt.savefig(savepath)


def fb_v_spikerate_fr_plot(spikerate_path, fb_path, bsr_path, savepath):
    spikerates = np.load(spikerate_path)    # each row a diff fb, each col a diff fi
    baserates = np.load(fb_path)
    bl_ratios = np.load(bsr_path)           # each row a diff fb, each col a diff fi

    res_freqs_idx = np.argmax(bl_ratios, axis=1)
    spikerate_at_fr = spikerates[range(len(baserates)), res_freqs_idx]  # dimensions of bl_ratios and spikerates are the same

    fig = plt.figure()
    ax = fig.add_subplot(15, 1, (1,14))
    ax.scatter(baserates, spikerate_at_fr)
    ax.plot([0, 36], [0, 36], linestyle='--')

    ax.set_xlabel(r"$f_B$")
    ax.set_ylabel(r"spiking rate at $f_r$")
    plt.tight_layout()
    plt.savefig(savepath)


def fr_hist_plotter(fr_array, all_fr, savepath):
    for i in range(len(fr_array)):
        all_fr[all_fr == fr_array[i]] = -i
    all_fr = all_fr * -1
    labels = [str(i) for i in fr_array]
    fig = plt.figure()
    ax = fig.add_subplot()
    bins = np.array(range(len(fr_array) + 1))
    ax.hist(all_fr, bins=bins, rwidth=0.5, align='left')
    ax.set_xticks(bins[:-1], labels, rotation=45)
    plt.savefig(savepath)
    plt.close()


if __name__ == '__main__':
    fi_array = np.array([0, 0.5, 1, 2, 3, 4, 5, 8, 9, 10, 12, 15, 16, 20, 25, 30])
