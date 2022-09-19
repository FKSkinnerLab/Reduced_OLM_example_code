import numpy as np
from neuron import h, gui
import record_1comp as r1
import efel
import scipy.signal as signal
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from currents_visualization import *

### Instantiate Model ###
h.load_file("init_1comp.hoc")
h.cvode_active(0)
h.dt = 0.1
h.steps_per_ms = 10
DNQX_Hold = 0.004
TTX_Hold = -0.0290514
# cell = h.cell

fi_vitro = np.array([0, 0.5, 1, 2, 3, 4, 5, 8, 9, 10, 12, 15, 16, 20, 25, 30])
ampb_vitro = currsteps = np.arange(0.025425, 0.025425+0.0026275*50, 0.0026275)    # 1hz is 25.425pA, 36 hz is 156.8pA


def base_rates_finder(ampb, runtime, filepath):
    base_rates = np.zeros(len(ampb))
    apc = h.APCount(h.soma(0.5))
    apc.thresh = -20
    h.ic_hold.delay = 0
    h.ic_hold.dur = runtime
    h.tstop = runtime

    for i in range(len(ampb)):
        h.ic_hold.amp = ampb[i]

        spike_times = h.Vector()
        apc.record(spike_times)

        h.run()

        spike_times = np.array(spike_times)
        base_rates[i] = len(spike_times) / (runtime / 1000)  # covert runtime from ms to s
    np.save(filepath, base_rates)


def spike_res_test(runtime, amp):
    h.ic_hold.delay = 0
    h.ic_hold.dur = runtime
    h.tstop = runtime
    h.ic_hold.amp = amp

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
    apc.thresh = -20
    spike_times = h.Vector()
    apc.record(spike_times)
    h.run()
    spike_times = np.array(spike_times, dtype=int) # in ms
    spike_train = np.zeros(runtime)
    spike_train[spike_times] = 1

    psd = signal.welch(spike_train, fs=1000, scaling='density', nperseg=int(len(spike_train)/2))
    # fs=1000 b/c spike train is based on runtime, which is in ms

    # plt.plot(v_trace)
    # plt.show()
    plt.plot(spike_train)
    plt.show()
    # plt.plot(psd[0], psd[1])
    # plt.show()
    # plt.plot(psd[0][0:100], psd[1][0:100])
    # plt.show()
    ide = np.argmax(psd[1])
    print("psd estimate:", psd[0][ide])
    print("true freq:", len(spike_times)/runtime*1000)


def psd_test(freq):
    trace = np.sin(np.linspace(0, 2*np.pi, 1000)*freq)
    psd = signal.welch(trace, fs=1000, scaling='density', nperseg=20000)

    plt.plot(trace)
    plt.show()
    plt.plot(psd[0], psd[1])
    plt.show()
    print(psd[0][np.argmax(psd[1])])


def baseline_ratio_test(runtime, amp, fi):
    h.ic_hold.delay = 0
    h.ic_hold.dur = runtime
    h.tstop = runtime
    h.ic_hold.amp = amp

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
    apc.thresh = -20
    spike_times = h.Vector()
    apc.record(spike_times)

    v_vec = r1.set_up_full_recording()[0]
    for i in range(1, len(fi)):
        stim.interval = 1000 / fi[i]  # interval in ms, but fi are in Hz
        stim.number = fi[i] * (runtime / 1000)

        h.run()
        spike_times = np.array(spike_times, dtype=int)  # in ms
        spike_train = np.zeros(runtime)
        spike_train[spike_times] = 1

        psd = signal.welch(spike_train, fs=1000, scaling='density', nperseg=int(len(spike_train) / 2))
        # fs=1000 b/c spike train is based on runtime, which is in ms

        v_trace = np.array(v_vec)

        ide = np.argmax(psd[1])
        print("psd estimate:", psd[0][ide])
        print("true freq:", len(spike_times) / runtime * 1000)

        plt.plot(v_trace)
        plt.show()
        plt.plot(spike_train)
        plt.show()
        plt.plot(psd[0], psd[1])
        plt.plot(psd[0][0:100], psd[1][0:100])
        plt.show()
        plt.close()


def resonance_exp_vitro(runtime, fi, ampb, blr_savepath, fr_savepath):
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
    apc.thresh = -20
    spike_timing = h.Vector()
    apc.record(spike_timing)

    # each row a diff fb, each col a diff fi
    results = np.zeros((len(ampb), len(fi)))
    firing_rates = np.zeros((len(ampb), len(fi)))

    for k in range(len(ampb)):
        h.ic_hold.amp = ampb[k]

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

    # np.save('data/baseline_ratio_vitro', results)
    # np.save('data/firing_rate_vitro', firing_rates)
    np.save(blr_savepath, results)
    np.save(fr_savepath, firing_rates)


def baseline_ratio_plotter(fi, fb_path, bsr_path, savepath):
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

    ax.scatter(res_freqs, np.zeros(len(res_freqs))+3.5e8, c=colors, cmap=cmap, norm=norm)

    # ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax.set_xlim(0, 16)
    ax.set_xticks(fi_str)
    ax.set_xticklabels(fi_str, rotation=45)
    ax.set_xlabel(r'$f_i$')
    ax.set_ylabel('Baseline Ratio')

    cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label(r'$f_B$')
    plt.savefig(savepath)


def resonant_freq_histogram_plotter(fi, fb_path, bsr_path, savepath):
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

    plt.savefig(savepath)


def fb_v_fr_plot(fi, fb_path, bsr_path, savepath):
    # bsr_path is the file that holds the baseline ratios
    bl_ratio = np.load(bsr_path)  # each row a diff fb, each col a diff fi
    baserates = np.load(fb_path)

    fig = plt.figure()
    ax = fig.add_subplot(15, 1, (1, 14))
    res_freqs = fi[np.argmax(bl_ratio, axis=1)]
    ax.scatter(baserates, res_freqs, color='b')
    ax.plot([0, baserates[-1]], [0, baserates[-1]], linestyle='--')
    ax.set_xlabel(r"$f_B$")
    ax.set_ylabel(r"$f_r$")
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
    ax.plot([0, baserates[-1]], [0, baserates[-1]], linestyle='--')

    ax.set_xlabel(r"$f_B$")
    ax.set_ylabel(r"spiking rate at $f_r$")

    plt.savefig(savepath)


if __name__ == '__main__':
    base_rates_finder(ampb_vitro, 10000, "data/fb_vitro.npy")
    resonance_exp_vitro(10000, fi_vitro, ampb_vitro, blr_savepath='data/baseline_ratio_vitro_inh',
                        fr_savepath='data/firing_rate_vitro_inh')

    baseline_ratio_plotter(fi_vitro, "data/fb_vitro.npy", "data/baseline_ratio_vitro_inh.npy",
                           "1_comp_plots/spike_res/baseline_ratio_vitro_inh.png")

    resonant_freq_histogram_plotter(fi_vitro, "data/fb_vitro.npy", "data/baseline_ratio_vitro_inh.npy",
                                    "1_comp_plots/spike_res/res_freq_vitro_hist_inh.png")
    fb_v_fr_plot(fi_vitro, "data/fb_vitro.npy", "data/baseline_ratio_vitro_inh.npy",
                 '1_comp_plots/spike_res/fb_v_fr_vitro_inh.png')
    fb_v_spikerate_fr_plot("data/firing_rate_vitro_inh.npy", "data/fb_vitro.npy", "data/baseline_ratio_vitro_inh.npy",
                           "1_comp_plots/spike_res/fb_v_spierate_fr_inh.png")
