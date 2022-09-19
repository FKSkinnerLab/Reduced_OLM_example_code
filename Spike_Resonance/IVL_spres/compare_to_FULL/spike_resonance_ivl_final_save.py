from spike_resonance_ivl_single import *
import numpy as np

def paras_sorter(para_path):
    para = np.load(para_path)
    values = []
    dtype = [('minr', float), ('maxr', float), ('ge', float), ('gi', float), ('stde', float), ('stdi', float)]
    for i in range(para.shape[0]):
        values.append((para[i][0], para[i][1], para[i][2], para[i][3], para[i][4], para[i][5]))
    para_sorted = np.array(values, dtype=dtype)
    para_sorted = np.sort(para_sorted, order=['minr', 'maxr'])
    return para_sorted


def spres_producer(savedir, lst, para_sorted, fi_array):
    for num in range(len(lst)):
        para = para_sorted[lst[num]]
        ivl_paras = np.array([para[2], para[3], para[4], para[5]])

        base_rates_finder(10000, f"{savedir}/{num}/fb_vivo.npy", ivl_paras=ivl_paras)
        p1 = mp.Process(target=resonance_exp_vivo,
                        args=(10000, fi_array, ivl_paras, [f'{savedir}/{num}/ext_baseline_ratio_vivo_single.npy',
                                                           f'{savedir}/{num}/ext_firing_rate_vivo_single.npy'],
                              np.array(range(50)) + 2021))
        p2 = mp.Process(target=resonance_exp_vivo_inh,
                        args=(10000, fi_array, ivl_paras, [f'{savedir}/{num}/inh_baseline_ratio_vivo_single.npy',
                                                           f'{savedir}/{num}/inh_firing_rate_vivo_single.npy'],
                              np.array(range(50)) + 2021))
        p1.start()
        p2.start()
        p1.join()
        p2.join()


def spres_plotter(savedir, img_dir, lst, fi_array):
    for term in ["ext", "inh"]:
        if term == "ext":
            height = 105
        else:
            height = 10.5
        for num in range(len(lst)):
            baseline_ratio_plotter(fi_array, f"{savedir}/{num}/fb_vivo.npy", height,
                                   f"{savedir}/{num}/{term}_baseline_ratio_vivo_single.npy",
                                   f"{img_dir}/{num}/{term}_baseline_ratio_vivo_single.png")

            resonant_freq_histogram_plotter(fi_array, f"{savedir}/{num}/{term}_baseline_ratio_vivo_single.npy",
                                            f"{img_dir}/{num}/{term}_res_freq_vitro_hist_single.png")

            fb_v_fr_plot(fi_array, f"{savedir}/{num}/fb_vivo.npy", f"{savedir}/{num}/{term}_baseline_ratio_vivo_single.npy",
                         f'{img_dir}/{num}/{term}_fb_v_fr_vivo_single.png')

            fb_v_spikerate_fr_plot(f"{savedir}/{num}/{term}_firing_rate_vivo_single.npy",
                                   f"{savedir}/{num}/fb_vivo.npy",
                                   f"{savedir}/{num}/{term}_baseline_ratio_vivo_single.npy",
                                   f"{img_dir}//{num}/{term}_fb_v_spikerate_fr_single.png")


def collected_process(fi_array, lst, para_path, savedir, img_dir, print_out=True):
    para_sorted = paras_sorter(para_path)
    if print_out:
        print(para_sorted[lst])
    spres_producer(savedir, lst, para_sorted, fi_array)
    spres_plotter(savedir, img_dir, lst, fi_array)


def exstats_plotter(num_states, exstats_path, save_dir, img_dir):
    exstats = np.load(exstats_path)
    for term in ['ext', 'inh']:
        seed_to_fr = None
        fb_vivo = None
        for i in range(num_states):
            temp = np.load(f"{save_dir}/{i}/{term}_seed_to_fr.npy")
            temp2 = np.load(f"{save_dir}/{i}/fb_vivo.npy")
            if seed_to_fr is None:
                seed_to_fr = temp
            else:
                seed_to_fr = np.vstack((seed_to_fr, temp))
            if fb_vivo is None:
                fb_vivo = temp2
            else:
                fb_vivo = np.vstack((fb_vivo, temp2))
        exstat_v_fr_plot(seed_to_fr, exstats[:num_states], fb_vivo,
                         "Subthreshold V mean", f"{img_dir}/{term}_vmean_v_fr.png")
        exstat_v_fr_plot(seed_to_fr, exstats[num_states:], fb_vivo,
                         "Subthreshold V variance", f"{img_dir}/{term}_vvar_v_fr.png")


def fb_v_fr_combiner(savedir, img_dir, lst, fi_array):
    for term in ["ext", "inh"]:
        fb_path_lst = []
        bsr_path_lst = []
        savepath = f'{img_dir}/{term}_fb_v_fr_vivo_comb.png'
        for num in range(len(lst)):
            fb_path_lst.append(f"{savedir}/{num}/fb_vivo.npy")
            bsr_path_lst.append(f"{savedir}/{num}/{term}_baseline_ratio_vivo_single.npy")
        combined_fb_v_fr_plot(fi_array, fb_path_lst, bsr_path_lst, savepath)


def all_state_fr_hist(num_states, fi_array, save_dir, img_dir):
    seed_to_fr = None
    for i in range(num_states):
        temp = np.load(f"{save_dir}/{i}/seed_to_fr.npy")[:, 1]

        if seed_to_fr is None:
            seed_to_fr = temp
        else:
            seed_to_fr = np.hstack((seed_to_fr, temp))
    fr_hist_plotter(fi_array, seed_to_fr, f'{img_dir}/comb_hist.png')


if __name__ == "__main__":
    mp.freeze_support()
    fi_array = np.array([0, 0.5, 1, 2, 3, 4, 5, 8, 9, 10, 12, 15, 16, 20, 25, 30])
    lst = [10, 40, 70, 120, 170, 220, 260, 300, 330, 350]
    savedir = "data/spres_ps5"                              # where to save the spike resonance results
    img_dir = "1_comp_plots/IVL/spike_res/spres_ps5"        # where the plots go
    para_path = "data/IVL_select/ps5_sl3_stats.npy"         # see ReadMe on IVL_para_search and "_sl3_stats"
    exstats_path = "data/IVL_select/ps5_sl3_exstats.npy"    # see ReadMe on IVL_para_search and "_sl3_exstats"

    collected_process(fi_array, lst, para_path, savedir, img_dir)

    exstats_plotter(num_states=10, exstats_path=exstats_path, save_dir=savedir, img_dir=img_dir)

    fb_v_fr_combiner(savedir, img_dir, lst, fi_array)

    all_state_fr_hist(num_states=10, fi_array=fi_array, save_dir=savedir, img_dir=img_dir)
