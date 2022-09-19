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


def spres_plotter(savedir, img_dir, lst, fi_array):
    for term in ["inh"]:
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


# def collected_process(fi_array, lst, para_path, savedir, img_dir, seeds, print_out=True, plot_out=True):
#     para_sorted = paras_sorter(para_path)
#     if print_out:
#         print(para_sorted[lst])
#     spres_producer(savedir, lst, para_sorted, fi_array, seeds)
#     if plot_out:
#         spres_plotter(savedir, img_dir, lst, fi_array)


def exstats_plotter(num_states, exstats_path, save_dir, img_dir):
    seed_to_fr = None
    fb_vivo = None
    exstats = np.load(exstats_path)
    for i in range(num_states):
        temp = np.load(f"{save_dir}/{i}/seed_to_fr.npy")
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
                     "Subthreshold V mean", f"{img_dir}/inh_vmean_v_fr.png")
    exstat_v_fr_plot(seed_to_fr, exstats[num_states:], fb_vivo,
                     "Subthreshold V variance", f"{img_dir}/inh_vvar_v_fr.png")


def fb_v_fr_combiner(savedir, img_dir, lst, fi_array):
    for term in ["ext", "inh"]:
        fb_path_lst = []
        bsr_path_lst = []
        savepath = f'{img_dir}/{term}_fb_v_fr_vivo_comb.png'
        for num in range(len(lst)):
            fb_path_lst.append(f"{savedir}/{num}/fb_vivo.npy")
            bsr_path_lst.append(f"{savedir}/{num}/{term}_baseline_ratio_vivo_single.npy")
        combined_fb_v_fr_plot(fi_array, fb_path_lst, bsr_path_lst, savepath)


def layer2_dic_sorter(dic, save_dir, save_name):
    for key in dic:
        results = None
        sorted_list = dic[key]
        sorted_list.sort()
        for seg_tuple in sorted_list:
            if results is None:
                results = seg_tuple[1]
            else:
                if len(results.shape) == 1:
                    results = np.hstack((results, seg_tuple[1]))
                else:
                    results = np.vstack((results, seg_tuple[1]))
        np.save(f"{save_dir}/{key}/{save_name}", results)


def func_wrap(para_path, all_seeds, seg_size, fi_array, lst, manager):
    proc_list = []
    para_sorted = paras_sorter(para_path)
    dic_lst = [manager.dict() for _ in range(3)]
    for dic in dic_lst:
        for i in range(len(lst)):
            dic[i] = manager.list()
    for seg_num in range(0, len(all_seeds), seg_size):
        for ivl_num in range(len(lst)):
            para = para_sorted[lst[ivl_num]]
            ivl_paras = np.array([para[2], para[3], para[4], para[5]])
            p1 = mp.Process(target=base_rates_finder,
                            args=(10000, dic_lst[0], ivl_paras, all_seeds[seg_num: seg_num + seg_size], ivl_num, seg_num))
            proc_list.append(p1)
            p2 = mp.Process(target=resonance_exp_vivo_inh,
                            args=(10000, fi_array, ivl_paras, all_seeds[seg_num: seg_num + seg_size],
                                  dic_lst[1], dic_lst[2], ivl_num, seg_num))
            proc_list.append(p2)
    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()

    return dic_lst


def func_wrap2(para_path, all_seeds, seg_size, fi_array, lst, manager):
    proc_list = []
    para_sorted = paras_sorter(para_path)
    dic_lst = [manager.dict() for _ in range(3)]
    for dic in dic_lst:
        for i in range(len(lst)):
            dic[i] = manager.list()
    for seg_num in range(0, len(all_seeds), seg_size):
        for ivl_num in range(len(lst)):
            para = para_sorted[lst[ivl_num]]
            ivl_paras = np.array([para[2], para[3], para[4], para[5]])
            p1 = mp.Process(target=base_rates_finder,
                            args=(10000, dic_lst[0], ivl_paras, all_seeds[seg_num: seg_num + seg_size], ivl_num, seg_num))
            proc_list.append(p1)
            p2 = mp.Process(target=resonance_exp_vivo,
                            args=(10000, fi_array, ivl_paras, all_seeds[seg_num: seg_num + seg_size],
                                  dic_lst[1], dic_lst[2], ivl_num, seg_num))
            proc_list.append(p2)
    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()

    return dic_lst


def all_state_fr_hist(states, fi_array, save_dir, img_dir, norm=False):
    seed_to_fr = None
    for i in states:
        temp = np.load(f"{save_dir}/{i}/seed_to_fr.npy")[:, 1]

        if seed_to_fr is None:
            seed_to_fr = temp
        else:
            seed_to_fr = np.hstack((seed_to_fr, temp))
    fr_hist_plotter(fi_array, seed_to_fr, f'{img_dir}/comb_hist.png', norm)


if __name__ == "__main__":
    mp.freeze_support()

    fi_array = np.array([0, 0.5, 1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21.0, 22.0,
                         23.0, 24.0, 25, 30])
    lst = [10, 40, 70, 120, 170, 220, 260, 300, 330, 350]

    all_seeds = np.array(range(5000)) + 2021
    seg_size = 50
    para_path = "data/ps5_sl3_stats.npy"
    savedir = "data/spres_ps5"
    img_dir = "plots/spres_ps5"
    exstats_path = "data/ps5_sl3_exstats.npy"

    save_names = ["fb_vivo.npy", "inh_baseline_ratio_vivo_single.npy", "inh_firing_rate_vivo_single.npy"]

    manager = mp.Manager()
    dic_lst = func_wrap(para_path, all_seeds, seg_size, fi_array, lst, manager)
    for i in range(3):
        layer2_dic_sorter(dic_lst[i], savedir, save_names[i])

    save_names = ["fb_vivo.npy", "ext_baseline_ratio_vivo_single.npy", "ext_firing_rate_vivo_single.npy"]

    manager = mp.Manager()
    dic_lst = func_wrap2(para_path, all_seeds, seg_size, fi_array, lst, manager)
    for i in range(5):
        layer2_dic_sorter(dic_lst[i], savedir, save_names[i])
