import numpy as np
import multiprocessing as mp


def seed_to_fr(bsr_path, seeds_array, fi_array):
    basline_ratio = np.load(bsr_path)
    results = np.zeros((seeds_array.shape[0], 2))
    for i in range(seeds_array.shape[0]):
        fr = fi_array[np.nanargmax(basline_ratio[i])]
        results[i][0] = seeds_array[i]
        results[i][1] = fr
    return results


# comment and uncommon parts to work for 20ks, 20-40ks, 40-120ks

"""20ks"""
# bsr_dirs = ["1comp_opt/re_all_20ks/spres_vivo_0_to_30_5000s", "1comp_opt/re_all_20ks/re_all_2",
#             "1comp_opt/re_all_20ks/re_all_3", "1comp_opt/re_all_20ks/re_all_4"]
# seeds_arrays = [np.array(range(5000)) + 2021, np.array(range(5000, 10000)) + 2021,
#                 np.array(range(10000, 15000)) + 2021, np.array(range(15000, 20000)) + 2021]
# final_path = f"1comp_opt/re_all_20ks/stfr_inh_20_40ks.npy"

"""20-40ks"""
# interm = "re_all_20_40ks"
# bsr_dirs = [f"1comp_opt/{interm}/re_all_1", f"1comp_opt/{interm}/re_all_2",
#             f"1comp_opt/{interm}/re_all_3", f"1comp_opt/{interm}/re_all_4"]
#
# seeds_arrays = [np.array(range(20000, 25000)) + 2021, np.array(range(25000, 30000)) + 2021,
#                 np.array(range(30000, 35000)) + 2021, np.array(range(35000, 40000)) + 2021]
#
# final_path = f"1comp_opt/{interm}/stfr_inh_20_40ks.npy"
# final_seed_to_fr = np.zeros((4, 5000, 10))

"""40-120ks"""
interm = "re_all_40_120ks"
bsr_dirs = [f"1comp_opt/{interm}/re_all_1", f"1comp_opt/{interm}/re_all_2",
            f"1comp_opt/{interm}/re_all_3", f"1comp_opt/{interm}/re_all_4"]

seeds_arrays = [np.array(range(40000, 60000)) + 2021, np.array(range(60000, 80000)) + 2021,
                np.array(range(80000, 100000)) + 2021, np.array(range(100000, 120000)) + 2021]

final_path = f"1comp_opt/{interm}/stfr_inh_40_120ks.npy"
final_seed_to_fr = np.zeros((4, 20000, 10))

"""Common parts shared between different number of seeds"""
num_states = 10
fi_list = np.array([0, 0.5, 1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21.0, 22.0, 23.0,
                    24.0, 25, 30])
fr_array = np.array([0.5, 1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21.0, 22.0, 23.0,
                     24.0, 25, 30])

for j in range(4):
    for i in range(num_states):
        local_stfr = seed_to_fr(bsr_path=f"{bsr_dirs[j]}/{i}/inh_baseline_ratio_vivo_single.npy",
                                seeds_array=seeds_arrays[j], fi_array=fi_list)
        final_seed_to_fr[j, :, i] = local_stfr[:, 1]

final_seed_to_fr = np.reshape(final_seed_to_fr, (-1, 10))
np.save(final_path, final_seed_to_fr)


