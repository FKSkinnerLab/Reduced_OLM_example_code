import numpy as np
from neuron import h, gui
import record_1comp as r1
import efel
import multiprocessing as mp
import IVL_select as ivs

# Instantiate Model
h.load_file("init_1comp.hoc")
h.cvode_active(0)
h.dt = 0.1
h.steps_per_ms = 10
DNQX_Hold = 0.004
TTX_Hold = -0.0290514


def paras_sorter(para_path):
    para = np.load(para_path)
    values = []
    dtype = [('minr', float), ('maxr', float), ('ge', float), ('gi', float), ('stde', float), ('stdi', float)]
    for i in range(para.shape[0]):
        values.append((para[i][0], para[i][1], para[i][2], para[i][3], para[i][4], para[i][5]))
    para_sorted = np.array(values, dtype=dtype)
    para_sorted = np.sort(para_sorted, order=['minr', 'maxr'])
    return para_sorted


def state_extra_stats(runtime, lst, loadpath, savepath):
    """Want to record subthreshold mean and variance of each seed of each selected parameter set"""
    paras = paras_sorter(loadpath)[lst]
    proc_list = []
    shared_dic = mp.Manager().dict()
    for i in range(len(lst)):
        p = mp.Process(target=state_extra_stats_single, args=(paras[i], runtime, shared_dic, i))
        p.start()
        proc_list.append(p)
    for p in proc_list:
        p.join()

    v_means = None
    v_vars = None
    for key in range(len(lst)):
        if shared_dic[key] is None:
            sing_result = (np.zeros(50) - 1, np.zeros(50) - 1)
        else:
            sing_result = shared_dic[key]

        if v_means is None:
            v_means = sing_result[0]
        else:
            v_means = np.vstack((v_means, sing_result[0]))

        if v_vars is None:
            v_vars = sing_result[1]
        else:
            v_vars = np.vstack((v_vars, sing_result[1]))

    if v_means is None or v_vars is None:
        print("Mission failed, we will get'em next time.")
    else:
        results = np.vstack((v_means, v_vars))
        np.save(savepath, results)


def state_extra_stats_single(paras, runtime, shared_dic, proc_num):
    h.tstop = runtime
    seeds = np.array(range(50)) + 2021
    OU_in = h.Gfluct2(h.soma(0.5))
    OU_in.E_i = -87.1
    OU_in.tau_e = 12.7
    OU_in.tau_i = 12.05
    t_vec = np.array((range(int(runtime / h.dt) + 1))) * h.dt
    v_hvec = r1.set_up_full_recording()[0]

    efel.api.setThreshold(-30)
    efel.api.setDerivativeThreshold(1)
    para_vmeans = np.zeros(50)
    para_vvars = np.zeros(50)
    OU_in.g_e0 = paras[2]
    OU_in.g_i0 = paras[3]
    OU_in.std_e = paras[4]
    OU_in.std_i = paras[5]

    for j in range(50):
        OU_in.new_seed(seeds[j])
        h.run()
        v_vec = np.array(v_hvec)
        v_mean, v_var = ivs.sub_thresh_calc(t_vec, v_vec)
        para_vmeans[j] = v_mean
        para_vvars[j] = v_var
    shared_dic[proc_num] = (para_vmeans, para_vvars)


if __name__ == "__main__":
    lst = [10, 40, 70, 120, 170, 220, 260, 300, 330, 350]   # the 10 reoresentative sets
    runtime = 10000
    loadpath = "ps5_sl3_stats.npy"
    savepath = "ps5_sl3_exstats.npy"
    state_extra_stats(runtime, lst, loadpath, savepath)
