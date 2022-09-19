from neuron import h
import record_1comp as r1
from IVL_helper import *
import multiprocessing as mp
import IVL_select
import numpy as np

# Instantiate Model
h.load_file("init_1comp.hoc")
h.cvode_active(0)
h.dt = 0.1
h.steps_per_ms = 10
DNQX_Hold = 0.004
TTX_Hold = -0.0290514
# cell = h.cell

# The objective is to find enough IVL states using the newly optimized gbars


def partial_func2(i1, i2, ge_array, gi_array, stde_array, stdi_array, dic, dest):
    runtime = 10000
    h.tstop = runtime

    OU_in = h.Gfluct2(h.soma(0.5))
    OU_in.E_i = -87.1

    if not dest:
        OU_in.tau_e = 12.7
        OU_in.tau_i = 12.05

    v_hvec = r1.set_up_full_recording()[0]

    t_vec = np.array((range(int(runtime / h.dt) + 1))) * h.dt

    results = None

    for j1 in range(stde_array.shape[0]):
        for j2 in range(stdi_array.shape[0]):
            OU_in.g_e0 = ge_array[i1]
            OU_in.g_i0 = gi_array[i2]
            OU_in.std_e = stde_array[j1]
            OU_in.std_i = stdi_array[j2]

            h.run()

            v_vec = np.array(v_hvec)
            score = evaluator3(t_vec, v_vec)
            if score == 4:
                value_array = np.array([ge_array[i1], gi_array[i2], stde_array[j1], stdi_array[j2]])
                if results is None:
                    results = value_array
                else:
                    results = np.vstack((results, value_array))

    dic[(i1, i2)] = results


def partial_run2(save_path, num, dest):
    """ge0 = 0.01; gi0=0.03-0.04, stde=0.0005-0.001, stdi=0.009-0.01"""
    proc_list = []
    result_dict = mp.Manager().dict()
    ge_value = np.linspace(0.01, 0.02, num)
    gi_value = np.linspace(0.02, 0.03, num)
    stde_value = np.linspace(0.001, 0.004, num)
    stdi_value = np.linspace(0.009, 0.01, num)
    for i1 in range(num):
        for i2 in range(num):
            p1 = mp.Process(target=partial_func2, args=(i1, i2, ge_value, gi_value, stde_value, stdi_value, result_dict, dest))
            proc_list.append(p1)

    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()

    results = None
    for key in list(result_dict):
        if result_dict[key] is None:
            sing_result = np.array([-1, -1, -1, -1])  # ge, gi, stde, stdi
        else:
            sing_result = result_dict[key]

        if results is None:
            results = sing_result
        else:
            results = np.vstack((results, sing_result))

    results = results[results[:, -1] != -1]

    np.save(save_path, results)


def opt_run1_1node(save_path, num, dest, max_proc=125):
    """Optimized full search 1 (OFS 1). Taken from ps4"""
    proc_list = []
    result_dict = mp.Manager().dict()
    g_value = np.linspace(0, 0.5, 50)
    # std_value = np.linspace(0, 0.01, 50)
    std_value = np.linspace(0, 1, 50)
    results = None
    index1 = num * num / max_proc
    for i1 in range(num):
        for i2 in range(num):
            p1 = mp.Process(target=opt_func1, args=(i1, i2, g_value, std_value, result_dict, dest))
            proc_list.append(p1)

            for p in proc_list:
                p.start()
            for p in proc_list:
                p.join()

        for key in list(result_dict):
            if result_dict[key] is None:
                sing_result = np.array([-1, -1, -1, -1, -1, -1])
            else:
                sing_result = result_dict[key]

            if results is None:
                results = sing_result
            else:
                results = np.vstack((results, sing_result))

        results = results[results[:, -1] != -1]

    np.save(save_path, results)


def opt_func1(i1, i2, g_array, std_array, dic, dest):
    g_value = g_array
    std_value = std_array

    runtime = 10000
    h.tstop = runtime

    OU_in = h.Gfluct2(h.soma(0.5))
    OU_in.E_i = -87.1

    if not dest:
        OU_in.tau_e = 12.7
        OU_in.tau_i = 12.05

    v_hvec = r1.set_up_full_recording()[0]

    t_vec = np.array((range(int(runtime / h.dt) + 1))) * h.dt

    results = None

    for j1 in range(50):
        for j2 in range(50):
            OU_in.g_e0 = g_value[i1]
            OU_in.g_i0 = g_value[i2]
            OU_in.std_e = std_value[j1]
            OU_in.std_i = std_value[j2]

            h.run()

            v_vec = np.array(v_hvec)
            score = evaluator3(t_vec, v_vec)

            if score == 4:
                v_mean, v_var = IVL_select.sub_thresh_calc(t_vec, v_vec)
                value_array = np.array([g_value[i1], g_value[i2], std_value[j1], std_value[j2], v_mean, v_var])
                if results is None:
                    results = value_array
                else:
                    results = np.vstack((results, value_array))

    dic[(i1, i2)] = results


def opt_run2_1node(save_path, num, dest, max_proc=125):
    """Optimized full search 1 (OFS 1). Taken from ps4"""
    proc_list = []
    result_dict = mp.Manager().dict()
    ge_value = np.linspace(0, 0.07, 50)
    gi_value = np.linspace(0, 0.5, 50)
    stde_value = np.linspace(0, 0.05, 50)
    stdi_value = np.linspace(0, 0.25, 50)
    results = None
    for i1 in range(num):
        for i2 in range(num):
            p1 = mp.Process(target=opt_func2, args=(i1, i2, ge_value, gi_value, stde_value,
                                                    stdi_value, result_dict, dest))
            proc_list.append(p1)
    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()
    for key in list(result_dict):
        if result_dict[key] is None:
            sing_result = np.array([-1, -1, -1, -1, -1, -1])
        else:
            sing_result = result_dict[key]

        if results is None:
            results = sing_result
        else:
            results = np.vstack((results, sing_result))
    results = results[results[:, -1] != -1]
    np.save(save_path, results)


def opt_func2(i1, i2, ge_array, gi_array, stde_array, stdi_array, dic, dest):
    ge_value = ge_array
    gi_value = gi_array
    stde_value = stde_array
    stdi_value = stdi_array

    runtime = 10000
    h.tstop = runtime

    OU_in = h.Gfluct2(h.soma(0.5))
    OU_in.E_i = -87.1

    if not dest:
        OU_in.tau_e = 12.7
        OU_in.tau_i = 12.05

    v_hvec = r1.set_up_full_recording()[0]

    t_vec = np.array((range(int(runtime / h.dt) + 1))) * h.dt

    results = None
    num3 = stde_array.size
    num4 = stdi_array.size
    for j1 in range(num3):
        for j2 in range(num4):
            OU_in.g_e0 = ge_value[i1]
            OU_in.g_i0 = gi_value[i2]
            OU_in.std_e = stde_value[j1]
            OU_in.std_i = stdi_value[j2]
            h.run()
            v_vec = np.array(v_hvec)
            score = evaluator3(t_vec, v_vec)

            if score == 4:
                v_mean, v_var = IVL_select.sub_thresh_calc(t_vec, v_vec)
                value_array = np.array([ge_value[i1], gi_value[i2], stde_value[j1], stdi_value[j2], v_mean, v_var])
                if results is None:
                    results = value_array
                else:
                    results = np.vstack((results, value_array))

    dic[(i1, i2)] = results


def partial_run3(save_path, num, dest):
    proc_list = []
    result_dict = mp.Manager().dict()
    ge_value = np.linspace(0, 0.03, num)    # ge and gi of the first 20-30 pairs of ofs2 (see ofs-100b)
    gi_value = np.linspace(0, 0.06, num)
    stde_value = np.linspace(0, 0.001, num)
    stdi_value = np.linspace(0, 0.03, num)
    for i1 in range(num):
        for i2 in range(num):
            p1 = mp.Process(target=opt_func2,
                            args=(i1, i2, ge_value, gi_value, stde_value, stdi_value, result_dict, dest))
            proc_list.append(p1)

    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()

    results = None
    for key in list(result_dict):
        if result_dict[key] is None:
            sing_result = np.array([-1, -1, -1, -1, -1, -1])  # ge, gi, stde, stdi
        else:
            sing_result = result_dict[key]

        if results is None:
            results = sing_result
        else:
            results = np.vstack((results, sing_result))

    results = results[results[:, -1] != -1]

    np.save(save_path, results)


def universal_run(paras, save_path, dest):
    proc_list = []
    result_dict = mp.Manager().dict()
    ge_value = paras[0]    # ge and gi of the first 20-30 pairs of ofs2 (see ofs-100b)
    gi_value = paras[1]
    stde_value = paras[2]
    stdi_value = paras[3]
    num1 = ge_value.size
    num2 = gi_value.size
    for i1 in range(num1):
        for i2 in range(num2):
            p1 = mp.Process(target=opt_func2,
                            args=(i1, i2, ge_value, gi_value, stde_value, stdi_value, result_dict, dest))
            proc_list.append(p1)

    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()

    results = None
    for key in list(result_dict):
        if result_dict[key] is None:
            sing_result = np.array([-1, -1, -1, -1, -1, -1])  # ge, gi, stde, stdi
        else:
            sing_result = result_dict[key]

        if results is None:
            results = sing_result
        else:
            results = np.vstack((results, sing_result))

    results = results[results[:, -1] != -1]

    np.save(save_path, results)


if __name__ == "__main__":
    mp.freeze_support()
    num = 15
    paras = np.array([np.linspace(0.003, 0.006, num), np.linspace(0.008, 0.012, num),
                      np.linspace(0, 0.0016, num), np.linspace(0, 0.01, num)])
    universal_run(paras, "ps5.npy", dest=False)
