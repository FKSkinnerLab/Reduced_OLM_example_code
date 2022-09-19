import numpy as np


def spread_calc_comb(ud_norm_sep, curr_lst, startend):
    """raw_prop is (len(fr_array), 7, # of sets, 2000). Calculate spread after summing currents indicated
    in curr_lst. startend is a tuple indicates portion of the sta being considered"""
    combined_curr = ud_norm_sep[:, curr_lst, :, :]    # (len(fr_array), k, # of sets, 2000)
    combined_curr = np.sum(combined_curr, axis=1)       # (len(fr_array), # of sets, 2000)
    spreads = 2 * np.sqrt(np.var(combined_curr, axis=1))    # (len(fr_array), 2000)
    spreads = np.sum(spreads, axis=1)                  # (len(fr_array), ))
    return spreads


def slope_calc_comb(ud_norm_sep, curr_lst, startend):
    """ud_norm is (len(fr_array), 7, # of sets, 2000). Calculate slopes and sum for currents indicated
    in curr_lst."""
    combined_curr = ud_norm_sep[:, curr_lst, :, :]  # (len(fr_array), k, # of sets, 2000)
    slopes = np.zeros(combined_curr.shape[0])          # (len(fr_array), )
    for i in range(combined_curr.shape[1]):
        curr = combined_curr[:, i, :, :]     # (len(fr_array), # of sets, 2000)
        curr = np.mean(curr, axis=1)      # (len(fr_array), 2000)
        x = np.linspace(start=startend[0]*0.1, stop=startend[1]*0.1, num=curr.shape[-1])
        slopes += np.polyfit(x, curr.T, deg=1)[0]
    return slopes


def undirected_norm_comb(raw_prop):
    """raw_prop is of shape (# of sets, len(fr_array), 7, 2000). We take absolute value and normalize for each fr,
    for each current type."""
    raw_prop = np.abs(raw_prop)
    raw_prop = np.transpose(raw_prop, axes=(1, 2, 0, 3))     # (len(fr_array), 7, # of sets, 2000)
    raw_prop = np.reshape(raw_prop, (*raw_prop.shape[:-2], -1))     # (len(fr_array), 7, 20000)
    maxima = np.amax(raw_prop, axis=2)                  # (len(fr_array), 7)
    minima = np.amin(raw_prop, axis=2)
    result = (raw_prop-np.reshape(minima, (27, 6, 1))) / np.reshape(maxima - minima, (27, 6, 1))
    return np.reshape(result, (27, 6, 10, -1))    # (len(fr_array), 7, # of sets, 2000)


def raw_proportion_comb(all_sta, startend):
    """all_sta is of shape (# of sets, len(fr_array), 50, 7, 2000). We want raw proportion within each of the
    (1, 1, 1, 7, 2000) (excluding voltage). The returned array will be without voltage"""
    all_sta = np.mean(all_sta[:, :, :, 1:, startend[0]:startend[1]], axis=2)      # (# of sets, len(fr_array), 7, 2000)
    all_total = np.sum(np.abs(all_sta), axis=2)             # (# of sets, len(fr_array), 2000)
    return all_sta / np.reshape(all_total, (all_sta.shape[0], all_sta.shape[1], 1, -1))


if __name__ == "__main__":
    shared_startend = (-1950, -250)
    all_sta = np.load("1comp_opt/STA_re/sta_50pts_final.npy", allow_pickle=True)

    raw_prop = raw_proportion_comb(all_sta, startend=shared_startend)
    all_sta = 0
    np.save("1comp_opt/STA_re/comb/sta_50pts_final_raw_prop_comb.npy", raw_prop)

    raw_prop = np.load("1comp_opt/STA_re/comb/sta_50pts_final_raw_prop_comb.npy")
    ud_norm_comb = undirected_norm_comb(raw_prop)
    raw_prop = 0
    np.save("1comp_opt/STA_re/comb/sta_50pts_final_udnorm_comb.npy", ud_norm_comb)

    ud_norm_comb = np.load("1comp_opt/STA_re/comb/sta_50pts_final_udnorm_comb.npy")
    for curr_lst, curr_name in zip([[5], [2], [2, 5]], ['h', 'm', 'hm']):           # IKa, IKdrf, Im, Il, INa, Ih
        spread_comb = spread_calc_comb(ud_norm_comb, curr_lst, startend=shared_startend)
        np.save(f"1comp_opt/STA_re/comb/sta_50pts_final_{curr_name}_spread_comb.npy", spread_comb)

        slope_comb = slope_calc_comb(ud_norm_comb, curr_lst, startend=shared_startend)
        np.save(f"1comp_opt/STA_re/comb/sta_50pts_final_{curr_name}_slope_comb.npy", slope_comb)


