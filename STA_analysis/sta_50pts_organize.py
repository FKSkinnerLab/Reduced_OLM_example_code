import numpy as np

filepaths = ["1comp_opt/re_all_20ks/sta_50pts.npy", "1comp_opt/re_all_20_40ks/sta_50pts_20_40ks.npy",
             "1comp_opt/re_all_40_120ks/sta_50pts_40_120ks.npy"]

stas = [np.load(name, allow_pickle=True) for name in filepaths]

final_path = "1comp_opt/STA_re/sta_50pts_final.npy"

# shape of files are (# of sets, len(fr_array), x, 7, 2000)
# desired final shape is (10, 27, 50, 7, 2000)

sta_final = np.empty((10, 27, 50, 7, 2000))
for i in range(10):
    for j in range(27):
        count = 0
        k = 0
        while count < 50:
            fut_count = stas[k][i, j].shape[0]
            sta_final[i, j][count: min(fut_count+count, 50)] = stas[k][i, j][: min(fut_count, 50-count)]
            count += fut_count
            k += 1

assert sta_final.shape == (10, 27, 50, 7, 2000)
np.save(final_path, sta_final)
