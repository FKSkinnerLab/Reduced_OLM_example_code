The folders in this repo contains cleaned up example code that represent the original code used to produce the data in
the paper. The example code is not directly run-able since many of them requires access to a computer cluster. In
addition, the example code serves as a demonstration of what was done to produce the results. Thus, it may not be the
most optimal way of coding, and is by no means the optimal way of using computation resources, since readability is
prioritized.


The following is a guide to the contents of the folders.


Models:
This folder contains the NMODL and HOC files of FULL and SINGLE. Compilation of the NMODL files into dll file by nmkdll
or equivalent is needed prior to running with the NURON simulator.

FULL_data_extraction:
This folder contains the python code used to extract current, voltage and e_feature values from the FULL model. They
are only usable when they are in the same folder as the model HOC files and dll file.

SINGLE_optimization:
This folder contains the python code used to optimize SINGLE and the e_features json file from the FULL. NSG was used to
run the code. NMODL files of SINGLE also need to be in the same folder as the content of SINGLE_optimization to run
optimization.

Spike_Resonance:
    In_vitro_spres:
    spike_resonance_ext and spike_resonance_inh conduct in vitro spike resonance analysis with the same perturbation
    frequencies and number of seeds as the in vitro conditions in Guet-McCreight and Skinner (2021). Both file needs to
    be in the same folder as SINGLE dll and HOC files.

    IVL_parasearch:
    NSG_run2 searches for parameters of the stochastic synapse that falls within the constraint of IVL states, requiring
    NMODL and HOC files of SINGLE to run in NSG.

    IVL_process_pipeline processes the data produced by NSG_run2, including imposing the extra subthreshold voltage mean
    and variance constraints, calculating max and min firing rates. It produces the files "_sl3.npy" and
    "_sl3_stats.npy".

    IVL_spres:
        compare_to_FULL:
        spike_resonance_ivl_final_save conducts spiking resonance analysis with the same perturbation frequencies and
        number of seeds as the IVL conditions in Guet-McCreight and Skinner (2021). This file needs to be in the same
        folder as SINGLE dll and HOC files.

        SINGLE_explore:
            0-20ks:
            spike_resonance_ivl_final_save_1 to 4 collectively run spike resonance analysis on the 10 representative
            sets with perturbation frequencies from 0 to 30 and seeds from 2021 to 22020. The save folder "data/spres_ps5"
            is an empty folder in the zipfile uploaded to NSG. spike_resonance_ivl_single, dll and hoc files of SINGLE,
            and "data/spres_ps5" needs to be in the same folder as spike_resonance_ivl_final_save_1 to 4.

            20-40ks:
            Similar to above but with seeds 22021 to 42020
            40-120ks:
            Similar to above but with seeds 42021 to 122020

            stfr_organize organizes the results from above into numpy files with shapes that are easier to work with.
            The folder names in stfr_organize refers to the local folder names that the content of "data/spres_ps5" was
            transferred to after downloading from NSG.

STA_analysis:
    STA_re_kelowna:
    STA_re_kelowna files run STA analysis based on the resonance frequencies obtained from SINGLE_explore, specifically
    the output files of stfr_organize. In addition to the output files of stfr_organize, dll and hoc files of SINGLE
    also need to be present.

    sta_50pts_organize organize output files from STA_re_kelowna into one file containing the required number of
    ISIs.

    sta_50pts_organized_proc processes the output file of sta_50pts_organize for undirected normalization, spread,
    and slope.