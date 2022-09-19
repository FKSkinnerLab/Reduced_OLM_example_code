The folders in this repo contains cleaned up example code that represent the original code used to produce the data in
the paper. The example code is not directly run-able since many of them requires access to a computer cluster. In
addition, the example code serves as a readable, and hopefully concise, demonstration of what was done to produce the
published result. Thus, it may not be the most optimal way of coding, and is by no means the optimal way of using
computation resources.

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

    IVL_parasearch:
    NSG_run2 searches for parameters of the stochastic synapse that falls within the constraint of IVL states, requiring
    NMODL and HOC files of SINGLE to run in NSG.

    IVL_process_pipeline processes the data produced by NSG_run2, including imposing the extra subthreshold voltage mean
    and variance constraints, calculating max and min firing rates. It produces the files "_sl3.npy" and
    "_sl3_stats.npy".

    IVL_spres:
        compare_to_FULL:
        spike_resonance_ivl_final_save conducts spiking resonance analysis with the same resonance frequencies and
        number of seeds as Guet-McCreight and Skinner (2021)

        SINGLE_explore:



STA_analysis: