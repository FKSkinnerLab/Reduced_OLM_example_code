import numpy as np
import efel
from neuron import h
import multiprocessing as mlt
import matplotlib.pyplot as plt
import bluepyopt.deapext.optimisations as bpop
import bluepyopt.ephys as ephys
import json

h.load_file("nrngui.hoc")
h.celsius = 34
h.v_init = -74

somatic_loc = ephys.locations.NrnSeclistLocation('somatic', seclist_name='somatic')

gbar_na_param = ephys.parameters.NrnSectionParameter(name='gna_Nasoma', param_name='gna_Nasoma',
                                                     value=70.986184915201491e-4, locations=[somatic_loc],
                                                     frozen=True)
gbar_kdrf_param = ephys.parameters.NrnSectionParameter(name='gbar_Ikdrf', param_name='gbar_Ikdrf',
                                                       value=115.46932938891074e-4, locations=[somatic_loc],
                                                       frozen=True)
gbar_m_param = ephys.parameters.NrnSectionParameter(name='gbar_IM', param_name='gbar_IM',
                                                    value=0.13738940328219354e-4, locations=[somatic_loc],
                                                    frozen=True)
gbar_ka_param = ephys.parameters.NrnSectionParameter(name='gbar_Ika', param_name='gbar_Ika',
                                                     value=76.077776610698493e-4, locations=[somatic_loc],
                                                     frozen=True)
gbar_h_param = ephys.parameters.NrnSectionParameter(name='gkhbar_Ih', param_name='gkhbar_Ih', value=1.06309e-5,
                                                    locations=[somatic_loc], frozen=True)

FULL_GBAR_PARAMS = [gbar_na_param, gbar_kdrf_param, gbar_m_param, gbar_ka_param, gbar_h_param]
FULL_GBAR_Bounds = {'gna_Nasoma': [10e-4, 100e-4], 'gbar_Ikdrf': [3e-4, 200e-4], 'gbar_IM': [0.05e-4, 1.5e-4], 'gbar_Ika': [1.25e-4, 120e-4], 'gkhbar_Ih': [0.05e-4, 1.5e-4]}

def generate_partial_cell():
    """Create cell without gbars being set"""
    # Set up 1 compartment morphology
    morph = ephys.morphologies.NrnFileMorphology("m1053s.swc")

    # Set up somatic section_list for mechanisms to be inserted
    somatic_loc = ephys.locations.NrnSeclistLocation('somatic', seclist_name='somatic')

    # Create Mechanisms with mod files
    kdrf_mech = ephys.mechanisms.NrnMODMechanism(name='kdrf', mod_path='Ikdrf.mod', suffix='Ikdrf',
                                                 locations=[somatic_loc])
    na_mech = ephys.mechanisms.NrnMODMechanism(name='na', mod_path='Nasoma.mod', suffix='Nasoma',
                                               locations=[somatic_loc])
    m_mech = ephys.mechanisms.NrnMODMechanism(name='m', mod_path='IMmintau.mod', suffix='IM', locations=[somatic_loc])
    h_mech = ephys.mechanisms.NrnMODMechanism(name='h', mod_path='Ih.mod', suffix='Ih', locations=[somatic_loc])
    ka_mech = ephys.mechanisms.NrnMODMechanism(name='ka', mod_path='IKa.mod', suffix='Ika', locations=[somatic_loc])
    l_mech = ephys.mechanisms.NrnMODMechanism(name='l', mod_path='Ipasssd.mod', suffix='passsd',
                                              locations=[somatic_loc])

    # NaS Parameters
    vshift_na_param = ephys.parameters.NrnSectionParameter(name='vshift_na', param_name='vshift_Nasoma',
                                                           value=-4.830346371483079, locations=[somatic_loc],
                                                           frozen=True)

    # Leak Parameters
    gbar_l_param = ephys.parameters.NrnSectionParameter(name='gbar_l', param_name='g_passsd', value=7.5833e-06,
                                                        locations=[somatic_loc], frozen=True)
    erev_l_param = ephys.parameters.NrnSectionParameter(name='erev_passsd', param_name='erev_passsd', value=-64.640,
                                                        locations=[somatic_loc], frozen=True)

    # Ih parameters
    param_eh_h = ephys.parameters.NrnSectionParameter(name='eh_h', param_name='eh', value=-34.0056,
                                                      locations=[somatic_loc], frozen=True)
    param_v_half_h = ephys.parameters.NrnSectionParameter(name='v_half_h', param_name='v_half_Ih', value=-103.69,
                                                          locations=[somatic_loc], frozen=True)
    param_k_h = ephys.parameters.NrnSectionParameter(name='k_h', param_name='k_Ih', value=9.9995804,
                                                     locations=[somatic_loc], frozen=True)

    param_t1_h = ephys.parameters.NrnSectionParameter(name='t1_h', param_name='t1_Ih', value=8.5657797,
                                                      locations=[somatic_loc], frozen=True)
    param_t2_h = ephys.parameters.NrnSectionParameter(name='t2_h', param_name='t2_Ih', value=0.0296317,
                                                      locations=[somatic_loc], frozen=True)
    param_t3_h = ephys.parameters.NrnSectionParameter(name='t3_h', param_name='t3_Ih', value=-6.9145,
                                                      locations=[somatic_loc], frozen=True)
    param_t4_h = ephys.parameters.NrnSectionParameter(name='t4_h', param_name='t4_Ih', value=0.1803,
                                                      locations=[somatic_loc], frozen=True)
    param_t5_h = ephys.parameters.NrnSectionParameter(name='t5_h', param_name='t5_Ih', value=4.3566601e-05,
                                                      locations=[somatic_loc], frozen=True)

    # Global properties
    param_ena = ephys.parameters.NrnSectionParameter(name='ena', param_name='ena', value=90, locations=[somatic_loc],
                                                     frozen=True)
    param_ek = ephys.parameters.NrnSectionParameter(name='ek', param_name='ek', value=-95, locations=[somatic_loc],
                                                    frozen=True)
    param_cm = ephys.parameters.NrnSectionParameter(name='cm', param_name='cm', value=0.27008, locations=[somatic_loc],
                                                    frozen=True)

    # Instantiate cell
    m1053s = ephys.models.CellModel(name='m1053s', morph=morph,
                                    mechs=[kdrf_mech, na_mech, m_mech, h_mech, ka_mech, l_mech],
                                    params=[vshift_na_param, gbar_l_param, erev_l_param, param_eh_h,
                                            param_v_half_h, param_k_h, param_t1_h, param_t2_h, param_t3_h, param_t4_h,
                                            param_t5_h, param_ena, param_ek, param_cm])

    return m1053s


def set_up_protocol():
    """Instantiate step and hold protocol"""
    soma_loc = ephys.locations.NrnSeclistCompLocation(name='soma', seclist_name='somatic', sec_index=0, comp_x=0.5)

    protocols = []
    for amp in [0.03, 0.06, 0.09]:
        sim_name = 'Step' + str(amp)
        step = ephys.stimuli.NrnSquarePulse(step_amplitude=amp, step_delay=1000, step_duration=2000,
                                            total_duration=4000, location=soma_loc)
        hold = ephys.stimuli.NrnSquarePulse(step_amplitude=0.004, step_delay=0, step_duration=4000, total_duration=4000,
                                            location=soma_loc)
        rec = ephys.recordings.CompRecording(name=sim_name + '.soma.v', location=soma_loc, variable='v')
        step_protocol = ephys.protocols.SweepProtocol(name=sim_name, stimuli=[step, hold], recordings=[rec],
                                                      cvode_active=True)
        protocols.append(step_protocol)

    full_protocol = ephys.protocols.SequenceProtocol('three_steps', protocols=protocols)
    return full_protocol


def create_evaluator(efeature_path, full_cell_model, full_protocol, opt_param_names):
    """Take a set of eafatures from json file and generate objective"""
    with open(efeature_path, "r") as outfile:
        efel_dict = json.load(outfile)
    objectives = []
    for trace_name in efel_dict:
        content = efel_dict[trace_name]
        features = []
        weights = []
        for feature_type in content:
            name = feature_type + "." + trace_name
            feature = ephys.efeatures.eFELFeature(name, efel_feature_name=feature_type,
                                                  recording_names={'': trace_name + ".soma.v"},
                                                  stim_start=1000, stim_end=3000,
                                                  exp_mean=content[feature_type]['Mean'],
                                                  exp_std=content[feature_type]['Weight'])
            # features.append(feature)
            # weights.append(content[feature_type]['Weight'])
            # objective = ephys.objectives.WeightedSumObjective(trace_name, features, weights)
            objective = ephys.objectives.SingletonObjective(name, feature)
            objectives.append(objective)
    calculator = ephys.objectivescalculators.ObjectivesCalculator(objectives)
    evaluator = ephys.evaluators.CellEvaluator(cell_model=full_cell_model, param_names=opt_param_names,
                                               fitness_protocols={full_protocol.name: full_protocol},
                                               fitness_calculator=calculator, sim=ephys.simulators.NrnSimulator())

    return evaluator


def single_optimize(full_gbar_list: list, gbar_bounds: dict, target_gbar_names: list):
    cell_model = generate_partial_cell()
    opt_param_names = []
    for gbar_param in full_gbar_list:
        if gbar_param.name not in target_gbar_names:
            cell_model.params[gbar_param.name] = gbar_param
        else:
            variable_param = ephys.parameters.NrnSectionParameter(gbar_param.name, param_name=gbar_param.param_name,
                                                                  bounds=gbar_bounds[gbar_param.name], locations=gbar_param.locations,
                                                                  frozen=False)
            cell_model.params[gbar_param.name] = variable_param
            opt_param_names.append(variable_param.name)

    full_protocol = set_up_protocol()
    cell_evaluator = create_evaluator("soma_ef_strat15.json", cell_model, full_protocol, opt_param_names)
    # with mlt.Manager() as manager:
    #     map = manager.dict
    optimisation = bpop.DEAPOptimisation(evaluator=cell_evaluator, offspring_size=100, mutpb=0.15, cxpb=0.85, selector_name='IBEA', eta=0.5)
    final_pop, hall_of_fame, logs, hist = optimisation.run(max_ngen=200)
    for i in range(len(opt_param_names)):
        print([opt_param_names[i], str(hall_of_fame[0][i])])


def multi_optimize():
    processes = []
    for gbar_param in FULL_GBAR_PARAMS:
        opt_run = mlt.Process(target=single_optimize, args=(FULL_GBAR_PARAMS, FULL_GBAR_Bounds, [gbar_param.name]))
        processes.append(opt_run)
    for process in processes:
        process.start()
    for process in processes:
        process.join()


if __name__ == "__main__":
    mlt.freeze_support()
    # multi_optimize()
    single_optimize(FULL_GBAR_PARAMS, FULL_GBAR_Bounds, [gbar.name for gbar in FULL_GBAR_PARAMS[:-1]])
