import numpy
from neuron import h
from record import *
import efel
import json

from currents_visualization import *

### Instantiate Model ###
h.load_file("init_model_tweaked.hoc")
h.cvode_active(0)
h.dt = 0.1
h.steps_per_ms = 10
DNQX_Hold = 0.004
TTX_Hold = -0.0290514
cell = h.cell
h.tstop = 4000


def make_soma_traces(target_dir, amps=(0.03, 0.06, 0.09)):
    for amp in amps:
        h.ic_hold.amp = DNQX_Hold
        if amp == -0.12 or amp == -0.09:
            h.ic_step.amp = amp + (TTX_Hold - DNQX_Hold)
        else:
            h.ic_step.amp = amp
        vvec = h.Vector()
        vvec.record(h.soma[1](0.5)._ref_v)
        h.run()
        vvec = numpy.array(vvec)
        file_name = "soma_" + str(amp * 1000) + "pA.npy"
        numpy.save(target_dir + file_name, vvec)


def make_dend_traces(target_dir):
    for location in [0, 31, 43, 46, 50]:
        for amp in [-0.12, 0.03, 0.06]:
            h.ic_hold.amp = DNQX_Hold
            if amp == -0.12 or amp == -0.09:
                h.ic_step.amp = amp + (TTX_Hold - DNQX_Hold)
            else:
                h.ic_step.amp = amp
            vvec = h.Vector()
            vvec.record(h.dend[location](0.5)._ref_v)
            h.run()
            vvec = numpy.array(vvec)
            file_name = "dend_" + str(location) + "_" + str(amp * 1000) + "pA.npy"
            numpy.save(target_dir + file_name, vvec)


def extract_ef(trace_path, target_dir, trace_amps=(0.03, 0.06, 0.09), location='soma'):
    efel.api.setThreshold(-20)
    efel.api.setDerivativeThreshold(1)
    locs = [""]
    if location == 'dend':
        locs = ['_0', '_31', '_43', '_46', '_50']
    for loc in locs:
        traces = []
        for amp in trace_amps:
                v_trace = numpy.load(trace_path + "{}{}_{}pA.npy".format(location, loc, amp*1000))
                trace = {}
                trace['V'] = v_trace
                trace['T'] = numpy.arange(0, int(len(v_trace) * h.dt) + h.dt, h.dt)
                trace['stim_start'] = [1000]
                trace['stim_end'] = [3000]
                traces.append(trace)
        traces_result = efel.getFeatureValues(traces, ['AP_amplitude_from_voltagebase', 'AP_width', 'AP_amplitude',
                                                   'AHP_time_from_peak', 'time_to_first_spike', 'voltage_base',
                                                   'AP_amplitude_change', 'AP_duration_half_width', 'AHP_depth',
                                                   'mean_frequency', 'AHP_slow_time', 'adaptation_index'])

        target_path = target_dir + "{}{}_ef_strat15.json".format(location, loc)
        save_efeatures(target_path, trace_amps, traces_result)


def save_efeatures(target_path, trace_amps, features, location='soma'):
    std_list = [0.1, 0.01, 0.1, 0.1, 1, 0.1, 0.001, 0.01, 0.1, 0.1, 0.001, 0.001]
    final_form = {}
    for i in range(len(trace_amps)):
        stim = 'Step' + str(trace_amps[i])
        final_form[stim] = {}
        k = 0
        for feature in features[i]:
            if features[i][feature] is not None:
                feature_dict = {}
                feature_data = features[i][feature]
                feature_dict['Std'] = std_list[k]
                feature_dict['Mean'] = numpy.mean(feature_data)
                feature_dict['Type'] = feature
                feature_dict['Weight'] = std_list[k]
                feature_dict['Stimulus'] = stim
                final_form[stim][feature] = feature_dict

                k += 1

    with open(target_path, 'w') as outfile:
        json.dump(final_form, outfile, indent=4)


if __name__ == '__main__':
    make_dend_traces("efeatures/traces/")
    make_soma_traces('e_features/traces/')
    extract_ef('e_features/traces/', 'e_features/', location='soma')

    make_soma_traces("e_features/traces/", (-0.12,))
