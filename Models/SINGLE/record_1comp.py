from neuron import h, gui
import numpy

"""This functino sets up recording for soma and D4 location and return them in a tuple of list of vectors"""


def set_up_full_recording():

    v_vecD4 = h.Vector()
    v_vecD4.record(h.soma(0.5)._ref_v)

    IhD4 = h.Vector()
    IhD4.record(h.soma(0.5)._ref_ih_Ih)

    INaD4 = h.Vector()
    INaD4.record(h.soma(0.5)._ref_ina_Nasoma)

    IKaD4 = h.Vector()
    IKaD4.record(h.soma(0.5)._ref_ik_Ika)

    IKdrfD4= h.Vector()
    IKdrfD4.record(h.soma(0.5)._ref_ik_Ikdrf)

    ImD4 = h.Vector()
    ImD4.record(h.soma(0.5)._ref_ik_IM)

    IlD4 = h.Vector()
    IlD4.record(h.soma(0.5)._ref_i_passsd)

    recording = [v_vecD4, IKaD4, IKdrfD4, ImD4, IlD4, INaD4, IhD4]

    return recording


def set_up_spike_count():
    apc = h.APCount(h.soma(0.5))
    apc.thresh = -20
    spike_timing = h.Vector()
    apc.record(spike_timing)
    return spike_timing

# def record_soma():
#     v_vecD4 = h.Vector()
#     v_vecD4.record(h.soma(0.5)._ref_v)
#
#     IhD4 = h.Vector()
#     IhD4.record(h.soma(0.5)._ref_ih_Ih)
#
#     INaD4 = h.Vector()
#     INaD4.record(h.soma(0.5)._ref_ina_Nasoma)
#
#     IKaD4 = h.Vector()
#     IKaD4.record(h.soma(0.5)._ref_ik_Ika)
#
#     IKdrfD4= h.Vector()
#     IKdrfD4.record(h.soma(0.5)._ref_ik_Ikdrf)
#
#     ImD4 = h.Vector()
#     ImD4.record(h.soma(0.5)._ref_ik_IM)
#
#     IlD4 = h.Vector()
#     IlD4.record(h.soma(0.5)._ref_i_passsd)
#
#     recording = [v_vecD4, IKaD4, IKdrfD4, ImD4, IlD4, INaD4, IhD4]
#
#     return recording


def record_V_cai_ica_dend():
    ica = h.Vector()
    ica.record(h.soma(0.5)._ref_ica)
    v = h.Vector()
    v.record(h.soma(0.5)._ref_v)
    cai = h.Vector()
    cai.record(h.soma(0.5)._ref_cai)

    return v, cai, ica
