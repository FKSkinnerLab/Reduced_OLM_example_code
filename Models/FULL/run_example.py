import numpy as np
from neuron import h, gui
from record import *
import matplotlib.pyplot as plt

from currents_visualization import *

### Instantiate Model ###
h.load_file("init_model_tweaked.hoc")
h.cvode_active(0)
h.dt = 0.1
h.steps_per_ms = 10
# cell = h.cell
h.ic_hold.amp = 0

h.ic_step.amp = 0.03        # in nA
h.ic_step.delay = 1000      # delay of step injection in ms
h.ic_step.dur = 2000        # duration of step injection

h.tstop = 4000              # duration of entire simulation in ms

recordings = set_up_full_recording(0)  # see file record.py for detail

h.run()

recordings = [np.array(i) for i in recordings]  # convert all HOC arrays to numpy arrays

v_vec = recordings[0]
t_vec = np.array(range(len(v_vec))) * h.dt

plt.plot(t_vec, v_vec)
plt.show()

fig = plotCurrentscape(recordings[0], recordings[1:])
plt.show()
