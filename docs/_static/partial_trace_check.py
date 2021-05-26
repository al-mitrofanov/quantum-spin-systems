import numpy as np
from partial_trace import partial_trace_a, partial_trace_b, partial_trace_c

a = np.array([[0.707], [0.707]])  # Polarization of the first spin
b = np.array([[1], [0]])  # Polarization of the second spin
c = np.array([[0.707], [-0.707]])  # Polarization of the third spin
state = np.kron(np.kron(a, b), c)  # State with 3 spins

state_a = partial_trace_a(state)
state_b = partial_trace_b(state)
state_c = partial_trace_c(state)
