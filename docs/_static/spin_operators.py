import numpy as np

UP = np.array([1, 0])
DOWN = np.array([0, 1])
LEFT = 1/np.sqrt(2) * np.array([1, 1])
RIGHT = 1/np.sqrt(2) * np.array([1, -1])
INSIDE = 1/np.sqrt(2) * np.array([1, 1j])
OUTSIDE = 1/np.sqrt(2) * np.array([1, -1j])

SX = 0.5 * np.array([(0, 1), (1, 0)])
SY = 0.5 * np.array([(0, -1j), (1j, 0)])
SZ = 0.5 * np.array([(1, 0), (0, -1)])

print(SX @ UP, SZ @ DOWN)
