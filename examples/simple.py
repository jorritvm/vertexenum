"""The README example: vertices of the unit square 0 <= x, y <= 1."""

import numpy as np

from vertexenum import enumerate_vertices

# unit square: x <= 1, -x <= 0, y <= 1, -y <= 0
A = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])
b = np.array([1, 0, 1, 0])

vertices = enumerate_vertices(A, b)
print(vertices)
