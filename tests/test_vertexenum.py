import warnings
from fractions import Fraction

import numpy as np
import pytest

from vertexenum import LrsError, UnboundedWarning, enumerate_vertices, lrs_version


def assert_vertex_sets_equal(actual, expected, tol=1e-9):
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    assert actual.shape == expected.shape, f"got {actual.shape}, want {expected.shape}"
    remaining = list(range(expected.shape[0]))
    for row in actual:
        hit = next(
            (k for k in remaining if np.allclose(row, expected[k], atol=tol)), None
        )
        assert hit is not None, f"unexpected vertex {row}; expected one of {expected}"
        remaining.remove(hit)
    assert not remaining


def test_lrs_version():
    assert lrs_version() == "lrslib v.5.1 2015.1.28"


def test_unit_square():
    A = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    b = [1, 0, 1, 0]
    verts = enumerate_vertices(A, b)
    assert_vertex_sets_equal(verts, [[0, 0], [1, 0], [0, 1], [1, 1]])


def test_unit_cube_3d():
    A = np.vstack([np.eye(3), -np.eye(3)])
    b = [1, 1, 1, 0, 0, 0]
    verts = enumerate_vertices(A, b)
    expected = [[x, y, z] for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    assert_vertex_sets_equal(verts, expected)


def test_simplex():
    # x, y, z >= 0, x + y + z <= 1
    A = np.vstack([-np.eye(3), np.ones((1, 3))])
    b = [0, 0, 0, 1]
    verts = enumerate_vertices(A, b)
    assert_vertex_sets_equal(verts, [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])


def test_fractional_vertices():
    # triangle with vertices (0,0), (1/2, 0), (0, 1/3): 2x + 3y <= 1, x,y >= 0
    A = [[2, 3], [-1, 0], [0, -1]]
    b = [1, 0, 0]
    verts = enumerate_vertices(A, b)
    assert_vertex_sets_equal(verts, [[0, 0], [0.5, 0], [0, 1 / 3]])


def test_float_and_fraction_inputs_agree():
    A_frac = [[Fraction(1, 2), Fraction(1, 3)], [-1, 0], [0, -1]]
    A_float = [[0.5, 1 / 3], [-1, 0], [0, -1]]
    b = [1, 0, 0]
    v1 = enumerate_vertices(A_frac, b)
    v2 = enumerate_vertices(A_float, b, max_denominator=10**6)
    assert v1.shape == v2.shape
    assert_vertex_sets_equal(v2, v1, tol=1e-5)


def test_r_package_documentation_example():
    # From the R package docs: system with vertices (0,0), (0,2), (2,0), (4/3, 4/3)
    A = [[-1, 0], [0, -1], [1, 2], [2, 1]]
    b = [0, 0, 4, 4]
    verts = enumerate_vertices(A, b)
    assert_vertex_sets_equal(verts, [[0, 0], [0, 2], [2, 0], [4 / 3, 4 / 3]])


def test_unbounded_warns_and_drops_rays():
    # first quadrant: only vertex is the origin, two extreme rays
    A = [[-1, 0], [0, -1]]
    b = [0, 0]
    with pytest.warns(UnboundedWarning):
        verts = enumerate_vertices(A, b)
    assert_vertex_sets_equal(verts, [[0, 0]])


def test_bounded_does_not_warn():
    A = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    b = [1, 0, 1, 0]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        enumerate_vertices(A, b)


def test_infeasible_raises():
    # x <= 0 and x >= 1
    A = [[1], [-1]]
    b = [0, -1]
    with pytest.raises(LrsError, match="infeasible"):
        enumerate_vertices(A, b)


def test_one_dimensional_interval():
    # 2 <= x <= 5
    A = [[1], [-1]]
    b = [5, -2]
    verts = enumerate_vertices(A, b)
    assert_vertex_sets_equal(verts, [[2], [5]])


def test_input_validation():
    with pytest.raises(ValueError, match="2-D"):
        enumerate_vertices([1, 2, 3], [1])
    with pytest.raises(ValueError, match="rows"):
        enumerate_vertices([[1, 0], [0, 1]], [1, 1, 1])
    with pytest.raises(ValueError, match="max_denominator"):
        enumerate_vertices([[1.0]], [1.0], max_denominator=0)
    with pytest.raises(ValueError, match="C long"):
        enumerate_vertices([[2**80]], [1])


def test_returns_float64_ndarray():
    verts = enumerate_vertices([[1], [-1]], [1, 0])
    assert isinstance(verts, np.ndarray)
    assert verts.dtype == np.float64
    assert verts.shape == (2, 1)
