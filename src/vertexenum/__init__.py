"""vertexenum: vertex enumeration of convex polyhedra via lrslib.

Given a system of linear inequalities ``A x <= b``, :func:`enumerate_vertices`
returns the vertices of the polyhedron ``{ x : A x <= b }``.  The heavy
lifting is done by David Avis's lrslib (reverse search), bundled and compiled
into this package; the interface mirrors the R package ``vertexenum`` by
Robert Robere.
"""

from __future__ import annotations

import ctypes
import warnings
from fractions import Fraction
from numbers import Integral, Rational

import numpy as np

from ._bindings import (
    C_LONG_MAX,
    C_LONG_MIN,
    VE_ERR_BADINPUT,
    VE_ERR_INIT,
    VE_ERR_NOBASIS,
    VE_ERR_NOMEM,
    get_lib,
    lrs_version,
)

__all__ = ["enumerate_vertices", "lrs_version", "LrsError", "UnboundedWarning"]

__version__ = "0.1.0"


class LrsError(RuntimeError):
    """Raised when the underlying lrslib computation fails."""


class UnboundedWarning(UserWarning):
    """Issued when the polyhedron is unbounded (extreme rays were dropped)."""


def _to_fraction(x, max_denominator: int) -> Fraction:
    """Convert a matrix entry to an exact rational.

    Integers and rationals (``Fraction``, ``Decimal``) convert exactly; floats
    are approximated with denominator at most ``max_denominator``.
    """
    if isinstance(x, Integral):
        return Fraction(int(x))
    if isinstance(x, Rational):
        return Fraction(x)
    f = Fraction(float(x))
    if f.denominator > max_denominator:
        f = f.limit_denominator(max_denominator)
    return f


def _check_range(f: Fraction, what: str) -> Fraction:
    if not (C_LONG_MIN < f.numerator <= C_LONG_MAX) or f.denominator > C_LONG_MAX:
        raise ValueError(
            f"{what} = {f} does not fit in the C long type used by lrslib "
            f"(range [{C_LONG_MIN}, {C_LONG_MAX}]); use smaller values or a "
            "smaller max_denominator"
        )
    return f


def enumerate_vertices(A, b, *, max_denominator: int = 1_000_000) -> np.ndarray:
    """Enumerate the vertices of the polyhedron ``{ x : A x <= b }``.

    Parameters
    ----------
    A : (m, n) array-like
        Coefficient matrix of the inequality system.  Entries may be ints,
        floats or :class:`fractions.Fraction`.
    b : (m,) array-like
        Right-hand side vector.
    max_denominator : int, optional
        lrslib computes in exact rational arithmetic, so float entries are
        approximated by the closest rational with denominator at most this
        value (exact inputs — ints and Fractions — are never altered).

    Returns
    -------
    numpy.ndarray of shape (k, n), dtype float64
        One row per vertex.  If the polyhedron is unbounded, its extreme rays
        are dropped from the output and an :class:`UnboundedWarning` is
        issued.  An empty (0, n) array means the system has no vertices.

    Raises
    ------
    LrsError
        If lrslib cannot find a starting basis — typically the system is
        infeasible (the polyhedron is empty).
    ValueError
        For malformed input or entries too large for lrslib's C long type.

    Examples
    --------
    Unit square ``0 <= x, y <= 1``:

    >>> A = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    >>> b = [1, 0, 1, 0]
    >>> enumerate_vertices(A, b).shape
    (4, 2)
    """
    A = np.asarray(A, dtype=object)
    b = np.asarray(b, dtype=object).reshape(-1)
    if A.ndim != 2:
        raise ValueError(f"A must be a 2-D matrix, got shape {A.shape}")
    m, n = A.shape
    if m < 1 or n < 1:
        raise ValueError(f"A must have at least one row and one column, got shape {A.shape}")
    if b.shape[0] != m:
        raise ValueError(f"b has {b.shape[0]} entries but A has {m} rows")
    if max_denominator < 1:
        raise ValueError("max_denominator must be >= 1")

    A_num = (ctypes.c_long * (m * n))()
    A_den = (ctypes.c_long * (m * n))()
    b_num = (ctypes.c_long * m)()
    b_den = (ctypes.c_long * m)()
    for i in range(m):
        f = _check_range(_to_fraction(b[i], max_denominator), f"b[{i}]")
        b_num[i] = f.numerator
        b_den[i] = f.denominator
        for j in range(n):
            f = _check_range(_to_fraction(A[i, j], max_denominator), f"A[{i},{j}]")
            A_num[i * n + j] = f.numerator
            A_den[i * n + j] = f.denominator

    lib = get_lib()
    out_rows = ctypes.c_long(0)
    out_data = ctypes.POINTER(ctypes.c_double)()
    status = lib.ve_enumerate(m, n, A_num, A_den, b_num, b_den,
                              ctypes.byref(out_rows), ctypes.byref(out_data))

    if status == VE_ERR_NOBASIS:
        raise LrsError(
            "lrslib could not find a starting basis; the system is likely "
            "infeasible (the polyhedron is empty)"
        )
    if status in (VE_ERR_INIT, VE_ERR_NOMEM):
        raise LrsError(f"lrslib initialization or allocation failed (status {status})")
    if status == VE_ERR_BADINPUT:
        raise ValueError("invalid input passed to lrslib (zero denominator or bad dimensions)")
    if status != 0:
        raise LrsError(f"unexpected lrslib status {status}")

    rows = out_rows.value
    if rows == 0:
        return np.empty((0, n), dtype=np.float64)
    try:
        raw = np.ctypeslib.as_array(out_data, shape=(rows, n + 1)).copy()
    finally:
        lib.ve_free(out_data)

    # column 0 flags each output row: 1.0 = vertex, 0.0 = extreme ray
    is_vertex = raw[:, 0] != 0.0
    if not is_vertex.all():
        warnings.warn(
            "the polyhedron is unbounded; its extreme rays were dropped from "
            "the output",
            UnboundedWarning,
            stacklevel=2,
        )
    return raw[is_vertex, 1:]
