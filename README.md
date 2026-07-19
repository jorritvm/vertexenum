# vertexenum

Vertex enumeration of convex polyhedra from Python: given a system of linear
inequalities `A x <= b`, compute the vertices of the polyhedron
`{ x : A x <= b }`.

This is a Python port of the R package
[vertexenum](https://cran.r-project.org/package=vertexenum) by Robert Robere.
The actual enumeration is done by [lrslib](https://cgm.cs.mcgill.ca/~avis/C/lrs.html)
(David Avis), which uses exact rational arithmetic and reverse search. The
same lrslib sources bundled with the R package (v5.1, using lrs's built-in
`lrsmp` multiprecision arithmetic, no GMP needed) are vendored in
`src/vertexenum/csrc/`, compiled into a small shared library, and accessed
via `ctypes`.

## Usage

```python
import numpy as np
from vertexenum import enumerate_vertices

# unit square: x <= 1, -x <= 0, y <= 1, -y <= 0
A = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])
b = np.array([1, 0, 1, 0])

enumerate_vertices(A, b)
# array([[0., 0.],
#        [1., 0.],
#        [0., 1.],
#        [1., 1.]])
```

- Returns a `(k, n)` float64 array with one vertex per row.
- Entries of `A` and `b` may be ints, floats, or `fractions.Fraction`.
  lrs computes in exact rational arithmetic, so floats are approximated by
  rationals with denominator at most `max_denominator` (default 1,000,000);
  ints and Fractions are used exactly.
- If the polyhedron is unbounded, its extreme rays are dropped and a
  `vertexenum.UnboundedWarning` is issued.
- An infeasible system (empty polyhedron) raises `vertexenum.LrsError`.
- `vertexenum.lrs_version()` reports the bundled lrslib version.

## Installation / building

The package builds the C sources at install time, so a C compiler is
required. On Windows use MinGW-w64 gcc (e.g. `winget install
BrechtSanders.WinLibs.POSIX.UCRT`); on Linux/macOS any gcc/clang works.

```
uv sync          # or: pip install .
uv run pytest    # run the tests
```

## Limitations (v1)

- Inputs must fit in the C `long` type lrslib uses internally (32-bit on
  Windows); out-of-range values raise `ValueError`.
- Not thread-safe: lrslib keeps global state, so don't call
  `enumerate_vertices` from multiple threads concurrently.
- Output is float64 (matching the R package); exact rational output is a
  possible future addition.
- The bundled lrslib is the 2015 version shipped with the R package;
  migrating to current lrslib (see `lrslib-073/`) is future work.

## License

GPL-2.0-or-later, inherited from lrslib and the R package's C interface.

## Author 
Jorrit Vander Mynsbrugge  
Written with Fable