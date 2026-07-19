# Developer documentation

This document explains how the `vertexenum` package is put together: what
lives in the C layer, what lives in the Python layer, how the build works,
and what to watch out for when changing any of it. For usage, see the
top-level `README.md`.

## Overview

```
                Python                              C (compiled to _vertexenum.dll / .so)
 ┌────────────────────────────────┐      ┌─────────────────────────────────────────┐
 │ vertexenum.enumerate_vertices  │      │ pyvertexenum.c        ve_enumerate()    │
 │   validate input               │      │   build lrs H-representation            │
 │   entries -> exact Fractions   │─────▶│   lrs_getfirstbasis / lrs_getsolution / │
 │   Fractions -> c_long num/den  │ctypes│   lrs_getnextbasis (reverse search)     │
 │   result -> numpy, drop rays   │◀─────│   collect rows, convert to doubles      │
 └────────────────────────────────┘      │ lrslib.c + lrsmp.c (vendored, v5.1)     │
                                         └─────────────────────────────────────────┘
```

The package answers one question: given `A x <= b`, what are the vertices of
the polyhedron? The algorithm (lexicographic reverse search, exact rational
arithmetic) is entirely David Avis's lrslib. Everything in this repository is
plumbing that gets matrices into lrslib and vertices back out.

The design follows the R package `vertexenum` (Robert Robere, in
`r_package/`), which this project ports: a thin C shim over lrslib plus a
thin high-level wrapper. Where the R original and this port differ, it is
called out below.

## Repository layout

| Path | Contents |
|---|---|
| `src/vertexenum/csrc/` | All C sources: vendored lrslib + our shim |
| `src/vertexenum/_bindings.py` | ctypes: locates and loads the shared library, declares signatures |
| `src/vertexenum/__init__.py` | Public API: `enumerate_vertices`, `lrs_version`, exceptions |
| `setup.py` | Compiles the C sources into a plain shared library |
| `pyproject.toml` | Package metadata, setuptools config, dependencies |
| `tests/test_vertexenum.py` | pytest suite (known polytopes, edge cases, validation) |

## The C layer

### Vendored lrslib (`lrslib.c/h`, `lrsmp.c/h`)

These four files should be treated as read-only third-party code. Facts that matter:

- **Version**: lrslib v5.1 (2015.1.28) — the exact sources the 2018 R package
  shipped, per the project's v1 goal of matching the R package's behavior.
- **Arithmetic**: lrs supports three arithmetic backends selected by
  preprocessor defines (`LONG`, `GMP`, or the default `lrsmp`). We compile
  with no defines, which selects **lrsmp**, lrs's own arbitrary-precision
  integer package (`lrsmp.c`). That is what makes the build self-contained —
  no GMP dependency.
- **Already modified by the R author.** These are not pristine upstream
  sources. Robert Robere commented out everything unsuitable for running
  inside a host process: all `printf`/`fprintf` on the computation paths, all
  `exit()` calls, and the signal handlers. `lrs_init` calls
  `lrs_mp_init(ZERO, NULL, NULL)`, so the global `lrs_ifp`/`lrs_ofp` file
  pointers end up **NULL** after init (see the null-device note below).
- **Silent overflow caveat**: in pristine lrslib, running out of
  multiprecision digits calls `digits_overflow()`, which prints a diagnostic
  and exits. In these sources that function is an empty stub, so a
  computation that exceeds `MAX_DIGITS` would proceed **silently with wrong
  results** rather than fail. This is inherited from the R package; in
  practice lrsmp's default digit budget is far beyond what typical inputs
  need, but it is the known sharp edge of v1.
- **`long` everywhere**: lrslib's whole API is in terms of C `long` — 32-bit
  on Windows, 64-bit on most other platforms. This is why the Python layer
  range-checks every numerator and denominator against the platform's
  `c_long` limits.
- **Global state**: lrslib keeps globals (`lrs_global_list`, file pointers,
  the arithmetic package's digit counters). Consequently the library is
  **not thread-safe** and `ve_enumerate` must not run concurrently. The
  globals are all defined inside `lrslib.c`/`lrsmp.c`; the driver does not
  need to define any (the globals at the top of the R shim were vestigial and
  were dropped from ours).

### The shim (`pyvertexenum.c`)

Everything we wrote in C is in this one file. It exists because lrslib's API
is hostile to direct ctypes binding: `lrs_mp` is a typedef for
`long[MAX_DIGITS+1]`, vectors/matrices are nested pointer types, and core
operations (`copy`, `zero`, ...) are C **macros**, which don't exist in a
compiled DLL at all. So, like the R package before us, we keep the lrslib
interaction in C and expose one flat, ctypes-friendly entry point.

Exported functions (`VE_EXPORT` = `__declspec(dllexport)` on Windows):

```c
int  ve_enumerate(long m, long n,
                  const long *A_num, const long *A_den,   /* row-major m*n */
                  const long *b_num, const long *b_den,   /* length m */
                  long *out_rows, double **out_data);
void ve_free(double *data);
const char *ve_lrs_version(void);
```

Semantics of `ve_enumerate`:

1. **Input** is the system `A x <= b` as exact rationals: entry (i, j) of A
   is `A_num[i*n + j] / A_den[i*n + j]`. The shim validates dimensions and
   nonzero denominators (`VE_ERR_BADINPUT`).
2. **Conversion to lrs's H-representation.** lrs wants each constraint as
   `b_i + (-A_i)·x >= 0`, i.e. row `[b_i, -A_i1, ..., -A_in]`. The shim does
   the negation itself, so its API keeps the natural `A x <= b` orientation.
   (In the R package this negation happened in R code instead.) Rows are fed
   in with `lrs_set_row(P, Q, row, num, den, GE)`.
3. **Reverse search loop** — the same loop as `Rvertexenum.c` and lrs's own
   demo drivers: `lrs_getfirstbasis`, then repeatedly `lrs_getsolution` for
   each column `0..P->d` and `lrs_getnextbasis` until the tree is exhausted.
   Output rows accumulate in an `lrs_mp_matrix` that starts at 8 rows and
   doubles (allocate new, `copy()` each entry, free old — lrs matrices can't
   be realloc'd).
4. **Output conversion.** Each collected row is an lrs output line of length
   `n+1`. lrs's encoding: if entry 0 is zero, the row is an **extreme ray**
   (entries are integers, converted with `mptodouble`); otherwise it is a
   **vertex** and entry 0 is the common denominator (entries converted with
   `rattodouble(entry, row[0])`, which makes column 0 itself come out as
   exactly 1.0). So in the returned matrix, column 0 is a 1.0/0.0
   vertex-vs-ray flag and columns 1..n are coordinates. This encoding is
   preserved as-is for the Python layer to interpret.
5. **Result ownership**: the doubles matrix is `malloc`'d in the shim and
   handed to the caller, who must release it with `ve_free`. (Python and the
   DLL may use different C runtimes, so the buffer must be freed by the same
   side that allocated it — never `free()` it from anywhere else.)

Status codes (mirrored as constants in `_bindings.py` — keep in sync):

| Code | Meaning |
|---|---|
| `VE_OK` (0) | success |
| `VE_ERR_INIT` (1) | `lrs_init`/`lrs_alloc_dat`/`lrs_alloc_dic` failed |
| `VE_ERR_NOMEM` (2) | allocation failed inside the shim |
| `VE_ERR_NOBASIS` (3) | `lrs_getfirstbasis` failed — typically the system is infeasible |
| `VE_ERR_BADINPUT` (4) | bad dimensions, NULL pointer, or zero denominator |

Two defensive details worth knowing:

- **Null device**: right after `lrs_init`, the shim opens `NUL` (Windows) /
  `/dev/null` and assigns it to `lrs_ifp`/`lrs_ofp`. The vendored sources
  leave those NULL, and while the hot paths never print, a few cold paths
  still contain `fprintf(lrs_ofp, ...)`; this guarantees none of them can
  dereference NULL inside the host Python process.
- **Cleanup discipline**: all error paths funnel through one `done:` block
  that frees whatever was allocated, in the right order (the lrs matrix is
  cleared before `Q` is freed, because clearing needs `Q->n`). Unlike the R
  original, no path leaks and no path can `error()`/abort the host process.

Improvements over `Rvertexenum.c`, for anyone diffing the two: real status
codes instead of `error()`/NULL returns, complete cleanup on every path, the
null-device guard, negation moved into C, and removal of dead code (the R
shim's unused globals and commented-out blocks).

## The Python layer

### `_bindings.py` — loading and declaring the DLL

- Finds the shared library by globbing `_vertexenum*.dll` / `.so` / `.dylib`
  **next to the module itself** (works for both editable installs, where the
  build drops it in `src/vertexenum/`, and wheels, where it sits in
  site-packages). A missing library raises an `OSError` explaining that the
  extension wasn't built.
- Loads it lazily via `get_lib()` (module import never triggers a DLL load —
  useful for tooling that imports the package without calling it).
- Declares ctypes signatures for the three exports and defines
  `C_LONG_MIN`/`C_LONG_MAX` from `ctypes.sizeof(c_long)`, so range checks are
  automatically correct for 32-bit-long platforms (Windows) and 64-bit ones.
- Re-declares the `VE_*` status codes; these must match `pyvertexenum.c`.

### `__init__.py` — the public API

`enumerate_vertices(A, b, *, max_denominator=1_000_000)` in stages:

1. **Validation**: `A` must be 2-D with at least one row and one column, `b`
   must match. (The R package required at least 2×2; that restriction was
   dropped — 1-D systems like `2 <= x <= 5` work and are tested.)
2. **Exact rationalization** (`_to_fraction`): integers and `Rational`s
   (`Fraction`, `Decimal`) pass through exactly; anything else goes through
   `float()` → `Fraction` → `limit_denominator(max_denominator)`. This
   replaces the R package's Farey-sequence approximation (`numbers::ratFarey`
   with a hard-coded limit of 10000) with the stdlib equivalent and a
   configurable, larger default.
3. **Range check** (`_check_range`): every numerator/denominator must fit in
   the platform C `long`, else `ValueError` with a pointer at
   `max_denominator`. Without this, values would be silently truncated at
   the ctypes boundary — the single easiest way to get wrong answers.
4. **Call** `ve_enumerate` with flat `c_long` arrays; map status codes to
   exceptions (`LrsError` for lrs failures, with an "infeasible" message for
   `VE_ERR_NOBASIS`; `ValueError` for `VE_ERR_BADINPUT`).
5. **Result handling**: wrap the returned buffer with
   `np.ctypeslib.as_array(...).copy()`, then `ve_free` it — the `.copy()` is
   load-bearing, since `as_array` only *views* C-owned memory. Rows with
   flag 0 (extreme rays, meaning the polyhedron is unbounded) are dropped
   with an `UnboundedWarning`; the flag column is stripped and a `(k, n)`
   float64 array of vertices is returned. Zero output rows yield an empty
   `(0, n)` array.

Public names: `enumerate_vertices`, `lrs_version`, `LrsError`,
`UnboundedWarning`.

## Build & tooling

### How the shared library is built

The extension is **not** a CPython extension module — there is no
`PyInit__vertexenum`, no `#include <Python.h>`; it is a plain C shared
library loaded with `ctypes.CDLL`. We still build it through setuptools so
that `pip install` / `uv sync` compile it automatically. That requires the
two standard overrides in `setup.py` (`build_ctypes_ext`):

- `get_export_symbols` — don't force-export a `PyInit_*` symbol (it doesn't
  exist; the default would make linking fail on Windows).
- `get_ext_filename` — name the artifact plain `_vertexenum.dll` / `.so`
  instead of tagging it with the CPython ABI (`.cp313-win_amd64.pyd`), since
  the binary is Python-version-independent.

Compilation covers `lrsmp.c`, `lrslib.c`, `pyvertexenum.c` with
`src/vertexenum/csrc` on the include path and **no preprocessor defines**
(that's what selects lrsmp arithmetic and keeps `PLRS`/`SIGNALS`/`TIMES`
code excluded).

Windows specifics, also in `setup.py`:

- `finalize_options` forces `compiler = "mingw32"` — we build with MinGW-w64
  gcc, the toolchain the R package used on CRAN, rather than MSVC (which
  would additionally choke on setuptools' default `/export:PyInit_...`).
  Install it with `winget install BrechtSanders.WinLibs.POSIX.UCRT` if
  needed; gcc must be on `PATH` when the build backend runs.
- `extra_link_args=["-static"]` — the DLL statically links the MinGW
  runtime, so importing the package never depends on MinGW being installed
  or on `PATH` for end users of a built wheel.

### Packaging config (`pyproject.toml`)

- Build backend: `setuptools.build_meta`, with `setup.py` present only for
  the extension logic — all metadata lives in `pyproject.toml`.
- src layout: `package-dir = {"" = "src"}`, single package `vertexenum`.
- The C sources are declared as package data so sdists are buildable.
- Runtime dependency: numpy only. Dev dependency group: pytest.
- License: `GPL-2.0-or-later`. This is not a choice — it is inherited from
  lrslib and the R shim the C code derives from, and it propagates to
  anything that bundles this package.

### Day-to-day commands

```
uv sync                # resolve deps, compile the extension, install editable
uv run pytest          # run the test suite
uv build               # build sdist + wheel (wheel is platform-specific)
```

`uv sync` installs the project in editable mode; setuptools' editable build
runs `build_ext` and drops `_vertexenum.dll` directly into `src/vertexenum/`
(gitignored via `*.dll`). **After editing any `.c` file you must re-run
`uv sync`** (or `uv pip install -e . --reinstall-package vertexenum`) —
Python-side edits are picked up immediately, C-side edits are not.

For fast C iteration you can bypass setuptools entirely; the module loader
only globs for the DLL, so a hand-built one works:

```
gcc -O2 -shared -static -Isrc/vertexenum/csrc -o src/vertexenum/_vertexenum.dll \
    src/vertexenum/csrc/lrsmp.c src/vertexenum/csrc/lrslib.c src/vertexenum/csrc/pyvertexenum.c
```

(Drop `-static` on Linux/macOS and use `-o ..._vertexenum.so`; add `-fPIC`.)

### Tests

`tests/test_vertexenum.py` checks polytopes with known vertex sets (unit
square, 3-D cube, simplex, triangles with fractional vertices, the example
from the R package's documentation), the unbounded / infeasible / 1-D edge
cases, float-vs-Fraction agreement, input validation, and the output
contract (float64 ndarray, correct shape). `assert_vertex_sets_equal`
compares vertex sets order-independently, since lrs's output order is an
implementation detail. When adding features, prefer adding a polytope whose
vertex set you can state exactly.

## Known limitations and future directions

Tracked from the original project plan:

- **Exact rational output**: the shim already holds exact `lrs_mp` rationals
  right up until the final doubles conversion; a v2 could return
  numerator/denominator pairs (e.g. as strings, since values can exceed any
  fixed-width integer) and surface `Fraction` vertices in Python.
- **Upgrade to current lrslib** (`lrslib-073/`): newer upstream has a
  reentrant API (`lrs_main`-style drivers, per-call state), 64-bit/GMP/hybrid
  arithmetic with overflow *detection*, and `mplrs` for parallelism. An
  upgrade would fix both the silent-overflow caveat and the thread-safety
  limitation, at the cost of redoing the shim against a changed API.
- **Thread safety / multiprocessing**: until then, if you need parallelism,
  use processes, not threads.
- **Input magnitude**: bounded by C `long` (32-bit on Windows). Documented
  in the README; enforced in `_check_range`.
