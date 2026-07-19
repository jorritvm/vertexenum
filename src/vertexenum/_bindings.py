"""ctypes bindings for the compiled lrslib shim (_vertexenum shared library)."""

from __future__ import annotations

import ctypes
from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().parent

# Status codes returned by ve_enumerate (keep in sync with pyvertexenum.c)
VE_OK = 0
VE_ERR_INIT = 1
VE_ERR_NOMEM = 2
VE_ERR_NOBASIS = 3
VE_ERR_BADINPUT = 4

# lrslib works in C ``long`` throughout: 32-bit on Windows, usually 64-bit
# elsewhere.  Inputs passed to the library must fit in this type.
C_LONG_MIN = -(2 ** (8 * ctypes.sizeof(ctypes.c_long) - 1))
C_LONG_MAX = -C_LONG_MIN - 1


def _find_library() -> Path:
    for pattern in ("_vertexenum*.dll", "_vertexenum*.so", "_vertexenum*.dylib"):
        matches = sorted(_PACKAGE_DIR.glob(pattern))
        if matches:
            return matches[0]
    raise OSError(
        f"could not find the compiled _vertexenum shared library in {_PACKAGE_DIR}; "
        "the package's C extension was not built. Reinstall the package "
        "(a C compiler such as MinGW-w64 gcc is required to build it)."
    )


def _load() -> ctypes.CDLL:
    lib = ctypes.CDLL(str(_find_library()))

    lib.ve_enumerate.argtypes = [
        ctypes.c_long,                    # m
        ctypes.c_long,                    # n
        ctypes.POINTER(ctypes.c_long),    # A_num, row-major m*n
        ctypes.POINTER(ctypes.c_long),    # A_den
        ctypes.POINTER(ctypes.c_long),    # b_num, length m
        ctypes.POINTER(ctypes.c_long),    # b_den
        ctypes.POINTER(ctypes.c_long),    # out_rows
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),  # out_data
    ]
    lib.ve_enumerate.restype = ctypes.c_int

    lib.ve_free.argtypes = [ctypes.POINTER(ctypes.c_double)]
    lib.ve_free.restype = None

    lib.ve_lrs_version.argtypes = []
    lib.ve_lrs_version.restype = ctypes.c_char_p

    return lib


_lib: ctypes.CDLL | None = None


def get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is None:
        _lib = _load()
    return _lib


def lrs_version() -> str:
    """Version string of the bundled lrslib, e.g. ``'lrslib v.5.1 2015.1.28'``."""
    return get_lib().ve_lrs_version().decode("ascii")
