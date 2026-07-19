"""Build the bundled lrslib + shim into a plain shared library loaded via ctypes.

The extension is not a CPython extension module (no PyInit_*), so we override
the parts of build_ext that assume one: no PyInit export symbol, and a plain
.dll/.so filename without the CPython ABI tag.
"""

import os

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

CSRC = os.path.join("src", "vertexenum", "csrc")


class CTypesExtension(Extension):
    pass


class build_ctypes_ext(build_ext):
    def finalize_options(self):
        super().finalize_options()
        # On Windows we build with MinGW-w64 gcc (the toolchain lrslib and the
        # original R package were developed against), not MSVC.
        if os.name == "nt" and self.compiler is None:
            self.compiler = "mingw32"

    def get_export_symbols(self, ext):
        if isinstance(ext, CTypesExtension):
            return ext.export_symbols
        return super().get_export_symbols(ext)

    def get_ext_filename(self, ext_name):
        if ext_name.split(".")[-1] == "_vertexenum":
            suffix = ".dll" if os.name == "nt" else ".so"
            return os.path.join(*ext_name.split(".")) + suffix
        return super().get_ext_filename(ext_name)


extension = CTypesExtension(
    "vertexenum._vertexenum",
    sources=[
        os.path.join(CSRC, "lrsmp.c"),
        os.path.join(CSRC, "lrslib.c"),
        os.path.join(CSRC, "pyvertexenum.c"),
    ],
    include_dirs=[CSRC],
    # self-contained DLL: don't depend on the MinGW runtime being on PATH
    extra_link_args=["-static"] if os.name == "nt" else [],
)

setup(
    ext_modules=[extension],
    cmdclass={"build_ext": build_ctypes_ext},
)
