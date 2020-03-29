from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os
import sysconfig

# courtesy to
# https://github.com/pybind/python_example/blob/master/setup.py

__version__ = '0.0.1'


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag): return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')

# import pdb
# pdb.set_trace()
class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args += opts
            ext.extra_link_args += link_opts
        build_ext.build_extensions(self)

suitesparse_path = os.path.abspath("./SuiteSparse")
eigen_path = os.path.abspath(".")

# library to link against (require the "make" command to have run)
LIBS = ["{}/KLU/Lib/libklu.a",
       "{}/BTF/Lib/libbtf.a",
       "{}/AMD/Lib/libamd.a",
       "{}/COLAMD/Lib/libcolamd.a",
       "{}/CXSparse/Lib/libcxsparse.a",
       "{}/SuiteSparse_config/libsuitesparseconfig.a"
       ]

LIBS = [el.format(suitesparse_path) for el in LIBS]

# include directory
INCLUDE_suitesparse = ["{}/SuiteSparse_config",
                       "{}/CXSparse/Include",
                       "{}/AMD/Include",
                       "{}/BTF/Include",
                       "{}/COLAMD/Include",
                       "{}/KLU/Include"
                       ]

INCLUDE_suitesparse = [el.format(suitesparse_path) for el in INCLUDE_suitesparse]
INCLUDE = INCLUDE_suitesparse
INCLUDE.append("{}/eigen".format(eigen_path))
# INCLUDE.append(os.path.abspath("."))

include_dirs = [
                # Path to pybind11 headers
                get_pybind_include(),
                get_pybind_include(user=True)
]
include_dirs += INCLUDE

# compiler options
extra_compile_args = ["-march=native", "-fext-numeric-literals"]
# -fext-numeric-literals is used for definition of complex number by some version of gcc
# -march=native is here to use the vectorization of the code offered by Eigen
ext_modules = [
    Extension(
        'pyklu_cpp',
        ['src/main.cpp', "src/KLUSolver.cpp", "src/GridModel.cpp", "src/DataConverter.cpp",
         "src/DataLine.cpp", "src/DataGeneric.cpp", "src/DataShunt.cpp", "src/DataTrafo.cpp",
         "src/DataLoad.cpp"],
        include_dirs=include_dirs,
        language='c++',
        extra_objects=LIBS,
        extra_compile_args=extra_compile_args
    )
]

setup(name='pyklu',
      version=__version__,
      author='Benjamin Donnot',
      author_email='benjamin.donnot@rte-france.com',
      url='TODO',
      description='A powerflow solver for pandapower written in c++ (using Eigen and klu) for pandapower',
      long_description='This fast powerflow can be used instead of the "pp.runpp()" function of pandapower under some '
                         'circumstances TODO',
      ext_modules=ext_modules,
      install_requires=['pybind11>=2.4', "pandapower"],
      setup_requires=['pybind11>=2.4'],
      cmdclass={'build_ext': BuildExt},
      zip_safe=False,
      packages=['pyklu'],
      keywords='powergrid optmization power-systems KLU sparse-matrix',
      classifiers=[
            'Development Status :: 4 - Beta',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "Natural Language :: English"
      ]
)