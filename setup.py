from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os
import warnings

__version__ = "0.5.1"
KLU_SOLVER_AVAILABLE = False

# courtesy to
# https://github.com/pybind/python_example/blob/master/setup.py
class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked.

    @author: Sylvain Corlay
    """

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

    @author: Sylvain Corlay
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
    The newer version is preferred over c++11 (when it is available).

    @author: Sylvain Corlay
    """
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag):
            return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')


class BuildExt(build_ext):
    """
    A custom build extension for adding compiler-specific options.
    @author: Sylvain Corlay
    """
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

        # for debug option
        if hasattr(self.compiler, "compiler"):
            print()
            print("Compiler options used:")
            print(self.compiler.compiler)
            print()

        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            opts.append("-DVERSION_INFO=\"%s\"" % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args += opts
            ext.extra_link_args += link_opts
        build_ext.build_extensions(self)

# Try to link against SuiteSparse (if available)
# check that they exist
suitesparse_path = os.path.abspath("./SuiteSparse")
eigen_path = os.path.abspath(".")
LIBS = ["{}/KLU/Lib/libklu.a",
        "{}/BTF/Lib/libbtf.a",
        "{}/AMD/Lib/libamd.a",
        "{}/COLAMD/Lib/libcolamd.a",
        "{}/CXSparse/Lib/libcxsparse.a",
        "{}/SuiteSparse_config/libsuitesparseconfig.a"
        ]
LIBS = [el.format(suitesparse_path) for el in LIBS]
exists_libs = True
for el in LIBS:
    if not os.path.exists(el):
        exists_libs = False
if exists_libs:
    # you will be able to use "SuiteSparse" and the faster "KLU" linear solver
    KLU_SOLVER_AVAILABLE = True

    # include directory
    INCLUDE_suitesparse = ["{}/SuiteSparse_config",
                           "{}/CXSparse/Include",
                           "{}/AMD/Include",
                           "{}/BTF/Include",
                           "{}/COLAMD/Include",
                           "{}/KLU/Include"
                           ]
    INCLUDE_suitesparse = [el.format(suitesparse_path) for el in INCLUDE_suitesparse]
else:
    # suitesparse, and in particular the KLU linear solver is not available.
    # we'll use a default solver (a bit slower)
    LIBS = []
    INCLUDE_suitesparse = []
    warnings.warn("SuiteSparse is not available on your system, or has not been compiled. The faster "
                  "\"KLU\" linear algebra solver will not be available. The \"SparseLU\" solver will however "
                  "be available, which is maybe ~30% slower than \"KLU\". If you are using grid2op there "
                  "will still be a huge benefit.")

INCLUDE = INCLUDE_suitesparse
INCLUDE.append("{}/eigen".format(eigen_path))

include_dirs = [
                # Path to pybind11 headers
                get_pybind_include(),
                get_pybind_include(user=True)
]
include_dirs += INCLUDE

# compiler options
extra_compile_args_tmp = ["-DNDEBUG"]
if sys.platform.startswith('linux'):
    # extra_compile_args_tmp = ["-fext-numeric-literals"]
    # -fext-numeric-literals is used for definition of complex number by some version of gcc
    extra_compile_args_tmp += []
elif sys.platform.startswith("darwin"):
    # extra_compile_args_tmp = ["-fsized-deallocation"]
    extra_compile_args_tmp += []
    # fix a bug in pybind11
    # https://github.com/pybind/pybind11/issues/1604
elif sys.platform.startswith("win32"):
    extra_compile_args_tmp += ["-D_USE_MATH_DEFINES"]
    # otherwise windows compiler does not import "M_PI" from the math header


# for even greater speed, you can add the "-march=native" flag. It does not work on all platform, that is
# why we deactivated it by default
# extra_compile_args_tmp += ["-march=native"]

# if you have installed some BLAS or LAPACKE libraries (on ubuntu sudo apt-get install libblas-dev liblapacke-dev)
# you can also trigger their use when using eigen.
# extra_compile_args_tmp += ["-DEIGEN_USE_BLAS", "-DEIGEN_USE_LAPACKE"]

extra_compile_args = extra_compile_args_tmp
src_files = ['src/main.cpp',
             "src/SparseLUSolver.cpp",
             "src/BaseConstants.cpp",
             "src/GridModel.cpp",
             "src/DataConverter.cpp",
             "src/DataLine.cpp",
             "src/DataGeneric.cpp",
             "src/DataShunt.cpp",
             "src/DataTrafo.cpp",
             "src/DataLoad.cpp",
             "src/DataGen.cpp",
             "src/DataSGen.cpp",
             "src/BaseNRSolver.cpp",
             "src/ChooseSolver.cpp",
             "src/GaussSeidelSolver.cpp",
             "src/GaussSeidelSynchSolver.cpp",
             "src/BaseSolver.cpp",
             "src/DCSolver.cpp"]

if KLU_SOLVER_AVAILABLE:
    src_files.append("src/KLUSolver.cpp")
    extra_compile_args_tmp.append("-DKLU_SOLVER_AVAILABLE")

ext_modules = [
    Extension(
        'lightsim2grid_cpp',
        src_files,
        include_dirs=include_dirs,
        language='c++',
        extra_objects=LIBS,
        extra_compile_args=extra_compile_args
    )
]

pkgs = {
    "required": [
        'pybind11>=2.4',
        "pandapower",
        "numpy",
        "scipy",
        "grid2op"
    ],
    "extras": {
        "docs": [
            "numpydoc>=0.9.2",
            "sphinx>=2.4.4",
            "sphinx-rtd-theme>=0.4.3",
            "sphinxcontrib-trio>=1.1.0",
            "autodocsumm>=0.1.13",
            # "m2r"
            "recommonmark",
        ],
        "benchmark": [
            "tabulate"
        ],
        "recommended": [
            "grid2op>=1.5.0"
        ],
        "test": [
            "grid2op>=1.5.0"
        ]
    }
}

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='LightSim2Grid',
      version=__version__,
      author='Benjamin Donnot',
      author_email='benjamin.donnot@rte-france.com',
      url='https://github.com/BDonnot/lightsim2grid/',
      description='LightSim2Grid implements a c++ backend targeting the Grid2Op platform.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      ext_modules=ext_modules,
      install_requires=pkgs["required"],
      extras_require=pkgs["extras"],
      setup_requires=['pybind11>=2.4'],
      cmdclass={'build_ext': BuildExt},
      zip_safe=False,
      packages=setuptools.find_packages(),
      keywords='pandapower powergrid simulator KLU Eigen c++',
      classifiers=[
            'Development Status :: 4 - Beta',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "Natural Language :: English"
      ]
     )
