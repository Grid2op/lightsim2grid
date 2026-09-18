export __COMPILE_MARCHNATIVE=1
export __O3_OPTIM=1
export PATH_CKTSO=`pwd`/cktso
export PATH_NICSLU=`pwd`/nicslu/nicslu202110

# NICSLU/CKTSO ship as shared libraries (unlike KLU/CXSparse, which are static
# and fully embedded), so the dynamic loader needs to find them at import time.
export LD_LIBRARY_PATH="`pwd`/nicslu/nicslu202110/linux/lib_centos6_x64_gcc482_fma/int32:`pwd`/cktso/ubuntu2004_x64_gcc940:$LD_LIBRARY_PATH"
