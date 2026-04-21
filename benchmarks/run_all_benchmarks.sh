export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:../nicslu/nicslu202110/linux/lib_centos6_x64_gcc482_fma/int32/:../cktso/ubuntu2004_x64_gcc940/"
echo "AC Solver..."
./benchmarks_ac.sh > res_benchmarks/benchmarks_ac.txt
echo "... done \n"
sleep 5
echo "DC solver..."
./benchmarks_dc.sh > res_benchmarks/benchmarks_dc.txt
echo "... done \n"
echo "Grid size..."
./benchmarks_grid_size.sh > res_benchmarks/benchmarks_grid_size.txt
echo "... done \n"
echo "Comparison with pypowsybl..."
./benchmark_pypowysbl.sh > res_benchmarks/benchmarks_pypowsybl.txt
echo "... done \n"
echo "Time serie and Contingency analysis..."
./benchmark_ts_ca.sh > res_benchmarks/benchmark_ts_ca.txt
echo "... done \n"
