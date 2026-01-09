echo "Solver"
python benchmark_solvers.py --env_name l2rpn_case14_sandbox --no_test --number 1000
sleep 1
echo "\n"
python benchmark_solvers.py --env_name l2rpn_neurips_2020_track2_small --no_test --number 1000
sleep 5
echo "\n\n DC solver"
python benchmark_dc_solvers.py --env_name l2rpn_case14_sandbox --no_test --number 8000
echo "\n"
sleep 1
python benchmark_dc_solvers.py --env_name l2rpn_neurips_2020_track2_small --no_test --number 8000
echo "\n\n grid size"
sleep 5
python benchmark_grid_size.py
echo "\n\n Vs pypowsybl"
sleep 5
python compare_lightsim2grid_pypowsybl.py --case ieee9
python compare_lightsim2grid_pypowsybl.py --case ieee14
python compare_lightsim2grid_pypowsybl.py --case ieee30
python compare_lightsim2grid_pypowsybl.py --case ieee57
python compare_lightsim2grid_pypowsybl.py --case ieee118
python compare_lightsim2grid_pypowsybl.py --case ieee300
