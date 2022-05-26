#!/bin/bash

rewienski_a=(2 2 10 10 6 8 8 4 4)

rewienski_b=(0.01 0.1 0.01 0.1 0.055 0.075 0.025 0.075 0.025)

for ((i = 0 ; i < ${#rewienski_a[@]} ; i++)); do

file="${rewienski_a[i]}_${rewienski_b[i]}_1d_burgers_rewienski_snapshot.prm"

echo "# Listing of Parameters"                                                                      >> $file
echo "# ---------------------"                                                                      >> $file
echo " "                                                                                            >> $file
echo "set dimension = 1 "                                                                           >> $file
echo "set test_type = flow_solver"                                                                  >> $file
echo "set pde_type = burgers_rewienski"                                                             >> $file
echo " "                                                                                            >> $file
echo "set use_weak_form = true"                                                                     >> $file
echo "set use_collocated_nodes = false"                                                             >> $file
echo " "                                                                                            >> $file
echo "subsection grid refinement study"                                                             >> $file
echo " set num_refinements = 10"                                                                    >> $file
echo " set poly_degree = 0"                                                                         >> $file
echo " set grid_left = 0.0"                                                                         >> $file
echo " set grid_right = 100.0"                                                                      >> $file
echo "end"                                                                                          >> $file
echo " "                                                                                            >> $file
echo "#Burgers parameters"                                                                          >> $file
echo "subsection burgers"                                                                           >> $file
echo " set rewienski_a = ${rewienski_a[i]}"                                                         >> $file
echo " set rewienski_b = ${rewienski_b[i]}"                                                         >> $file
echo "end"                                                                                          >> $file
echo " "                                                                                            >> $file
echo "subsection flow_solver"                                                                       >> $file
echo " set flow_case_type = burgers_rewienski_snapshot"                                             >> $file
echo " set final_time = 0.5"                                                                        >> $file
echo " set steady_state = true"                                                                     >> $file
echo "end"                                                                                          >> $file
echo " "                                                                                            >> $file
echo "subsection ODE solver "                                                                       >> $file
echo " set nonlinear_max_iterations            = 50"                                               >> $file
echo " set nonlinear_steady_residual_tolerance = 1e-15"                                             >> $file
echo " set print_iteration_modulo              = 1"                                                 >> $file
echo " set output_solution_vector_modulo       = 1"                                                 >> $file
echo " set solutions_table_filename = ${rewienski_a[i]}_${rewienski_b[i]}_snapshot">> $file
echo " set ode_solver_type                     = implicit"                                          >> $file
echo " end"                                                                                         >> $file
echo " "                                                                                            >> $file
echo "subsection manufactured solution convergence study"                                           >> $file
echo " set use_manufactured_source_term = true"                                                     >> $file
echo "end"                                                                                          >> $file


dir=$(pwd)
mpirun "-n" "1" "$1/PHiLiP_1D" "-i" "${dir}/${file}"

rm ${file}
done

echo Done!