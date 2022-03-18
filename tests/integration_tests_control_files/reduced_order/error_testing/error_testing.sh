#!/bin/bash

rewienski_a=( 3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000)

rewienski_b=( 0.01  0.0115   0.0127 0.0156 0.0178 0.0190 0.0206 0.0225 0.0236 0.0259 0.0280 0.0312
              0.0340 0.0365 0.0400 0.0418 0.0485 0.0524 0.0600 0.0650 0.07   0.0750 0.0800 0.0850 0.0900 0.0950 0.1)

for ((i = 0 ; i < ${#rewienski_a[@]} ; i++)); do

file="${rewienski_a[i]}_${rewienski_b[i]}_1d_burgers_rewienski_test.prm"

echo "# Listing of Parameters"                                                                      >> $file
echo "# ---------------------"                                                                      >> $file
echo " "                                                                                            >> $file
echo "set dimension = 1 "                                                                           >> $file
echo "set test_type = POD_adaptation"                                                                  >> $file
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
echo "subsection flow_solver"  >> $file
echo "  set flow_case_type = burgers_rewienski_snapshot"  >> $file
echo "end"  >> $file
echo "#Burgers parameters"                                                                          >> $file
echo "subsection burgers"                                                                           >> $file
echo " set rewienski_a = ${rewienski_a[i]}"                                                         >> $file
echo " set rewienski_b = ${rewienski_b[i]}"                                                         >> $file
echo "end"                                                                                          >> $file
echo " "                                                                                            >> $file
echo "#Reduced order parameters">> $file
echo "subsection reduced order">> $file
echo "  set adaptation_tolerance = 1E-16">> $file
echo "  set adapt_coarse_basis_constant = 2">> $file
echo "  set coarse_basis_dimension = 12">> $file
echo "  set fine_basis_dimension =12">> $file
echo "  set path_to_search = ./burgers_rewienski_1param_consistent">> $file
echo "  set method_of_snapshots = true">> $file
echo "  set consider_error_sign = false">> $file
echo "end">> $file
echo " ">> $file
echo "subsection linear solver">> $file
echo "  set linear_solver_type = direct">> $file
echo "end">> $file
echo " ">> $file
echo "subsection ODE solver">> $file
echo "  set nonlinear_max_iterations = 50">> $file
echo "  set nonlinear_steady_residual_tolerance = 5e-16">> $file
echo "  set print_iteration_modulo  = 1">> $file
echo "  set ode_solver_type = implicit">> $file
echo "end">> $file
echo " ">> $file
echo "subsection manufactured solution convergence study">> $file
echo "  set use_manufactured_source_term = true">> $file
echo "end">> $file



dir=$(pwd)
mpirun "-n" "1" "$HOME/Codes/PHiLiP/cmake-build-release/bin/PHiLiP_1D" "-i" "${dir}/${file}"

rm ${file}
done

echo Done!