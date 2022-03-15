#!/bin/bash

rewienski_a=( 3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000
              3.0000   3.0000     3.0000 3.000 3.000)

rewienski_b=(  0.055 0.0325 0.0775 0.02125 0.06625 0.04375 0.08875 0.015625
                0.060625 0.038125 0.0985938 0.0971875 0.09578 0.094375 0.0929688
                0.0915625 0.0901563 0.0873438 0.0859375 0.0845313 0.083125 0.0817188)

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
echo "  set coarse_basis_dimension = 17">> $file
echo "  set fine_basis_dimension =17">> $file
echo "  set path_to_search = ./burgers_rewienski_gapupper2">> $file
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
echo "  set nonlinear_steady_residual_tolerance = 1e-14">> $file
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