#!/bin/bash

rewienski_a=( 3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3)

rewienski_b=( 0.0100000000000000	0.0118367346938776	0.0136734693877551	0.0155102040816327	0.0173469387755102	0.0191836734693878
              0.0210204081632653	0.0228571428571429	0.0246938775510204	0.0265306122448980	0.0283673469387755	0.0302040816326531
              0.0320408163265306	0.0338775510204082	0.0357142857142857	0.0375510204081633	0.0393877551020408	0.0412244897959184
              0.0430612244897959	0.0448979591836735	0.0467346938775510	0.0485714285714286	0.0504081632653061	0.0522448979591837
              0.0540816326530612	0.0559183673469388	0.0577551020408163	0.0595918367346939	0.0614285714285714	0.0632653061224490
              0.0651020408163265	0.0669387755102041	0.0687755102040816	0.0706122448979592	0.0724489795918367	0.0742857142857143
              0.0761224489795918	0.0779591836734694	0.0797959183673469	0.0816326530612245	0.0834693877551020	0.0853061224489796
              0.0871428571428572	0.0889795918367347	0.0908163265306123	0.0926530612244898	0.0944897959183674	0.0963265306122449
              0.0981632653061225	0.100000000000000)

for ((i = 0 ; i < ${#rewienski_a[@]} ; i++)); do

file="${rewienski_a[i]}_${rewienski_b[i]}_1d_burgers_rewienski_solution_snapshots_steady.prm"

echo "# Listing of Parameters"                                                                      >> $file   
echo "# ---------------------"                                                                      >> $file   
echo " "                                                                                            >> $file   
echo "set dimension = 1 "                                                                           >> $file   
echo "set test_type = finite_difference_sensitivity"                                                                  >> $file
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
echo " set sensitivity_table_filename = ${rewienski_b[i]}_sensitivity_snapshots_steady"                 >> $file
echo " set steady_state = true"                                                                     >> $file   
echo "end"                                                                                          >> $file   
echo " "                                                                                            >> $file   
echo "subsection ODE solver "                                                                       >> $file   
echo " set initial_time_step = 0.1"                                                                 >> $file   
echo " set nonlinear_max_iterations            = 500"                                               >> $file   
echo " set nonlinear_steady_residual_tolerance = 1e-12"                                             >> $file   
echo " set print_iteration_modulo              = 1"                                                 >> $file   
echo " set output_solution_vector_modulo       = 1"                                                 >> $file   
echo " set solutions_table_filename = ${rewienski_a[i]}_${rewienski_b[i]}_solution_snapshots_steady">> $file
echo " set ode_solver_type                     = implicit"                                          >> $file   
echo " end"                                                                                         >> $file   
echo " "                                                                                            >> $file   
echo "subsection manufactured solution convergence study"                                           >> $file   
echo " set use_manufactured_source_term = true"                                                     >> $file   
echo "end"                                                                                          >> $file   


dir=$(pwd)
mpirun "-n" "1" "$HOME/Codes/PHiLiP/cmake-build-release/bin/PHiLiP_1D" "-i" "${dir}/${file}"

rm ${file}
done

echo Done!