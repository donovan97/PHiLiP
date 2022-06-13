#!/bin/bash

mach=(0.5
      0.9
      0.9
      0.7
      0.5
      0.7
      0.6
      0.8
      0.55
      0.75
      0.675606423723048
      0.644669190880549
      0.827514215729673
      0.868460675499731
      0.9
      0.821407462389499)

alpha=(4
       4
       0
       2
       0
       1.33333333333333
       2.66666666666667
       0.444444444444446
       1.77777777777778
       3.11111111111111
       4
       0
       3.6673627841082
       1.31601972817496
       0.834803142735038
       0.97359735186394)

for ((i = 0 ; i < ${#mach[@]} ; i++)); do

file="${mach[i]}_${alpha[i]}_naca0012.prm"

echo "set test_type = flow_solver" >> $file
echo "set dimension = 2" >> $file
echo "set pde_type  = euler" >> $file
echo "" >> $file
echo "set conv_num_flux = roe" >> $file
echo "set diss_num_flux = bassi_rebay_2" >> $file
echo "" >> $file
echo "set use_split_form = false" >> $file
echo "" >> $file
echo "subsection artificial dissipation" >> $file
echo "	set add_artificial_dissipation = true" >> $file
echo "end" >> $file
echo "" >> $file
echo "set overintegration = 0" >> $file
echo "" >> $file
echo "subsection euler" >> $file
echo "  set reference_length = 1.0" >> $file
echo "  set mach_infinity = ${mach[i]}" >> $file
echo "  set angle_of_attack = ${alpha[i]}" >> $file
echo "end" >> $file
echo "" >> $file
echo "subsection linear solver" >> $file
echo "  subsection gmres options" >> $file
echo "    set ilut_atol                 = 1e-4" >> $file
echo "    set ilut_rtol                 = 1.00001" >> $file
echo "    set ilut_drop                 = 0.0" >> $file
echo "    set ilut_fill                 = 10" >> $file
echo "    set linear_residual_tolerance = 1e-13" >> $file
echo "    set max_iterations            = 2000" >> $file
echo "    set restart_number            = 200" >> $file
echo "  end" >> $file
echo "end" >> $file
echo "" >> $file
echo "subsection ODE solver" >> $file
echo "  set output_solution_every_x_steps = 1" >> $file
echo "  set nonlinear_max_iterations            = 2000" >> $file
echo "  set nonlinear_steady_residual_tolerance = 1e-15" >> $file
echo "  set ode_solver_type  = implicit" >> $file
echo "  set initial_time_step = 1e3" >> $file
echo "  set time_step_factor_residual = 15.0" >> $file
echo "  set time_step_factor_residual_exp = 2" >> $file
echo "  #set print_iteration_modulo              = 1" >> $file
echo "  set output_solution_vector_modulo       = 1" >> $file
echo "  set solutions_table_filename = ${mach[i]}_${alpha[i]}_solution_snapshot" >> $file
echo "end" >> $file
echo "" >> $file
echo "subsection grid refinement study" >> $file
echo " set poly_degree = 0" >> $file
echo " set num_refinements = 0" >> $file
echo "end" >> $file
echo "" >> $file
echo "subsection flow_solver" >> $file
echo "  set flow_case_type = naca0012" >> $file
echo "  set input_mesh_filename = naca0012_hopw_ref1" >> $file
echo "  set steady_state = true" >> $file
echo "end" >> $file

dir=$(pwd)
mpirun "-n" "4" "$1/PHiLiP_2D" "-i" "${dir}/${file}"

rm ${file}
done

echo Done!