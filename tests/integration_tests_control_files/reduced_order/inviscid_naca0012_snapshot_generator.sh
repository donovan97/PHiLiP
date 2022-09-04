#!/bin/bash

mach=(0.53141
      0.76931
      0.69289
      0.64643
      0.72853
      0.84855
      0.57817
      0.51415
      0.66708
      0.60021
      0.83636
      0.78634
      0.55223
      0.54099
      0.86372
      0.75664
      0.50061
      0.82930
      0.73053
      0.74410
      0.77878
      0.58359
      0.63266
      0.87873
      0.89243
      0.62242
      0.68449
      0.88421
      0.80635
      0.67647
      0.52927
      0.85769
      0.65920
      0.79251
      0.81307
      0.61777
      0.59600
      0.71600
      0.70347
      0.56437)

alpha=(2.65734
       3.15271
       2.76382
       1.28005
       3.86530
       1.09970
       4.30085
       4.13370
       4.01247
       3.37283
       1.70110
       2.94577
       3.39192
       0.54214
       3.07837
       4.67450
       0.08562
       1.61688
       4.80579
       2.23066
       1.78834
       1.94259
       1.43101
       0.36499
       2.37561
       3.62337
       1.20967
       2.59361
       0.14454
       0.85104
       4.52497
       4.38670
       4.88508
       2.08307
       3.99344
       2.32748
       0.62716
       0.88308
       3.70783
       0.42840)

for ((i = 0 ; i < ${#mach[@]} ; i++)); do

file="${mach[i]}_${alpha[i]}_naca0012.prm"

echo "set run_type = flow_simulation" >> $file
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
echo "    set linear_residual_tolerance = 1e-04" >> $file
echo "    set max_iterations            = 2000" >> $file
echo "    set restart_number            = 300" >> $file
echo "  end" >> $file
echo "end" >> $file
echo "" >> $file
echo "subsection ODE solver" >> $file
echo "  #set output_solution_every_x_steps = 1" >> $file
echo "  set nonlinear_max_iterations            = 50" >> $file
echo "  set nonlinear_steady_residual_tolerance = 1e-14" >> $file
echo "  set ode_solver_type  = implicit" >> $file
echo "  set initial_time_step = 1e3" >> $file
echo "  set time_step_factor_residual = 15.0" >> $file
echo "  set time_step_factor_residual_exp = 2" >> $file
echo "  set print_iteration_modulo              = 1" >> $file
echo "  set output_final_steady_state_solution_to_file       = true" >> $file
echo "  set steady_state_final_solution_filename = ${mach[i]}_${alpha[i]}_solution_snapshot" >> $file
echo "end" >> $file
echo "" >> $file
echo "subsection grid refinement study" >> $file
echo " set num_refinements = 0" >> $file
echo "end" >> $file
echo "" >> $file
echo "subsection flow_solver" >> $file
echo "  set flow_case_type = naca0012" >> $file
echo "  set poly_degree = 0" >> $file
echo "  set steady_state = true" >> $file
echo "  set steady_state_polynomial_ramping = true" >> $file
echo "  subsection grid" >> $file
echo "    set input_mesh_filename = ../../meshes/naca0012_hopw_ref1" >> $file
echo "  end" >> $file
echo "end" >> $file

dir=$(pwd)
mpirun "-n" "4" "$1/PHiLiP_2D" "-i" "${dir}/${file}"

rm ${file}
done

echo Done!