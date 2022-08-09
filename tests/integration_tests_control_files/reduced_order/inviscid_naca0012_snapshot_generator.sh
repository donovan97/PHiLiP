#!/bin/bash

mach=(0.3
     0.425
     0.3625
     0.4875
     0.33125
     0.45625
     0.39375
     0.51875
     0.31562
     0.44063
     0.37812
     0.50313
     0.34687
     0.47188
     0.40938
     0.53438
     0.30781
     0.43281
     0.37031
     0.49531
     0.33906
     0.46406
     0.40156
     0.52656
     0.32344
     0.44844
     0.38594
     0.51094
     0.35469
     0.47969
     0.41719
     0.54219
     0.30391
     0.42891
     0.36641
     0.49141)

alpha=(0
       1
       2
       0.33333
       1.33333
       2.33333
       0.66667
       1.66667
       2.66667
       0.11111
       1.11111
       2.11111
       0.44444
       1.44444
       2.44444
       0.77778
       1.77778
       2.77778
       0.22222
       1.22222
       2.22222
       0.55556
       1.55556
       2.55556
       0.88889
       1.88889
       2.88889
       0.03704
       1.03704
       2.03704
       0.37037
       1.37037
       2.37037
       0.7037
       1.7037
       2.7037)

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
echo "  set poly_degree = 3" >> $file
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