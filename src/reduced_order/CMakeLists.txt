set(SOURCE
    pod_basis_state.cpp
    pod_basis_types.cpp
    pod_adaptation_base.cpp
    pod_basis_sensitivity.cpp
    pod_basis_sensitivity_types.cpp
    pod_fine_basis_adaptation.cpp
    pod_full_dim_adaptation.cpp full_order_solution.cpp full_order_solution.h reduced_order_solution.cpp reduced_order_solution.h rom_test_location.cpp rom_test_location.h pod_basis_online.cpp pod_basis_online.h rbf_interpolation.cpp rbf_interpolation.h delaunay.cpp delaunay.h)

foreach(dim RANGE 1 3)
    # Output library
    string(CONCAT PODLib POD_${dim}D)
    add_library(${PODLib} STATIC ${SOURCE})

    # Library dependency
    string(CONCAT ParameterLib ParametersLibrary)
    target_link_libraries(${PODLib} ${ParameterLib})

    # Setup target with deal.II
    if(NOT DOC_ONLY)
        DEAL_II_SETUP_TARGET(${PODLib})
    endif()

    target_compile_definitions(${PODLib} PRIVATE PHILIP_DIM=${dim})

    # Setup target with deal.II
    unset(PODLib)

endforeach()