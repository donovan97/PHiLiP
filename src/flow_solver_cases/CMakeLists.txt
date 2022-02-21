set(FLOW_CASES_SOURCE
    flow_solver_case_base.cpp
    periodic_cube_flow.cpp
    1D_burgers_rewienski_snapshot.cpp
    1d_burgers_viscous_snapshot.cpp)

foreach(dim RANGE 1 3)
    # Output library
    string(CONCAT FlowCasesLib FlowCases_${dim}D)
    add_library(${FlowCasesLib} STATIC ${FLOW_CASES_SOURCE})

    # Library dependency
    string(CONCAT ParameterLib ParametersLibrary)
    target_link_libraries(${FlowCasesLib} ${ParameterLib})

    # Setup target with deal.II
    if(NOT DOC_ONLY)
        DEAL_II_SETUP_TARGET(${FlowCasesLib})
    endif()

    target_compile_definitions(${FlowCasesLib} PRIVATE PHILIP_DIM=${dim})

    # Setup target with deal.II
    unset(FlowCasesLib)

endforeach()
