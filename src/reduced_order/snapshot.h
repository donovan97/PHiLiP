#ifndef __SNAPSHOT__
#define __SNAPSHOT__

#include <fstream>
#include <iostream>
#include <filesystem>
#include "functional/functional.h"
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector_operation.h>
#include "parameters/all_parameters.h"
#include "dg/dg.h"
#include "pod_interfaces.h"
#include "reduced_order_solution.h"
#include "full_order_solution.h"
#include "linear_solver/linear_solver.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

/// Class to hold information about the reduced-order solution
template <int dim, int nstate>
class Snapshot
{
public:
    /// Constructor
    Snapshot(double parameter, std::shared_ptr<ROMSolution<dim, nstate>> rom_solution);

    /// Constructor
    Snapshot(double parameter, std::shared_ptr<FOMSolution<dim, nstate>> fom_solution);

    /// Destructor
    virtual ~Snapshot() {};

    void add_FOM(std::shared_ptr<FOMSolution<dim, nstate>> fom_solution);

    void add_ROM(std::shared_ptr<ROMSolution<dim, nstate>> rom_solution);

    void compute_FOM_to_initial_ROM_error();

    void compute_initial_rom_to_final_rom_error(std::shared_ptr<ProperOrthogonalDecomposition::POD<dim>> pod_updated);

    void compute_initial_rom_to_final_rom_error_estimate();

    void compute_total_error();

    double parameter;

    std::shared_ptr<ROMSolution<dim, nstate>> rom_solution;

    std::shared_ptr<FOMSolution<dim, nstate>> fom_solution;

    double fom_to_initial_rom_error;

    double initial_rom_to_final_rom_error;

    double total_error;

};

}
}


#endif
