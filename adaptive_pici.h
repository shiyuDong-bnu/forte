/*
 *@BEGIN LICENSE
 *
 * Libadaptive: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef _pi_fci_h_
#define _pi_fci_h_

#include <fstream>

#include <libmints/wavefunction.h>
#include <liboptions/liboptions.h>
#include <physconst.h>

#include "integrals.h"
#include "string_determinant.h"
#include "bitset_determinant.h"

namespace psi{ namespace libadaptive{

enum PropagatorType {LinearPropagator,QuadraticPropagator};

/**
 * @brief The SparsePathIntegralCI class
 * This class implements an a sparse path-integral FCI algorithm
 */
class AdaptivePathIntegralCI : public Wavefunction
{
public:
    // ==> Class Constructor and Destructor <==

    /**
     * Constructor
     * @param wfn The main wavefunction object
     * @param options The main options object
     * @param ints A pointer to an allocated integral object
     */
    AdaptivePathIntegralCI(boost::shared_ptr<Wavefunction> wfn, Options &options, ExplorerIntegrals* ints);

    /// Destructor
    ~AdaptivePathIntegralCI();

    // ==> Class Interface <==

    /// Compute the energy
    double compute_energy();

    double compute_energy_parallel();

private:

    // ==> Class data <==

    /// A reference to the options object
    Options& options_;
    /// The molecular integrals required by Explorer
    ExplorerIntegrals* ints_;
    /// The maximum number of threads
    int num_threads_;
    /// The type of propagator used
    PropagatorType propagator_;
    /// A string that describes the propagator type
    std::string propagator_description_;
    /// The wave function symmetry
    int wavefunction_symmetry_;
    /// The symmetry of each orbital in Pitzer ordering
    std::vector<int> mo_symmetry_;
    /// The number of correlated molecular orbitals
    int ncmo_;
    /// The number of correlated molecular orbitals per irrep
    Dimension ncmopi_;
    /// The nuclear repulsion energy
    double nuclear_repulsion_energy_;
    /// The reference determinant
    StringDeterminant reference_determinant_;
    /// The PT2 energy correction
    std::vector<double> multistate_pt2_energy_correction_;

    /// The threshold applied to the primary space
    double spawning_threshold_;
    /// The threshold applied during the initial guess
    double initial_guess_spawning_threshold_;
    /// The size of the time step (TAU)
    double time_step_;
    /// The number of roots computed
    int nroot_;
    /// The maximum number of iterations
    int maxiter_;

    /// Estimate the energy via a variational procedure?
    bool variational_estimate_;
    /// The frequency of variational estimation of the energy
    int energy_estimate_freq_;

    /// Use an adaptive time step?
    bool adaptive_beta_;
    /// Shift the Hamiltonian?
    bool do_shift_;
    ///
    bool do_adaptive_initial_guess_;



    /// Prescreen the spawning step of single excitations?
    bool do_prescreen_spawning_;
    /// The tollerance factor applied when prescreening singles
    double prescreening_tollerance_factor_;

    /// The energy convergence criterium
    double e_convergence_;

    /// Maximum value of the one-electron coupling
    double new_max_one_HJI_;
    double old_max_one_HJI_;

    /// Maximum value of the two-electron coupling
    double new_max_two_HJI_;
    double old_max_two_HJI_;

    /// Number of determinants v
    /// Number of determinants visited during a time step
    size_t ndet_visited_;
    /// Number of determinants accepted during a time step
    size_t ndet_accepted_;
    /// Number of determinants spawned during a time step
    size_t nspawned_;
    /// Number of determinants that don't spawn
    size_t nzerospawn_;

    ///
    bool do_dynamic_prescreening_;
    std::map<BitsetDeterminant,std::pair<double,double>> dets_max_couplings_;


    // ==> Class functions <==

    /// All that happens before we compute the energy
    void startup();

    /// Print information about this calculation
    void print_info();

    /// Print a wave function
    void print_wfn(std::vector<BitsetDeterminant> space, std::vector<double> C);

    /// Initial wave function guess
    double initial_guess(std::vector<BitsetDeterminant>& dets,std::vector<double>& C);

    ///
    /**
    * Propagate the wave function by a step of length tau
    * @param dets The set of determinants that form the wave function at time n
    * @param C The wave function coefficients at time n
    * @param order The order of the integration algorithm
    * @param tau The time step in a.u.
    * @param spawning_threshold The threshold used to accept or reject spawning events
    * @param S An energy shift subtracted from the Hamiltonian
    */
    void propagate(PropagatorType propagator,std::vector<BitsetDeterminant>& dets,std::vector<double>& C,double tau,double spawning_threshold_,double S);

    void propagate_first_order(std::vector<BitsetDeterminant>& dets,std::vector<double>& C,double tau,double spawning_threshold,double S);

    void propagate_second_order(std::vector<BitsetDeterminant>& dets,std::vector<double>& C,double tau,double spawning_threshold,double S);

    std::map<std::string, double> estimate_energy(std::vector<BitsetDeterminant>& dets,std::vector<double>& C);

    double estimate_proj_energy(std::vector<BitsetDeterminant>& dets,std::vector<double>& C);

    double estimate_var_energy(std::vector<BitsetDeterminant>& dets,std::vector<double>& C);

    /// Perform a time step
    double time_step_optimized(double spawning_threshold,BitsetDeterminant& detI, double CI, std::map<BitsetDeterminant,double>& new_space_C, double E0);

    /// Apply tau x H to a determinant
    size_t apply_tau_H(double tau,double spawning_threshold,BitsetDeterminant& detI, double CI, std::map<BitsetDeterminant,double>& new_space_C, double E0);

    size_t apply_tau_H_spawning(double tau,double spawning_threshold,BitsetDeterminant& detI, double CI, std::map<BitsetDeterminant,double>& new_space_C, double E0,std::pair<double,double>& max_coupling);
};

}} // End Namespaces

#endif // _pi_fci_h_