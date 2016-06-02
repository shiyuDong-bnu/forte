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

#ifndef _dsrg_mrpt3_h_
#define _dsrg_mrpt3_h_

#include <fstream>
#include <boost/assign.hpp>

#include <liboptions/liboptions.h>
#include <libmints/wavefunction.h>

#include "integrals.h"
#include "ambit/blocked_tensor.h"
#include "reference.h"
#include "helpers.h"
#include "blockedtensorfactory.h"
#include "dsrg_time.h"
#include "dsrg_source.h"

using namespace ambit;
namespace psi{ namespace forte{

class DSRG_MRPT3 : public Wavefunction
{
public:
    /**
     * DSRG_MRPT3 Constructor
     * @param ref_wfn The reference wavefunction object
     * @param options The main options object
     * @param ints A pointer to an allocated integral object
     * @param mo_space_info The MOSpaceInfo object
     */
    DSRG_MRPT3(Reference reference, SharedWavefunction ref_wfn, Options &options,
               std::shared_ptr<ForteIntegrals> ints, std::shared_ptr<MOSpaceInfo> mo_space_info);

    /// Destructor
    ~DSRG_MRPT3();

    /// Compute the DSRG-MRPT3 energy
    double compute_energy();

    /// Compute the DSRG-MRPT3 energy with relaxed reference (once)
    double compute_energy_relaxed();

protected:
    // => Class initialization and termination <= //

    /// Called in the constructor
    void startup();
    /// Called in the destructor
    void cleanup();
    /// Print a summary of the options
    void print_summary();
    /// Print levels
    int print_;

    /// The reference object
    Reference reference_;
    /// The energy of the reference
    double Eref_;
    /// The frozen-core energy
    double frozen_core_energy_;

    /// The molecular integrals required by MethodBase
    std::shared_ptr<ForteIntegrals>  ints_;

    /// MO space info
    std::shared_ptr<MOSpaceInfo> mo_space_info_;

    /// List of alpha core MOs
    std::vector<size_t> acore_mos_;
    /// List of alpha active MOs
    std::vector<size_t> aactv_mos_;
    /// List of alpha virtual MOs
    std::vector<size_t> avirt_mos_;
    /// List of beta core MOs
    std::vector<size_t> bcore_mos_;
    /// List of beta active MOs
    std::vector<size_t> bactv_mos_;
    /// List of beta virtual MOs
    std::vector<size_t> bvirt_mos_;

    /// Alpha core label
    std::string acore_label_;
    /// Alpha active label
    std::string aactv_label_;
    /// Alpha virtual label
    std::string avirt_label_;
    /// Beta core label
    std::string bcore_label_;
    /// Beta active label
    std::string bactv_label_;
    /// Beta virtual label
    std::string bvirt_label_;

    /// Map from space label to list of MOs
    std::map<char, std::vector<size_t>> label_to_spacemo_;

    /// If ERI density fitted or Cholesky decomposed
    bool eri_df_;
    /// Auxiliary MOs
    std::vector<size_t> aux_mos_;
    /// Auxiliary space label
    std::string aux_label_;

    /// Fill up two-electron integrals
    void build_tei(BlockedTensor& V);
    /// Fill up density matrix and density cumulants
    void build_density();
    /// Build Fock matrix and diagonal Fock matrix elements
    void build_fock_half();
    /// Build Fock matrix when two-electron integrals are fully stored
    void build_fock_full();

    /// Are orbitals semi-canonicalized?
    bool semi_canonical_;
    /// Check if orbitals are semi-canonicalized
    bool check_semicanonical();
    /// Diagonal elements of Fock matrices
    std::vector<double> Fa_;
    std::vector<double> Fb_;

    // => DSRG related <= //

    /// The flow parameter
    double s_;
    /// Source operator
    std::string source_;
    /// The dsrg source operator
    std::shared_ptr<DSRG_SOURCE> dsrg_source_;
    /// Threshold for the Taylor expansion of f(z) = (1-exp(-z^2))/z
    double taylor_threshold_;

    /// Effective alpha one-electron integrals (used in denormal ordering)
    std::vector<double> aone_eff_;
    /// Effective beta one-electron integrals (used in denormal ordering)
    std::vector<double> bone_eff_;

    /// Renormalize Fock matrix
    void renormalize_F(const bool& plusone = true);
    /// Renormalize two-electron integrals
    void renormalize_V(const bool& plusone = true);

    // => Tensors <= //

    /// Kevin's Tensor Wrapper
    std::shared_ptr<BlockedTensorFactory> BTF_;
    /// Tensor type for AMBIT
    TensorType tensor_type_;

    /// One-particle density matrix
    ambit::BlockedTensor Gamma1_;
    /// One-hole density matrix
    ambit::BlockedTensor Eta1_;
    /// Two-body denisty cumulant
    ambit::BlockedTensor Lambda2_;
    /// Three-body density cumulant
    ambit::BlockedTensor Lambda3_;

    /// Generalized Fock matrix (bare or renormalized)
    ambit::BlockedTensor F_;
    /// Zeroth-order Hamiltonian (bare diagonal blocks of Fock)
    ambit::BlockedTensor F0th_;
    /// Generalized Fock matrix (bare off-diagonal blocks)
    ambit::BlockedTensor F1st_;
    /// Two-electron integrals (bare or renormalized)
    ambit::BlockedTensor V_;
    /// Three-index integrals
    ambit::BlockedTensor B_;
    /// Single excitation amplitude
    ambit::BlockedTensor T1_;
    /// Double excitation amplitude
    ambit::BlockedTensor T2_;
    /// One-body transformed Hamiltonian (active only)
    ambit::BlockedTensor Hbar1_;
    /// Two-body transformed Hamiltonian (active only)
    ambit::BlockedTensor Hbar2_;
    /// One-body temp tensor ([[H0th,A1st],A1st] or 1st-order amplitudes)
    ambit::BlockedTensor O1_;
    /// Two-body temp tensor ([[H0th,A1st],A1st] or 1st-order amplitudes)
    ambit::BlockedTensor O2_;

    /// Diagonal blocks of Fock matrix
    ambit::BlockedTensor Fdiag_;
    /// Unitary matrix to block diagonal Fock
    ambit::BlockedTensor U_;


    // => Amplitude <= //

    /// Compute T2 amplitudes
    void compute_t2();
    /// Check T2 and store large amplitudes
    void check_t2();
    /// Norm of T2
    double T2norm_;
    /// Max (with sign) of T2
    double T2max_;

    /// Compute T1 amplitudes
    void compute_t1();
    /// Check T1 and store large amplitudes
    void check_t1();
    /// Norm of T1
    double T1norm_;
    /// Max (with sign) of T1
    double T1max_;

    /// Number of amplitudes will be printed in amplitude summary
    int ntamp_;
    /// Print amplitudes summary
    void print_amp_summary(const std::string& name,
                           const std::vector<std::pair<std::vector<size_t>, double>>& list, const double &norm,
                           const size_t& number_nonzero);

    /// Threshold for amplitudes considered as intruders
    double intruder_tamp_;
    /// List of large amplitudes
    std::vector<std::pair<std::vector<size_t>, double>> lt1a_;
    std::vector<std::pair<std::vector<size_t>, double>> lt1b_;
    std::vector<std::pair<std::vector<size_t>, double>> lt2aa_;
    std::vector<std::pair<std::vector<size_t>, double>> lt2ab_;
    std::vector<std::pair<std::vector<size_t>, double>> lt2bb_;
    /// Print intruder analysis
    void print_intruder(const std::string& name,
                        const std::vector<std::pair<std::vector<size_t>, double>>& list);


    // => Energy terms <= //

    /// compute second-order energy and transformed Hamiltonian
    double compute_energy_pt2();
    /// compute third-order energy contribution 1.0 / 12.0 * [[[H0th,A1st],A1st],A1st], computed before pt2
    double compute_energy_pt3_1();
    /// compute third-order energy contribution 1.0 / 2.0  * [H1st + Hbar1st,A2nd]
    double compute_energy_pt3_2();
    /// compute third-order energy contribution 1.0 / 2.0  * [Hbar2nd,A1st]
    double compute_energy_pt3_3();

    /// Timings for computing the commutators
    DSRG_TIME dsrg_time_;

    /// Compute zero-body Hbar truncated to 2nd-order
    double Hbar0_;

    /// Compute zero-body term of commutator [H1, T1]
    void H1_T1_C0(BlockedTensor& H1, BlockedTensor& T1, const double& alpha, double& C0);
    /// Compute zero-body term of commutator [H1, T2]
    void H1_T2_C0(BlockedTensor& H1, BlockedTensor& T2, const double& alpha, double& C0);
    /// Compute zero-body term of commutator [H2, T1]
    void H2_T1_C0(BlockedTensor& H2, BlockedTensor& T1, const double& alpha, double& C0);
    /// Compute zero-body term of commutator [H2, T2]
    void H2_T2_C0(BlockedTensor& H2, BlockedTensor& T2, const double& alpha, double& C0);

    /// Compute one-body term of commutator [H1, T1]
    void H1_T1_C1(BlockedTensor& H1, BlockedTensor& T1, const double& alpha, BlockedTensor& C1);
    /// Compute one-body term of commutator [H1, T2]
    void H1_T2_C1(BlockedTensor& H1, BlockedTensor& T2, const double& alpha, BlockedTensor& C1);
    /// Compute one-body term of commutator [H2, T1]
    void H2_T1_C1(BlockedTensor& H2, BlockedTensor& T1, const double& alpha, BlockedTensor& C1);
    /// Compute one-body term of commutator [H2, T2]
    void H2_T2_C1(BlockedTensor& H2, BlockedTensor& T2, const double& alpha, BlockedTensor& C1);

    /// Compute two-body term of commutator [H2, T1]
    void H2_T1_C2(BlockedTensor& H2, BlockedTensor& T1, const double& alpha, BlockedTensor& C2);
    /// Compute two-body term of commutator [H1, T2]
    void H1_T2_C2(BlockedTensor& H1, BlockedTensor& T2, const double& alpha, BlockedTensor& C2);
    /// Compute two-body term of commutator [H2, T2]
    void H2_T2_C2(BlockedTensor& H2, BlockedTensor& T2, const double& alpha, BlockedTensor& C2);

    /// Compute one-body term of commutator [V, T1], V is constructed from B (DF / CD)
    void V_T1_C1_DF(BlockedTensor& B, BlockedTensor& T1, const double& alpha, BlockedTensor& C1);
    /// Compute one-body term of commutator [V, T2], V is constructed from B (DF / CD)
    void V_T2_C1_DF(BlockedTensor& B, BlockedTensor& T2, const double& alpha, BlockedTensor& C1);
    /// Compute two-body term of commutator [V, T1], V is constructed from B (DF / CD)
    void V_T1_C2_DF(BlockedTensor& B, BlockedTensor& T1, const double& alpha, BlockedTensor& C2);
    /// Compute two-body term of commutator [V, T2], V is constructed from B (DF / CD)
    void V_T2_C2_DF(BlockedTensor& B, BlockedTensor& T2, const double& alpha, BlockedTensor& C2);

    /// Compute two-body term of commutator [V, T2] in batching, particle-particle contraction
    void V_T2_C2_DF_PP(BlockedTensor& B, BlockedTensor& T2, const double& alpha, BlockedTensor& C2, const std::vector<std::vector<size_t>>& sub_virt_mos);
    /// Compute two-body term of commutator [V, T2] in batching, particle-hole contraction (exchange part)
    void V_T2_C2_DF_PH_EX(BlockedTensor& B, BlockedTensor& T2, const double& alpha, BlockedTensor& C2, const std::vector<std::vector<size_t>>& sub_virt_mos);

    // => Reference relaxation <= //

    /// Relaxation type
    std::string relax_ref_;

    /// Transfer integrals for FCI
    void transfer_integrals();
    /// Diagonalize the diagonal blocks of the Fock matrix
    std::vector<std::vector<double>> diagonalize_Fock_diagblocks(BlockedTensor& U);
    /// Separate an 2D ambit::Tensor according to its irrep
    ambit::Tensor separate_tensor(ambit::Tensor& tens, const Dimension& irrep, const int& h);
    /// Combine a separated 2D ambit::Tensor
    void combine_tensor(ambit::Tensor& tens, ambit::Tensor& tens_h, const Dimension& irrep, const int& h);
};

}} // End Namespaces

#endif // _dsrg_mrpt3_h_
