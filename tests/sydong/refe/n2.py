# MR-LDSRG(2) test on reading previous amplitudes files
# The number of commutators in the BCH expansion
# is truncated to speed up the test.
import psi4
import forte
refrhf  = -108.867618373021401
refdsrg = -109.100877299246122 # default dsrg_rsc_ncomm
refdsrg = -109.100837616506638 # keep only 4-nested commutator
psi4.set_output_file("testn2.out",False)
psi4.set_num_threads(10)
psi4.set_memory("3 gb")
mol=psi4.geometry(
"""
  0 1
  N
  N  1 R
  R = 1.1
"""
)
psi4.set_options(
   {
        "basis":"6-31g",
        "scf_type":"pk",
        "reference" : "rhf",
        "d_convergence"  :   8,
        "e_convergence"   : 12,
   }
)
forte_options={
  "active_space_solver"  :  "fci",
  "correlation_solver"   :  "mrdsrg",
  "corr_level"           : "ldsrg2",
  "restricted_docc"      : [2,0,0,0,0,2,0,0],
  "active"               : [1,0,1,1,0,1,1,1],
  "dsrg_s"               : 1.0,
  "e_convergence"        :  1e-8 ,
  "r_convergence"         : 1e-7,
  "dsrg_read_amps"       :  True,
  "dsrg_diis_start"       : 1,
  "dsrg_rsc_ncomm"        : 4,
}

Escf, wfn = psi4.energy('scf', return_wfn=True)

# fix orbital phase
Ca = wfn.Ca().clone()
nirrep = wfn.nirrep()
coldim = Ca.coldim()
rowdim = Ca.rowdim()
for h in range(nirrep):
    for i in range(coldim[h]):
        v = Ca.get(h, 0, i)
        if v < 0:
            for j in range(rowdim[h]):
                Ca.set(h, j, i, -1.0 * Ca.get(h, j, i))
wfn.Ca().copy(Ca)

psi4.energy('forte', forte_options=forte_options,ref_wfn=wfn)
# compare_values(refdsrg, variable("CURRENT ENERGY"), 8, "MR-LDSRG(2) unrelaxed energy")
