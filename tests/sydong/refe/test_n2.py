import psi4
import forte
psi4.set_output_file("testn2.out",False)
psi4.set_num_threads(10)
psi4.set_memory("3 gb")
psi4.set_options(
   {
        "basis":"cc-pVDZ",
        "scf_type":"pk",
        "frozen_docc":[1,0,0,0],
        "active":[3,0,1,1],
   }
)
forte_options={
    "reference":"cascsf",
    "ACTIVE_SPACE_SOLVER":"detci",
    "CORRELATION_SOLVER" : "sa-mrdsrg",
    "restricted_docc":         [2,0,0,0,0,2,0,0],
    "active":                  [1,0,1,1,0,1,1,1],
    "DSRG_S":0.5,
    "CORR_LEVEL":"ldsrg2",
}
mol=psi4.geometry(
"""
  0 1
  N
  N  1 R
  R = 1.1
"""
)
e,wfn=psi4.energy("casscf",return_wfn=True)

e_forte=psi4.energy("forte",forte_options=forte_options,ref_wfn=wfn)

