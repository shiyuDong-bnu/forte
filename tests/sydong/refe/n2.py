import psi4
import forte
psi4.set_output_file("n2.out",False)
psi4.set_num_threads(10)
psi4.set_memory("3 gb")
psi4.set_options(
   {
  "basis"            :     "6-31g" ,
  "reference"        :     "rhf",
  "scf_type"         :     "pk",
  "d_convergence"    :     "8",
  "e_convergence"    :     "12",
   }
)
forte_options={
    "ACTIVE_SPACE_SOLVER":"fci",
    "CORRELATION_SOLVER" : "sa-mrdsrg",
    "CORR_LEVEL"         :"ldsrg2",
    "restricted_docc":         [2,0,0,0,0,2,0,0],
    "active":                  [1,0,1,1,0,1,1,1],
    "DSRG_S":1.0,
}
mol=psi4.geometry(
"""
  0 1
  N
  N  1 R
  R = 1.1
"""
)
e,wfn=psi4.energy("scf",return_wfn=True)

e_forte=psi4.energy("forte",forte_options=forte_options,ref_wfn=wfn)
