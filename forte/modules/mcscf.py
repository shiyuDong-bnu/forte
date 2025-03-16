# import forte
from typing import List
from .module import Module
from forte.data import ForteData
from forte._forte import (
    to_state_nroots_map,
    make_active_space_solver,
    make_mcscf_two_step,
    MOSpaceInfo,
    make_mo_space_info_from_map,
)


class MCSCF(Module):
    """
    A module to perform MCSCF calculations.
    """

    def __init__(self, solver_type: str = "FCI"):
        """
        Parameters
        ----------
        solver_type: str
            The type of the active space solver.
        """
        super().__init__()
        self.solver_type = solver_type

    def _run(self, data: ForteData) -> ForteData:
        state_map = to_state_nroots_map(data.state_weights_map)

        data.active_space_solver = make_active_space_solver(
            self.solver_type, state_map, data.scf_info, data.mo_space_info, data.options
        )
        mcscf = make_mcscf_two_step(
            data.active_space_solver, data.state_weights_map, data.scf_info, data.options, data.mo_space_info, data.ints
        )
        energy = mcscf.compute_energy()

        rdm1=mcscf.get_rdm1()
        rdm2=mcscf.get_rdm2()

        data.results.add("mcscf energy", energy, "MCSCF energy", "hartree")
        data.results.add("mcscf rdm1", rdm1, "MCSCF rdm", "None")
        data.results.add("mcscf rdm2", rdm2, "MCSCF rdm", "None")

        return data
