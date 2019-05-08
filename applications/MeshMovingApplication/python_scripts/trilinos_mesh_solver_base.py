from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import baseclass
from KratosMultiphysics.MeshMovingApplication.mesh_solver_base import MeshSolverBase


class TrilinosMeshSolverBase(MeshSolverBase):
    def __init__(self, model, custom_settings):
        super(TrilinosMeshSolverBase, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMeshSolverBase]:: Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "linear_solver_settings" : {
                "solver_type" : "amgcl",
                "smoother_type":"ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation",
                "max_iteration": 200,
                "provide_coordinates": false,
                "gmres_krylov_space_dimension": 100,
                "verbosity" : 0,
                "tolerance": 1e-7,
                "scaling": false,
                "block_size": 1,
                "use_block_matrices_if_possible" : true,
                "coarse_enough" : 5000
            }
        }""")
        this_defaults.AddMissingParameters(super(TrilinosMeshSolverBase, cls).GetDefaultSettings())
        return this_defaults

    #### Public user interface functions ####

    def AddVariables(self):
        super(TrilinosMeshSolverBase, self).AddVariables()
        self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMeshSolverBase]:: Variables ADDED.")

    def ImportModelPart(self):
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMeshSolverBase]:: ", "Importing model part.")
        from KratosMultiphysics.TrilinosApplication.trilinos_import_model_part_utility import TrilinosImportModelPartUtility
        self.trilinos_model_part_importer = TrilinosImportModelPartUtility(self.mesh_model_part, self.settings)
        self.trilinos_model_part_importer.ImportModelPart()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMeshSolverBase]:: ", "Finished importing model part.")

    def PrepareModelPart(self):
        super(TrilinosMeshSolverBase, self).PrepareModelPart()
        # Construct the mpi-communicator
        self.trilinos_model_part_importer.CreateCommunicators()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMeshSolverBase]::", "ModelPart prepared for Solver.")

    def Finalize(self):
        super(TrilinosMeshSolverBase, self).Finalize()
        self.get_mesh_motion_solving_strategy().Clear() # needed for proper finalization of MPI

    #### Specific internal functions ####

    def get_communicator(self):
        if not hasattr(self, '_communicator'):
            self._communicator = TrilinosApplication.CreateCommunicator()
        return self._communicator

    #### Private functions ####

    def _create_linear_solver(self):
        from KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory import ConstructSolver
        return ConstructSolver(self.settings["mesh_motion_linear_solver_settings"])

    def _create_mesh_motion_solving_strategy(self):
        raise Exception("Mesh motion solver must be created by the derived class.")