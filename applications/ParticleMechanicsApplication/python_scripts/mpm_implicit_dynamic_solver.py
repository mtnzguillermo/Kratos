from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from KratosMultiphysics.ParticleMechanicsApplication.mpm_solver import MPMSolver

def CreateSolver(model, custom_settings):
    return MPMImplicitDynamicSolver(model, custom_settings)

class MPMImplicitDynamicSolver(MPMSolver):

    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(MPMImplicitDynamicSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MPMImplicitDynamicSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"   : "newmark",
            "damp_factor_m" : 0.0,
            "newmark_beta"  : 0.25
        }""")
        this_defaults.AddMissingParameters(super(MPMImplicitDynamicSolver, cls).GetDefaultSettings())
        return this_defaults

    def AddVariables(self):
        super(MPMImplicitDynamicSolver, self).AddVariables()
        self._add_dynamic_variables(self.grid_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[MPMImplicitDynamicSolver]:: ", "Variables ADDED")

    #### Private functions ####

    def _create_solution_scheme(self):
        grid_model_part = self.GetGridModelPart()
        domain_size = self._get_domain_size()
        block_size  = domain_size
        if (self.settings["pressure_dofs"].GetBool()):
            block_size += 1

        # Setting the time integration schemes
        scheme_type = self.settings["scheme_type"].GetString()
        if(scheme_type == "newmark"):
            damp_factor_m = 0.0
            newmark_beta = self.settings["newmark_beta"].GetDouble()
        elif(scheme_type == "bossak"):
            damp_factor_m = self.settings["damp_factor_m"].GetDouble()
            newmark_beta = self.settings["newmark_beta"].GetDouble()
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\""
            raise Exception(err_msg)

        is_dynamic = self._is_dynamic()

        return KratosParticle.MPMResidualBasedBossakScheme( grid_model_part,
                                                            domain_size,
                                                            block_size,
                                                            damp_factor_m,
                                                            newmark_beta,
                                                            is_dynamic)
    def _is_dynamic(self):
        return True