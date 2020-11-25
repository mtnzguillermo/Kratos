from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from sys import argv

import KratosMultiphysics as Kratos
import KratosMultiphysics.FreeSurfaceApplication as KratosFreeSurface

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.process_factory import KratosProcessFactory

import problem_settings

class FreeSurfaceAnalysis(AnalysisStage):
    '''Main script for free surface fluid simulations using the navier_stokes family of python solvers.'''

    def __init__(self,model,project_parameters):
        if (type(model) != Kratos.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        self.model = model
        self.project_parameters = project_parameters

        self._GetSolver().AddVariables() # this creates the solver and adds the variables
        print("fluid solver created")

    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        self._GetSolver().Initialize()

    def RunSolutionLoop(self):
        time = 0.0
        step = 0
        next_output_time = self._GetSolver().output_dt
        while(time < self._GetSolver().final_time):

            if(step < self._GetSolver().number_of_inital_steps):
                max_Dt = self._GetSolver().initial_time_step
            else:
                max_Dt = self._GetSolver().original_max_dt
                # progressively increment the safety factor
                # in the steps that follow a reduction of it
                self._GetSolver().safety_factor = self._GetSolver().safety_factor * 1.2
                if(self._GetSolver().safety_factor > self._GetSolver().max_safety_factor):
                    self._GetSolver().safety_factor = self._GetSolver().max_safety_factor

            Dt = self._GetSolver().EstimateTimeStep(self._GetSolver().safety_factor, max_Dt)

            time = time + Dt
            self._GetSolver().main_model_part.CloneTimeStep(time)

            print("******** CURRENT TIME = ", time)

            if(step >= 3):
                self._GetSolver().Solve()

                check_dt = self._GetSolver().EstimateTimeStep(0.95, max_Dt)

                if(check_dt < Dt):
                    print("***********************************************************")
                    print("***********************************************************")
                    print("***********************************************************")
                    print("            *** REDUCING THE TIME STEP ***")
                    print("***********************************************************")
                    print("***********************************************************")
                    print("***********************************************************")

                    # we found a velocity too large! we need to reduce the time step
                    self._GetSolver().fluid_solver.ReduceTimeStep(self._GetSolver().main_model_part, time)  # this is to set the database to the value at the beginning of the step

                    safety_factor *= problem_settings.reduction_on_failure
                    reduced_dt = fluid_solver.EstimateTimeStep(safety_factor, max_Dt)

                    print("time before reduction= ", time)
                    time = time - Dt + reduced_dt
                    print("reduced time = ", time)
                    print("Dt = ", Dt)
                    print("reduced_dt = ", reduced_dt)

                    fluid_solver.fluid_solver.ReduceTimeStep(fluid_model_part, time)  # this is to set the database to the value at the beginning of the step

                    fluid_solver.Solve()

            if(time >= next_output_time):
                if not self.project_parameters.single_output_file:
                    # writing mesh
                    self._GetSolver().gid_io.InitializeMesh(time)
                    self._GetSolver().gid_io.WriteMesh((self._GetSolver().main_model_part).GetMesh())
                    self._GetSolver().gid_io.FinalizeMesh()
                    self._GetSolver().gid_io.InitializeResults(time, (self._GetSolver().main_model_part).GetMesh())

                self._GetSolver().gid_io.WriteNodalResults(Kratos.PRESSURE, self._GetSolver().main_model_part.Nodes, time, 0)
                self._GetSolver().gid_io.WriteNodalResults(Kratos.POROSITY, self._GetSolver().main_model_part.Nodes, time, 0)
                self._GetSolver().gid_io.WriteNodalResults(Kratos.VELOCITY, self._GetSolver().main_model_part.Nodes, time, 0)
                self._GetSolver().gid_io.WriteNodalResults(Kratos.DISTANCE, self._GetSolver().main_model_part.Nodes, time, 0)
                self._GetSolver().gid_io.WriteNodalResults(Kratos.PRESS_PROJ, self._GetSolver().main_model_part.Nodes, time, 0)
                self._GetSolver().gid_io.WriteNodalResults(Kratos.LIN_DARCY_COEF, self._GetSolver().main_model_part.Nodes, time, 0)
                self._GetSolver().gid_io.WriteNodalResults(Kratos.NONLIN_DARCY_COEF, self._GetSolver().main_model_part.Nodes, time, 0)
                self._GetSolver().gid_io.Flush()

                if not self.project_parameters.single_output_file:
                    self._GetSolver().gid_io.FinalizeResults()

                next_output_time = time + self._GetSolver().output_dt

                self._GetSolver().out = 0

            self._GetSolver().out = self._GetSolver().out + 1
            step = step + 1

        if self.project_parameters.single_output_file:
            self._GetSolver().gid_io.FinalizeResults()

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        #for process in self._GetListOfProcesses():
            #process.ExecuteFinalize()

        self._GetSolver().Finalize()

        Kratos.Logger.PrintInfo(self._GetSimulationName(), "Analysis -END- ")

    def _CreateSolver(self):
        return self.CreateSolver(self.model, self.project_parameters)

    def CreateSolver(self, model, solver_settings):
        from importlib import import_module
        module_full = 'KratosMultiphysics.FreeSurfaceApplication.edgebased_levelset_solver'
        solver = import_module(module_full).CreateSolver(model, solver_settings)
        return solver

if __name__ == '__main__':
    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    # if len(argv) == 2: # ProjectParameters is being passed from outside
    #     parameter_file_name = argv[1]
    # else: # using default name
    #     parameter_file_name = "ProjectParameters.json"

    # with open(parameter_file_name,'r') as parameter_file:
    #     parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = FreeSurfaceAnalysis(model,problem_settings)
    simulation.Run()
