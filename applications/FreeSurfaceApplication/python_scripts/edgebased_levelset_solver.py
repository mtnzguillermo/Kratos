from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FreeSurfaceApplication as KratosFreeSurface

def CreateSolver(model, custom_settings):
    return EdgeBasedLevelSetSolver(model, custom_settings)

class EdgeBasedLevelSetSolver(PythonSolver):

    def __init__(self, model, settings):

        self.model = model
        self.settings = settings

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = None
        self.condition_name = None
        self.min_buffer_size = 3

        # Either retrieve the model part from the model or create a new one
        self.model_part_name = str(self.settings.problem_name)

        if self.model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(self.model_part_name):
            self.main_model_part = self.model.GetModelPart(self.model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(self.model_part_name)

        self.domain_size = self.settings.domain_size

        if self.domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.AUX_INDEX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESS_PROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POROSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIAMETER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LIN_DARCY_COEF)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NONLIN_DARCY_COEF)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRUCTURE_VELOCITY)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Adjoint fluid solver variables added correctly.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE,self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("EdgeBasedLevelSetSolver", "EdgeBasedLevelSet solver DOFs added correctly.")

    def Initialize(self):
        print("entered in EdgeBasedLevelSetSolver python constructor")

        # reading the fluid part
        gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
        multifile = KratosMultiphysics.MultiFileFlag.MultipleFiles
        deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

        # selecting output format
        if self.settings.print_layers:
            self.gid_io = KratosFreeSurface.EdgebasedGidIO(self.settings.problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
        else:
            self.gid_io = KratosMultiphysics.GidIO(self.settings.problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

        model_part_io_fluid = KratosMultiphysics.ModelPartIO(self.settings.problem_name)
        model_part_io_fluid.ReadModelPart(self.main_model_part)

        # setting up the buffer size: SHOULD BE DONE AFTER READING!!!
        self.main_model_part.SetBufferSize(2)

        # adding dofs
        self.AddDofs()

        # we assume here that all of the internal nodes are marked with a negative distance
        # set the distance of all of the internal nodes to a small value
        small_value = 0.0001
        n_active = 0
        for node in self.main_model_part.Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if(dist < 0.0):
                n_active = n_active + 1
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -small_value)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, small_value)

        if(n_active == 0):
            raise "ERROR. At least one node has to be initialized with a distance lesser than 0"

        # make sure that the porosity is not zero on any node (set by default to fluid only)
        for node in self.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.POROSITY) == 0.0):
                node.SetSolutionStepValue(KratosMultiphysics.POROSITY, 0, 1.0)
            if(node.GetSolutionStepValue(KratosMultiphysics.DIAMETER) == 0.0):
                node.SetSolutionStepValue(KratosMultiphysics.DIAMETER, 0, 1.0)

        # data of the problem
        # constructing the solver
        self.body_force = KratosMultiphysics.Vector(3)
        self.body_force[0] = self.settings.body_force_x
        self.body_force[1] = self.settings.body_force_y
        self.body_force[2] = self.settings.body_force_z
        if(self.body_force[0] == 0.0 and self.body_force[1] == 0.0 and self.body_force[2] == 0.0):
            raise "ERROR. Body Force cannot be a ZERO VECTOR"

        self.density = self.settings.density
        self.viscosity = self.settings.viscosity
        self.stabdt_pressure_factor = 1.0
        self.stabdt_convection_factor = 0.01
        self.use_mass_correction = True
        self.redistance_frequency = 5
        self.step = 0
        self.extrapolation_layers = 5
        self.tau2_factor = 0.0
        self.edge_detection_angle = 45.0
        self.assume_constant_pressure = True
        self.timer = KratosMultiphysics.Timer()

        self.use_parallel_distance_calculation = False
        # 0 = None; 1 = Ergun; 2 = Custom A y B;
        self.compute_porous_resistance_law = 0

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(
            self.main_model_part, number_of_avg_elems, number_of_avg_nodes)
        (self.neighbour_search).Execute()

        # erase isolated notes
        eliminate_isolated = KratosMultiphysics.EliminateIsolatedNodesProcess(self.main_model_part)
        eliminate_isolated.Execute()

        # definition of the solvers
#        pDiagPrecond = DiagonalPreconditioner()
#        self.pressure_linear_solver =  CGSolver(1e-3, 5000,pDiagPrecond)
 #       self.pressure_linear_solver =  CGSolver(1e-3, 5000)

#        pDiagPrecond = DiagonalPreconditioner()
#        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
        self.pressure_linear_solver = KratosMultiphysics.BICGSTABSolver(1e-3, 5000)

        # initializing the press proj to -body_force
        press_proj_init = KratosMultiphysics.Vector(3)
        press_proj_init[0] = self.body_force[0] * self.density
        press_proj_init[1] = self.body_force[1] * self.density
        press_proj_init[2] = self.body_force[2] * self.density
        for node in self.main_model_part.Nodes:
            eps = node.GetSolutionStepValue(KratosMultiphysics.POROSITY)
            node.SetSolutionStepValue(KratosMultiphysics.PRESS_PROJ, 0, press_proj_init * eps)

        self.redistance_frequency = self.settings.redistance_frequency
        self.extrapolation_layers = int(self.settings.extrapolation_layers)
        self.stabdt_pressure_factor = self.settings.stabdt_pressure_factor
        self.stabdt_convection_factor = self.settings.stabdt_convection_factor
        self.use_mass_correction = self.settings.use_mass_correction
        self.tau2_factor = self.settings.tau2_factor
        self.edge_detection_angle = self.settings.edge_detection_angle
        self.assume_constant_pressure = self.settings.assume_constant_pressure
        self.compute_porous_resistance_law = int(self.settings.compute_porous_resistance_law)  # 0 = None; 1 = Ergun; 2 = Custom;
        # print "compute_porous_resistance_law   ", fluid_solver.compute_porous_resistance_law
        # using MKLPardisosolver ----> it has to be compiled in kratos!!
        # fluid_solver.pressure_linear_solver = MKLPardisoSolver()

        # build the edge data structure
        if(self.domain_size == 2):
            self.matrix_container = KratosFreeSurface.MatrixContainer2D()
        else:
            self.matrix_container = KratosFreeSurface.MatrixContainer3D()
        self.matrix_container.ConstructCSRVector(self.main_model_part)
        self.matrix_container.BuildCSRData(self.main_model_part)
        # for 3D problems we need to evaluate the condition's neighbours
        if(self.domain_size == 3):
            self.condition_neighbours_finder = KratosMultiphysics.FindConditionsNeighboursProcess(
                self.main_model_part, self.domain_size, 10)
            self.condition_neighbours_finder.Execute()
        # constructing the solver
        if(self.domain_size == 2):
            if(self.use_parallel_distance_calculation == False):
                self.distance_utils = KratosMultiphysics.SignedDistanceCalculationUtils2D()
            else:
                self.distance_utils = KratosMultiphysics.ParallelDistanceCalculator2D()

            self.fluid_solver = KratosFreeSurface.EdgeBasedLevelSet2D(
                self.matrix_container,
                self.main_model_part,
                self.viscosity,
                self.density,
                self.body_force,
                self.use_mass_correction,
                self.edge_detection_angle,
                self.stabdt_pressure_factor,
                self.stabdt_convection_factor,
                self.edge_detection_angle,
                self.assume_constant_pressure)
        else:
            if(self.use_parallel_distance_calculation == False):
                self.distance_utils = KratosMultiphysics.SignedDistanceCalculationUtils3D()
            else:
                self.distance_utils = KratosMultiphysics.ParallelDistanceCalculator3D()

            self.fluid_solver = KratosFreeSurface.EdgeBasedLevelSet3D(
                self.matrix_container,
                self.main_model_part,
                self.viscosity,
                self.density,
                self.body_force,
                self.use_mass_correction,
                self.edge_detection_angle,
                self.stabdt_pressure_factor,
                self.stabdt_convection_factor,
                self.edge_detection_angle,
                self.assume_constant_pressure)
#
        self.max_edge_size = self.distance_utils.FindMaximumEdgeSize(
            self.main_model_part)
        self.distance_size = self.max_edge_size * 3.0
        print("###################### max distance = ", self.distance_size)

        self.fluid_solver.SetShockCapturingCoefficient(0.0)

#        self.reorder = True
#        self.distance_tools = BodyDistanceCalculationUtils()
        nneg = 0
        npos = 0
        for node in self.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0):
                nneg = nneg + 1
            else:
                npos = npos + 1

        print("nneg=", nneg)
        print("npos=", npos)

        self.fluid_solver.Initialize()
        nneg = 0
        npos = 0
        for node in self.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0):
                nneg = nneg + 1
            else:
                npos = npos + 1

        print("nneg=", nneg)
        print("npos=", npos)
        self.Redistance()

#        for node in self.main_model_part.Nodes:
#            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
#            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,1,dist)
#        self.Redistance()

        if(self.settings.wall_law_y > 1e-10):
            self.fluid_solver.ActivateWallResistance(problem_settings.wall_law_y)

        # settings to be changed
        self.max_Dt = self.settings.max_time_step
        self.initial_Dt = 0.001 * self.max_Dt
        self.final_time = self.settings.max_time
        self.output_dt = self.settings.output_dt
        self.safety_factor = self.settings.safety_factor

        self.number_of_inital_steps = self.settings.number_of_inital_steps
        self.initial_time_step = self.settings.initial_time_step
        self.out = 0

        self.original_max_dt = self.max_Dt

        # mesh to be printed
        if self.settings.single_output_file:
            self.mesh_name = 0.0
            self.gid_io.InitializeMesh(self.mesh_name)
            self.gid_io.WriteMesh(self.main_model_part.GetMesh())
            self.gid_io.FinalizeMesh()
            self.gid_io.Flush()

            self.gid_io.InitializeResults(self.mesh_name, self.main_model_part.GetMesh())

        self.max_safety_factor = self.safety_factor

        print("**********************************************")
        print("finished EdgeBasedLevelSetSolver initialize")

    #
    def Redistance(self):
        if(self.use_parallel_distance_calculation == False):
            self.distance_utils.CalculateDistances(
                self.main_model_part,
                KratosMultiphysics.DISTANCE,
                self.distance_size)
        else:
            print("max distance", self.distance_size)
            print("max extrapolation layers", self.extrapolation_layers)
            self.distance_utils.CalculateDistances(
                self.main_model_part,
                KratosMultiphysics.DISTANCE,
                KratosMultiphysics.NODAL_AREA,
                self.extrapolation_layers,
                self.distance_size)

    #
    #
    def FluidOnlySolve(self):
        if (self.extrapolation_layers < 3):
            print("insufficient number of extrapolation layers. Minimum is 3")
            raise ValueError

        print("entered in EdgeBasedLevelSetSolver fluid only solve")
        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)

        (self.fluid_solver).SolveStep1()
        (self.fluid_solver).SolveStep2(self.pressure_linear_solver)
        (self.fluid_solver).SolveStep3()

        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)
        print("finished EdgeBasedLevelSetSolver fluid only solve")

    #
    #
    def Solve(self):
        if (self.extrapolation_layers < 3):
            print("insufficient number of extrapolation layers. Minimum is 3")
            raise ValueError

        self.timer.Start("Calculate Porous Resistance Law")
        (self.fluid_solver).CalculatePorousResistanceLaw(
            self.compute_porous_resistance_law)
        self.timer.Stop("Calculate Porous Resistance Law")

        self.timer.Start("Update Fixed Velocity Values")
        (self.fluid_solver).UpdateFixedVelocityValues()
        self.timer.Stop("Update Fixed Velocity Values")

        self.timer.Start("Extrapolate Values")
        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)
        self.timer.Stop("Extrapolate Values")

        # convect levelset function
       # self.convection_solver.Solve();
        self.timer.Start("Convect Distance")
        (self.fluid_solver).ConvectDistance()
        self.timer.Stop("Convect Distance")

        if(self.step == self.redistance_frequency):
            self.timer.Start("Redistance")
            self.Redistance()
            self.timer.Stop("Redistance")
            self.step = 0
            # print "redistance was executed"
        self.step += 1

        # solve fluid
        self.timer.Start("Solve Step 1")
        (self.fluid_solver).SolveStep1()
        self.timer.Stop("Solve Step 1")
        self.timer.Start("Solve Step 2")
        (self.fluid_solver).SolveStep2(self.pressure_linear_solver)
        self.timer.Stop("Solve Step 2")
        self.timer.Start("Solve Step 3")
        (self.fluid_solver).SolveStep3()
        self.timer.Stop("Solve Step 3")

# if(self.step == self.redistance_frequency):
# self.Redistance()
# self.step = 0
# print "redistance was executed"
# self.step += 1

    #
#
#    def Solve(self):
#        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)
#
# convect levelset function
#        (self.fluid_solver).ConvectDistance()
#
#        convection_success = (self.fluid_solver).CheckDistanceConvection()
#        if(convection_success == False):
# time step reduction is needed
# print "############### distance convection failed!! ###############"
#            return False
#            errrrrrrrrr
#
# solve fluid
#        (self.fluid_solver).SolveStep1();
#        (self.fluid_solver).SolveStep2(self.pressure_linear_solver);
#        (self.fluid_solver).SolveStep3();
#
#        if(self.step == self.redistance_frequency):
#            self.Redistance()
#            self.step = 0
#            print "redistance was executed"
#        self.step += 1
#
#        return True

    #
    #
    def EstimateTimeStep(self, safety_factor, max_Dt):
        dt = (self.fluid_solver).ComputeTimeStep(safety_factor, max_Dt)

        if(dt > max_Dt):
            dt = max_Dt

        # print dt

        return dt

    #
    #
    def EstimateBoundedTimeStep(self, safety_factor, max_Dt):
        dt = (self.fluid_solver).ComputeBoundedTimeStep(safety_factor, max_Dt)

        if(dt > max_Dt):
            dt = max_Dt

        # print dt

        return dt

    #
    #
    def CalculateInitialPressureDistribution(self):
        # prova!
        dt_aux = 1e-6
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, dt_aux)
        (self.fluid_solver).SolveStep1()
        aaa = KratosMultiphysics.Vector(3)
        aaa[0] = self.body_force[0] * dt_aux
        aaa[1] = self.body_force[1] * dt_aux
        aaa[2] = self.body_force[2] * dt_aux
        for node in self.main_model_part.Nodes:
            if(node.IsFixed(KratosMultiphysics.VELOCITY_X) == False):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, aaa)

#        for node in self.main_model_part.Nodes:
#            print node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)

        (self.fluid_solver).SolveStep2(self.pressure_linear_solver)

        zero = Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0
        for node in self.main_model_part.Nodes:
            if(node.IsFixed(KratosMultiphysics.VELOCITY_X) == False):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, zero)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.0)
