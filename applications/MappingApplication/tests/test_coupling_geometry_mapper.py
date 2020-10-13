import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
from KratosMultiphysics import KratosUnittest
data_comm = KM.DataCommunicator.GetDefault()
import os

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "mdpa_files", file_name)

class TestCouplingGeometryMapper(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(self):
        SetDefaultMappingParameters(self)
        SetupModelParts(self)
        CreateMapper(self)

    def test_map_displacements(self):
        reference_displacement = 1.0
        SetConstantVariable(self.interface_model_part_origin,KM.DISPLACEMENT,reference_displacement)
        self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
        mapped_results = GetInterfaceResult(self.interface_model_part_destination,KM.DISPLACEMENT)
        reference_result = [1.0, 1.0, 1.0, 0.9999999999999999, 0.9999999999999999, 0.9999999999999999, 1.0000000000000004, 1.0000000000000004, 1.0000000000000004, 0.9999999999999998, 0.9999999999999998, 0.9999999999999998, 1.0000000000000004, 1.0000000000000004, 1.0000000000000004]
        self.assertVectorAlmostEqual(mapped_results,reference_result)

    def test_inverse_map_forces(self):
        reference_force = 1.0
        SetConstantVariable(self.interface_model_part_destination,KM.FORCE,reference_force)
        self.mapper.InverseMap(KM.FORCE, KM.FORCE,KratosMapping.Mapper.USE_TRANSPOSE)
        mapped_results = GetInterfaceResult(self.interface_model_part_origin,KM.FORCE)
        reference_result = [0.2501913843336516, 0.2501913843336516, 0.2501913843336516, 1.292449164738119, 1.292449164738119, 1.292449164738119, 0.7004426959773076, 0.7004426959773076, 0.7004426959773076, 0.8902247629970814, 0.8902247629970814, 0.8902247629970814, 0.9474688054530651, 0.9474688054530651, 0.9474688054530651, 0.919223186500775, 0.919223186500775, 0.919223186500775]
        self.assertVectorAlmostEqual(mapped_results,reference_result)


class TestDualMortarCouplingGeometryMapper(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.mapper_parameters = KM.Parameters("""{
            "mapper_type": "coupling_geometry",
            "echo_level" : 0,
            "precompute_mapping_matrix" : false,
			"dual_mortar": true,
			"consistency_scaling" : true,
			"modeler_name" : "MappingGeometriesModeler",
            "modeler_parameters":{
						"origin_model_part_name" : "origin",
						"destination_model_part_name" : "destination",
						"is_interface_sub_model_parts_specified" : true,
						"origin_interface_sub_model_part_name" : "origin.line_tri",
						"destination_interface_sub_model_part_name" : "destination.line_quad"
					}
        }""")

        SetupModelParts(self)
        CreateMapper(self)

    def test_dual_mortar(self):
        reference_displacement = 1.0
        SetConstantVariable(self.interface_model_part_origin,KM.DISPLACEMENT,reference_displacement)
        self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
        mapped_results = GetInterfaceResult(self.interface_model_part_destination,KM.DISPLACEMENT)
        reference_result = [1.0, 1.0, 1.0, 0.9999999999999999, 0.9999999999999999, 0.9999999999999999, 1.0000000000000004, 1.0000000000000004, 1.0000000000000004, 0.9999999999999998, 0.9999999999999998, 0.9999999999999998, 1.0000000000000004, 1.0000000000000004, 1.0000000000000004]
        self.assertVectorAlmostEqual(mapped_results,reference_result)


class TestComputeMappingMatrixCouplingGeometryMapper(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.mapper_parameters = KM.Parameters("""{
            "mapper_type": "coupling_geometry",
            "echo_level" : 0,
            "precompute_mapping_matrix" : true,
			"dual_mortar": false,
			"consistency_scaling" : true,
			"modeler_name" : "MappingGeometriesModeler",
            "modeler_parameters":{
						"origin_model_part_name" : "origin",
						"destination_model_part_name" : "destination",
						"is_interface_sub_model_parts_specified" : true,
						"origin_interface_sub_model_part_name" : "origin.line_tri",
						"destination_interface_sub_model_part_name" : "destination.line_quad"
					}
        }""")

        SetupModelParts(self)
        CreateMapper(self)

    def test_precompute_mapping_matrix(self):
        reference_displacement = 1.0
        SetConstantVariable(self.interface_model_part_origin,KM.DISPLACEMENT,reference_displacement)
        self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
        mapped_results = GetInterfaceResult(self.interface_model_part_destination,KM.DISPLACEMENT)
        reference_result = [1.0, 1.0, 1.0, 0.9999999999999999, 0.9999999999999999, 0.9999999999999999, 1.0000000000000004, 1.0000000000000004, 1.0000000000000004, 0.9999999999999998, 0.9999999999999998, 0.9999999999999998, 1.0000000000000004, 1.0000000000000004, 1.0000000000000004]
        self.assertVectorAlmostEqual(mapped_results,reference_result)



def SetupModelParts(self):
    self.model = KM.Model()
    self.model_part_origin = self.model.CreateModelPart("origin")
    self.model_part_destination = self.model.CreateModelPart("destination")

    self.model_part_origin.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    self.model_part_origin.AddNodalSolutionStepVariable(KM.FORCE)

    self.model_part_destination.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    self.model_part_destination.AddNodalSolutionStepVariable(KM.FORCE)

    origin_mdpa_file_name = "cube_tri"
    destination_mdpa_file_name = "cube_quad"

    ReadModelPart(self.model_part_origin, origin_mdpa_file_name)
    ReadModelPart(self.model_part_destination, destination_mdpa_file_name)


def ReadModelPart(model_part, mdpa_file_name):
    import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER
    KM.ModelPartIO(GetFilePath(mdpa_file_name), import_flags).ReadModelPart(model_part)

def SetDefaultMappingParameters(self):
    self.mapper_parameters = KM.Parameters("""{
            "mapper_type": "coupling_geometry",
            "echo_level" : 0,
            "precompute_mapping_matrix" : false,
			"dual_mortar": false,
			"consistency_scaling" : true,
			"modeler_name" : "MappingGeometriesModeler",
            "modeler_parameters":{
						"origin_model_part_name" : "origin",
						"destination_model_part_name" : "destination",
						"is_interface_sub_model_parts_specified" : true,
						"origin_interface_sub_model_part_name" : "origin.line_tri",
						"destination_interface_sub_model_part_name" : "destination.line_quad"
					}
        }""")

def CreateMapper(self):
    self.mapper_type = self.mapper_parameters["mapper_type"].GetString()

    origin_interface_string = self.mapper_parameters["modeler_parameters"]["origin_interface_sub_model_part_name"].GetString()
    self.interface_model_part_origin =self.model.GetModelPart(origin_interface_string)

    dest_interface_string = self.mapper_parameters["modeler_parameters"]["destination_interface_sub_model_part_name"].GetString()
    self.interface_model_part_destination =self.model.GetModelPart(dest_interface_string)

    if data_comm.IsDistributed():
        self.mapper = KratosMapping.MapperFactory.CreateMPIMapper(
            self.model_part_origin, self.model_part_destination, self.mapper_parameters)
    else:
        self.mapper = KratosMapping.MapperFactory.CreateMapper(
            self.model_part_origin, self.model_part_destination, self.mapper_parameters)

def SetConstantVariable(model_part,variable,reference_value):
    for node in model_part.Nodes:
        var_x = reference_value
        var_y = reference_value
        var_z = reference_value
        node.SetSolutionStepValue(variable, KM.Vector([var_x, var_y, var_z]))

def GetInterfaceResult(model_part,variable):
    var_vector = []
    for node in model_part.Nodes:
        nodal_result = node.GetSolutionStepValue(variable)
        for dof in nodal_result:
            var_vector.append(dof)
    return var_vector


if __name__ == '__main__':
    KratosUnittest.main()
