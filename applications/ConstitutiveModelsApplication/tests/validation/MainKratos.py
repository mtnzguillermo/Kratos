import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication
import MainMaterial

MainMaterial.Solution("shear_parameters.json","shear_materials.json").Run()

