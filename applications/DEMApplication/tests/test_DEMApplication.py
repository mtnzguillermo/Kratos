import KratosMultiphysics
import KratosMultiphysics.DEMApplication as DEMApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import Logger

import test_guis
import test_kinematic_constraints
import test_particle_creator_destructor
import test_wall_creator_destructor
import test_analytics
import test_glued_particles
import test_restart
import test_DEM_2D
import test_DEM_3D_contact
import test_DEM_2D_contact
import test_DEM_3D_restitution
import test_DEM_2D_restitution
import test_DEM_3D_continuum
import test_DEM_2D_inlet
import test_DEM_3D_inlet
import test_DEM_2D_control_module
import test_post_process
import test_friction_decay
import test_forces_and_moments
import test_random_variable
import test_DEM_schemes

def AssembleTestSuites():

    ''' Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    #smallSuite.addTest(test_guis.TestGUIs("test_GUIs_1"))
    #smallSuite.addTest(test_guis.TestGUIs("test_GUIs_2"))
    # smallSuite.addTest(test_kinematic_constraints.TestKinematicConstraints("test_KinematicConstraints_1"))
    # smallSuite.addTest(test_particle_creator_destructor.TestParticleCreatorDestructor("test_CreateSphericParticle1"))
    # smallSuite.addTest(test_particle_creator_destructor.TestParticleCreatorDestructor("test_CreateSphericParticle2"))
    # smallSuite.addTest(test_wall_creator_destructor.TestWallCreatorDestructor("test_CreateWallTriangle"))
    # smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_1"))
    # #smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_2"))
    # #smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_3"))
    # smallSuite.addTest(test_glued_particles.TestGluedParticles("test_Glued_Particles_1"))
    # smallSuite.addTest(test_DEM_2D.TestDEM2D("test_DEM2D_1"))
    # smallSuite.addTest(test_DEM_3D_contact.TestDEM3DContact("test_DEM3D_contact"))
    # smallSuite.addTest(test_DEM_2D_contact.TestDEM2DContact("test_DEM2D_contact"))

    smallSuite.addTest(test_DEM_2D_inlet.TestDEM2DInlet("test_DEM2D_inlet"))
    smallSuite.addTest(test_DEM_3D_inlet.TestDEM3DInlet("test_DEM3D_inlet"))
    smallSuite.addTest(test_random_variable.TestRandomVariable("test_random_variable"))
    smallSuite.addTest(test_random_variable.TestRandomVariable("test_random_variable"))

    smallSuite.addTest(test_DEM_schemes.TestDEMSchemes("test_ForwardEuler"))
    smallSuite.addTest(test_DEM_schemes.TestDEMSchemes("test_Taylor"))
    smallSuite.addTest(test_DEM_schemes.TestDEMSchemes("test_Symplectic"))
    smallSuite.addTest(test_DEM_schemes.TestDEMSchemes("test_Verlet"))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']

    nightSuite.addTest(test_restart.TestRestartOneBall("test_execution"))
    nightSuite.addTest(test_restart.TestRestartTwoBalls("test_execution"))
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nightly and you can use to validate
    validationSuite = suites['validation']

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        smallSuite
        #KratosUnittest.TestLoader().loadTestsFromTestCases([])
    )

    return suites

if __name__ == '__main__':
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
    KratosUnittest.runTests(AssembleTestSuites())
