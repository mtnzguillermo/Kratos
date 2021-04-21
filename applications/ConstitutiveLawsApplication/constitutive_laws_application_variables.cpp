// KRATOS   ___                _   _ _         _   _             __                       _
//        / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//

#include "constitutive_laws_application_variables.h"

namespace Kratos
{
// Constitutive laws variables
KRATOS_CREATE_VARIABLE(Vector, HIGH_CYCLE_FATIGUE_COEFFICIENTS)
KRATOS_CREATE_VARIABLE(double, FATIGUE_REDUCTION_FACTOR)
KRATOS_CREATE_VARIABLE(int, LOCAL_NUMBER_OF_CYCLES)
KRATOS_CREATE_VARIABLE(double, WOHLER_STRESS)
KRATOS_CREATE_VARIABLE(double, REVERSION_FACTOR_RELATIVE_ERROR)
KRATOS_CREATE_VARIABLE(double, MAX_STRESS_RELATIVE_ERROR)
KRATOS_CREATE_VARIABLE(double, MAX_STRESS)
KRATOS_CREATE_VARIABLE(double, THRESHOLD_STRESS)
KRATOS_CREATE_VARIABLE(bool, CYCLE_INDICATOR)
KRATOS_CREATE_VARIABLE(double, CYCLES_TO_FAILURE)
KRATOS_CREATE_VARIABLE(double, TIME_INCREMENT)
KRATOS_CREATE_VARIABLE(bool, DAMAGE_ACTIVATION)
KRATOS_CREATE_VARIABLE(double, PREVIOUS_CYCLE);
KRATOS_CREATE_VARIABLE(double, CYCLE_PERIOD)
KRATOS_CREATE_VARIABLE(bool, ADVANCE_STRATEGY_APPLIED);

// Constitutive laws variables
KRATOS_CREATE_VARIABLE(bool, INELASTIC_FLAG)
KRATOS_CREATE_VARIABLE(double, INFINITY_YIELD_STRESS)
KRATOS_CREATE_VARIABLE(double, YIELD_STRESS_TENSION)
KRATOS_CREATE_VARIABLE(Vector, PLASTIC_STRAIN_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, BACK_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Matrix, PLASTIC_DEFORMATION_GRADIENT)
KRATOS_CREATE_VARIABLE(double, YIELD_STRESS_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, DILATANCY_ANGLE)
KRATOS_CREATE_VARIABLE(int, SOFTENING_TYPE)
KRATOS_CREATE_VARIABLE(int, SOFTENING_TYPE_COMPRESSION)
KRATOS_CREATE_VARIABLE(int, HARDENING_CURVE)
KRATOS_CREATE_VARIABLE(int, MAX_NUMBER_NL_CL_ITERATIONS)
KRATOS_CREATE_VARIABLE(double, VISCOUS_PARAMETER)
KRATOS_CREATE_VARIABLE(double, DELAY_TIME)
KRATOS_CREATE_VARIABLE(double, MAXIMUM_STRESS)
KRATOS_CREATE_VARIABLE(double, MAXIMUM_STRESS_POSITION)
KRATOS_CREATE_VARIABLE(double, PLASTIC_DISSIPATION_LIMIT_LINEAR_SOFTENING)
KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRESS)
KRATOS_CREATE_VARIABLE(double, FRICTION_ANGLE)
KRATOS_CREATE_VARIABLE(double, COHESION)
KRATOS_CREATE_VARIABLE(double, DAMAGE)
KRATOS_CREATE_VARIABLE(double, DAMAGE_MATRIX)
KRATOS_CREATE_VARIABLE(double, DAMAGE_FIBER)
KRATOS_CREATE_VARIABLE(double, THRESHOLD)
KRATOS_CREATE_VARIABLE(Matrix, INTEGRATED_STRESS_TENSOR)
KRATOS_CREATE_VARIABLE(Matrix, PLASTIC_STRAIN_TENSOR)
KRATOS_CREATE_VARIABLE(Matrix, BACK_STRESS_TENSOR)
KRATOS_CREATE_VARIABLE(Vector, CURVE_FITTING_PARAMETERS)
KRATOS_CREATE_VARIABLE(bool, TANGENCY_REGION2)
KRATOS_CREATE_VARIABLE(Vector, PLASTIC_STRAIN_INDICATORS)
KRATOS_CREATE_VARIABLE(double, EQUIVALENT_PLASTIC_STRAIN)
KRATOS_CREATE_VARIABLE(Vector, KINEMATIC_PLASTICITY_PARAMETERS)
KRATOS_CREATE_VARIABLE(int, KINEMATIC_HARDENING_TYPE)
KRATOS_CREATE_VARIABLE(bool, CONSIDER_PERTURBATION_THRESHOLD)
KRATOS_CREATE_VARIABLE(int, TANGENT_OPERATOR_ESTIMATION)
KRATOS_CREATE_VARIABLE(Matrix, TENSION_STRESS_TENSOR)
KRATOS_CREATE_VARIABLE(Matrix, COMPRESSION_STRESS_TENSOR)
KRATOS_CREATE_VARIABLE(Vector, TENSION_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, COMPRESSION_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, EFFECTIVE_TENSION_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, EFFECTIVE_COMPRESSION_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Matrix, CAUCHY_STRESS_TENSOR_FIBER)
KRATOS_CREATE_VARIABLE(Matrix, CAUCHY_STRESS_TENSOR_MATRIX)
KRATOS_CREATE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR_MATRIX)
KRATOS_CREATE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR_FIBER)
KRATOS_CREATE_VARIABLE(double, EXPONENTIAL_SATURATION_YIELD_STRESS)
KRATOS_CREATE_VARIABLE(double, ACCUMULATED_PLASTIC_STRAIN)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EULER_ANGLES)
KRATOS_CREATE_VARIABLE(Vector, HARDENING_MODULI_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, LAYER_EULER_ANGLES)
KRATOS_CREATE_VARIABLE(double, OGDEN_BETA_1)
KRATOS_CREATE_VARIABLE(double, OGDEN_BETA_2)
KRATOS_CREATE_VARIABLE(Vector, MULTI_LINEAR_ELASTICITY_MODULI)
KRATOS_CREATE_VARIABLE(Vector, MULTI_LINEAR_ELASTICITY_STRAINS)

// D+D- Damage Constitutive laws variables
KRATOS_CREATE_VARIABLE(double, DAMAGE_TENSION)
KRATOS_CREATE_VARIABLE(double, DAMAGE_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, THRESHOLD_TENSION)
KRATOS_CREATE_VARIABLE(double, THRESHOLD_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRESS_TENSION)
KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRESS_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY_DAMAGE_PROCESS)
KRATOS_CREATE_VARIABLE(double, PLASTIC_DAMAGE_PROPORTION)

// D+D- Damage Constitutive laws variables, additional Masonry 2D & 3D
KRATOS_CREATE_VARIABLE(double, DAMAGE_ONSET_STRESS_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, BIAXIAL_COMPRESSION_MULTIPLIER)
KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY_TENSION)
KRATOS_CREATE_VARIABLE(double, RESIDUAL_STRESS_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_C1)
KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_C2)
KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_C3)
KRATOS_CREATE_VARIABLE(double, YIELD_STRAIN_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, SHEAR_COMPRESSION_REDUCTOR)
KRATOS_CREATE_VARIABLE(double, TRIAXIAL_COMPRESSION_COEFFICIENT)
KRATOS_CREATE_VARIABLE(int, INTEGRATION_IMPLEX)
KRATOS_CREATE_VARIABLE(int, TENSION_YIELD_MODEL)

// For anisotropy + orthotrophy
// The ratios between the yield strength in the isotropic space and the anisotropic space
// at each direction in local coordinates ratio_x = ft / ft,x
KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS( ISOTROPIC_ANISOTROPIC_YIELD_RATIO );
KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS( ORTHOTROPIC_ELASTIC_CONSTANTS );
KRATOS_CREATE_VARIABLE( double, SERIAL_PARALLEL_EQUILIBRIUM_TOLERANCE );



}
