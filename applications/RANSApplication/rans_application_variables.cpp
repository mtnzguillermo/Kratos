//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
#include "rans_application_variables.h"

namespace Kratos
{
    // incompressible potential flow specific variables
    KRATOS_CREATE_VARIABLE( double, VELOCITY_POTENTIAL )
    KRATOS_CREATE_VARIABLE( double, PRESSURE_POTENTIAL )
    KRATOS_CREATE_VARIABLE( int, RANS_IS_INLET )
    KRATOS_CREATE_VARIABLE( int, RANS_IS_OUTLET )
    KRATOS_CREATE_VARIABLE( int, RANS_IS_STRUCTURE )

    // residual based flux corrected stabilization variables
    KRATOS_CREATE_VARIABLE( double, RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT )
    KRATOS_CREATE_VARIABLE( double, RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT )

    // algebraic flux corrected stabilization variables
    KRATOS_CREATE_VARIABLE( double, AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX )
    KRATOS_CREATE_VARIABLE( double, AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX )
    KRATOS_CREATE_VARIABLE( double, AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT )
    KRATOS_CREATE_VARIABLE( double, AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT )

    // common auxiliary variables
    KRATOS_CREATE_VARIABLE( double, RANS_AUXILIARY_VARIABLE_1 )
    KRATOS_CREATE_VARIABLE( double, RANS_AUXILIARY_VARIABLE_2 )

    // k-epsilon-high-re turbulence modelling variables
    KRATOS_CREATE_VARIABLE_WITH_TIME_DERIVATIVE(double, TURBULENT_KINETIC_ENERGY_RATE, RANS_AUXILIARY_VARIABLE_1)
    KRATOS_CREATE_VARIABLE_WITH_TIME_DERIVATIVE(double, TURBULENT_KINETIC_ENERGY, TURBULENT_KINETIC_ENERGY_RATE)

    KRATOS_CREATE_VARIABLE_WITH_TIME_DERIVATIVE(double, TURBULENT_ENERGY_DISSIPATION_RATE_2, RANS_AUXILIARY_VARIABLE_2)
    KRATOS_CREATE_VARIABLE_WITH_TIME_DERIVATIVE(double, TURBULENT_ENERGY_DISSIPATION_RATE, TURBULENT_ENERGY_DISSIPATION_RATE_2)

    KRATOS_CREATE_VARIABLE( double, TURBULENT_KINETIC_ENERGY_SIGMA )
    KRATOS_CREATE_VARIABLE( double, TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_CREATE_VARIABLE( double, TURBULENCE_RANS_C_MU )
    KRATOS_CREATE_VARIABLE( double, TURBULENCE_RANS_C1 )
    KRATOS_CREATE_VARIABLE( double, TURBULENCE_RANS_C2 )

    // k-omega turbulence modelling specific additional variables
    KRATOS_CREATE_VARIABLE_WITH_TIME_DERIVATIVE(double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, RANS_AUXILIARY_VARIABLE_2)
    KRATOS_CREATE_VARIABLE_WITH_TIME_DERIVATIVE(double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2)

    KRATOS_CREATE_VARIABLE( double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_CREATE_VARIABLE( double, TURBULENCE_RANS_BETA )
    KRATOS_CREATE_VARIABLE( double, TURBULENCE_RANS_GAMMA )

    // k-omega-sst turbulence modelling specific additional variables
    KRATOS_CREATE_VARIABLE( double, TURBULENT_KINETIC_ENERGY_SIGMA_1 )
    KRATOS_CREATE_VARIABLE( double, TURBULENT_KINETIC_ENERGY_SIGMA_2 )
    KRATOS_CREATE_VARIABLE( double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1 )
    KRATOS_CREATE_VARIABLE( double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 )
    KRATOS_CREATE_VARIABLE( double, TURBULENCE_RANS_A1 )
    KRATOS_CREATE_VARIABLE( double, TURBULENCE_RANS_BETA_1 )
    KRATOS_CREATE_VARIABLE( double, TURBULENCE_RANS_BETA_2 )
    KRATOS_CREATE_VARIABLE( double, VON_KARMAN )

    // wall function condition specific additional variables
    KRATOS_CREATE_VARIABLE( double, RANS_Y_PLUS )
    KRATOS_CREATE_VARIABLE( double, RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT )
    KRATOS_CREATE_VARIABLE( double, WALL_SMOOTHNESS_BETA )
    KRATOS_CREATE_VARIABLE( int, RANS_IS_WALL_FUNCTION_ACTIVE )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FRICTION_VELOCITY )

    // formulation specific variables
    KRATOS_CREATE_VARIABLE( std::vector<std::string>, ANALYSIS_STEPS )
    KRATOS_CREATE_VARIABLE( std::string, WALL_MODEL_PART_NAME )
    KRATOS_CREATE_VARIABLE( double, NUMBER_OF_NEIGHBOUR_CONDITIONS )

}
