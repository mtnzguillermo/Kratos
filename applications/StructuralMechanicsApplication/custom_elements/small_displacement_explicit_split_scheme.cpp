// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// System includes

// External includes


// Project includes
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/small_displacement_explicit_split_scheme.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SmallDisplacementExplicitSplitScheme::SmallDisplacementExplicitSplitScheme( IndexType NewId, GeometryType::Pointer pGeometry )
    : SmallDisplacement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementExplicitSplitScheme::SmallDisplacementExplicitSplitScheme( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : SmallDisplacement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementExplicitSplitScheme::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementExplicitSplitScheme>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementExplicitSplitScheme::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementExplicitSplitScheme>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementExplicitSplitScheme::~SmallDisplacementExplicitSplitScheme()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementExplicitSplitScheme::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    SmallDisplacementExplicitSplitScheme::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementExplicitSplitScheme>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * dimension;

    // Compiting the nodal mass
    if (rDestinationVariable == NODAL_MASS ) {

        VectorType element_mass_vector(mat_size);
        this->CalculateLumpedMassVector(element_mass_vector);

        VectorType element_stiffness_vector(mat_size);
        CalculateLumpedStiffnessVector(element_stiffness_vector, rCurrentProcessInfo);

        VectorType element_damping_vector(mat_size);
        CalculateLumpedDampingVector(element_damping_vector, rCurrentProcessInfo);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            array_1d<double, 3>& r_nodal_stiffness = r_geom[i].GetValue(NODAL_DIAGONAL_STIFFNESS);
            array_1d<double, 3>& r_nodal_damping = r_geom[i].GetValue(NODAL_DIAGONAL_DAMPING);
            const IndexType index = i * dimension;

            #pragma omp atomic
            r_geom[i].GetValue(NODAL_MASS) += element_mass_vector[index];

            for (IndexType j = 0; j < dimension; ++j) {
                #pragma omp atomic
                r_nodal_stiffness[j] += element_stiffness_vector[index+j];

                #pragma omp atomic
                r_nodal_damping[j] += element_damping_vector[index+j];
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    // const auto& r_prop = this->GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType element_size = dimension * number_of_nodes;

    Vector internal_forces = ZeroVector(element_size);
    this->CalculateInternalForces(internal_forces,rCurrentProcessInfo);

    // VectorType displacement_vector(element_size);
    // GetValuesVector(displacement_vector);

    // VectorType k_a(element_size);
    // MatrixType stiffness_matrix( element_size, element_size );
    // noalias(stiffness_matrix) = ZeroMatrix(element_size,element_size);
    // this->CalculateLeftHandSide(stiffness_matrix, rCurrentProcessInfo);
    // noalias(k_a) = prod(stiffness_matrix,displacement_vector);

    // VectorType k_hat_a(element_size);
    // MatrixType non_diagonal_stiffness_matrix( element_size, element_size );
    // noalias(non_diagonal_stiffness_matrix) = ZeroMatrix(element_size,element_size);
    // this->CalculateLeftHandSide(non_diagonal_stiffness_matrix, rCurrentProcessInfo);
    // for (IndexType i = 0; i < element_size; ++i)
    //     non_diagonal_stiffness_matrix(i,i) = 0.0;
    // noalias(k_hat_a) = prod(non_diagonal_stiffness_matrix,displacement_vector);

    // VectorType element_mass_vector(element_size);
    // this->CalculateLumpedMassVector(element_mass_vector);
    // VectorType M_d_a(element_size);
    // MatrixType diagonal_mass_matrix( element_size, element_size );
    // noalias(diagonal_mass_matrix) = ZeroMatrix(element_size,element_size);
    // for (IndexType i = 0; i < element_size; ++i)
    //     diagonal_mass_matrix(i,i) = element_mass_vector[i];
    // noalias(M_d_a) = prod(diagonal_mass_matrix,displacement_vector);

    // VectorType k_d_a(element_size);
    // MatrixType diagonal_stiffness_matrix( element_size, element_size );
    // noalias(diagonal_stiffness_matrix) = ZeroMatrix(element_size,element_size);
    // for (IndexType i = 0; i < element_size; ++i)
    //     diagonal_stiffness_matrix(i,i) = stiffness_matrix(i,i);
    // noalias(k_d_a) = prod(diagonal_stiffness_matrix,displacement_vector);

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = dimension * i;
            array_1d<double, 3>& r_external_forces = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            array_1d<double, 3>& r_internal_forces = GetGeometry()[i].FastGetSolutionStepValue(NODAL_INERTIA);
            // array_1d<double, 3>& r_k_hat_a = GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_ACCELERATION);

            for (IndexType j = 0; j < dimension; ++j) {
                #pragma omp atomic
                r_external_forces[j] += rRHSVector[index + j] + internal_forces[index + j];

                #pragma omp atomic
                r_internal_forces[j] += internal_forces[index + j];

                // #pragma omp atomic
                // r_k_hat_a[j] += k_hat_a[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::CalculateLumpedDampingVector(
    VectorType& rDampingVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * dimension;

    // Clear Vector
    if (rDampingVector.size() != mat_size) {
        rDampingVector.resize(mat_size, false);
    }
    noalias(rDampingVector) = ZeroVector(mat_size);

    if (rCurrentProcessInfo[USE_CONSISTENT_MASS_MATRIX] == true) {

        // Rayleigh Damping Vector (C= alpha*M + beta*K)

        // Get Damping Coefficients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
        double alpha = 0.0;
        if( GetProperties().Has(RAYLEIGH_ALPHA) )
            alpha = GetProperties()[RAYLEIGH_ALPHA];
        else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
            alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
        double beta  = 0.0;
        if( GetProperties().Has(RAYLEIGH_BETA) )
            beta = GetProperties()[RAYLEIGH_BETA];
        else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
            beta = rCurrentProcessInfo[RAYLEIGH_BETA];

        // 1.-Calculate mass Vector:
        if (alpha > std::numeric_limits<double>::epsilon()) {
            VectorType mass_vector(mat_size);
            CalculateLumpedMassVector(mass_vector);
            for (IndexType i = 0; i < mat_size; ++i)
                rDampingVector[i] += alpha * mass_vector[i];
        }

        // 2.-Calculate Stiffness Vector:
        if (beta > std::numeric_limits<double>::epsilon()) {
            VectorType stiffness_vector(mat_size);
            CalculateLumpedStiffnessVector(stiffness_vector,rCurrentProcessInfo);
            for (IndexType i = 0; i < mat_size; ++i)
                rDampingVector[i] += beta * stiffness_vector[i];
        }

    } else {

        // Critical damping vector (C=2*xi*sqrt(K*M))

        double xi_damping  = 0.0;
        if( GetProperties().Has(XI_DAMPING) )
            xi_damping = GetProperties()[XI_DAMPING];
        else if( rCurrentProcessInfo.Has(XI_DAMPING) )
            xi_damping = rCurrentProcessInfo[XI_DAMPING];

        VectorType mass_vector(mat_size);
        CalculateLumpedMassVector(mass_vector);
        VectorType stiffness_vector(mat_size);
        CalculateLumpedStiffnessVector(stiffness_vector,rCurrentProcessInfo);
        for (IndexType i = 0; i < mat_size; ++i)
            rDampingVector[i] = 2.0 * xi_damping * sqrt(stiffness_vector[i] * mass_vector[i]);

    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::CalculateLumpedMassVector(VectorType& rMassVector) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassVector.size() != mat_size)
        rMassVector.resize( mat_size, false );

    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    // LUMPED MASS MATRIX
    const double total_mass = GetGeometry().DomainSize() * density * thickness;

    Vector lumping_factors;
    lumping_factors = GetGeometry().LumpingFactors( lumping_factors );

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const double temp = lumping_factors[i] * total_mass;
        for ( IndexType j = 0; j < dimension; ++j ) {
            IndexType index = i * dimension + j;
            rMassVector[index] = temp;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::CalculateLumpedStiffnessVector(
    VectorType& rStiffnessVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    auto& r_geom = this->GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * dimension;

    // Clear Vector
    if (rStiffnessVector.size() != mat_size) {
        rStiffnessVector.resize(mat_size, false);
    }

    MatrixType stiffness_matrix( mat_size, mat_size );
    noalias(stiffness_matrix) = ZeroMatrix(mat_size,mat_size);
    // ProcessInfo temp_process_information = rCurrentProcessInfo;
    this->CalculateLeftHandSide(stiffness_matrix, rCurrentProcessInfo);
    for (IndexType i = 0; i < mat_size; ++i)
        rStiffnessVector[i] = stiffness_matrix(i,i);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::CalculateInternalForces(
    VectorType& rInternalForces,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    // Resizing as needed the RHS
    if ( rInternalForces.size() != mat_size )
        rInternalForces.resize( mat_size, false );

    rInternalForces = ZeroVector( mat_size ); //resetting RHS

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    double int_to_reference_weight;

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        noalias( rInternalForces ) += int_to_reference_weight * prod( trans( this_kinematic_variables.B ), this_constitutive_variables.StressVector );
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

// void SmallDisplacementExplicitSplitScheme::CalculateNoDiagonalDampingMatrix(
//     MatrixType& rNoDiagonalDampingMatrix,
//     const ProcessInfo& rCurrentProcessInfo
//     )
// {
//     KRATOS_TRY

//     auto& r_geom = this->GetGeometry();
//     const SizeType dimension = r_geom.WorkingSpaceDimension();
//     const SizeType number_of_nodes = r_geom.size();
//     const SizeType mat_size = number_of_nodes * dimension;

//     // Clear matrix
//     if (rNoDiagonalDampingMatrix.size1() != mat_size || rNoDiagonalDampingMatrix.size2() != mat_size) {
//         rNoDiagonalDampingMatrix.resize(mat_size, mat_size, false);
//     }
//     rNoDiagonalDampingMatrix = ZeroMatrix(mat_size, mat_size);

//     double beta  = 0.0;
//     if( GetProperties().Has(RAYLEIGH_BETA) )
//         beta = GetProperties()[RAYLEIGH_BETA];
//     else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
//         beta = rCurrentProcessInfo[RAYLEIGH_BETA];

//     if (beta > std::numeric_limits<double>::epsilon()) {
//         MatrixType non_diagonal_stiffness_matrix( mat_size, mat_size );
//         noalias(non_diagonal_stiffness_matrix) = ZeroMatrix(mat_size,mat_size);
//         this->CalculateLeftHandSide(non_diagonal_stiffness_matrix, rCurrentProcessInfo);
//         for (IndexType i = 0; i < mat_size; ++i) {
//             non_diagonal_stiffness_matrix(i,i) = 0.0;
//         }
//         noalias(rNoDiagonalDampingMatrix) = beta * non_diagonal_stiffness_matrix;
//     }

//     KRATOS_CATCH("")
// }

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacement );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacement );
}

} // Namespace Kratos


