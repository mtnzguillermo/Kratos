#include "includes/kratos_flags.h"

#include "custom_elements/embedded_fluid_element_discontinuous_edge.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/symbolic_navier_stokes.h"

#include "utilities/element_size_calculator.h"
#include "custom_utilities/embedded_discontinuous_edge_data.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

namespace Kratos {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::EmbeddedFluidElementDiscontinuousEdge(IndexType NewId)
    : EmbeddedFluidElementDiscontinuous<TBaseElement>(NewId)
{}

template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, const NodesArrayType& ThisNodes)
    : EmbeddedFluidElementDiscontinuous<TBaseElement>(NewId,ThisNodes)
{}

template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, typename Geometry<NodeType>::Pointer pGeometry)
    : EmbeddedFluidElementDiscontinuous<TBaseElement>(NewId,pGeometry)
{}

template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, typename Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties)
    : EmbeddedFluidElementDiscontinuous<TBaseElement>(NewId,pGeometry,pProperties)
{}


template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::~EmbeddedFluidElementDiscontinuousEdge()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TBaseElement >
Element::Pointer EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<EmbeddedFluidElementDiscontinuousEdge>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TBaseElement >
Element::Pointer EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::Create(
    IndexType NewId,
    typename Geometry<NodeType>::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<EmbeddedFluidElementDiscontinuousEdge>(NewId, pGeom, pProperties);
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Call the base element initialize method to set the constitutive law
    TBaseElement::Initialize(rCurrentProcessInfo);

    // Initialize the ELEMENTAL_DISTANCES variable (make it threadsafe)
    if (!this->Has(ELEMENTAL_DISTANCES)) {
        Vector zero_vector(NumNodes, 0.0);
        this->SetValue(ELEMENTAL_DISTANCES, zero_vector);
    }

    // NEW: Initialize the ELEMENTAL_EDGE_DISTANCES variable
    if (!this->Has(ELEMENTAL_EDGE_DISTANCES)) {
        Vector zero_vector(NumEdges, 0.0);
        this->SetValue(ELEMENTAL_EDGE_DISTANCES, zero_vector);
    }

    // Initialize the nodal EMBEDDED_VELOCITY variable (make it threadsafe)
    const array_1d<double,3> zero_vel = ZeroVector(3);
    for (auto &r_node : this->GetGeometry()) {
        r_node.SetLock();
        if (!r_node.Has(EMBEDDED_VELOCITY)) {
            r_node.SetValue(EMBEDDED_VELOCITY, zero_vel);
        }
        r_node.UnSetLock();
    }

    KRATOS_CATCH("");
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize){
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    if (rRightHandSideVector.size() != LocalSize){
        rRightHandSideVector.resize(LocalSize, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    EmbeddedDiscontinuousEdgeElementData data;
    data.Initialize(*this, rCurrentProcessInfo);
    this->InitializeGeometryData(data);

    // Iterate over the positive side volume integration points
    const unsigned int number_of_positive_gauss_points = data.PositiveSideWeights.size();
    for (unsigned int g = 0; g < number_of_positive_gauss_points; ++g){
        const size_t gauss_pt_index = g;
        this->UpdateIntegrationPointData(data, gauss_pt_index, data.PositiveSideWeights[g], row(data.PositiveSideN, g), data.PositiveSideDNDX[g]);
        this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
    }

    // Iterate over the negative side volume integration points
    const unsigned int number_of_negative_gauss_points = data.NegativeSideWeights.size();
    for (unsigned int g = 0; g < number_of_negative_gauss_points; ++g){
        const size_t gauss_pt_index = g + number_of_positive_gauss_points;
        this->UpdateIntegrationPointData(data, gauss_pt_index, data.NegativeSideWeights[g], row(data.NegativeSideN, g), data.NegativeSideDNDX[g]);
        this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
    }

    // If the element is cut, add the interface contributions
    if (data.IsCut()) {
        // Add the base element boundary contribution on the positive interface
        const size_t volume_gauss_points = number_of_positive_gauss_points + number_of_negative_gauss_points;
        const unsigned int number_of_positive_interface_gauss_points = data.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < number_of_positive_interface_gauss_points; ++g){
            const size_t gauss_pt_index = g + volume_gauss_points;
            this->UpdateIntegrationPointData(data, gauss_pt_index, data.PositiveInterfaceWeights[g], row(data.PositiveInterfaceN, g), data.PositiveInterfaceDNDX[g]);
            this->AddBoundaryTraction(data, data.PositiveInterfaceUnitNormals[g], rLeftHandSideMatrix, rRightHandSideVector);
        }

        // Add the base element boundary contribution on the negative interface
        const unsigned int number_of_negative_interface_gauss_points = data.NegativeInterfaceWeights.size();
        for (unsigned int g = 0; g < number_of_negative_interface_gauss_points; ++g){
            const size_t gauss_pt_index = g + volume_gauss_points + number_of_positive_interface_gauss_points;
            this->UpdateIntegrationPointData(data, gauss_pt_index, data.NegativeInterfaceWeights[g], row(data.NegativeInterfaceN, g), data.NegativeInterfaceDNDX[g]);
            this->AddBoundaryTraction(data, data.NegativeInterfaceUnitNormals[g], rLeftHandSideMatrix, rRightHandSideVector);
        }

        // Add the Nitsche Navier boundary condition implementation (Winter, 2018)
        data.InitializeBoundaryConditionData(rCurrentProcessInfo);
        this->AddNormalPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
        this->AddNormalSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, data); // NOTE: IMPLEMENT THE SKEW-SYMMETRIC ADJOINT IF IT IS NEEDED IN THE FUTURE. CREATE A IS_SKEW_SYMMETRIC ELEMENTAL FLAG.
        this->AddTangentialPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
        this->AddTangentialSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, data); // NOTE: IMPLEMENT THE SKEW-SYMMETRIC ADJOINT IF IT IS NEEDED IN THE FUTURE. CREATE A IS_SKEW_SYMMETRIC ELEMENTAL FLAG.
    } else if (data.IsIncised()) {
        // TODO: do extra stuff
        data.InitializeBoundaryConditionData(rCurrentProcessInfo);
        this->AddVelocityGradientPenalty(rLeftHandSideMatrix, rRightHandSideVector, data);
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);

    // If the element is split, integrate sigma.n over the interface
    // Note that in the ausas formulation, both interface sides need to be integrated
    if (rVariable == DRAG_FORCE) {
        EmbeddedDiscontinuousEdgeElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        // Calculate the drag force
        this->CalculateDragForce(data, rOutput);
    } else if (rVariable == DRAG_FORCE_CENTER) {
        EmbeddedDiscontinuousEdgeElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        // Calculate the drag force location
        this->CalculateDragForceCenter(data, rOutput);
    } else {
        TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <class TBaseElement>
int EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    int out = EmbeddedDiscontinuousEdgeElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    return TBaseElement::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TBaseElement>
std::string EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedFluidElementDiscontinuousEdge #" << this->Id();
    return buffer.str();
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedFluidElementDiscontinuousEdge" << Dim << "D" << NumNodes << "N"
             << std::endl
             << "on top of ";
    TBaseElement::PrintInfo(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// Operations

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::InitializeGeometryData(EmbeddedDiscontinuousEdgeElementData& rData) const
{
    rData.PositiveIndices.clear();
    rData.NegativeIndices.clear();
    rData.PositiveEdgeIndices.clear();

    // Number of positive and negative distance function values
    for (uint8_t i = 0; i < EmbeddedDiscontinuousEdgeElementData::NumNodes; ++i){
        if (rData.ElementalDistances[i] > 0.0){
            rData.NumPositiveNodes++;
            rData.PositiveIndices.push_back(i);
        } else {
            rData.NumNegativeNodes++;
            rData.NegativeIndices.push_back(i);
        }
    }

    // Number of positive edge distance ratios
    for (uint8_t i = 0; i < EmbeddedDiscontinuousEdgeElementData::NumEdges; ++i){
        if (!(rData.EdgeDistances[i] < 0.0)){
            rData.NumPositiveEdges++;
            rData.PositiveEdgeIndices.push_back(i);
        }
    }
    KRATOS_WATCH(rData.NumPositiveEdges);

    if (rData.IsCut()){
        this->DefineCutGeometryData(rData);
    } else if (rData.IsIncised()) {
        this->DefineIncisedGeometryData(rData);
    } else {
        this->DefineStandardGeometryData(rData);
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::DefineIncisedGeometryData(EmbeddedDiscontinuousEdgeElementData& rData) const
{
    // call standard geometry definition to make element work like EmbeddedFluidElementDiscontinuous
    this->DefineStandardGeometryData(rData);
    KRATOS_WATCH("is incised!");

    /* TODO
    For intersected element:
    // Auxiliary distance vector for the element subdivision utility
    Vector elemental_distances = rData.ElementalDistances;

    ModifiedShapeFunctions::Pointer p_calculator =
        EmbeddedDiscontinuousInternals::GetShapeFunctionCalculator<EmbeddedDiscontinuousElementData::Dim, EmbeddedDiscontinuousElementData::NumNodes>(
            *this,
            elemental_distances);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN,
        rData.PositiveSideDNDX,
        rData.PositiveSideWeights,
        GeometryData::GI_GAUSS_2);

    // Negative side volume
    p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeSideN,
        rData.NegativeSideDNDX,
        rData.NegativeSideWeights,
        GeometryData::GI_GAUSS_2);

    // Positive side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN,
        rData.PositiveInterfaceDNDX,
        rData.PositiveInterfaceWeights,
        GeometryData::GI_GAUSS_2);

    // Negative side interface
    p_calculator->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeInterfaceN,
        rData.NegativeInterfaceDNDX,
        rData.NegativeInterfaceWeights,
        GeometryData::GI_GAUSS_2);

    // Positive side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals,
        GeometryData::GI_GAUSS_2);

    // Negative side interface normals
    p_calculator->ComputeNegativeSideInterfaceAreaNormals(
        rData.NegativeInterfaceUnitNormals,
        GeometryData::GI_GAUSS_2);

    // Normalize the normals
    // Note: we calculate h here (and we don't use the value in rData.ElementSize)
    // because rData.ElementSize might still be uninitialized: some data classes define it at the Gauss point.
    double h = ElementSizeCalculator<Dim,NumNodes>::MinimumElementSize(this->GetGeometry());
    const double tolerance = std::pow(1e-3 * h, Dim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
    this->NormalizeInterfaceNormals(rData.NegativeInterfaceUnitNormals, tolerance);
    */
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::AddVelocityGradientPenalty(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedDiscontinuousEdgeElementData& rData) const
{
    // Obtain the previous iteration velocity solution
    /*array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);

    // Substract the embedded nodal velocity to the previous iteration solution
    const auto &r_geom = this->GetGeometry();
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        const auto &r_i_emb_vel = r_geom[i_node].GetValue(EMBEDDED_VELOCITY);
        for (unsigned int d = 0; d < Dim; ++d) {
            values(i_node * BlockSize + d) -= r_i_emb_vel(d);
        }
    }*/

    // Compute penalty coefficient for velocity gradient penalization
    const double pen_coef = this->ComputeVelGradPenaltyCoefficient(rData);
    KRATOS_WATCH(pen_coef);

    // penalty_term = integral-areaOrvol (penalty_const * SF *vel_grad)

}

template <class TBaseElement>
double EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::ComputeVelGradPenaltyCoefficient(const EmbeddedDiscontinuousEdgeElementData& rData) const
{
    // Get user-defined penalization constant
    const double c_p = rData.VelGradPenaltyCoefficient;

    // Compute the element average velocity norm // TODO: include mesh velocity?
    double v_norm = 0.0;
    for (uint8_t comp = 0; comp < Dim; ++comp){
        double aux_vel = 0.0;
        for (unsigned int j = 0; j < NumNodes; ++j){
            aux_vel += rData.Velocity(j,comp);
        }
        aux_vel /= NumNodes;
        v_norm += aux_vel*aux_vel;
    }
    v_norm = std::sqrt(v_norm);

    // Calculate penality coefficient // TODO: use effective viscosity?
    const double h = rData.ElementSize;
    const double d = rData.Dim;
    const double rho = rData.Density;
    const double mu = rData.DynamicViscosity;
    double vol = 0.0;
    const auto &r_geom = this->GetGeometry();
    if (Dim == 2){
        vol = rData.CalculateElementArea(r_geom);
    } else {
        vol = rData.CalculateElementVolume(r_geom);
    }
    const double penalty_coef = c_p * ( mu * std::pow(h, d-1.0) + rho * v_norm * std::pow(h, d) ) / vol;

    return penalty_coef;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::CalculateDragForce(
    EmbeddedDiscontinuousEdgeElementData& rData,
    array_1d<double,3>& rDragForce) const
{
    // Initialize the embedded element data
    const unsigned int number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    const unsigned int number_of_negative_gauss_points = rData.NegativeSideWeights.size();
    const size_t volume_gauss_points = number_of_positive_gauss_points + number_of_negative_gauss_points;

    if (rData.IsCut()){
        // Integrate positive interface side drag
        const unsigned int n_int_pos_gauss = rData.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {
            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points,
                rData.PositiveInterfaceWeights[g],
                row(rData.PositiveInterfaceN, g),
                rData.PositiveInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i){
                rDragForce(i) -= shear_proj(i);
            }
            rDragForce += rData.Weight * p_gauss * aux_unit_normal;
        }

        // Integrate negative interface side drag
        const unsigned int n_int_neg_gauss = rData.NegativeInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_neg_gauss; ++g) {
            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points + n_int_pos_gauss,
                rData.NegativeInterfaceWeights[g],
                row(rData.NegativeInterfaceN, g),
                rData.NegativeInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i){
                rDragForce(i) -= shear_proj(i);
            }
            rDragForce += rData.Weight * p_gauss * aux_unit_normal;
        }
    } else if (rData.IsIncised()) {
        // TODO: do stuff
        // KRATOS_WATCH("in incised CalculateDragForce");
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::CalculateDragForceCenter(
    EmbeddedDiscontinuousEdgeElementData& rData,
    array_1d<double,3>& rDragForceLocation) const
{
    const auto &r_geometry = this->GetGeometry();
    array_1d<double,3> tot_drag = ZeroVector(3);
    const unsigned int number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    const unsigned int number_of_negative_gauss_points = rData.NegativeSideWeights.size();
    const size_t volume_gauss_points = number_of_positive_gauss_points + number_of_negative_gauss_points;

    if (rData.IsCut()){
        // Get the positive interface continuous shape functions
        // We use these ones to interpolate the position of the intersection Gauss pt.
        // Note that we take advantage of the fact that the positive and negative interface Gauss pt. coincide
        Vector pos_int_continuous_weights;
        Matrix pos_int_continuous_N;
        typename EmbeddedDiscontinuousEdgeElementData::ShapeFunctionsGradientsType pos_int_continuous_DN_DX;
        auto p_continuous_sh_func_calculator = EmbeddedDiscontinuousInternals::GetContinuousShapeFunctionCalculator<Dim, NumNodes>(*this, rData.ElementalDistances);
        p_continuous_sh_func_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
            pos_int_continuous_N,
            pos_int_continuous_DN_DX,
            pos_int_continuous_weights,
            GeometryData::GI_GAUSS_2);

        // Integrate positive interface side drag
        const unsigned int n_int_pos_gauss = rData.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {
            // Obtain the Gauss pt. coordinates using the standard shape functions
            array_1d<double,3> g_coords = ZeroVector(3);
            const auto g_shape_functions = row(pos_int_continuous_N, g);
            for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
                g_coords += g_shape_functions[i_node] * r_geometry[i_node].Coordinates();
            }

            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points,
                rData.PositiveInterfaceWeights[g],
                row(rData.PositiveInterfaceN, g),
                rData.PositiveInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, 3> p_proj = rData.Weight * p_gauss * aux_unit_normal;
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i){
                tot_drag(i) -= shear_proj(i);
                rDragForceLocation(i) += g_coords(i) * p_proj(i);
                rDragForceLocation(i) -= g_coords(i) * shear_proj(i);
            }
            tot_drag += p_proj;
        }

        // Integrate negative interface side drag
        const unsigned int n_int_neg_gauss = rData.NegativeInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_neg_gauss; ++g) {
            // Obtain the Gauss pt. coordinates using the standard shape functions
            array_1d<double,3> g_coords = ZeroVector(3);
            const auto g_shape_functions = row(pos_int_continuous_N, g);
            for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
                g_coords += g_shape_functions[i_node] * r_geometry[i_node].Coordinates();
            }

            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points + n_int_pos_gauss,
                rData.NegativeInterfaceWeights[g],
                row(rData.NegativeInterfaceN, g),
                rData.NegativeInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, 3> p_proj = rData.Weight * p_gauss * aux_unit_normal;
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i){
                tot_drag(i) -= shear_proj(i);
                rDragForceLocation(i) += g_coords(i) * p_proj(i);
                rDragForceLocation(i) -= g_coords(i) * shear_proj(i);
            }
            tot_drag += p_proj;
        }

        // Divide the obtained result by the total drag
        rDragForceLocation(0) /= tot_drag(0);
        rDragForceLocation(1) /= tot_drag(1);
        if (Dim == 3) {
            rDragForceLocation(2) /= tot_drag(2);
        }
    } else if (rData.IsIncised()) {
        // TODO: do stuff
    }
}

// serializer

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement);
    // TODO save member variables? e.g. rSerializer.save("NU",mNU); - none right now
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement);
    // TODO load member variables? e.g. rSerializer.load("NU",mNU); - none right now
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace EmbeddedDiscontinuousEdgeInternals {

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<2, 3>(const Element& rElement, const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Triangle2D3AusasModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<3, 4>(const Element& rElement, const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4AusasModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<2, 3>(
    const Element& rElement,
    const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Triangle2D3ModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<3, 4>(
    const Element& rElement,
    const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4ModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class EmbeddedFluidElementDiscontinuousEdge< QSVMS< TimeIntegratedQSVMSData<2,3> > >;
template class EmbeddedFluidElementDiscontinuousEdge< QSVMS< TimeIntegratedQSVMSData<3,4> > >;

template class EmbeddedFluidElementDiscontinuousEdge< SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > >;
template class EmbeddedFluidElementDiscontinuousEdge< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}