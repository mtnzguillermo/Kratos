//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_EMBEDDED_DISCONTINUOUS_EDGE_DATA_H)
#define KRATOS_EMBEDDED_DISCONTINUOUS_EDGE_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"
#include "custom_utilities/embedded_discontinuous_data.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< class TFluidData >
class EmbeddedDiscontinuousEdgeData : public EmbeddedDiscontinuousData<TFluidData>
{
public:

///@name Type Definitions
///@{

using typename TFluidData::NodalScalarData;
using typename TFluidData::NodalVectorData;

typedef array_1d<double, (TFluidData::Dim - 1) * 3> EdgeScalarData;

using typename EmbeddedDiscontinuousData<TFluidData>::ShapeFunctionsGradientsType;
using typename EmbeddedDiscontinuousData<TFluidData>::InterfaceNormalsType;

///@}
///@name Public Members
///@{

double VelGradPenaltyCoefficient;

constexpr static unsigned int NumEdges = (TFluidData::Dim - 1) * 3;

EdgeScalarData EdgeDistances;
size_t NumPositiveEdges;
std::vector< size_t > PositiveEdgeIndices;

///@}
///@name Public Operations
///@{

/**
 * @brief Discontinuous embedded formulation data container initialization
 * This method initializes the discontinuous embedded formulation data container. This implies to intialize
 * the base formulation data container as well as to get the elemental distances from the elemental variable
 * ELEMENTAL_DISTANCES and to get the ELEMENTAL_EDGE_DISTANCES
 * (note that this requires the ELEMENTAL_DISTANCES and ELEMENTAL_EDGE_DISTANCES to be set before this operation).
 * @see CalculateDiscontinuousEdgeDistanceToSkinProcess
 * @param rElement reference to the element that owns the data container
 * @param rProcessInfo reference to the current ProcessInfo container
 */
void Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) override
{
    EmbeddedDiscontinuousData<TFluidData>::Initialize(rElement, rProcessInfo);
    //this->FillFromElementData(EdgeDistances, ELEMENTAL_EDGE_DISTANCES, rElement);  //INFO: no matching method tin fluid_element_data
    EdgeDistances = rElement.GetValue(ELEMENTAL_EDGE_DISTANCES);

    NumPositiveEdges = 0;
}

/**
 * @brief Fills the boundary condition data fields
 * This method needs to be called in cut elements. It fills the data structure fields related to the boundary
 * condition imposition (slip length and penalty coefficient) by retrieving their value from the ProcessInfo.
 * @param rProcessInfo
 */
virtual void InitializeBoundaryConditionData(const ProcessInfo& rProcessInfo) override
{
    EmbeddedDiscontinuousData<TFluidData>::InitializeBoundaryConditionData(rProcessInfo);
    this->FillFromProcessInfo(VelGradPenaltyCoefficient, V_GRAD_PENALTY_COEFFICIENT, rProcessInfo);
}

/**
 * @brief Discontinous embedded formulation data container check
 * Simple discontinuous embedded formulation data container check. The base formulation data container is
 * checked as well. Returns 0 if the check process succeeds.
 * @param rElement reference to the element that owns the data container
 * @param rProcessInfo reference to the current ProcessInfo container
 * @return int returns 0 if the check process succeeds
 */
static int Check(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_CHECK_VARIABLE_KEY(SLIP_LENGTH);
    KRATOS_CHECK_VARIABLE_KEY(PENALTY_COEFFICIENT);
    KRATOS_CHECK_VARIABLE_KEY(V_GRAD_PENALTY_COEFFICIENT);
    KRATOS_CHECK_VARIABLE_KEY(ELEMENTAL_DISTANCES);
    KRATOS_CHECK_VARIABLE_KEY(ELEMENTAL_EDGE_DISTANCES);

    int out = TFluidData::Check(rElement,rProcessInfo);
    return out;
}

/** TODO: what is the definition of incised we consider?? min 2 edges in 3D?
 * @brief Checks if the current element is incised (not intersected completely)
 * Checks if the current element is incised by checking the number of positive distance edges.
 * @return true if the element is incised
 * @return false if the element is not incised
 */
const bool IsIncised()
{
    return (0 < NumPositiveEdges && NumPositiveEdges < TFluidData::Dim);
}

/** TODO: does a function like this already exist?
 * @brief Calculates area of a triangular element
 * @param rGeometry the geometry of the element
 * @return area of triangular element
 */
const double CalculateTriangleArea(const Geometry<Node<3> >& rGeometry) const
{
    // calculate edge vectors
    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();

    //calculate components of cross product
    const double c3 = x10 * y20 - y10 * x20;

    // return area
    return ( 0.5 * std::abs(c3) );
}

/** TODO: does a function like this already exist?
 * @brief Calculates volume of a tetrahedral element
 * @param rGeometry the geometry of the element
 * @return volume of tetrahedral element
 */
const double CalculateTetrahedronVolume(const Geometry<Node<3> >& rGeometry) const
{
    // calculate edge vectors
    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double z20 = rGeometry[2].Z() - rGeometry[0].Z();

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();

    //calculate components of cross product
    const double c1 = y10 * z20 - z10 * y20;
    const double c2 = z10 * x20 - x10 * z20;
    const double c3 = x10 * y20 - y10 * x20;

    //calculate dot product
    const double d = c1 * x30 + c2 * y30 + c3 * z30;

    // return volume
    return ( 1.0/6.0 * std::abs(d) );
}


///@}

};

///@}

///@}

}

#endif