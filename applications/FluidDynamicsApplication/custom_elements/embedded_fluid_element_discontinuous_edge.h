//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla, Franziska Wahl
//  Co-authors:      Jordi Cotela
//

#ifndef KRATOS_EMBEDDED_FLUID_ELEMENT_DISCONTINUOUS_EDGE_H
#define KRATOS_EMBEDDED_FLUID_ELEMENT_DISCONTINUOUS_EDGE_H

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "modified_shape_functions/modified_shape_functions.h"

#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"
#include "custom_elements/embedded_fluid_element_discontinuous.h"

#include "custom_utilities/embedded_discontinuous_edge_data.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template< class TBaseElement >
class EmbeddedFluidElementDiscontinuousEdge : public EmbeddedFluidElementDiscontinuous<TBaseElement>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedFluidElementDiscontinuousEdge
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EmbeddedFluidElementDiscontinuousEdge);

    /// Node type (default is: Node<3>)
    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::NodeType;

    /// Definition of nodes container type, redefined from GeometryType
    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::NodesArrayType;

    /// Vector type for local contributions to the linear system
    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::VectorType;

    /// Matrix type for local contributions to the linear system
    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::MatrixType;

    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::IndexType;

    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::SizeType;

    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::EquationIdVectorType;

    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::DofsVectorType;

    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::DofsArrayType;

    /// Type for shape function values container
    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::ShapeFunctionDerivativesArrayType;

    using EmbeddedFluidElementDiscontinuous<TBaseElement>::Dim;
    using EmbeddedFluidElementDiscontinuous<TBaseElement>::NumNodes;
    using EmbeddedFluidElementDiscontinuous<TBaseElement>::BlockSize;
    using EmbeddedFluidElementDiscontinuous<TBaseElement>::LocalSize;
    using EmbeddedFluidElementDiscontinuous<TBaseElement>::StrainSize;

    using typename EmbeddedFluidElementDiscontinuous<TBaseElement>::BaseElementData;
    using EmbeddedDiscontinuousEdgeElementData = EmbeddedDiscontinuousEdgeData< BaseElementData >;

    static constexpr unsigned int NumEdges = (Dim -1) *3;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    EmbeddedFluidElementDiscontinuousEdge(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, typename Geometry<NodeType>::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, typename Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    ~EmbeddedFluidElementDiscontinuousEdge() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new EmbeddedFluidElementDiscontinuous element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            Properties::Pointer pProperties) const override;

    /// Create a new element of this type using given geometry
    /**
     * Returns a pointer to a new FluidElement element, created using given input
     * @param NewId the ID of the new element
     * @param pGeom a pointer to the geomerty to be used to create the element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            typename Geometry<NodeType>::Pointer pGeom,
                            Properties::Pointer pProperties) const override;

    /// Set up the element for solution.  NEW
    /** For EmbeddedFluidElementDiscontinuous, this initializes the discontinuous
     * level set (ELEMENTAL_DISTANCES) and the nodal imposed velocity (EMBEDDED_VELOCITY)
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;


    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Calculates both LHS and RHS contributions TODO: if Is_Incised
    /**
     * Computes the LHS and RHS elementar matrices. If the element is split
     * includes the contribution of the level set boundary condition imposition.
     * TODO NEW incised element
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Computes an elemental 3 components array value
    /**
     * Given a 3 components array variable, this function computes its value inside the element.
     * If the function has not implemented this variable computation, calls the base class one.
     * Drag force and drag center calculation implemented.
     * @param rVariable Variable to be computed
     * @param rOutput Reference to the output array
     * @param rCurrentProcessInfo Reference to the process info
     */
    void Calculate(
        const Variable<array_1d<double, 3>> &rVariable,
        array_1d<double, 3> &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /** INFO: no override, because use of different element data type
     * @brief Current element data structure initialization
     * This method checks if the element is intersected or incised and calls the elemental data filling methods accordingly.
     * @param rData reference to the element data structure
     */
    void InitializeGeometryData(EmbeddedDiscontinuousEdgeElementData& rData) const;

    /**
     * @brief Incised element geometry data fill
     * This method sets the data structure geometry fields (shape functions, gradients, interface normals, ...) for an
     * incised element. To do that, the modified shape functions utility is firstly created and then called to perform
     * all operations in both the positive and negative sides of the element.
     * @param rData reference to the element data structure
     *
     * TODO: fill method
     */
    void DefineIncisedGeometryData(EmbeddedDiscontinuousEdgeElementData& rData) const;

    /**
    * This method adds a penalization for velocity gradients in incised elements.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddVelocityGradientPenalty(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedDiscontinuousEdgeElementData& rData) const;

    /**
    * This method computes the penalty coefficient for the velocity gradient in incised elements
    * Works only for triangle and tetrahedron!
    * @param rData reference to element data structure
    */
    double ComputeVelGradPenaltyCoefficient(const EmbeddedDiscontinuousEdgeElementData& rData) const;


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /** TODO: if Is_Incised
     * @brief Calculates the drag force
     * For an intersected element, this method calculates the drag force.
     * Note that the drag force includes both the shear and the pressure contributions.
     * @param rData reference to the embedded elemental data
     * @param rDragForce reference to the computed drag force
     */
    void CalculateDragForce(
        EmbeddedDiscontinuousEdgeElementData& rData,
        array_1d<double,3>& rDragForce) const;

    /** TODO: if Is_Incised
     * @brief Calculates the location of the drag force
     * For an intersected element, this method calculates the drag force location.
     * Note that the drag force includes both the shear and the pressure contributions.
     * @param rData reference to the embedded elemental data
     * @param rDragForce reference to the computed drag force
     */
    void CalculateDragForceCenter(
        EmbeddedDiscontinuousEdgeElementData& rData,
        array_1d<double,3>& rDragForceLocation) const;


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    EmbeddedFluidElementDiscontinuousEdge& operator=(EmbeddedFluidElementDiscontinuousEdge const& rOther);

    /// Copy constructor.
    EmbeddedFluidElementDiscontinuousEdge(EmbeddedFluidElementDiscontinuousEdge const& rOther);

    ///@}


}; // Class EmbeddedFluidElementDiscontinuousEdge

namespace EmbeddedDiscontinuousEdgeInternals {

template <size_t TDim, size_t TNumNodes>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator(
    const Element &rElement,
    const Vector &rElementalDistances);

template <size_t TDim, size_t TNumNodes>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator(
    const Element &rElement,
    const Vector &rElementalDistances);
}

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >>(
    std::istream& rIStream,
    EmbeddedFluidElementDiscontinuousEdge<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const EmbeddedFluidElementDiscontinuousEdge<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_FLUID_ELEMENT_DISCONTINUOUS_EDGE_H
