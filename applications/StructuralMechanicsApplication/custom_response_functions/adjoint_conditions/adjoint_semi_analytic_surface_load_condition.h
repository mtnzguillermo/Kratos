// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Mahmoud Sesa
//

// System includes
#if !defined(ADJOINT_SEMI_ANALYTIC_SURFACE_LOAD_CONDITION )
#define  ADJOINT_SEMI_ANALYTIC_SURFACE_LOAD_CONDITION

// System includes

// External includes

// Project includes
#include "adjoint_semi_analytic_base_condition.h"

namespace Kratos
{

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

template <typename TPrimalCondition>
class AdjointSemiAnalyticSurfaceLoadCondition
    : public AdjointSemiAnalyticBaseCondition<TPrimalCondition>
{
public:
    ///@name Type Definitions
    ///@{

    // redefine the typedefs because of templated base class
    typedef AdjointSemiAnalyticBaseCondition<TPrimalCondition> BaseType;
    typedef Node<3>  NodeType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;
    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;
    typedef typename BaseType::DofsVectorType DofsVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::IntegrationMethod IntegrationMethod;
    typedef typename BaseType::GeometryDataType GeometryDataType;

    /// Counted pointer of AdjointSemiAnalyticSurfaceLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( AdjointSemiAnalyticSurfaceLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    AdjointSemiAnalyticSurfaceLoadCondition(IndexType NewId = 0)
    : AdjointSemiAnalyticBaseCondition<TPrimalCondition>(NewId)
    {
    }

    AdjointSemiAnalyticSurfaceLoadCondition(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : AdjointSemiAnalyticBaseCondition<TPrimalCondition>(NewId, pGeometry)
    {
    }

    AdjointSemiAnalyticSurfaceLoadCondition(IndexType NewId,
                        typename GeometryType::Pointer pGeometry,
                        typename PropertiesType::Pointer pProperties)
    : AdjointSemiAnalyticBaseCondition<TPrimalCondition>(NewId, pGeometry, pProperties)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointSemiAnalyticSurfaceLoadCondition<TPrimalCondition>>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(IndexType NewId,
                              typename GeometryType::Pointer pGeometry,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointSemiAnalyticSurfaceLoadCondition<TPrimalCondition>>(
            NewId, pGeometry, pProperties);
    }

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo ) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo ) override;

    void GetValuesVector(Vector& rValues, int Step = 0 ) override;

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;


    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;


    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                         std::vector< array_1d<double, 3 > >& rValues,
                         const ProcessInfo& rCurrentProcessInfo) override
    {
        this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                          std::vector< array_1d<double, 3 > >& rOutput,
                          const ProcessInfo& rCurrentProcessInfo) override;

    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


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
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, AdjointSemiAnalyticBaseCondition<TPrimalCondition>);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AdjointSemiAnalyticBaseCondition<TPrimalCondition> );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class AdjointSemiAnalyticSurfaceLoadCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // ADJOINT_SEMI_ANALYTIC_SURFACE_LOAD_CONDITION  defined


