// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_INTEGRATION_H
#define MAPPER_VERTEX_MORPHING_INTEGRATION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing.h"
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/normal_calculation_utils.h"
#include "spaces/ublas_space.h"
#include "processes/find_conditions_neighbours_process.h"
#include "utilities/math_utils.h"
#include "shape_optimization_application.h"
#include "filter_function.h"

// ==============================================================================

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

/// Short class definition.
/** Detail class definition.
*/

class MapperVertexMorphingIntegration : public MapperVertexMorphing
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperVertexMorphingIntegration
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingIntegration);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingIntegration( ModelPart& designSurface, Parameters& optimizationSettings )
        : MapperVertexMorphing(designSurface, optimizationSettings)
    {
        SetIntegrationMethod(optimizationSettings);
        FindNeighbourConditions();
    }

    /// Destructor.
    virtual ~MapperVertexMorphingIntegration()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void SetIntegrationMethod( Parameters& optimizationSettings )
    {
        int integrationMethod = optimizationSettings["design_variables"]["integration_method"].GetInt();
        mAreaWeightedNodeSum = false;
        if (integrationMethod == 0)
            mAreaWeightedNodeSum = true;
        else if (integrationMethod == 1)
            mIntegrationMethod = GeometryData::GI_GAUSS_1;
        else if (integrationMethod == 2)
            mIntegrationMethod = GeometryData::GI_GAUSS_2;
        else if (integrationMethod == 3)
            mIntegrationMethod = GeometryData::GI_GAUSS_3;
        else if (integrationMethod == 4)
            mIntegrationMethod = GeometryData::GI_GAUSS_4;
        else if (integrationMethod == 5)
            mIntegrationMethod = GeometryData::GI_GAUSS_5;
        else
        {
            std::cout << "\n> Integration method " << integrationMethod << " not valid! USING DEFAULT: 2 " << std::endl;
            mIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
    }

    // --------------------------------------------------------------------------
    void FindNeighbourConditions()
    {

            // store neighbouring information
        std::cout << "> Computing neighbour conditions ..." << std::endl;
        FindConditionsNeighboursProcess find_conditions_neighbours_process(mrDesignSurface,
                                                        mrDesignSurface.GetProcessInfo()[DOMAIN_SIZE]);
        find_conditions_neighbours_process.Execute();
    }

    // --------------------------------------------------------------------------
    void ComputeWeightForAllNeighbors(  ModelPart::NodeType& node_i,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
        {
            // Get node information
            ModelPart::NodeType& node_j = *neighbor_nodes[j_itr];

            // Get all neighbour conditions
            const WeakPointerVector<Condition>& rConditions = node_j.GetValue(NEIGHBOUR_CONDITIONS);

            // loop conditions
            for(unsigned int c_itr=0; c_itr<rConditions.size(); c_itr++)
            {
                // Get geometry of current condition
                Condition rCondition = rConditions[c_itr];
                Condition::GeometryType& geom_i = rCondition.GetGeometry();

                if (mAreaWeightedNodeSum){

                    // Computation of weight according specified weighting function
                    // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
                    double Aij = mpFilterFunction->compute_weight(node_j.Coordinates(),node_i.Coordinates());
                    Aij *= geom_i.DomainSize();

                    // Add values to list
                    list_of_weights[j_itr] += Aij;

                    // Computed for integration of weighting function later using post-scaling
                    sum_of_weights += Aij;
                }
                else
                {
                    // Get geometry information of current condition
                    unsigned int n_nodes = geom_i.size();
                    int localNodeIndex = -1;
                    for (unsigned int node_ctr=0; node_ctr<n_nodes; node_ctr++)
                    {
                        if (geom_i[node_ctr].Id() == node_j.Id())
                            localNodeIndex = node_ctr;
                    }

                    // Evaluate shape functions of design surface according specified integration method
                    const Condition::GeometryType::IntegrationPointsArrayType& integrationPoints = geom_i.IntegrationPoints(mIntegrationMethod);
                    const unsigned int numberOfIntegrationPoints = integrationPoints.size();
                    const Matrix& N_container = geom_i.ShapeFunctionsValues(mIntegrationMethod);

                    for ( unsigned int pointNumber = 0; pointNumber < numberOfIntegrationPoints; pointNumber++ )
                    {

                        // Get FEM-shape-function-value for current integration point
                        Vector N_FEM_GPi = row( N_container, pointNumber);

                        // Get gp coordinates
                        NodeType::CoordinatesArrayType gp_i_coord;
                        geom_i.GlobalCoordinates(gp_i_coord, integrationPoints[pointNumber].Coordinates());

                        // Computation of weight according specified weighting function
                        // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
                        double Aij = mpFilterFunction->compute_weight(gp_i_coord,node_i.Coordinates());

                        // multiply with evaluation of shape function at gauss point
                        Aij *= geom_i.ShapeFunctionValue(pointNumber,localNodeIndex,mIntegrationMethod);;

                        // Get weight for integration
                        Aij *= integrationPoints[pointNumber].Weight();

                        // consider jacobian
                        Aij *= geom_i.DeterminantOfJacobian(pointNumber,mIntegrationMethod);

                        // Add values to list
                        list_of_weights[j_itr] += Aij;

                        // Computed for integration of weighting function later using post-scaling
                        sum_of_weights += Aij;
                    }
                }
            }
        }
    }

    // ==============================================================================

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "MapperVertexMorphingIntegration";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperVertexMorphingIntegration";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    Element::IntegrationMethod mIntegrationMethod;
    bool mAreaWeightedNodeSum;

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      MapperVertexMorphingIntegration& operator=(MapperVertexMorphingIntegration const& rOther);

    /// Copy constructor.
//      MapperVertexMorphingIntegration(MapperVertexMorphingIntegration const& rOther);


    ///@}

}; // Class MapperVertexMorphingIntegration

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_INTEGRATION_H
