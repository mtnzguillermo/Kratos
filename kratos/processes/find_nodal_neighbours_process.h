//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/find_global_nodal_elemental_neighbours_process.h"

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

/** 
 * @class FindNodalNeighboursProcess
 * @ingroup KratosCore
 * @brief This method allows to look for neighbours in a triangular or tetrahedral mesh
 * @details It checks the connectivity of the elements
 * @author Riccardo Rossi
*/
class KRATOS_API(KRATOS_CORE) FindNodalNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// The nodes container type
    typedef ModelPart::NodesContainerType NodesContainerType;
    
    /// The elements container type
    typedef ModelPart::ElementsContainerType ElementsContainerType;
    
    /// The node definition
    typedef Node<3> NodeType;
    
    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;
    
    /// The index type
    typedef std::size_t IndexType;
    
    /// The size type
    typedef std::size_t SizeType;
    
    /// Pointer definition of FindNodalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindNodalNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The modelpart to be processed
    */
    FindNodalNeighboursProcess(ModelPart& rModelPart) 
        : mrModelPart(rModelPart)
    {
         auto& r_comm = mrModelPart.GetCommunicator().GetDataCommunicator();
        mpNodeNeighboursCalculator = Kratos::make_unique<FindGlobalNodalNeighboursProcess>(r_comm, mrModelPart);
        mpElemNeighboursCalculator = Kratos::make_unique<FindGlobalNodalElementalNeighboursProcess>(r_comm, mrModelPart);

        KRATOS_INFO("FindNodalNeighboursProcess") << 
            R"(please call separetely FindGlobalNodalNeighboursProcess 
            and FindGlobalNodalElementalNeighboursProcess. 
            The two calculations are currently independent,
             hence memory savings can be achieved)" << std::endl;
    }
    
    /**
     * @brief Default constructor (deprecated one)
     * @param rModelPart The modelpart to be processed
     * @param AverageElements Expected number of neighbour elements per node.
     * @param AverageNodes Expected number of neighbour Nodes
    */
    FindNodalNeighboursProcess(
      ModelPart& rModelPart, 
      const SizeType AverageElements, 
      const SizeType AverageNodes
      ) : FindNodalNeighboursProcess(rModelPart)
    {
        KRATOS_WARNING("FindNodalNeighboursProcess") << "Parameters AverageElements and AverageNodes are currently ignored. This constructor will be removed on the 2 of April 2020" << std::endl;
    }
    
    /// Destructor.
    ~FindNodalNeighboursProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method esxecutes the neighbour search
     */
    void Execute() override;

    /**
     * @brief This method clears the neighbours found
     */
    void ClearNeighbours();

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
    std::string Info() const override
    {
        return "FindNodalNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindNodalNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
    
    ModelPart& mrModelPart;                                                                          /// The modelpart to be processed
    std::unique_ptr<FindGlobalNodalElementalNeighboursProcess> mpElemNeighboursCalculator = nullptr; /// FindGlobalNodalElementalNeighboursProcess pointer 
    std::unique_ptr<FindGlobalNodalNeighboursProcess> mpNodeNeighboursCalculator = nullptr;          /// FindGlobalNodalNeighboursProcess pointer 
    
    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method allows to add an unique weak pointer to a vector of global pointers
     * @param rVector The vector containing the global pointers 
     * @param Candidate The candidate to add
     */
    template< class TDataType > 
    void  AddUniqueWeakPointer(
        GlobalPointersVector<TDataType>& rVector, 
        const GlobalPointer<TDataType> Candidate
        )
    {
        auto i = rVector.begin();
        auto endit = rVector.end();
        while ( i != endit && (i)->Id() != (Candidate)->Id()) {
            i++;
        }
        if( i == endit ) {
            rVector.push_back(Candidate);
        }

    }

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
    FindNodalNeighboursProcess& operator=(FindNodalNeighboursProcess const& rOther);

    /// Copy constructor.
    //FindNodalNeighboursProcess(FindNodalNeighboursProcess const& rOther);


    ///@}

}; // Class FindNodalNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindNodalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindNodalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined 


