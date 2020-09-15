//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//


#if !defined(KRATOS_PARALLEL_HELPERS_H_INCLUDED)
#define KRATOS_PARALLEL_HELPERS_H_INCLUDED


// System includes

// External includes

// Project includes


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class ParallelHelpers
{
public:
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParallelHelpers() = delete;

    ///@}
    ///@name Operations
    ///@{

    static int GetNumThreads();

    ///@}

}; // Class ParallelHelpers

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PARALLEL_HELPERS_H_INCLUDED defined