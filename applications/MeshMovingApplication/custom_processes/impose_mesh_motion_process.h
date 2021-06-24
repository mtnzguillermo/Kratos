//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//	                Kratos default license: kratos/license.txt
//
//  Main Authors:   Máté Kelemen
//

#ifndef KRATOS_IMPOSE_MESH_MOTION_PROCESS_INCLUDED
#define KRATOS_IMPOSE_MESH_MOTION_PROCESS_INCLUDED

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "utilities/quaternion.h"
#include "utilities/interval_utility.h"

namespace Kratos
{

///@name Kratos Classes
///@{


/** Impose a rotation followed by translation on a ModelPart
 *  @details The transformation is equivalent to:
 *  1) Translation to the reference frame (offset the origin)
 *  2) Specified rotation
 *  3) Reverse translation from the reference frame (undo origin offset)
 *  4) Specified translation
 *  @note angles in radians
 */
class KRATOS_API(MESH_MOVING_APPLICATION) ImposeMeshMotionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SetHMapProcess
    KRATOS_CLASS_POINTER_DEFINITION(ImposeMeshMotionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor
     *  @details The rotation can be defined by either "euler_angles"
     *  or a "rotation_axis" and "rotation_angle" pair. The following parameters can be
     *  defined parametrically (see @ref{GenericFunctionUtility}):
     *  "euler_angles", "rotation_axis", "reference_point", "rotation_angle", "translation_vector"
     * 
     *  Default parameters:
     *  {
     *      "model_part_name"       : "",
     *      "interval"              : [0.0, "End"],
     *      "rotation_definition"   : "rotation_axis",
     *      "euler_angles"          : [0.0, 0.0, 0.0],
     *      "rotation_axis"         : [0.0, 0.0, 1.0],
     *      "reference_point"       : [0.0, 0.0, 0.0]
     *      "rotation_angle"        : 0,
     *      "translation_vector"    : [0.0, 0.0, 0.0],
     *  }
     *  @note The euler angles follow the convention specified by @ref{Quaternion} (Z, -X', Z")
     */
    ImposeMeshMotionProcess(Model& rModel, Parameters parameters);

    /** Constructor
     *  @details The rotation can be defined by either "euler_angles"
     *  or a "rotation_axis" and "rotation_angle" pair. The following parameters can be
     *  defined parametrically (see @ref{GenericFunctionUtility}):
     *  "euler_angles", "rotation_axis", "reference_point", "rotation_angle", "translation_vector"
     * 
     *  Default parameters:
     *  {
     *      "model_part_name"       : "",
     *      "interval"              : [0.0, "End"],
     *      "rotation_definition"   : "rotation_axis",
     *      "euler_angles"          : [0.0, 0.0, 0.0],
     *      "rotation_axis"         : [0.0, 0.0, 1.0],
     *      "reference_point"       : [0.0, 0.0, 0.0]
     *      "rotation_angle"        : 0,
     *      "translation_vector"    : [0.0, 0.0, 0.0],
     *  }
     *  @note The euler angles follow the convention specified by @ref{Quaternion} (Z, -X', Z")
     */
    ImposeMeshMotionProcess(ModelPart& rModelPart, Parameters parameters);

    ///@}
    ///@name Operations
    ///@{

    virtual void ExecuteInitializeSolutionStep() override;

    virtual const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Class name as string
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rStream) const override;

    ///@}

private:
    ///@name Private operations
    ///@{

    void ParseAndSetConstantTransform(const std::string& rRotationDefinition,
                                      const Parameters& rEulerAngles,
                                      const Parameters& rRotationAxis,
                                      const Parameters& rRotationAngle,
                                      const Parameters& rReferencePoint,
                                      const Parameters& rTranslationVector);

    void ParseAndSetParametricTransform(const std::string& rRotationDefinition,
                                        const Parameters& rEulerAngles,
                                        const Parameters& rRotationAxis,
                                        const Parameters& rRotationAngle,
                                        const Parameters& rReferencePoint,
                                        const Parameters& rTranslationVector);

    ///@}
    ///@name Member Variables
    ///@{
        
    ModelPart& mrModelPart;

    IntervalUtility mIntervalUtility;

    std::function<array_1d<double,3>(const Node<3>&)> mTransformFunctor;

    ///@}
};

///@}

///@}
///@name Input and output
///@{

/// Dump info to output stream
std::ostream& operator<<(std::ostream& rStream,
                         const ImposeMeshMotionProcess& rThis);

///@}

} // namespace Kratos

#endif // KRATOS_IMPOSE_MESH_MOTION_PROCESS_INCLUDED