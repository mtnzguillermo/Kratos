//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    KratosAppGenerator
//
//

// System includes
// See the header file

// External includes
// See the header file

// Project includes
// See the header file

// Application includes
#include "surface_smoothing_process.h"
// See the header file

namespace Kratos
{

/* Public functions *******************************************************/
SurfaceSmoothingProcess::SurfaceSmoothingProcess(
    ModelPart& rModelPart,
    TLinearSolver::Pointer plinear_solver)
    : Process(),
      mrModelPart(rModelPart)    
{
    // Member variables initialization
    // Nothing!

    // Generate an auxilary model part and populate it by elements of type MySimpleElement
    CreateAuxModelPart();

    auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);
    //auto p_builder_solver = Kratos::make_unique<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);

    InitializeSolutionStrategy(plinear_solver, p_builder_solver);
}

void SurfaceSmoothingProcess::CreateAuxModelPart()
{
    ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& r_elems = mrModelPart.Elements();
    
    Model& current_model = mrModelPart.GetModel();
    if(current_model.HasModelPart( mAuxModelPartName ))
        current_model.DeleteModelPart( mAuxModelPartName );

    // Adding DISTANCE to the solution variables is not needed if it is already a solution variable of the problem
    mrModelPart.AddNodalSolutionStepVariable(DISTANCE_AUX);

    // Ensure that the nodes have distance as a DOF
    VariableUtils().AddDof<Variable<double> >(DISTANCE_AUX, mrModelPart);

    // Ensure that the nodes have DISTANCE_AUX as a DOF !!!! NOT NEEDED HERE IF IT IS DONE IN THE PYTHON SCRIPT !!!!
    /* #pragma omp parallel for
    for (int k = 0; k < static_cast<int>(mrModelPart.NumberOfNodes()); ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        it_node->AddDof(DISTANCE_AUX);
    } */

    // Generate AuxModelPart
    ModelPart& r_smoothing_model_part = current_model.CreateModelPart( mAuxModelPartName );

    Element::Pointer p_smoothing_element = Kratos::make_intrusive<SurfaceSmoothingElement>();

    ConnectivityPreserveModeler modeler;
    modeler.GenerateModelPart(mrModelPart, r_smoothing_model_part, *p_smoothing_element);

    const double delta_time = mrModelPart.pGetProcessInfo()->GetValue(DELTA_TIME);
    r_smoothing_model_part.pGetProcessInfo()->SetValue(DELTA_TIME, delta_time);
        
}

void SurfaceSmoothingProcess::Execute()
{
    KRATOS_TRY;
    
    const unsigned int NumNodes = mrModelPart.NumberOfNodes();
    const unsigned int NumElements = mrModelPart.NumberOfElements();

    std::vector<double> DistDiff(NumNodes);
    std::vector<double> DistDiffAvg(NumNodes);
    std::vector<double> NumNeighbors(NumNodes);

    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        it_node->Free(DISTANCE_AUX);
        const double distance = it_node->FastGetSolutionStepValue(DISTANCE);
        it_node->FastGetSolutionStepValue(DISTANCE_AUX) = distance;

        //if ( it_node->GetValue(IS_STRUCTURE) == 1.0 ){
        //    it_node->Fix(DISTANCE_AUX);
        //}
    }

    //Model& current_model = mrModelPart.GetModel();
    //ModelPart& r_smoothing_model_part = current_model.GetModelPart( mAuxModelPartName );
    //r_smoothing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,1);

    KRATOS_INFO("SurfaceSmoothingProcess") << "About to solve the LSE" << std::endl;
    mp_solving_strategy->Solve();
    KRATOS_INFO("SurfaceSmoothingProcess") << "LSE is solved" << std::endl;

    //#pragma omp parallel for
    //for (unsigned int k = 0; k < NumNodes; ++k) {
    //    auto it_node = mrModelPart.NodesBegin() + k;
    //    const double distance = it_node->FastGetSolutionStepValue(DISTANCE_AUX);
    //    it_node->FastGetSolutionStepValue(DISTANCE) = distance;
    //}

    //r_smoothing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,2);

    //KRATOS_INFO("SurfaceSmoothingProcess") << "About to solve the LSE" << std::endl;
    //mp_solving_strategy->Solve();
    //KRATOS_INFO("SurfaceSmoothingProcess") << "LSE is solved" << std::endl;

    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        DistDiff[it_node->Id()-1] = it_node->FastGetSolutionStepValue(DISTANCE_AUX) - it_node->FastGetSolutionStepValue(DISTANCE);
    } 

    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    for (unsigned int k = 0; k < NumElements; ++k) {
        auto it_elem = mrModelPart.ElementsBegin() + k;
        auto geom = it_elem->pGetGeometry();

        for (unsigned int i=0; i < num_nodes; i++){
            for (unsigned int j=i+1; j < num_nodes; j++){
                const int iId = (*geom)[i].Id() - 1;
                const int jId = (*geom)[j].Id() - 1;
                NumNeighbors[iId] += 1.0;
                NumNeighbors[jId] += 1.0;

                DistDiffAvg[iId] += DistDiff[jId];
                DistDiffAvg[jId] += DistDiff[iId];
            }
        }
    }

    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        //KRATOS_INFO("SurfaceSmoothingProcess, Nodes Id") << it_node->Id() << std::endl;
        if (NumNeighbors[it_node->Id()-1] != 0.0 /* && it_node->GetValue(IS_STRUCTURE) == 0.0 */){
             it_node->FastGetSolutionStepValue(DISTANCE_AUX) = 
                it_node->FastGetSolutionStepValue(DISTANCE_AUX) - 1.0/NumNeighbors[it_node->Id()-1]*DistDiffAvg[it_node->Id()-1];
        }
        
        it_node->Free(DISTANCE_AUX);
    }       

    KRATOS_CATCH("");
}

void SurfaceSmoothingProcess::ExecuteInitialize()
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void SurfaceSmoothingProcess::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void SurfaceSmoothingProcess::ExecuteInitializeSolutionStep() {
//Nothing
}

void SurfaceSmoothingProcess::ExecuteFinalizeSolutionStep() {
//Nothing
}

/* Protected functions ****************************************************/

//void SurfaceSmoothingProcess::CheckAndStoreVariablesList(const std::vector<std::string>& rVariableStringArray){
//Nothing
//}

}; // namespace Kratos