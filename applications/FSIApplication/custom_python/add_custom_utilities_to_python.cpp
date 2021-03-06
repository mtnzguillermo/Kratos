//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/FSI_utils.h"
#include "custom_utilities/partitioned_fsi_utilities.hpp"

namespace Kratos
{

namespace Python
{

void AddCustomUtilitiesToPython(pybind11::module &m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, Matrix, Vector > TSpace;

    py::class_<FSIUtils>(m,"FSIUtils")
        .def(py::init<>())
        .def("CheckPressureConvergence",&FSIUtils::CheckPressureConvergence)
        .def("StructuralPressurePrediction",&FSIUtils::StructuralPressurePrediction)
        ;

    py::class_<PartitionedFSIUtilities<TSpace,double,2>, PartitionedFSIUtilities<TSpace,double,2>::Pointer>(m,"PartitionedFSIUtilitiesDouble2D")
        .def(py::init<>())
        .def("GetInterfaceArea", &PartitionedFSIUtilities<TSpace, double, 2>::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &PartitionedFSIUtilities<TSpace, double, 2>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues", &PartitionedFSIUtilities<TSpace, double, 2>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualNorm", &PartitionedFSIUtilities<TSpace, double, 2>::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &PartitionedFSIUtilities<TSpace, double, 2>::ComputeInterfaceResidualVector)
        .def("ComputeAndPrintFluidInterfaceNorms", &PartitionedFSIUtilities<TSpace, double, 2>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &PartitionedFSIUtilities<TSpace, double, 2>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &PartitionedFSIUtilities<TSpace, double, 2>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &PartitionedFSIUtilities<TSpace, double, 2>::CheckCurrentCoordinatesStructure)
        .def("SetUpInterfaceVector", &PartitionedFSIUtilities<TSpace, double, 2>::SetUpInterfaceVector)
        .def("InitializeInterfaceVector", &PartitionedFSIUtilities<TSpace, double, 2>::InitializeInterfaceVector)
        .def("CreateCouplingSkin", &PartitionedFSIUtilities<TSpace, double, 2>::CreateCouplingSkin)
        .def("EmbeddedPressureToPositiveFacePressureInterpolator", &PartitionedFSIUtilities<TSpace, double, 2>::EmbeddedPressureToPositiveFacePressureInterpolator)
        .def("CalculateTractionFromPressureValues", [](PartitionedFSIUtilities<TSpace,double,2>& rPartitionedFSIUtilities, ModelPart& rModelPart, const Variable<double>& rPressureVariable, const Variable<array_1d<double,3>>& rTractionVariable, const bool SwapTractionSign){rPartitionedFSIUtilities.CalculateTractionFromPressureValues(rModelPart, rPressureVariable, rTractionVariable, SwapTractionSign);})
        .def("CalculateTractionFromPressureValues", [](PartitionedFSIUtilities<TSpace,double,2>& rPartitionedFSIUtilities, ModelPart& rModelPart, const Variable<double>& rPositivePressureVariable, const Variable<double>& rNegativePressureVariable, const Variable<array_1d<double,3>>& rTractionVariable, const bool SwapTractionSign){rPartitionedFSIUtilities.CalculateTractionFromPressureValues(rModelPart, rPositivePressureVariable, rNegativePressureVariable, rTractionVariable, SwapTractionSign);})
        ;

    py::class_<PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>, PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::Pointer>(m, "PartitionedFSIUtilitiesArray2D")
        .def(py::init<>())
        .def("GetInterfaceArea", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualNorm", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::ComputeInterfaceResidualVector)
        .def("ComputeAndPrintFluidInterfaceNorms", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::CheckCurrentCoordinatesStructure)
        .def("SetUpInterfaceVector", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::SetUpInterfaceVector)
        .def("InitializeInterfaceVector", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::InitializeInterfaceVector)
        .def("CreateCouplingSkin", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::CreateCouplingSkin)
        .def("EmbeddedPressureToPositiveFacePressureInterpolator", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 2>::EmbeddedPressureToPositiveFacePressureInterpolator)
        .def("CalculateTractionFromPressureValues", [](PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>& rPartitionedFSIUtilities, ModelPart& rModelPart, const Variable<double>& rPressureVariable, const Variable<array_1d<double,3>>& rTractionVariable, const bool SwapTractionSign){rPartitionedFSIUtilities.CalculateTractionFromPressureValues(rModelPart, rPressureVariable, rTractionVariable, SwapTractionSign);})
        .def("CalculateTractionFromPressureValues", [](PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>& rPartitionedFSIUtilities, ModelPart& rModelPart, const Variable<double>& rPositivePressureVariable, const Variable<double>& rNegativePressureVariable, const Variable<array_1d<double,3>>& rTractionVariable, const bool SwapTractionSign){rPartitionedFSIUtilities.CalculateTractionFromPressureValues(rModelPart, rPositivePressureVariable, rNegativePressureVariable, rTractionVariable, SwapTractionSign);})
        ;

    py::class_<PartitionedFSIUtilities<TSpace,double,3>, PartitionedFSIUtilities<TSpace,double,3>::Pointer>(m,"PartitionedFSIUtilitiesDouble3D")
        .def(py::init<>())
        .def("GetInterfaceArea", &PartitionedFSIUtilities<TSpace,double,3>::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &PartitionedFSIUtilities<TSpace,double,3>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues", &PartitionedFSIUtilities<TSpace,double,3>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualNorm", &PartitionedFSIUtilities<TSpace,double,3>::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &PartitionedFSIUtilities<TSpace,double,3>::ComputeInterfaceResidualVector)
        .def("ComputeAndPrintFluidInterfaceNorms", &PartitionedFSIUtilities<TSpace,double,3>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &PartitionedFSIUtilities<TSpace,double,3>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &PartitionedFSIUtilities<TSpace,double,3>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &PartitionedFSIUtilities<TSpace,double,3>::CheckCurrentCoordinatesStructure)
        .def("SetUpInterfaceVector", &PartitionedFSIUtilities<TSpace, double, 3>::SetUpInterfaceVector)
        .def("InitializeInterfaceVector", &PartitionedFSIUtilities<TSpace,double,3>::InitializeInterfaceVector)
        .def("CreateCouplingSkin", &PartitionedFSIUtilities<TSpace,double,3>::CreateCouplingSkin)
        .def("EmbeddedPressureToPositiveFacePressureInterpolator", &PartitionedFSIUtilities<TSpace,double,3>::EmbeddedPressureToPositiveFacePressureInterpolator)
        .def("CalculateTractionFromPressureValues", [](PartitionedFSIUtilities<TSpace,double,3>& rPartitionedFSIUtilities, ModelPart& rModelPart, const Variable<double>& rPressureVariable, const Variable<array_1d<double,3>>& rTractionVariable, const bool SwapTractionSign){rPartitionedFSIUtilities.CalculateTractionFromPressureValues(rModelPart, rPressureVariable, rTractionVariable, SwapTractionSign);})
        .def("CalculateTractionFromPressureValues", [](PartitionedFSIUtilities<TSpace,double,3>& rPartitionedFSIUtilities, ModelPart& rModelPart, const Variable<double>& rPositivePressureVariable, const Variable<double>& rNegativePressureVariable, const Variable<array_1d<double,3>>& rTractionVariable, const bool SwapTractionSign){rPartitionedFSIUtilities.CalculateTractionFromPressureValues(rModelPart, rPositivePressureVariable, rNegativePressureVariable, rTractionVariable, SwapTractionSign);})
        ;

    py::class_<PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>, PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::Pointer>(m,"PartitionedFSIUtilitiesArray3D")
        .def(py::init<>())
        .def("GetInterfaceArea", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualNorm", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::ComputeInterfaceResidualVector)
        .def("ComputeAndPrintFluidInterfaceNorms", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::CheckCurrentCoordinatesStructure)
        .def("SetUpInterfaceVector", &PartitionedFSIUtilities<TSpace, array_1d<double, 3>, 3>::SetUpInterfaceVector)
        .def("InitializeInterfaceVector", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::InitializeInterfaceVector)
        .def("CreateCouplingSkin", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::CreateCouplingSkin)
        .def("EmbeddedPressureToPositiveFacePressureInterpolator", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::EmbeddedPressureToPositiveFacePressureInterpolator)
        .def("CalculateTractionFromPressureValues", [](PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>& rPartitionedFSIUtilities, ModelPart& rModelPart, const Variable<double>& rPressureVariable, const Variable<array_1d<double,3>>& rTractionVariable, const bool SwapTractionSign){rPartitionedFSIUtilities.CalculateTractionFromPressureValues(rModelPart, rPressureVariable, rTractionVariable, SwapTractionSign);})
        .def("CalculateTractionFromPressureValues", [](PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>& rPartitionedFSIUtilities, ModelPart& rModelPart, const Variable<double>& rPositivePressureVariable, const Variable<double>& rNegativePressureVariable, const Variable<array_1d<double,3>>& rTractionVariable, const bool SwapTractionSign){rPartitionedFSIUtilities.CalculateTractionFromPressureValues(rModelPart, rPositivePressureVariable, rNegativePressureVariable, rTractionVariable, SwapTractionSign);})
        ;

}

}  // namespace Python.

} // Namespace Kratos
