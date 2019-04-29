//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Chung To Sang
//


// External includes

// Project includes
#include "add_custom_utilities_to_python.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/renumbering_nodes_utility_for_plasma_dynamics.h"


namespace Kratos
{

namespace Python
{




void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RenumberingNodesUtilityForPlasmaDynamics> (m, "RenumberingNodesUtilityForPlasmaDynamics")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&,ModelPart&,ModelPart&>())
        .def("Renumber", &RenumberingNodesUtilityForPlasmaDynamics::Renumber)
        .def("UndoRenumber", &RenumberingNodesUtilityForPlasmaDynamics::UndoRenumber)
        ;  





}

}  // namespace Python.
} // Namespace Kratos
