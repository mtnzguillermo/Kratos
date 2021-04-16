//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "containers/data_value_container.h"
#include "includes/variables.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"

namespace Kratos {

    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(TestPointerVectorSet, KratosCoreFastSuite) {
            // create model and model part
            Model model;
            ModelPart& r_model_part = model.CreateModelPart("Main");
            r_model_part.SetBufferSize(3);

            // create 2 elements
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        	r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        	r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        	r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
            r_model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
        	r_model_part.CreateNewNode(6, 1.0, 0.0, 1.0);
        	r_model_part.CreateNewNode(7, 0.0, 1.0, 1.0);
        	r_model_part.CreateNewNode(8, 0.0, 0.0, 2.0);
        	std::vector<ModelPart::IndexType> elNodes1 {1, 2, 3, 4};
        	std::vector<ModelPart::IndexType> elNodes2 {5, 6, 7, 8};
           	auto p_elem_prop = r_model_part.CreateNewProperties(1);
        	r_model_part.CreateNewElement("Element3D4N", 1, elNodes1, p_elem_prop);
        	r_model_part.CreateNewElement("Element3D4N", 2, elNodes2, p_elem_prop);
        	Element::Pointer p_element1 = r_model_part.pGetElement(1);
        	Element::Pointer p_element2 = r_model_part.pGetElement(2);

            // set TO_ERASE = true first element
            p_element1->Set(TO_ERASE, true);
            // remove element 1 from model part
            r_model_part.RemoveElements();

            // list with only element 2
            std::vector<int> ids = {2};

            for(auto& element : r_model_part.Elements()){
                KRATOS_CHECK(element.Id() == 2);
            }
        }
    }
} // namespace Kratos.