// Project includes
#include "split_forward_euler_scheme.h"

namespace Kratos {

    void SplitForwardEulerScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning SplitForwardEulerScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void SplitForwardEulerScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning SplitForwardEulerScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void SplitForwardEulerScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node < 3 >& i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            const array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3]) {

        if(StepFlag == 1) //PREDICT
        {
            const array_1d<double, 3 >& aux_impulse = i.GetValue(AUX_IMPULSE);
            const double theta_2 = i.GetValue(THETA_2);
            const double nodal_damping = i.GetValue(NODAL_DAMPING);

            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    delta_displ[k] = (delta_t*aux_impulse[k] + (mass - delta_t*(1.0-theta_2)*nodal_damping)*displ[k]) /
                                     (mass + delta_t*theta_2*nodal_damping) - displ[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                    vel[k] = delta_displ[k] / delta_t;
                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
        else if(StepFlag == 2) //CORRECT
        {
            array_1d<double, 3 >& aux_impulse = i.GetValue(AUX_IMPULSE);
            const array_1d<double, 3 >& external_forces = i.GetValue(EXTERNAL_FORCES);
            const array_1d<double, 3 >& internal_forces = i.GetValue(INTERNAL_FORCES);
            const array_1d<double, 3 >& internal_forces_old = i.GetValue(INTERNAL_FORCES_OLD);
            const double beta = i.GetValue(BETA_RAYLEIGH);
            const double theta_1 = i.GetValue(THETA_1);
            const double nodal_stiffness = i.GetValue(NODAL_STIFFNESS);

            for (int k = 0; k < 3; k++) {
                aux_impulse[k] += delta_t*external_forces[k] - (beta+delta_t*theta_1)*internal_forces[k] +
                                  (beta-delta_t*(1.0-theta_1))*internal_forces_old[k] +
                                  delta_t*beta*nodal_stiffness*vel[k];
            }
        }

    }

    void SplitForwardEulerScheme::CalculateNewRotationalVariablesOfSpheres(
                int StepFlag,
                Node < 3 >& i,
                const double moment_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        if(StepFlag == 1) //PREDICT
        {
            const array_1d<double, 3 >& rotational_aux_impulse = i.GetValue(ROTATIONAL_AUX_IMPULSE);
            const double theta_2 = i.GetValue(THETA_2);
            const double rotational_nodal_damping = i.GetValue(ROTATIONAL_NODAL_DAMPING);

            for (int k = 0; k < 3; k++) {
                if (Fix_Ang_vel[k] == false) {
                    delta_rotation[k] = (delta_t*rotational_aux_impulse[k] + (moment_of_inertia - delta_t*(1.0-theta_2)*rotational_nodal_damping)*rotated_angle[k]) /
                                        (moment_of_inertia + delta_t*theta_2*rotational_nodal_damping) - rotated_angle[k];
                    rotated_angle[k] += delta_rotation[k];
                    angular_velocity[k] = delta_rotation[k] / delta_t;
                } else {
                    delta_rotation[k] = delta_t * angular_velocity[k];
                    rotated_angle[k] += delta_rotation[k];
                }
            }
        }
        else if(StepFlag == 2) //CORRECT
        {
            array_1d<double, 3 >& rotational_aux_impulse = i.GetValue(ROTATIONAL_AUX_IMPULSE);
            const array_1d<double, 3 >& rotational_external_forces = i.GetValue(ROTATIONAL_EXTERNAL_FORCES);
            const array_1d<double, 3 >& rotational_internal_forces = i.GetValue(ROTATIONAL_INTERNAL_FORCES);
            const array_1d<double, 3 >& rotational_internal_forces_old = i.GetValue(ROTATIONAL_INTERNAL_FORCES_OLD);
            const double beta = i.GetValue(BETA_RAYLEIGH);
            const double theta_1 = i.GetValue(THETA_1);
            const double rotational_nodal_stiffness = i.GetValue(ROTATIONAL_NODAL_STIFFNESS);

            for (int k = 0; k < 3; k++) {
                rotational_aux_impulse[k] += delta_t*rotational_external_forces[k] - (beta+delta_t*theta_1)*rotational_internal_forces[k] +
                                             (beta-delta_t*(1.0-theta_1))*rotational_internal_forces_old[k] +
                                             delta_t*beta*rotational_nodal_stiffness*angular_velocity[k];
            }
        }
    }

    void SplitForwardEulerScheme::CalculateNewRotationalVariablesOfRigidBodyElements(
                int StepFlag,
                Node < 3 >& i,
                const array_1d<double, 3 > moments_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        array_1d<double, 3 >& local_angular_velocity = i.FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY);

        array_1d<double, 3 > local_angular_acceleration, local_torque, angular_acceleration;

        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, torque, local_torque);
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
        CalculateLocalAngularAccelerationByEulerEquations(local_angular_velocity, moments_of_inertia, local_torque, moment_reduction_factor, local_angular_acceleration);
        GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_acceleration, angular_acceleration);

        // TODO: this is from symplectic euler...
        UpdateRotationalVariables(StepFlag, i, rotated_angle, delta_rotation, angular_velocity, angular_acceleration, delta_t, Fix_Ang_vel);

        double ang = DEM_INNER_PRODUCT_3(delta_rotation, delta_rotation);

        if (ang) {
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        } //if ang
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
    }

    void SplitForwardEulerScheme::UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        for (int k = 0; k < 3; k++) {
            if (Fix_Ang_vel[k] == false) {
                angular_velocity[k] += delta_t * angular_acceleration[k];
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            } else {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            }
        }
    }

    void SplitForwardEulerScheme::CalculateLocalAngularAccelerationByEulerEquations(
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration) {

        for (int j = 0; j < 3; j++) {
            local_angular_acceleration[j] = (local_torque[j] - (local_angular_velocity[(j + 1) % 3] * moments_of_inertia[(j + 2) % 3] * local_angular_velocity[(j + 2) % 3] - local_angular_velocity[(j + 2) % 3] * moments_of_inertia[(j + 1) % 3] * local_angular_velocity[(j + 1) % 3])) / moments_of_inertia[j];
            local_angular_acceleration[j] = local_angular_acceleration[j] * moment_reduction_factor;
        }
    }
} //namespace Kratos
