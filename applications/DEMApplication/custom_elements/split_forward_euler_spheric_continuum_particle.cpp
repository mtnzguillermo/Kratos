// System includes
#include <string>
#include <iostream>
#include <cmath>

#include <fstream>

// External includes

// Project includes
#include "split_forward_euler_spheric_continuum_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/discrete_particle_configure.h"
#include "custom_strategies/schemes/glued_to_wall_scheme.h"


namespace Kratos
{
// using namespace GeometryFunctions;

SplitForwardEulerSphericContinuumParticle::SplitForwardEulerSphericContinuumParticle()
    : SphericContinuumParticle() {
}

SplitForwardEulerSphericContinuumParticle::SplitForwardEulerSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : SphericContinuumParticle(NewId, pGeometry) {
}

SplitForwardEulerSphericContinuumParticle::SplitForwardEulerSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : SphericContinuumParticle(NewId, pGeometry, pProperties) {
}

SplitForwardEulerSphericContinuumParticle::SplitForwardEulerSphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericContinuumParticle(NewId, ThisNodes) {
}

Element::Pointer SplitForwardEulerSphericContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SplitForwardEulerSphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

/// Destructor.
SplitForwardEulerSphericContinuumParticle::~SplitForwardEulerSphericContinuumParticle(){
}

// SplitForwardEulerSphericContinuumParticle& SplitForwardEulerSphericContinuumParticle::operator=(const SplitForwardEulerSphericContinuumParticle& rOther) {
//     SphericContinuumParticle::operator=(rOther);

//     return *this;
// }

void SplitForwardEulerSphericContinuumParticle::Initialize(const ProcessInfo& r_process_info)
{
    KRATOS_TRY

    SphericContinuumParticle::Initialize(r_process_info);

    NodeType& node = GetGeometry()[0];

    node.SetValue(NODAL_DAMPING,0.0);
    node.SetValue(ROTATIONAL_NODAL_DAMPING,0.0);
    node.SetValue(NODAL_STIFFNESS,0.0);
    node.SetValue(ROTATIONAL_NODAL_STIFFNESS,0.0);
    node.SetValue(ALPHA_RAYLEIGH,0.0);
    node.SetValue(BETA_RAYLEIGH,0.0);
    node.SetValue(THETA_1,0.0);
    node.SetValue(THETA_2,0.0);
    array_1d<double,3> zero_vector = ZeroVector(3);
    node.SetValue(EXTERNAL_FORCES,zero_vector);
    node.SetValue(ROTATIONAL_EXTERNAL_FORCES,zero_vector);
    node.SetValue(INTERNAL_FORCES,zero_vector);
    node.SetValue(ROTATIONAL_INTERNAL_FORCES,zero_vector);
    node.SetValue(INTERNAL_FORCES_OLD,zero_vector);
    node.SetValue(ROTATIONAL_INTERNAL_FORCES_OLD,zero_vector);
    node.SetValue(AUX_IMPULSE,zero_vector);
    node.SetValue(ROTATIONAL_AUX_IMPULSE,zero_vector);

    KRATOS_CATCH( "" )
}

void SplitForwardEulerSphericContinuumParticle::CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity)
{
    KRATOS_TRY

    // Creating a data buffer to store those variables that we want to reuse so that we can keep function parameter lists short

    SphericParticle::BufferPointerType p_buffer = CreateParticleDataBuffer(this); // all memory will be freed once this shared pointer goes out of scope
    SphericParticle::ParticleDataBuffer& data_buffer = *p_buffer;
    data_buffer.SetBoundingBox(r_process_info[DOMAIN_IS_PERIODIC], r_process_info[DOMAIN_MIN_CORNER], r_process_info[DOMAIN_MAX_CORNER]);

    NodeType& this_node = GetGeometry()[0];

    data_buffer.mDt = dt;
    data_buffer.mMultiStageRHS = false;

    array_1d<double, 3> additional_forces = ZeroVector(3);
    array_1d<double, 3> additionally_applied_moment = ZeroVector(3);
    array_1d<double, 3>& elastic_force       = this_node.FastGetSolutionStepValue(ELASTIC_FORCES);
    array_1d<double, 3>& contact_force       = this_node.FastGetSolutionStepValue(CONTACT_FORCES);
    array_1d<double, 3>& rigid_element_force = this_node.FastGetSolutionStepValue(RIGID_ELEMENT_FORCE);

    const double& mass = this_node.FastGetSolutionStepValue(NODAL_MASS);
    const double& moment_of_inertia = this_node.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
    double& nodal_stiffness = this_node.GetValue(NODAL_STIFFNESS);
    double& nodal_rotational_stiffness = this_node.GetValue(ROTATIONAL_NODAL_STIFFNESS);
    double& nodal_damping = this_node.GetValue(NODAL_DAMPING);
    double& nodal_rotational_damping = this_node.GetValue(ROTATIONAL_NODAL_DAMPING);

    double& alpha_rayleigh = this_node.GetValue(ALPHA_RAYLEIGH);
    double& beta_rayleigh = this_node.GetValue(BETA_RAYLEIGH);
    double& theta_1 = this_node.GetValue(THETA_1);
    double& theta_2 = this_node.GetValue(THETA_2);
    // NOTE: we define alpha and beta constant for all particles
    alpha_rayleigh = r_process_info[ALPHA_RAYLEIGH];
    beta_rayleigh = r_process_info[BETA_RAYLEIGH];
    theta_1 = r_process_info[THETA_1];
    theta_2 = r_process_info[THETA_2];

    array_1d<double, 3>& external_forces = this_node.GetValue(EXTERNAL_FORCES);
    array_1d<double, 3>& rotational_external_forces = this_node.GetValue(ROTATIONAL_EXTERNAL_FORCES);
    array_1d<double, 3>& internal_forces = this_node.GetValue(INTERNAL_FORCES);
    array_1d<double, 3>& rotational_internal_forces = this_node.GetValue(ROTATIONAL_INTERNAL_FORCES);
    array_1d<double, 3>& internal_forces_old = this_node.GetValue(INTERNAL_FORCES_OLD);
    array_1d<double, 3>& rotational_internal_forces_old = this_node.GetValue(ROTATIONAL_INTERNAL_FORCES_OLD);

    noalias(internal_forces_old) = contact_force;
    noalias(rotational_internal_forces_old) = mContactMoment;

    mContactMoment.clear();
    elastic_force.clear();
    contact_force.clear();
    rigid_element_force.clear();
    nodal_stiffness = 0.0;
    nodal_rotational_stiffness = 0.0;

    InitializeForceComputation(r_process_info);

    double RollingResistance = 0.0;

    ComputeBallToBallStiffness(data_buffer, nodal_stiffness, nodal_rotational_stiffness);

    ComputeBallToRigidFaceStiffness(data_buffer, nodal_stiffness, nodal_rotational_stiffness);

    // Rayleigh Damping
    // TODO: does it make sense to use the same alpha and beta both for translational and rotational DOFs?
    nodal_damping = alpha_rayleigh * mass + beta_rayleigh * nodal_stiffness;
    nodal_rotational_damping = alpha_rayleigh * moment_of_inertia + beta_rayleigh * nodal_rotational_stiffness;

    // KRATOS_WATCH(this->Id())
    // KRATOS_WATCH(mass)
    // KRATOS_WATCH(nodal_stiffness)
    // KRATOS_WATCH(nodal_damping)
    // KRATOS_WATCH(moment_of_inertia)
    // KRATOS_WATCH(nodal_rotational_stiffness)
    // KRATOS_WATCH(nodal_rotational_damping)

    ComputeBallToBallContactForce(data_buffer, r_process_info, elastic_force, contact_force, RollingResistance);

    ComputeBallToRigidFaceContactForce(data_buffer, elastic_force, contact_force, RollingResistance, rigid_element_force, r_process_info);

    if (this->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)){
        ComputeAdditionalForces(additional_forces, additionally_applied_moment, r_process_info, gravity);
        #ifdef KRATOS_DEBUG
        DemDebugFunctions::CheckIfNan(additional_forces, "NAN in Additional Force in RHS of Ball");
        DemDebugFunctions::CheckIfNan(additionally_applied_moment, "NAN in Additional Torque in RHS of Ball");
        #endif
    }

    // ROLLING FRICTION
    if (this->Is(DEMFlags::HAS_ROTATION) && !data_buffer.mMultiStageRHS) {
        if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !data_buffer.mMultiStageRHS) {
            array_1d<double, 3>& rolling_resistance_moment = this_node.FastGetSolutionStepValue(ROLLING_RESISTANCE_MOMENT);
            rolling_resistance_moment.clear();

            ComputeRollingFriction(rolling_resistance_moment, RollingResistance, data_buffer.mDt);
        }
    }

    noalias(external_forces) = additional_forces;
    noalias(rotational_external_forces) = additionally_applied_moment;
    // NOTE: here the contact forces have no damping forces because we set them to zero
    noalias(internal_forces) = contact_force;
    noalias(rotational_internal_forces) = mContactMoment;

    array_1d<double,3>& total_forces = this_node.FastGetSolutionStepValue(TOTAL_FORCES);
    array_1d<double,3>& total_moment = this_node.FastGetSolutionStepValue(PARTICLE_MOMENT);
    array_1d<double,3>& velocity = this_node.FastGetSolutionStepValue(VELOCITY);
    array_1d<double,3>& angular_velocity = this_node.FastGetSolutionStepValue(ANGULAR_VELOCITY);

    // TODO: check sign of diagonal damping force
    // TODO: can be diagonal damping force added globally ? I think it must be added globally. However,
    //       if the local damping force is used to calculate other local forces, there will be a small
    //       error because the diagonal damping force won't have been added until now.
    // total_forces[0] = contact_force[0] + additional_forces[0] + r_process_info[BETA_RAYLEIGH] * nodal_stiffness * velocity[0];
    // total_forces[1] = contact_force[1] + additional_forces[1] + r_process_info[BETA_RAYLEIGH] * nodal_stiffness * velocity[1];
    // total_forces[2] = contact_force[2] + additional_forces[2] + r_process_info[BETA_RAYLEIGH] * nodal_stiffness * velocity[2];

    // total_moment[0] = mContactMoment[0] + additionally_applied_moment[0] + r_process_info[BETA_RAYLEIGH] * nodal_rotational_stiffness * angular_velocity[0];
    // total_moment[1] = mContactMoment[1] + additionally_applied_moment[1] + r_process_info[BETA_RAYLEIGH] * nodal_rotational_stiffness * angular_velocity[1];
    // total_moment[2] = mContactMoment[2] + additionally_applied_moment[2] + r_process_info[BETA_RAYLEIGH] * nodal_rotational_stiffness * angular_velocity[2];

    total_forces[0] = contact_force[0] + additional_forces[0];
    total_forces[1] = contact_force[1] + additional_forces[1];
    total_forces[2] = contact_force[2] + additional_forces[2];

    total_moment[0] = mContactMoment[0] + additionally_applied_moment[0];
    total_moment[1] = mContactMoment[1] + additionally_applied_moment[1];
    total_moment[2] = mContactMoment[2] + additionally_applied_moment[2];

    ApplyGlobalDampingToContactForcesAndMoments(total_forces, total_moment);

    #ifdef KRATOS_DEBUG
    DemDebugFunctions::CheckIfNan(total_forces, "NAN in Total Forces in RHS of Ball");
    DemDebugFunctions::CheckIfNan(total_moment, "NAN in Total Torque in RHS of Ball");
    #endif

    FinalizeForceComputation(data_buffer);
    KRATOS_CATCH("")
}

void SplitForwardEulerSphericContinuumParticle::ComputeBallToBallContactForce(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                const ProcessInfo& r_process_info,
                                                                array_1d<double, 3>& rElasticForce,
                                                                array_1d<double, 3>& rContactForce,
                                                                double& RollingResistance)
{
    KRATOS_TRY

    NodeType& this_node = this->GetGeometry()[0];
    DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mMyCoors, this_node)

    const int time_steps = r_process_info[TIME_STEPS];

    const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
    const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    Vector& cont_ini_neigh_area            = this->GetValue(NEIGHBOURS_CONTACT_AREAS);
    int NeighbourSize = mNeighbourElements.size();
    GetGeometry()[0].GetSolutionStepValue(NEIGHBOUR_SIZE) = NeighbourSize;

    for (int i = 0; data_buffer.SetNextNeighbourOrExit(i); ++i) {

        if (mNeighbourElements[i] == NULL) continue;
        if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;
        SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
        data_buffer.mpOtherParticle = neighbour_iterator;

        unsigned int neighbour_iterator_id = data_buffer.mpOtherParticle->Id();

        noalias(data_buffer.mOtherToMeVector) = this->GetGeometry()[0].Coordinates() - data_buffer.mpOtherParticle->GetGeometry()[0].Coordinates();

        const double& other_radius = data_buffer.mpOtherParticle->GetRadius();

        data_buffer.mDistance = DEM_MODULUS_3(data_buffer.mOtherToMeVector);
        double radius_sum = GetRadius() + other_radius;

        double initial_delta = GetInitialDelta(i);

        double initial_dist = radius_sum - initial_delta;
        double indentation = initial_dist - data_buffer.mDistance;
        double myYoung = GetYoung();
        double myPoisson = GetPoisson();

        double kn_el = 0.0;
        double kt_el = 0.0;
        double DeltDisp[3] = {0.0};
        double RelVel[3] = {0.0};
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)
        bool sliding = false;

        double contact_tau = 0.0;
        double contact_sigma = 0.0;
        double failure_criterion_state = 0.0;
        double acumulated_damage = 0.0;

        // Getting neighbor properties
        double other_young = data_buffer.mpOtherParticle->GetYoung();
        double other_poisson = data_buffer.mpOtherParticle->GetPoisson();
        double equiv_poisson;
        if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
        else { equiv_poisson = 0.0; }

        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
        double calculation_area = 0.0;
        const double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));

        if (i < (int)mContinuumInitialNeighborsSize) {
            mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, cont_ini_neigh_area, i, calculation_area); //some Constitutive Laws get a value, some others calculate the value.
            mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area, this, neighbour_iterator);
        }

        EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel, data_buffer.mLocalCoordSystem, data_buffer.mOldLocalCoordSystem, vel, delta_displ);

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            RelativeDisplacementAndVelocityOfContactPointDueToRotationQuaternion(DeltDisp, RelVel, data_buffer.mOldLocalCoordSystem, other_radius, data_buffer.mDt, ang_vel, neighbour_iterator);
        }

        RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(r_process_info, DeltDisp, RelVel, data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, neighbour_iterator);

        double LocalDeltDisp[3] = {0.0};
        double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
        double LocalElasticExtraContactForce[3] = {0.0};
        double GlobalElasticContactForce[3] = {0.0};
        double GlobalElasticExtraContactForce[3] = {0.0};
        double TotalGlobalElasticContactForce[3] = {0.0};
        double OldLocalElasticContactForce[3] = {0.0};

        //FilterNonSignificantDisplacements(DeltDisp, RelVel, indentation);

        GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltDisp, LocalDeltDisp);

        RotateOldContactForces(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, mNeighbourElasticContactForces[i]);
        RotateOldContactForces(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, mNeighbourElasticExtraContactForces[i]);

        GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);

        GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i][0];
        GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i][1];
        GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i][2];

        GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //TODO: can we remove this? We should overwrite LocalElasticContactForce afterwards

        double ViscoDampingLocalContactForce[3] = {0.0};
        double ElasticLocalRotationalMoment[3] = {0.0};
        double ViscoLocalRotationalMoment[3] = {0.0};
        double cohesive_force =  0.0;
        double LocalRelVel[3] = {0.0};
        GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, RelVel, LocalRelVel);

        if (i < (int)mContinuumInitialNeighborsSize) {

            mContinuumConstitutiveLawArray[i]->CheckFailure(i, this, neighbour_iterator);

            mContinuumConstitutiveLawArray[i]->CalculateForcesRayleigh(r_process_info,
                                                                OldLocalElasticContactForce,
                                                                LocalElasticContactForce,
                                                                LocalElasticExtraContactForce,
                                                                data_buffer.mLocalCoordSystem,
                                                                LocalDeltDisp,
                                                                kn_el,
                                                                kt_el,
                                                                contact_sigma,
                                                                contact_tau,
                                                                failure_criterion_state,
                                                                equiv_young,
                                                                equiv_shear,
                                                                indentation,
                                                                calculation_area,
                                                                acumulated_damage,
                                                                this,
                                                                neighbour_iterator,
                                                                i,
                                                                r_process_info[TIME_STEPS],
                                                                sliding,
                                                                LocalRelVel,
                                                                ViscoDampingLocalContactForce);

        } else if (indentation > 0.0) {
            const double previous_indentation = indentation + LocalDeltDisp[2];
            mDiscontinuumConstitutiveLaw->CalculateForcesRayleigh(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce,
                    LocalDeltDisp, LocalRelVel, indentation, previous_indentation,
                    ViscoDampingLocalContactForce, cohesive_force, this, data_buffer.mpOtherParticle, sliding, data_buffer.mLocalCoordSystem);
        } else { //Not bonded and no idata_buffer.mpOtherParticlendentation
            LocalElasticContactForce[0] = 0.0;      LocalElasticContactForce[1] = 0.0;      LocalElasticContactForce[2] = 0.0;
            ViscoDampingLocalContactForce[0] = 0.0; ViscoDampingLocalContactForce[1] = 0.0; ViscoDampingLocalContactForce[2] = 0.0;
            cohesive_force= 0.0;
        }

        // Transforming to global forces and adding up
        double LocalContactForce[3] = {0.0};
        double GlobalContactForce[3] = {0.0};

        if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < (int)mContinuumInitialNeighborsSize)) { // We leave apart the discontinuum neighbors (the same for the walls). The neighbor would not be able to do the same if we activate it.
            mContinuumConstitutiveLawArray[i]->AddPoissonContribution(equiv_poisson, data_buffer.mLocalCoordSystem, LocalElasticContactForce[2], calculation_area, mSymmStressTensor, this, neighbour_iterator, r_process_info, i, indentation);
        }

        array_1d<double, 3> other_ball_to_ball_forces(3,0.0);
        ComputeOtherBallToBallForces(other_ball_to_ball_forces);

        AddUpForcesAndProject(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, LocalContactForce, LocalElasticContactForce, LocalElasticExtraContactForce, GlobalContactForce,
                                GlobalElasticContactForce, GlobalElasticExtraContactForce, TotalGlobalElasticContactForce,ViscoDampingLocalContactForce, 0.0, other_ball_to_ball_forces, rElasticForce, rContactForce, i, r_process_info); //TODO: replace the 0.0 with an actual cohesive force for discontinuum neighbours

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            ComputeMoments(LocalContactForce[2], TotalGlobalElasticContactForce, RollingResistance, data_buffer.mLocalCoordSystem[2], data_buffer.mpOtherParticle, indentation, false, i);
            if (i < (int)mContinuumInitialNeighborsSize && mIniNeighbourFailureId[i] == 0) {
                mContinuumConstitutiveLawArray[i]->ComputeParticleRotationalMoments(this, neighbour_iterator, equiv_young, data_buffer.mDistance, calculation_area,
                                                                                    data_buffer.mLocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment, equiv_poisson, indentation);
            }
            // TODO: why ViscoLocalRotationalMoment is 0.0 ?
            AddUpMomentsAndProject(data_buffer.mLocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
        }

        if (r_process_info[CONTACT_MESH_OPTION] == 1 && (i < (int)mContinuumInitialNeighborsSize) && this->Id() < neighbour_iterator_id) {
            double total_local_elastic_contact_force[3] = {0.0};
            total_local_elastic_contact_force[0] = LocalElasticContactForce[0] + LocalElasticExtraContactForce[0];
            total_local_elastic_contact_force[1] = LocalElasticContactForce[1] + LocalElasticExtraContactForce[1];
            total_local_elastic_contact_force[2] = LocalElasticContactForce[2] + LocalElasticExtraContactForce[2];
            CalculateOnContinuumContactElements(i, total_local_elastic_contact_force, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage, time_steps);
        }

        if (this->Is(DEMFlags::HAS_STRESS_TENSOR) /*&& (i < mContinuumInitialNeighborsSize)*/) {
            AddNeighbourContributionToStressTensor(r_process_info, TotalGlobalElasticContactForce, data_buffer.mLocalCoordSystem[2], data_buffer.mDistance, radius_sum, this);
        }

        AddContributionToRepresentativeVolume(data_buffer.mDistance, radius_sum, calculation_area);

        ComputeForceWithNeighbourFinalOperations();

        /*if (i < mContinuumInitialNeighborsSize) {
            DEM_COPY_SECOND_TO_FIRST_3(mArrayOfDeltaDisplacements[i], DeltDisp);
        }*/
    } // for each neighbor

    ComputeBrokenBondsRatio();

    KRATOS_CATCH("")
} //  ComputeBallToBallContactForce

void SplitForwardEulerSphericContinuumParticle::ComputeBallToRigidFaceContactForce(SphericParticle::ParticleDataBuffer & data_buffer,
                                                         array_1d<double, 3>& r_elastic_force,
                                                         array_1d<double, 3>& r_contact_force,
                                                         double& RollingResistance,
                                                         array_1d<double, 3>& rigid_element_force,
                                                         const ProcessInfo& r_process_info)
{
    KRATOS_TRY

    RenewData();

    std::vector<DEMWall*>& rNeighbours   = this->mNeighbourRigidFaces;
    array_1d<double, 3> velocity         = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& AngularVel = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

    for (unsigned int i = 0; i < rNeighbours.size(); i++) {
        DEMWall* wall = rNeighbours[i];
        if(wall == NULL) continue;
        if(this->Is(DEMFlags::STICKY)) {
            DEMIntegrationScheme& dem_scheme = this->GetTranslationalIntegrationScheme();
            GluedToWallScheme* p_glued_scheme = dynamic_cast<GluedToWallScheme*>(&dem_scheme);
            Condition* p_condition = p_glued_scheme->pGetCondition();
            if(p_condition == wall) continue;
        }
        if(wall->IsPhantom()){
            wall->CheckSide(this);
            continue;
        }

        double LocalElasticContactForce[3]       = {0.0};
        double GlobalElasticContactForce[3]      = {0.0};
        double ViscoDampingLocalContactForce[3]  = {0.0};
        double cohesive_force                    =  0.0;
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)
        array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
        array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);
        bool sliding = false;

        double ini_delta = GetInitialDeltaWithFEM(i);
        double DistPToB = 0.0;

        int ContactType = -1;
        array_1d<double, 4>& Weight = this->mContactConditionWeights[i];

        rNeighbours[i]->ComputeConditionRelativeData(i, this, data_buffer.mLocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

        if (ContactType == 1 || ContactType == 2 || ContactType == 3) {

            double indentation = -(DistPToB - GetInteractionRadius()) - ini_delta;
            double DeltDisp[3] = {0.0};
            double DeltVel [3] = {0.0};

            DeltVel[0] = velocity[0] - wall_velocity_at_contact_point[0];
            DeltVel[1] = velocity[1] - wall_velocity_at_contact_point[1];
            DeltVel[2] = velocity[2] - wall_velocity_at_contact_point[2];

            // For translation movement delta displacement
            const array_1d<double, 3>& delta_displ  = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            DeltDisp[0] = delta_displ[0] - wall_delta_disp_at_contact_point[0];
            DeltDisp[1] = delta_displ[1] - wall_delta_disp_at_contact_point[1];
            DeltDisp[2] = delta_displ[2] - wall_delta_disp_at_contact_point[2];

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                const array_1d<double,3> negative_delta_rotation = -1.0*GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
                array_1d<double, 3> actual_arm_vector, old_arm_vector;
                actual_arm_vector[0] = -data_buffer.mLocalCoordSystem[2][0] * DistPToB;
                actual_arm_vector[1] = -data_buffer.mLocalCoordSystem[2][1] * DistPToB;
                actual_arm_vector[2] = -data_buffer.mLocalCoordSystem[2][2] * DistPToB;

                double tangential_vel[3]           = {0.0};
                double tangential_displacement_due_to_rotation[3]  = {0.0};
                GeometryFunctions::CrossProduct(AngularVel, actual_arm_vector, tangential_vel);

                Quaternion<double> NegativeDeltaOrientation = Quaternion<double>::Identity();
                GeometryFunctions::OrientationFromRotationAngle(NegativeDeltaOrientation, negative_delta_rotation);

                NegativeDeltaOrientation.RotateVector3(actual_arm_vector, old_arm_vector);

                // Contribution of the rotation
                tangential_displacement_due_to_rotation[0] = (actual_arm_vector[0] - old_arm_vector[0]);
                tangential_displacement_due_to_rotation[1] = (actual_arm_vector[1] - old_arm_vector[1]);
                tangential_displacement_due_to_rotation[2] = (actual_arm_vector[2] - old_arm_vector[2]);

                DEM_ADD_SECOND_TO_FIRST(DeltVel, tangential_vel)
                DEM_ADD_SECOND_TO_FIRST(DeltDisp, tangential_displacement_due_to_rotation)
            }

            double LocalDeltDisp[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltDisp, LocalDeltDisp);

            double OldLocalElasticContactForce[3] = {0.0};

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, mNeighbourRigidFacesElasticContactForce[i], OldLocalElasticContactForce);
            const double previous_indentation = indentation + LocalDeltDisp[2];
            data_buffer.mLocalRelVel[0] = 0.0;
            data_buffer.mLocalRelVel[1] = 0.0;
            data_buffer.mLocalRelVel[2] = 0.0;

            if (indentation > 0.0) {

                GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltVel, data_buffer.mLocalRelVel);
                mDiscontinuumConstitutiveLaw->CalculateForcesRayleighWithFEM(r_process_info,
                                                                    OldLocalElasticContactForce,
                                                                    LocalElasticContactForce,
                                                                    LocalDeltDisp,
                                                                    data_buffer.mLocalRelVel,
                                                                    indentation,
                                                                    previous_indentation,
                                                                    ViscoDampingLocalContactForce,
                                                                    cohesive_force,
                                                                    this,
                                                                    wall,
                                                                    sliding);
            }

            double LocalContactForce[3]  = {0.0};
            double GlobalContactForce[3] = {0.0};

            AddUpFEMForcesAndProject(data_buffer.mLocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce,
                                     GlobalElasticContactForce, ViscoDampingLocalContactForce, cohesive_force, r_elastic_force,
                                     r_contact_force, mNeighbourRigidFacesElasticContactForce[i], mNeighbourRigidFacesTotalContactForce[i]);

            rigid_element_force[0] -= GlobalContactForce[0];
            rigid_element_force[1] -= GlobalContactForce[1];
            rigid_element_force[2] -= GlobalContactForce[2];

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                ComputeMoments(LocalContactForce[2], GlobalContactForce, RollingResistance, data_buffer.mLocalCoordSystem[2], this, indentation, true, i); //WARNING: sending itself as the neighbor!!
            }

            //WEAR
            if (wall->GetProperties()[COMPUTE_WEAR]) {
                const double area              = Globals::Pi * GetInteractionRadius() * GetInteractionRadius();
                const double inverse_of_volume = 1.0 / (4.0 * 0.333333333333333 * area * GetInteractionRadius());
                ComputeWear(data_buffer.mLocalRelVel,
                            data_buffer.mDt,
                            sliding,
                            inverse_of_volume,
                            LocalElasticContactForce[2],
                            wall);
            } //wall->GetProperties()[COMPUTE_WEAR] if

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
                AddWallContributionToStressTensor(GlobalElasticContactForce, data_buffer.mLocalCoordSystem[2], DistPToB, 0.0);
            }
        } //ContactType if
    } //rNeighbours.size loop

    auto& list_of_point_condition_pointers = this->GetValue(WALL_POINT_CONDITION_POINTERS);
    auto& neighbour_point_faces_elastic_contact_force = this->GetValue(WALL_POINT_CONDITION_ELASTIC_FORCES);
    auto& neighbour_point_faces_total_contact_force = this->GetValue(WALL_POINT_CONDITION_TOTAL_FORCES);

    for (unsigned int i = 0; i < list_of_point_condition_pointers.size(); i++) {
        Condition* wall_condition = list_of_point_condition_pointers[i];

        array_1d<double, 3> wall_coordinates = wall_condition->GetGeometry().Center();

        double RelVel[3] = { 0.0 };
        double LocalElasticContactForce[3] = { 0.0 };
        double GlobalElasticContactForce[3] = { 0.0 };
        double ViscoDampingLocalContactForce[3] = { 0.0 };
        double cohesive_force = 0.0;

        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)

        const Matrix& r_N = wall_condition->GetGeometry().ShapeFunctionsValues();

        array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
        for (IndexType i = 0; i < wall_condition->GetGeometry().size(); ++i) {
            wall_delta_disp_at_contact_point -= r_N(0, i)*wall_condition->GetGeometry()[i].GetSolutionStepValue(DELTA_DISPLACEMENT, 0);
        }

        array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);
        for (IndexType i = 0; i < wall_condition->GetGeometry().size(); ++i) {
            wall_velocity_at_contact_point += r_N(0, i) * wall_condition->GetGeometry()[i].GetSolutionStepValue(VELOCITY, 0);
        }

        bool sliding = false;
        double ini_delta = GetInitialDeltaWithFEM(i);

        array_1d<double, 3> cond_to_me_vect;
        noalias(cond_to_me_vect) = GetGeometry()[0].Coordinates() - wall_coordinates;

        double DistPToB = DEM_MODULUS_3(cond_to_me_vect);

        double indentation = -(DistPToB - GetInteractionRadius()) - ini_delta;

        double DeltDisp[3] = { 0.0 };
        double DeltVel[3] = { 0.0 };

        DeltVel[0] = velocity[0] - wall_velocity_at_contact_point[0];
        DeltVel[1] = velocity[1] - wall_velocity_at_contact_point[1];
        DeltVel[2] = velocity[2] - wall_velocity_at_contact_point[2];

        // For translation movement delta displacement
        const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        DeltDisp[0] = delta_displ[0] - wall_delta_disp_at_contact_point[0];
        DeltDisp[1] = delta_displ[1] - wall_delta_disp_at_contact_point[1];
        DeltDisp[2] = delta_displ[2] - wall_delta_disp_at_contact_point[2];

        Node<3>::Pointer other_particle_node = this->GetGeometry()[0].Clone();
        other_particle_node->GetSolutionStepValue(DELTA_DISPLACEMENT) = wall_delta_disp_at_contact_point;
        data_buffer.mpOtherParticleNode = &*other_particle_node;
        DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mOtherToMeVector, cond_to_me_vect)
        data_buffer.mDistance = DistPToB;
        data_buffer.mDomainIsPeriodic = false;
        EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel, data_buffer.mLocalCoordSystem, data_buffer.mOldLocalCoordSystem, velocity, delta_displ);

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            const array_1d<double, 3>& delta_rotation = GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);

            array_1d<double, 3> actual_arm_vector;
            actual_arm_vector[0] = -data_buffer.mLocalCoordSystem[2][0] * DistPToB;
            actual_arm_vector[1] = -data_buffer.mLocalCoordSystem[2][1] * DistPToB;
            actual_arm_vector[2] = -data_buffer.mLocalCoordSystem[2][2] * DistPToB;

            double tangential_vel[3] = { 0.0 };
            double tangential_displacement_due_to_rotation[3] = { 0.0 };
            GeometryFunctions::CrossProduct(AngularVel, actual_arm_vector, tangential_vel);
            GeometryFunctions::CrossProduct(delta_rotation, actual_arm_vector, tangential_displacement_due_to_rotation);

            DEM_ADD_SECOND_TO_FIRST(DeltVel, tangential_vel)
            DEM_ADD_SECOND_TO_FIRST(DeltDisp, tangential_displacement_due_to_rotation)
        }

        double LocalDeltDisp[3] = { 0.0 };
        GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltDisp, LocalDeltDisp);

        double OldLocalElasticContactForce[3] = { 0.0 };

        GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, neighbour_point_faces_elastic_contact_force[i], OldLocalElasticContactForce);
        const double previous_indentation = indentation + LocalDeltDisp[2];
        data_buffer.mLocalRelVel[0] = 0.0;
        data_buffer.mLocalRelVel[1] = 0.0;
        data_buffer.mLocalRelVel[2] = 0.0;

        if (indentation > 0.0)
        {
            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltVel, data_buffer.mLocalRelVel);

            mDiscontinuumConstitutiveLaw->CalculateForcesRayleighWithFEM(
                r_process_info,
                OldLocalElasticContactForce,
                LocalElasticContactForce,
                LocalDeltDisp,
                data_buffer.mLocalRelVel,
                indentation,
                previous_indentation,
                ViscoDampingLocalContactForce,
                cohesive_force,
                this,
                wall_condition,
                sliding);
        }

        double LocalContactForce[3] = { 0.0 };
        double GlobalContactForce[3] = { 0.0 };

        AddUpFEMForcesAndProject(
            data_buffer.mLocalCoordSystem,
            LocalContactForce,
            LocalElasticContactForce,
            GlobalContactForce,
            GlobalElasticContactForce,
            ViscoDampingLocalContactForce,
            cohesive_force,
            r_elastic_force,
            r_contact_force,
            neighbour_point_faces_elastic_contact_force[i],
            neighbour_point_faces_total_contact_force[i]);

        rigid_element_force[0] -= GlobalContactForce[0];
        rigid_element_force[1] -= GlobalContactForce[1];
        rigid_element_force[2] -= GlobalContactForce[2];

        Vector GlobalContactForceVector(3);
        DEM_COPY_SECOND_TO_FIRST_3(GlobalContactForceVector, GlobalContactForce)
        wall_condition->SetValue(EXTERNAL_FORCES_VECTOR, GlobalContactForceVector);

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            ComputeMoments(LocalContactForce[2], GlobalContactForce, RollingResistance, data_buffer.mLocalCoordSystem[2], this, indentation, true, i); //WARNING: sending itself as the neighbor!!
        }

        if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
            AddWallContributionToStressTensor(GlobalElasticContactForce, data_buffer.mLocalCoordSystem[2], DistPToB, 0.0);
        }
    }

    KRATOS_CATCH("")
}// ComputeBallToRigidFaceContactForce

void SplitForwardEulerSphericContinuumParticle::ComputeBallToBallStiffness(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                            double& r_nodal_stiffness,
                                                                            double& r_nodal_rotational_stiffness)
{
    KRATOS_TRY

    Vector& cont_ini_neigh_area = this->GetValue(NEIGHBOURS_CONTACT_AREAS);

    for (int i = 0; data_buffer.SetNextNeighbourOrExit(i); ++i) {

        if (mNeighbourElements[i] == NULL) continue;
        if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;
        SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
        data_buffer.mpOtherParticle = neighbour_iterator;

        noalias(data_buffer.mOtherToMeVector) = this->GetGeometry()[0].Coordinates() - data_buffer.mpOtherParticle->GetGeometry()[0].Coordinates();
        const double& other_radius = data_buffer.mpOtherParticle->GetRadius();
        data_buffer.mDistance = DEM_MODULUS_3(data_buffer.mOtherToMeVector);
        double radius_sum = GetRadius() + other_radius;
        double initial_delta = GetInitialDelta(i);
        double initial_dist = radius_sum - initial_delta;
        double indentation = initial_dist - data_buffer.mDistance;
        double myYoung = GetYoung();
        double myPoisson = GetPoisson();

        // Getting neighbor properties
        double other_young = data_buffer.mpOtherParticle->GetYoung();
        double other_poisson = data_buffer.mpOtherParticle->GetPoisson();
        double equiv_poisson;
        if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
        else { equiv_poisson = 0.0; }
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
        double calculation_area = 0.0;

        double normal_stiffness = 0.0;
        double tangential_stiffness = 0.0;

        if (i < (int)mContinuumInitialNeighborsSize) {
            mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, cont_ini_neigh_area, i, calculation_area);
            mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(normal_stiffness, tangential_stiffness, initial_dist, equiv_young, equiv_poisson, calculation_area, this, neighbour_iterator);
        }
        r_nodal_stiffness += CalculateStiffnessNorm(normal_stiffness,tangential_stiffness);

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            if (i < (int)mContinuumInitialNeighborsSize && mIniNeighbourFailureId[i] == 0) {
                // TODO: is this right ?
                double arm_length = GetInteractionRadius() - indentation;
                double normal_rotational_stiffness = normal_stiffness * arm_length*1.0;
                double tangential_rotational_stiffness = tangential_stiffness * arm_length*1.0;

                r_nodal_rotational_stiffness += CalculateStiffnessNorm(normal_rotational_stiffness,tangential_rotational_stiffness);
            }
        }
    }// for each neighbor

    KRATOS_CATCH("")
}

void SplitForwardEulerSphericContinuumParticle::ComputeBallToRigidFaceStiffness(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                                double& r_nodal_stiffness,
                                                                                double& r_nodal_rotational_stiffness)
{
    KRATOS_TRY

    std::vector<DEMWall*>& rNeighbours = this->mNeighbourRigidFaces;

    for (unsigned int i = 0; i < rNeighbours.size(); i++) {
        DEMWall* wall = rNeighbours[i];
        if(wall == NULL) continue;
        if(this->Is(DEMFlags::STICKY)) {
            DEMIntegrationScheme& dem_scheme = this->GetTranslationalIntegrationScheme();
            GluedToWallScheme* p_glued_scheme = dynamic_cast<GluedToWallScheme*>(&dem_scheme);
            Condition* p_condition = p_glued_scheme->pGetCondition();
            if(p_condition == wall) continue;
        }
        if(wall->IsPhantom()){
            wall->CheckSide(this);
            continue;
        }

        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
        array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
        array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);

        double ini_delta = GetInitialDeltaWithFEM(i);
        double DistPToB = 0.0;

        int ContactType = -1;
        array_1d<double, 4>& Weight = this->mContactConditionWeights[i];

        rNeighbours[i]->ComputeConditionRelativeData(i, this, data_buffer.mLocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

        if (ContactType == 1 || ContactType == 2 || ContactType == 3) {

            double indentation = -(DistPToB - GetInteractionRadius()) - ini_delta;

            if (indentation > 0.0) {

                double normal_stiffness, tangential_stiffness;
                mDiscontinuumConstitutiveLaw->InitializeContactWithFEM(this,wall,indentation);
                normal_stiffness = mDiscontinuumConstitutiveLaw->mKn;
                tangential_stiffness = mDiscontinuumConstitutiveLaw->mKt;

                r_nodal_stiffness += CalculateStiffnessNorm(normal_stiffness,tangential_stiffness);

                if (this->Is(DEMFlags::HAS_ROTATION)) {
                    // TODO: is this right ?
                    double arm_length = GetInteractionRadius() - indentation;
                    double normal_rotational_stiffness = normal_stiffness * arm_length*1.0;
                    double tangential_rotational_stiffness = tangential_stiffness * arm_length*1.0;

                    r_nodal_rotational_stiffness += CalculateStiffnessNorm(normal_rotational_stiffness,tangential_rotational_stiffness);
                }
            }
        } //ContactType if
    } //rNeighbours.size loop

    KRATOS_CATCH("")
}

double SplitForwardEulerSphericContinuumParticle::CalculateStiffnessNorm(const double& r_normal_stiffness, const double& r_tangential_stiffness) {
    return std::sqrt(r_normal_stiffness*r_normal_stiffness+2.0*r_tangential_stiffness*r_tangential_stiffness);
}

}  // namespace Kratos.
