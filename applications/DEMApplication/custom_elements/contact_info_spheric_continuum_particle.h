//   $Author: Joaquín Irazábal jirazabal@cimne.upc.edu

#if !defined(KRATOS_CONTACT_INFO_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED)
#define  KRATOS_CONTACT_INFO_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "spheric_continuum_particle.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) ContactInfoSphericContinuumParticle : public SphericContinuumParticle
    {
        public:

        /// Pointer definition of ContactInfoSphericContinuumParticle
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ContactInfoSphericContinuumParticle);

        ContactInfoSphericContinuumParticle();
        ContactInfoSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry);
        ContactInfoSphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes);
        ContactInfoSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
        ContactInfoSphericContinuumParticle(Element::Pointer p_continuum_spheric_particle);

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor.
        virtual ~ContactInfoSphericContinuumParticle(){}

        ContactInfoSphericContinuumParticle& operator=(const ContactInfoSphericContinuumParticle& rOther);

        virtual std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "ContactInfoSphericContinuumParticle" ;
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "ContactInfoSphericContinuumParticle";
        }

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const override
        {
	    //mpGeometry->PrintData(rOStream);
        }

        virtual double GetParticleInitialCohesion();
        void   SetParticleInitialCohesionFromProperties(double* particle_initial_cohesion);
        virtual double GetAmountOfCohesionFromStress();
        void   SetAmountOfCohesionFromStressFromProperties(double* amount_of_cohesion_from_stress);
        virtual double GetParticleConicalDamageContactRadius();
        void   SetParticleConicalDamageContactRadiusFromProperties(double* particle_contact_radius);
        virtual double GetParticleConicalDamageMaxStress();
        void   SetParticleConicalDamageMaxStressFromProperties(double* particle_max_stress);
        virtual double GetParticleConicalDamageGamma();
        void   SetParticleConicalDamageGammaFromProperties(double* particle_gamma);
        virtual double GetLevelOfFouling();
        void   SetLevelOfFoulingFromProperties(double* level_of_fouling);

        double SlowGetParticleInitialCohesion();
        double SlowGetAmountOfCohesionFromStress();
        double SlowGetParticleConicalDamageContactRadius();
        double SlowGetParticleConicalDamageMaxStress();
        double SlowGetParticleConicalDamageGamma();
        double SlowGetLevelOfFouling();

        std::vector<double> mNeighbourContactRadius;
        std::vector<double> mNeighbourRigidContactRadius;
        std::vector<double> mNeighbourIndentation;
        std::vector<double> mNeighbourRigidIndentation;
        std::vector<double> mNeighbourTgOfStatFriAng;
        std::vector<double> mNeighbourTgOfDynFriAng;
        std::vector<double> mNeighbourRigidTgOfStatFriAng;
        std::vector<double> mNeighbourRigidTgOfDynFriAng;
        std::vector<double> mNeighbourContactStress;
        std::vector<double> mNeighbourRigidContactStress;
        std::vector<double> mNeighbourCohesion;
        std::vector<double> mNeighbourRigidCohesion;

        protected:

        void ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids,
                                                std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces) override;

        void ComputeNewRigidFaceNeighboursHistoricalData() override;

        private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

    }; // Class ContactInfoSphericContinuumParticle

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, ContactInfoSphericContinuumParticle& rThis){return rIStream;}

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const ContactInfoSphericContinuumParticle& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

}  // namespace Kratos.

#endif // KRATOS_CONTACT_INFO_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED  defined
