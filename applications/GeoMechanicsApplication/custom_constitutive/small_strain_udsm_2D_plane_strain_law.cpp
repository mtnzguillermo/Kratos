// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// System includes

// External includes

#include "custom_constitutive/small_strain_udsm_2D_plane_strain_law.hpp"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainUDSM2DPlaneStrainLaw::SmallStrainUDSM2DPlaneStrainLaw()
   : SmallStrainUDSM3DLaw()
   {
    KRATOS_TRY;

    KRATOS_CATCH("")

   }

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
SmallStrainUDSM2DPlaneStrainLaw::
   SmallStrainUDSM2DPlaneStrainLaw(const SmallStrainUDSM2DPlaneStrainLaw &rOther)
   : SmallStrainUDSM3DLaw(rOther)
{
   KRATOS_TRY;

   KRATOS_CATCH("");
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainUDSM2DPlaneStrainLaw::Clone() const
{
   KRATOS_TRY;

   return Kratos::make_shared<SmallStrainUDSM2DPlaneStrainLaw>(*this);

   KRATOS_CATCH("");
}

//********************************ASSIGNMENT******************************************
//************************************************************************************
SmallStrainUDSM2DPlaneStrainLaw 
  &SmallStrainUDSM2DPlaneStrainLaw::operator=(SmallStrainUDSM2DPlaneStrainLaw const &rOther)
{
   KRATOS_TRY;

   SmallStrainUDSM3DLaw::operator=(rOther);

   return *this;

   KRATOS_CATCH("");
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallStrainUDSM2DPlaneStrainLaw::~SmallStrainUDSM2DPlaneStrainLaw() {}


//************************************************************************************
//************************************************************************************

void SmallStrainUDSM2DPlaneStrainLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues)
{
   const Vector& rStrainVector = rValues.GetStrainVector();

   mDeltaStrainVector[INDEX_3D_XX] = rStrainVector(INDEX_2D_PLANE_STRESS_XX) - mStrainVectorFinalized[INDEX_3D_XX];
   mDeltaStrainVector[INDEX_3D_YY] = rStrainVector(INDEX_2D_PLANE_STRESS_YY) - mStrainVectorFinalized[INDEX_3D_YY];
   mDeltaStrainVector[INDEX_3D_XY] = rStrainVector(INDEX_2D_PLANE_STRESS_XY) - mStrainVectorFinalized[INDEX_3D_XY];
}

void SmallStrainUDSM2DPlaneStrainLaw::SetExternalStressVector(Vector& rStressVector)
{
   KRATOS_TRY;

   rStressVector(INDEX_2D_PLANE_STRAIN_XX) = mStressVector[INDEX_3D_XX];
   rStressVector(INDEX_2D_PLANE_STRAIN_YY) = mStressVector[INDEX_3D_YY];
   rStressVector(INDEX_2D_PLANE_STRAIN_XY) = mStressVector[INDEX_3D_XY];
   rStressVector(INDEX_2D_PLANE_STRAIN_ZZ) = mStressVector[INDEX_3D_ZZ];

   KRATOS_CATCH("")
}


void SmallStrainUDSM2DPlaneStrainLaw::SetInternalStressVector(const Vector& rStressVector)
{
   KRATOS_TRY;
   std::fill(mStressVectorFinalized.begin(), mStressVectorFinalized.end(), 0.0);

   mStressVectorFinalized[INDEX_3D_XX] = rStressVector(INDEX_2D_PLANE_STRAIN_XX);
   mStressVectorFinalized[INDEX_3D_YY] = rStressVector(INDEX_2D_PLANE_STRAIN_YY);
   mStressVectorFinalized[INDEX_3D_XY] = rStressVector(INDEX_2D_PLANE_STRAIN_XY);
   mStressVectorFinalized[INDEX_3D_ZZ] = rStressVector(INDEX_2D_PLANE_STRAIN_ZZ);

   KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
   KRATOS_TRY;
   std::fill(mStrainVectorFinalized.begin(), mStrainVectorFinalized.end(), 0.0);

   mStrainVectorFinalized[INDEX_3D_XX] = rStrainVector(INDEX_2D_PLANE_STRESS_XX);
   mStrainVectorFinalized[INDEX_3D_YY] = rStrainVector(INDEX_2D_PLANE_STRESS_YY);
   mStrainVectorFinalized[INDEX_3D_XY] = rStrainVector(INDEX_2D_PLANE_STRESS_XY);

   KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::CopyConstitutiveMatrix( ConstitutiveLaw::Parameters &rValues,
                                                              Matrix& rConstitutiveMatrix )
{
   KRATOS_TRY;

   if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM])
   {
      // transfer fortran style matrix to C++ style
      for (unsigned int i = 0; i < VoigtSizePlaneStress; i++) {
         for (unsigned int j = 0; j < VoigtSizePlaneStress; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[getIndex3D(static_cast<indexStress2DPlaneStress>(j))][getIndex3D(static_cast<indexStress2DPlaneStress>(i))];
         }
      }
   }
   else
   {
      for (unsigned int i = 0; i < VoigtSizePlaneStress; i++) {
         for (unsigned int j = 0; j < VoigtSizePlaneStress; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[getIndex3D(static_cast<indexStress2DPlaneStress>(i))][getIndex3D(static_cast<indexStress2DPlaneStress>(j))];
         }
      }
   }

   KRATOS_CATCH("")
}

indexStress3D SmallStrainUDSM2DPlaneStrainLaw::getIndex3D(indexStress2DPlaneStress index2D)
{
   KRATOS_TRY;

   switch (index2D)
   {
      case INDEX_2D_PLANE_STRESS_XX:
        return INDEX_3D_XX;
      case INDEX_2D_PLANE_STRESS_YY:
        return INDEX_3D_YY;
      case INDEX_2D_PLANE_STRESS_XY:
        return INDEX_3D_XY;
      default:
        KRATOS_THROW_ERROR(std::invalid_argument, "invalid index: ", index2D);
   }

   KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters &rValues)
{
   const Vector& rStrainVector = rValues.GetStrainVector();
   std::fill(mStrainVectorFinalized.begin(), mStrainVectorFinalized.end(), 0.0);

   mStrainVectorFinalized[INDEX_3D_XX] = rStrainVector(INDEX_2D_PLANE_STRESS_XX);
   mStrainVectorFinalized[INDEX_3D_YY] = rStrainVector(INDEX_2D_PLANE_STRESS_YY);
   mStrainVectorFinalized[INDEX_3D_XY] = rStrainVector(INDEX_2D_PLANE_STRESS_XY);

}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainUDSM2DPlaneStrainLaw::
   CalculateCauchyGreenStrain( ConstitutiveLaw::Parameters& rValues,
                               Vector& rStrainVector )
{
   //1.-Compute total deformation gradient
   const Matrix& F = rValues.GetDeformationGradientF();

   // for shells/membranes in case the DeformationGradient is of size 3x3
   BoundedMatrix<double, 2, 2> F2x2;
   for (unsigned int i = 0; i<2; ++i)
      for (unsigned int j = 0; j<2; ++j)
         F2x2(i, j) = F(i, j);

   Matrix E_tensor = prod(trans(F2x2), F2x2);

   for (unsigned int i = 0; i<2; ++i)
      E_tensor(i, i) -= 1.0;

   E_tensor *= 0.5;
   noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

//----------------------------------------------------------------------------------------
Vector& SmallStrainUDSM2DPlaneStrainLaw::
   GetValue( const Variable<Vector> &rThisVariable,
             Vector &rValue )
{
   // KRATOS_INFO("0-SmallStrainUDSM2DPlaneStrainLaw::GetValue()") << std::endl;

   if (rThisVariable == STATE_VARIABLES)
   {
      SmallStrainUDSM3DLaw::GetValue(rThisVariable, rValue );
   }
   else if (rThisVariable == CAUCHY_STRESS_VECTOR)
   {
      if (rValue.size() != VoigtSize)
         rValue.resize(VoigtSize);

      rValue[INDEX_2D_PLANE_STRAIN_XX] = mStressVectorFinalized[INDEX_3D_XX];
      rValue[INDEX_2D_PLANE_STRAIN_YY] = mStressVectorFinalized[INDEX_3D_YY];
      rValue[INDEX_2D_PLANE_STRAIN_XY] = mStressVectorFinalized[INDEX_3D_XY];
      rValue[INDEX_2D_PLANE_STRAIN_ZZ] = mStressVectorFinalized[INDEX_3D_ZZ];

   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::GetValue()") << std::endl;

    return rValue;
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM2DPlaneStrainLaw::SetValue( const Variable<Vector>& rThisVariable,
                                                const Vector& rValue,
                                                const ProcessInfo& rCurrentProcessInfo )
{
   // KRATOS_INFO("02-SmallStrainUDSM2DPlaneStrainLaw::SetValue()") << std::endl;

   if (rThisVariable == STATE_VARIABLES)
   {
      SmallStrainUDSM3DLaw::SetValue(rThisVariable, rValue, rCurrentProcessInfo );
   }
   else if (rThisVariable == CAUCHY_STRESS_VECTOR)
   {
      if (rValue.size() == VoigtSize) 
      {
         this->SetInternalStressVector(rValue);
      }
   }

   // KRATOS_INFO("12-SmallStrainUDSM2DPlaneStrainLaw::SetValue()") << std::endl;
}

} // Namespace Kratos
