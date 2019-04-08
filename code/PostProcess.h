/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */

// Computes based on the system displacement solution the total strain and stress
// per element, Integrates the local forces and moments and transforms them in the
// Global Cartesian system. Thereupon it computes the effective axial and
// torsion related fascicle modulus 
#ifndef _POSTPROCESS_H_
#define _POSTPROCESS_H_ 1

#include "DefsMTUsed.h"

class PostProcessing
{
public:
  Vector3d Force_Z_Global; //=Vector3d::Zero();
  Vector3d Moment_Z_Global; //=Vector3d::Zero();
  double Force_Z;
  double Moment_Z;
  double E_eff, Mom_Eff;

public:
  PostProcessing();
  Vector3d Transform_Local_To_Global(const double a, const double b, const double g, const double k, const double t,  const double Xic, const double Etac, const double Force_L_1, const double Force_L_2, const double Force_L_3);
  void Compute_Axial_Force(InputParams& InputData, int N_Elements, vector<ElemIntegrCoord>& ElemCoords, Double2DVec& Connectivity,vector<Elem_Stiff_Data>& Stiff_Data, VectorXd& Elem_Ax_Strain, VectorXd& U_Global,
    double Ef, double Em, VectorXi& Identities);

  // Constructor used with arguments Input and meshing data classes
  PostProcessing(InputParams& InputData, int N_Elements, vector<ElemIntegrCoord>& ElemCoords, Double2DVec& Connectivity, vector<Elem_Stiff_Data>& Stiff_Data, VectorXd& Elem_Ax_Strain, VectorXd& U_Global, VectorXi& Identities)
  {
    Force_Z_Global=Vector3d::Zero();
    Moment_Z_Global=Vector3d::Zero();


    // Compute the Axial Loading Vector
    Compute_Axial_Force(InputData, N_Elements,ElemCoords, Connectivity, Stiff_Data, Elem_Ax_Strain, U_Global, InputData.Ef, InputData.Em, Identities);
  }
  ~PostProcessing() {};
};

Vector3d PostProcessing::Transform_Local_To_Global(const double a, const double b, const double g, const double k, const double t, const double Xic,const double Etac, const double Force_L_1, const double Force_L_2, const double Force_L_3 )
{
  Matrix3d Transf_Matrix=Matrix3d::Zero();
  Vector3d Local_Forces=Vector3d::Zero();

  Local_Forces(0)=Force_L_1;  Local_Forces(1)=Force_L_2; Local_Forces(2)=Force_L_3;
  Transf_Matrix<< -1, 0, t*Etac, 0, -b/g, (-Xic+a)/g, 0, a/g, b/g;

  return Transf_Matrix*Local_Forces;
}

void PostProcessing::Compute_Axial_Force(InputParams& InputData,int N_Elements, vector<ElemIntegrCoord>& ElemCoords, Double2DVec& Connectivity, vector<Elem_Stiff_Data>& Stiff_Data, VectorXd& Elem_Ax_Strain, VectorXd& U_Global,
  double Ef, double Em, VectorXi &Identities)
{
  int ILO, ISY;
  VectorXd Elem_DOfs, U_Elem;
  VectorXd Elem_T_Strain, Elem_T_Stress;
  double Moment_Local_1, Moment_Local_2, Moment_Local_3;
  double Force_Local_1, Force_Local_2, Force_Local_3;
  double Force_Y_Check, Force_Z_Check, Moment_Z_Check;

  // Allocations
  Force_Z_Check=0; Moment_Z_Check=0;
  Elem_DOfs.resize(9);
  U_Elem.resize(9);
  Elem_T_Strain.resize(6);
  Elem_T_Stress.resize(6);

  // Element loop and integrations
  for (int it=0; it<N_Elements; ++it)
  {
    // Loop to get the Dof IDS of each element
    for (unsigned int k=0; k<3; ++k)
    {
      ILO=k*3;
      ISY=3*(Connectivity[it][k]-1);
      for (unsigned int l=0; l<3; ++l)
      {
        Elem_DOfs(ILO+l)=ISY+l;
        U_Elem(ILO+l)=U_Global(Elem_DOfs(ILO+l));
      }
    }

    // Compute the strain and stress per element and  integrate
    Elem_T_Strain=Stiff_Data[it].Beta_M*U_Elem+Elem_Ax_Strain;
    Elem_T_Stress=Stiff_Data[it].C_Tensor*Elem_T_Strain;

    //Compute the Moments and forces in the local system
    Moment_Local_1=Elem_T_Stress(2)*0.5*Stiff_Data[it].Det_Jacob*ElemCoords[it].Etac; //Xic;
    Moment_Local_2=-Elem_T_Stress(2)*0.5*Stiff_Data[it].Det_Jacob*ElemCoords[it].Xic;
    Moment_Local_3=(Elem_T_Stress(1)*ElemCoords[it].Xic-Elem_T_Stress(0)*ElemCoords[it].Etac)*0.5*Stiff_Data[it].Det_Jacob;

    Force_Local_1=Elem_T_Stress(4)*0.5*Stiff_Data[it].Det_Jacob;
    Force_Local_2=Elem_T_Stress(3)*0.5*Stiff_Data[it].Det_Jacob;
    Force_Local_3=Elem_T_Stress(2)*0.5*Stiff_Data[it].Det_Jacob;

    Force_Z_Check=Force_Z_Check+Force_Local_3;
    Force_Y_Check= Force_Y_Check+Force_Local_2;
    Moment_Z_Check=Moment_Z_Check+Moment_Local_3;

    // Compute total force vector and total moment vector
    Force_Z_Global=Force_Z_Global+Transform_Local_To_Global(InputData.a, InputData.b, InputData.gam, InputData.k, InputData.t, ElemCoords[it].Xic, ElemCoords[it].Etac, Force_Local_1, Force_Local_2, Force_Local_3 );

    Moment_Z_Global=Moment_Z_Global+Transform_Local_To_Global(InputData.a, InputData.b, InputData.gam, InputData.k, InputData.t, ElemCoords[it].Xic, ElemCoords[it].Etac, Moment_Local_1, Moment_Local_2, Moment_Local_3 );

  } // Loop over elements
  Force_Z=Force_Z_Check;
  Moment_Z=Moment_Z_Global(2)+Force_Z_Global(1)*InputData.a;
  // Calculate the effective fascicle axial modulus and normalized effective moment
  E_eff=abs(Force_Z/InputData.AppliedExternalStrain)/(M_PI*pow(InputData.r,2));
  Mom_Eff=abs(Moment_Z/InputData.AppliedExternalStrain)/(M_PI*pow(InputData.r,2));
  
  
}

#endif
