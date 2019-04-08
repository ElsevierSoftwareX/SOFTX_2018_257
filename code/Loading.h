/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */


// Compute the Loading vector on the local Serret-Frenet axis frame

#ifndef _LOADING_H_
#define _LOADING_H_ 1

#include "DefsMTUsed.h"

class Loading
{
public:
  double Epsilon_Z;
  VectorXd Elem_Ax_Strain;
  VectorXd AxialLoad;

public:
  Loading();
  void Compute_Axial_Load(const int N_Elems, const int N_dofs, const Double2DVec Connectivity, InputParams& InputData, Elem_Stiff_Mat& Stiffness);

  // Constructor used with arguments Input and meshing data classes
  Loading(InputParams& InputData, MeshingData& MeshingData_C, Elem_Stiff_Mat& Stiff):Epsilon_Z(InputData.AppliedExternalStrain)
  {
    // Compute the Axial Loading Vector
    Compute_Axial_Load(MeshingData_C.N_Elements, MeshingData_C.N_dofs, MeshingData_C.Connectivity, InputData, Stiff);
  }
  ~Loading() {};
};

void Loading::Compute_Axial_Load(const int N_Elems, const int N_dofs, const Double2DVec Connectivity, InputParams& InputData, Elem_Stiff_Mat& Stiffness)
{
  vector<int> Node_Dofs;
  VectorXd  ElemLoad;
  MatrixXd CheckMult;

  // Inititialize the global Loading vector
  AxialLoad.resize(N_dofs);
  AxialLoad << VectorXd::Zero(N_dofs);

  // Inititialize the local element strain vector
  Elem_Ax_Strain.resize(6);
  Elem_Ax_Strain << VectorXd::Zero(6);
  Elem_Ax_Strain(0)=0;  Elem_Ax_Strain(1)=InputData.k*InputData.a*Epsilon_Z; Elem_Ax_Strain(2)=InputData.b*InputData.t*Epsilon_Z;
  Elem_Ax_Strain(3)=2*0.5*(InputData.k*InputData.b+InputData.t*InputData.a)*Epsilon_Z;  Elem_Ax_Strain(4)=0; Elem_Ax_Strain(5)=0;

  // Initialize the Node_Dof and ElemLoad vector
  Node_Dofs.resize(9);
  ElemLoad.resize(9);
  CheckMult.resize(6,6);

  // Loop over all the elements to create the load vector
  // N_Elems
  for (int it=0; it<N_Elems; ++it)
  {
    // Get the local DOFS out of the connectivity matrix
    Node_Dofs[0]=(Connectivity[it][0]-1)*3+1;
    Node_Dofs[1]=(Connectivity[it][0]-1)*3+2;
    Node_Dofs[2]=(Connectivity[it][0]-1)*3+3;
    Node_Dofs[3]=(Connectivity[it][1]-1)*3+1;
    Node_Dofs[4]=(Connectivity[it][1]-1)*3+2;
    Node_Dofs[5]=(Connectivity[it][1]-1)*3+3;
    Node_Dofs[6]=(Connectivity[it][2]-1)*3+1;
    Node_Dofs[7]=(Connectivity[it][2]-1)*3+2;
    Node_Dofs[8]=(Connectivity[it][2]-1)*3+3;

    //Compute element load vector
    ElemLoad << VectorXd::Zero(9);
    ElemLoad=-Stiffness.Stiff_Data[it].Beta_M.transpose()*Stiffness.Stiff_Data[it].C_Tensor*Elem_Ax_Strain*Stiffness.Stiff_Data[it].Det_Jacob;
    for (int DOF_ID=0; DOF_ID<9; ++DOF_ID)
    {
      AxialLoad(Node_Dofs[DOF_ID]-1)=AxialLoad(Node_Dofs[DOF_ID]-1)+ElemLoad(DOF_ID);
    }
  } // Loop over elements
}

#endif
