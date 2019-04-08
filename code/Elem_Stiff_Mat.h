/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */
// Computes the element stiffness matrices based on the input and meshing 
// specifications previously defined. The finite element linear matrix operator B
// and material tensor C pertinent to the helical curvilinear domain are used.
#ifndef _ELEM_STIFF_MAT_H_
#define _ELEM_STIFF_MAT_H_

#include "DefsMTUsed.h"

struct Christ_Symb
{
  double g131, g311, g132, g312, g133, g313, g231, g321, g331, g332, g333;
};

struct Elem_Stiff_Data
{
  Christ_Symb Chris;
  double Det_Jacob;
  MatrixXd Beta_M; //=Eigen::MatrixXd::Zero(6,9);
  MatrixXd C_Tensor;//=Eigen::MatrixXd::Zero(6,6);
  MatrixXd Ke; //=Eigen::MatrixXd::Zero(9,9);
};

class Elem_Stiff_Mat
{
public:
  vector<Elem_Stiff_Data> Stiff_Data;
  MatrixXd K_total;

public:
  Elem_Stiff_Mat();
  void ComputeChristophel(InputParams InputData, MeshingData MeshingData_C);
  void Compute_Beta_M(InputParams InputData, MeshingData MeshingData_C);
  void Compute_Tensor_C(InputParams InputData, MeshingData  MeshingData_C, double E1, double E2, VectorXi& Identities);
  void Compute_Mat_Ke(const int N_Elems);
  void Print_Check_Function(const MatrixXd MatToPrint, const int NumberRows, const int NumberOfCols);
  void Compute_Total_Stiff(const int N_Elems, const int N_dofs, const Double2DVec Connectivity);

  // Constructor used with arguments Input and meshing data classes
  Elem_Stiff_Mat(InputParams InputData, MeshingData MeshingData_C, VectorXi& Identities)
  {
    // Allocate memory only for the outer struct
    Stiff_Data.reserve(MeshingData_C.Connectivity.size());

    // Compute the Christoffel Symbols
    ComputeChristophel(InputData, MeshingData_C);

    // Compute the Beta matrix for each element
    Compute_Beta_M(InputData, MeshingData_C);

    // Compute the C matric for each element
    Compute_Tensor_C(InputData, MeshingData_C, InputData.Ef, InputData.Em, Identities);

    // Compute the element stiffness matrix Ke
    Compute_Mat_Ke(MeshingData_C.N_Elements);

    Compute_Total_Stiff(MeshingData_C.N_Elements,MeshingData_C.N_dofs, MeshingData_C.Connectivity);
  }
  ~Elem_Stiff_Mat() {};
};

void Elem_Stiff_Mat::ComputeChristophel(InputParams InputData, MeshingData MeshingData_C)
{
  int Checks;
  double Curv, Tort, Xi, Eta;

  Checks=Stiff_Data.size();
  Curv=InputData.k;
  Tort=InputData.t;
  for (size_t it=0; it<MeshingData_C.Connectivity.size(); ++it)
  {
    // Store locally the variables needed for space gain
    Xi=MeshingData_C.IntegrNodeData[it].Xic;
    Eta=MeshingData_C.IntegrNodeData[it].Etac;

    // Christophel Symbol Components
    Stiff_Data[it].Chris.g131=-(Curv*Tort*Eta)/(1.0-Curv*Xi); //MeshingData_C.IntegrNodeData[it].Xic;
    Stiff_Data[it].Chris.g311=-(Curv*Tort*Eta)/(1.0-Curv*Xi);
    Stiff_Data[it].Chris.g132=(Curv*Tort*Xi)/(1.0-Curv*Xi)+Tort;
    Stiff_Data[it].Chris.g312=(Curv*Tort*Xi)/(1.0-Curv*Xi)+Tort;
    Stiff_Data[it].Chris.g133=-Curv/(1.0-Curv*Xi);
    Stiff_Data[it].Chris.g313=-Curv/(1.0-Curv*Xi);
    Stiff_Data[it].Chris.g231=-Tort;
    Stiff_Data[it].Chris.g321=-Tort;
    Stiff_Data[it].Chris.g331=(Curv*pow(Tort*Eta,2))/(1.0-Curv*Xi)+Curv*(1-Curv*Xi)-(pow(Tort,2))*Xi;
    Stiff_Data[it].Chris.g332=-(Curv*(pow(Tort,2))*Xi*Eta)/(1.0-Curv*Xi)-(pow(Tort,2))*Eta;
    Stiff_Data[it].Chris.g333=(Curv*Tort*Eta)/(1.0-Curv*Xi);
  }
  Checks=Checks+1;
}


void Elem_Stiff_Mat::Compute_Beta_M(InputParams InputData, MeshingData MeshingData_C)
{
  int dummy,  NodeN;
  double Jac_Det, Phi_1, Phi_2;
  Vector3d Sh_F(1/3.0, 1.0/3.0, 1.0/3.0), Sh_FX(-1.0,1.0,0), Sh_FY(-1.0,0.0,1.0);
  Matrix2d JAC, JACIN=MatrixXd::Zero(2,2);

  for (size_t Elem_Number=0; Elem_Number<MeshingData_C.Connectivity.size(); ++Elem_Number)
  {
    // Compute Jacobian
    //
    // J=[du/dl1 dy/dl1]
    //   [du/dl2  dy/dl2]
    JAC(0,0)=0; JAC(0,1)=0; JAC(1,0)=0; JAC(1,1)=0;
    for (size_t it=0; it<3; ++it)
    {
      NodeN=MeshingData_C.Connectivity[Elem_Number][it];
      JAC(0,0)=JAC(0,0)+Sh_FX(it)*MeshingData_C.coords[NodeN-1][0];
      JAC(0,1)=JAC(0,1)+Sh_FX(it)*MeshingData_C.coords[NodeN-1][1];
      JAC(1,0)=JAC(1,0)+Sh_FY(it)*MeshingData_C.coords[NodeN-1][0];
      JAC(1,1)=JAC(1,1)+Sh_FY(it)*MeshingData_C.coords[NodeN-1][1];
    }

    // Determinant of the Jacobian
    Jac_Det=JAC(0,0)*JAC(1,1)-JAC(0,1)*JAC(1,0);
    Stiff_Data[Elem_Number].Det_Jacob=Jac_Det;

    // Inverse of the Jacobian
    JACIN(0,0)=0; JACIN(0,1)=0; JACIN(1,0)=0; JACIN(1,1)=0;
    if (Jac_Det>0.0000000001)  // 1e-10
    {
      JACIN(0,0)= JAC(1,1)/Jac_Det;
      JACIN(1,1)= JAC(0,0)/Jac_Det;
      JACIN(0,1)=-JAC(0,1)/Jac_Det;
      JACIN(1,0)=-JAC(1,0)/Jac_Det;
    } else {
      cout << "The Jacobian is 0!\n" << endl;
      break;
    }

    // Compute Beta-Matrix
    Stiff_Data[Elem_Number].Beta_M.resize(6,9);
    Stiff_Data[Elem_Number].Beta_M << MatrixXd::Zero(6,9);
    for (size_t it=0; it<3; ++it)
    {
      Phi_1=Sh_FX(it)*JACIN(0,0)+Sh_FY(it)*JACIN(0,1);
      Phi_2=Sh_FX(it)*JACIN(1,0)+Sh_FY(it)*JACIN(1,1);
      dummy=it*3;
      // Substitute the data to the Beta matrix
      Stiff_Data[Elem_Number].Beta_M(0,dummy)  =Phi_1;
      Stiff_Data[Elem_Number].Beta_M(1,dummy+1)=Phi_2;
      Stiff_Data[Elem_Number].Beta_M(2,dummy)  =-Stiff_Data[Elem_Number].Chris.g331*Sh_F(it);
      Stiff_Data[Elem_Number].Beta_M(2,dummy+1)=-Stiff_Data[Elem_Number].Chris.g332*Sh_F(it);
      Stiff_Data[Elem_Number].Beta_M(2,dummy+2)=-Stiff_Data[Elem_Number].Chris.g333*Sh_F(it);
      Stiff_Data[Elem_Number].Beta_M(3,dummy)  =-Stiff_Data[Elem_Number].Chris.g231*Sh_F(it);
      Stiff_Data[Elem_Number].Beta_M(3,dummy+2)=0.5*Phi_2;
      Stiff_Data[Elem_Number].Beta_M(4,dummy)  =-Stiff_Data[Elem_Number].Chris.g131*Sh_F(it);
      Stiff_Data[Elem_Number].Beta_M(4,dummy+1)=-Stiff_Data[Elem_Number].Chris.g132*Sh_F(it);
      Stiff_Data[Elem_Number].Beta_M(4,dummy+2)=0.5*Phi_1-Stiff_Data[Elem_Number].Chris.g133*Sh_F(it);
      Stiff_Data[Elem_Number].Beta_M(5,dummy ) =0.5*Phi_2;
      Stiff_Data[Elem_Number].Beta_M(5,dummy+1)=0.5*Phi_1;
    }
  }
} // End of Function B

void Elem_Stiff_Mat::Compute_Tensor_C(InputParams InputData, MeshingData  MeshingData_C,
  double E1, double E2, VectorXi& Identities)
{
  double Factor_1, Factor_2;
  double Xi, Eta, Curv, Tort;
  double G, C11, C12, C13, C21, C22, C23, C31, C32, C33;

  // Compute stiffness factors
  for (size_t Elem_Number=0; Elem_Number<MeshingData_C.Connectivity.size(); ++Elem_Number)
  {
    double Ef=E1; 
    double Em=E2; 

    if (Identities[Elem_Number]>0.5)
    {
      Factor_1=(InputData.nu*Ef)/((1.0+InputData.nu)*(1.0-2.0*InputData.nu));
      Factor_2=(             Ef)/(2.0*(1.0+InputData.nu));
    } else {
      Factor_1=(InputData.nu*Em)/((1.0+InputData.nu)*(1.0-2.0*InputData.nu));
      Factor_2=(             Em)/(2.0*(1.0+InputData.nu));
    }

    // Compute the contravariant basis components of the element
    Curv=InputData.k;
    Tort=InputData.t;
    Xi=MeshingData_C.IntegrNodeData[Elem_Number].Xic;
    Eta=MeshingData_C.IntegrNodeData[Elem_Number].Etac;

    G  = pow(1.0-Curv*Xi,2);
    C11= (G+pow(Tort*Eta,2))/G;
    C12=-(pow(Tort,2)*Xi*Eta)/G;
    C13= ((Tort*Eta))/G;
    C21=-((pow(Tort,2)*Xi*Eta))/G;
    C22= (G+pow(Tort*Xi,2))/G;
    C23=-((Tort*Xi))/G;
    C31= ((Tort*Eta))/G;
    C32=-((Tort*Xi))/G;
    C33= (1.0)/G;

    // Preallocate all elements and enter the stiffness matrix values
    Stiff_Data[Elem_Number].C_Tensor.resize(6,6);
    Stiff_Data[Elem_Number].C_Tensor<<MatrixXd::Zero(6,6);
    // Compute all stiffness components
    Stiff_Data[Elem_Number].C_Tensor(0,0)=Factor_1*C11*C11+Factor_2*(C11*C11+C11*C11);
    Stiff_Data[Elem_Number].C_Tensor(0,1)=Factor_1*C11*C22+Factor_2*(C12*C12+C12*C12);
    Stiff_Data[Elem_Number].C_Tensor(0,2)=Factor_1*C11*C33+Factor_2*(C13*C13+C13*C13);
    Stiff_Data[Elem_Number].C_Tensor(0,3)=Factor_1*C11*C23+Factor_2*(C12*C13+C13*C12);
    Stiff_Data[Elem_Number].C_Tensor(0,4)=Factor_1*C11*C13+Factor_2*(C11*C13+C13*C11);
    Stiff_Data[Elem_Number].C_Tensor(0,5)=Factor_1*C11*C12+Factor_2*(C11*C12+C12*C11);

    Stiff_Data[Elem_Number].C_Tensor(1,0)=Stiff_Data[Elem_Number].C_Tensor(0,1);
    Stiff_Data[Elem_Number].C_Tensor(1,1)=Factor_1*C22*C22+Factor_2*(C22*C22+C22*C22);
    Stiff_Data[Elem_Number].C_Tensor(1,2)=Factor_1*C22*C33+Factor_2*(C23*C23+C23*C23);
    Stiff_Data[Elem_Number].C_Tensor(1,3)=Factor_1*C22*C23+Factor_2*(C22*C23+C23*C22);
    Stiff_Data[Elem_Number].C_Tensor(1,4)=Factor_1*C22*C13+Factor_2*(C21*C23+C23*C21);
    Stiff_Data[Elem_Number].C_Tensor(1,5)=Factor_1*C22*C12+Factor_2*(C21*C22+C22*C21);

    Stiff_Data[Elem_Number].C_Tensor(2,0)=Stiff_Data[Elem_Number].C_Tensor(0,2);
    Stiff_Data[Elem_Number].C_Tensor(2,1)=Stiff_Data[Elem_Number].C_Tensor(1,2);
    Stiff_Data[Elem_Number].C_Tensor(2,2)=Factor_1*C33*C33+Factor_2*(C33*C33+C33*C33);
    Stiff_Data[Elem_Number].C_Tensor(2,3)=Factor_1*C33*C23+Factor_2*(C32*C33+C33*C32);
    Stiff_Data[Elem_Number].C_Tensor(2,4)=Factor_1*C33*C13+Factor_2*(C31*C33+C33*C31);
    Stiff_Data[Elem_Number].C_Tensor(2,5)=Factor_1*C33*C12+Factor_2*(C31*C32+C32*C31);

    Stiff_Data[Elem_Number].C_Tensor(3,0)=Stiff_Data[Elem_Number].C_Tensor(0,3);
    Stiff_Data[Elem_Number].C_Tensor(3,1)=Stiff_Data[Elem_Number].C_Tensor(1,3);
    Stiff_Data[Elem_Number].C_Tensor(3,2)=Stiff_Data[Elem_Number].C_Tensor(2,3);
    Stiff_Data[Elem_Number].C_Tensor(3,3)=Factor_1*C23*C23+Factor_2*(C22*C33+C23*C32);
    Stiff_Data[Elem_Number].C_Tensor(3,4)=Factor_1*C23*C13+Factor_2*(C21*C33+C23*C31);
    Stiff_Data[Elem_Number].C_Tensor(3,5)=Factor_1*C23*C12+Factor_2*(C21*C32+C22*C31);

    Stiff_Data[Elem_Number].C_Tensor(4,0)=Stiff_Data[Elem_Number].C_Tensor(0,4);
    Stiff_Data[Elem_Number].C_Tensor(4,1)=Stiff_Data[Elem_Number].C_Tensor(1,4);
    Stiff_Data[Elem_Number].C_Tensor(4,2)=Stiff_Data[Elem_Number].C_Tensor(2,4);
    Stiff_Data[Elem_Number].C_Tensor(4,3)=Stiff_Data[Elem_Number].C_Tensor(3,4);
    Stiff_Data[Elem_Number].C_Tensor(4,4)=Factor_1*C13*C13+Factor_2*(C11*C33+C13*C31);
    Stiff_Data[Elem_Number].C_Tensor(4,5)=Factor_1*C13*C12+Factor_2*(C11*C32+C12*C31);

    Stiff_Data[Elem_Number].C_Tensor(5,0)=Stiff_Data[Elem_Number].C_Tensor(0,5);
    Stiff_Data[Elem_Number].C_Tensor(5,1)=Stiff_Data[Elem_Number].C_Tensor(1,5);
    Stiff_Data[Elem_Number].C_Tensor(5,2)=Stiff_Data[Elem_Number].C_Tensor(2,5);
    Stiff_Data[Elem_Number].C_Tensor(5,3)=Stiff_Data[Elem_Number].C_Tensor(3,5);
    Stiff_Data[Elem_Number].C_Tensor(5,4)=Stiff_Data[Elem_Number].C_Tensor(4,5);
    Stiff_Data[Elem_Number].C_Tensor(5,5)=Factor_1*C12*C12+Factor_2*(C11*C22+C12*C21);  
    }  // End of the loop over the elements
}  // End of Compute_Tensor_C

void Elem_Stiff_Mat::Compute_Mat_Ke(const int N_Elems)
{
  MatrixXd Ke_Trial=MatrixXd::Zero(9,9);

  // Preallocate the size of all stiffness element matrices
  for (int it=0; it<N_Elems; ++it)
  {
    Stiff_Data[it].Ke=MatrixXd::Zero(9,9);
  }

  // Loop over the elements and compute their stiffness matrix
  for (int it=0; it<N_Elems; ++it)
  {
    // Compute the matrix to be added
    Ke_Trial=Stiff_Data[it].Beta_M.transpose()*Stiff_Data[it].C_Tensor*Stiff_Data[it].Beta_M*Stiff_Data[it].Det_Jacob;

    // Single integration point Coordweight=1
    Stiff_Data[it].Ke=Ke_Trial;
  
  }
}

void Elem_Stiff_Mat::Compute_Total_Stiff(const int N_Elems, const int N_dofs, const Double2DVec Connectivity)
{
  int Local_Add, Local_ID, Row_Glo, Col_Glob;
  vector<int> MappingElem;

  // preallocate the total stiffness matrix
  K_total.resize(N_dofs,N_dofs);
  K_total<< MatrixXd::Zero(N_dofs,N_dofs);

  printf("sizeof(K_total)=%dx%d -> %.2f Mbytes\n", N_dofs, N_dofs, N_dofs*N_dofs*sizeof(double)/(1024.0*1024.0));

  MappingElem.resize(9);

  // Assemble the total stiffness
  for (int i=0; i<N_Elems; ++i)
  {
    // Compute mapping vector for each element
    for (int i3=0; i3<3; ++i3)
    {
      Local_Add=(i3)*3;
      Local_ID=(Connectivity[i][i3]-1)*3;
      for (int j3=0; j3<3; ++j3){
        MappingElem[Local_Add+j3]=Local_ID+j3;
      }
    }

    // loop over the elements of the stiffness matrix of each element and position it in the global
    for (int i_r=0; i_r<9; ++i_r)
    {
      Row_Glo=MappingElem[i_r];
      for (int i_c=0; i_c<9; ++i_c)
      {
        Col_Glob=MappingElem[i_c];
        K_total(Row_Glo,Col_Glob)=K_total(Row_Glo,Col_Glob)+Stiff_Data[i].Ke(i_r,i_c);
      }
    }
  }  // N_Elems
}

void Elem_Stiff_Mat::Print_Check_Function(const MatrixXd MatToPrint, const int NumberRows, const int NumberOfCols)
{
  int NumRowsUsed, NumColsUsed;

  // Check the size of the matrix and appoint loop sizes
  if (MatToPrint.size()>NumberRows)
  {
    NumRowsUsed=NumberRows;
  } else {
    NumRowsUsed=MatToPrint.size();
  }

  NumColsUsed=NumberOfCols;

  for (int it_R=0; it_R<NumRowsUsed; ++it_R)
  {
    for (int it_C=0; it_C<NumColsUsed; ++it_C)
    {
      cout << MatToPrint(it_R,it_C) << ",  ";
      if (it_C==(NumColsUsed-1))
          cout << ",  " << endl;
      }
  }
}

#endif
