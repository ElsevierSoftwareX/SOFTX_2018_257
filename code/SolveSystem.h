/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */

// Computes the solution of the K*u=F system in the local helix Serret-Frenet axis,
// after applying the boundary conditions. The boundary conditions are
// implemented so that there is no radial (Xn) displacement for all nodes
// on the local Etac axis (Xn=0) and there is no circumferential displacement
// for all nodes belonging to Xn axis (Etac=0)

#ifndef _SOLVESYSTEM_H_
#define _SOLVESYSTEM_H_ 1

#include "DefsMTUsed.h"

class SolveSystem
{
public:
  VectorXd Axial_Load_Sol;
  double Effective_Poisson_Ratio;
public:
  SolveSystem();
  void Apply_Boundary_Solve(Double2DVec& coords,MatrixXd & K_total, VectorXd & AxialLoad, int N_Slices, int N_Radians, double AppliedExternalStrain);

  // Constructor used with arguments Input and meshing data classes
  SolveSystem(InputParams& InputData, Double2DVec& coords, Double2DVec& Connectivity, MatrixXd& K_total, VectorXd& AxialLoad)
  {
    // Apply the boundary conditions and solve the system
    Apply_Boundary_Solve(coords, K_total, AxialLoad, InputData.MeshNSlices, InputData.MeshNRadians, InputData.AppliedExternalStrain);
  }
  ~SolveSystem() {};
};


void SolveSystem::Apply_Boundary_Solve(Double2DVec& coords, MatrixXd & K_total, VectorXd & AxialLoad, int N_Slices, int N_Radians, double AppliedExternalStrain)
{
  int column, N_Nodes;
  MatrixXd Total_Stiff_Transf;
  VectorXd AxialLoad_Transf;
  vector<int> ConstDofs;
  vector<double> D_EBC;

  //Constrain the central node of the helix in all directions
  ConstDofs.push_back(0); ConstDofs.push_back(1); ConstDofs.push_back(2);
  D_EBC.push_back(0.0); D_EBC.push_back(0.0);  D_EBC.push_back(0.0);
  N_Nodes=K_total.rows()/3;

  for (int iter=0; iter<N_Nodes; ++iter)
  {
    if (abs(coords[iter][1])<0.0001)
    {   // abs -> fabs
      ConstDofs.push_back((iter)*3+1);
      D_EBC.push_back(0.0);
    }
    if (abs(coords[iter][0])<0.0001)
    {   // abs -> fabs
      ConstDofs.push_back((iter)*3+0);
      D_EBC.push_back(0.0);
    }
  }

  // Apply the boundary condition matrix transformation
  for (unsigned int i_con=0; i_con<ConstDofs.size(); ++i_con)
  {
    column=ConstDofs[i_con ];
    for (unsigned int i_r=0; i_r<K_total.rows(); ++i_r)
    {
      AxialLoad(i_r)=AxialLoad(i_r)-D_EBC[i_con]*K_total(i_r,i_r);
      K_total(column,i_r)=0.0;
      K_total(i_r,column)=0.0;
    }
    K_total(column,column)=1.0;
    AxialLoad(column)=D_EBC[i_con];
  } // end of for loops

  // System Solution
  Axial_Load_Sol.resize(K_total.rows());

  MyTimer timer;
  timer.start();
  Axial_Load_Sol=K_total.llt().solve(AxialLoad);   //K_total.colPivHouseholderQr().solve(AxialLoad);
  timer.stop("Linear system solver");

  // Compute the mean radial displacements at the domain outer boundary
  int NodesToCheck=N_Slices*N_Radians;
  int Node_Of_Interest=0;
  vector<double> U_Radial_Boundary;

  // Loop over the boundary nodes
  for (int i=0; i<NodesToCheck; ++i)
  {
    Node_Of_Interest=N_Nodes-NodesToCheck+i;
    U_Radial_Boundary.push_back(sqrt(pow(Axial_Load_Sol[Node_Of_Interest*3],2)+pow(Axial_Load_Sol[Node_Of_Interest*3+1],2)));
  }
 
  double sum=0;
  for (size_t i=0; i<U_Radial_Boundary.size(); ++i){
    sum +=U_Radial_Boundary[i];
  }

  Effective_Poisson_Ratio=-sum/U_Radial_Boundary.size()/AppliedExternalStrain;
 
}// End of the system solution

#endif
