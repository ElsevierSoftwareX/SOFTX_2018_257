/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */

#include "DefsMTUsed.h"
#include "MyTimer.h"
#include "InputData.h"
#include "MeshingData.h"
#include "FibrilIdentifier.h"
#include "Elem_Stiff_Mat.h"
#include "Loading.h"
#include "SolveSystem.h"
#include "PostProcess.h"


int main(int argc, char *argv[])
{
  //Main Program Run
  //Read from file basic parameters
  char *input_param_file = (char *)"InputHelixData.txt";
  if (argc == 2)
    input_param_file = argv[1];

  printf(">> Input Parameters File: %s\n", input_param_file);

  InputParams InputData(input_param_file);
  MeshingData MeshingDataC(InputData);    // Uses Input in an Automated manner
 
  printf(">> Parameters:\n");
  printf(">> Angle=%.3f, Centerline=%.3f\n",InputData.theta, InputData.a);
  printf(">> Cf=%.3f\n", InputData.Cf);
  printf(">> Ef=%.3f, Em=%.3f\n", InputData.Ef, InputData.Em);
  printf(">> #Slices=%d, #Radians=%d\n", InputData.MeshNSlices, InputData.MeshNRadians);

  MyTimer timer;

  timer.start();
  VectorXi Identities=FibrilIdentifier(MeshingDataC, InputData.Cf);
  timer.stop("Loading Fibrils");

  // Stiffness Matrix Data
  timer.start();
  Elem_Stiff_Mat Stiff(InputData, MeshingDataC, Identities);
  timer.stop("Stiffness matrix data");

  // Normal loading data
  timer.start();
  Loading AxialLoading(InputData, MeshingDataC, Stiff);
  timer.stop("Normal loading data");

  // Solve system
  timer.start();
  SolveSystem Solution(InputData, MeshingDataC.coords, MeshingDataC.Connectivity, Stiff.K_total, AxialLoading.AxialLoad);
  timer.stop("Solve system");

  // Compute Forces and moments
  timer.start();
  PostProcessing Results(InputData, MeshingDataC.N_Elements,MeshingDataC.IntegrNodeData,MeshingDataC.Connectivity,Stiff.Stiff_Data,
                         AxialLoading.Elem_Ax_Strain,Solution.Axial_Load_Sol, Identities);
  timer.stop("Compute forces and moments");

  cout << "The effective fascicle modulus is: " << Results.E_eff << endl;
  cout << "The effective Poisson's ratio is: " << Solution.Effective_Poisson_Ratio << endl;
  cout << "The effective normalized moment is: " << Results.Mom_Eff << endl;
  

  return 0;
}
