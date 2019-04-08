/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */
//  Reads the input geometric and material data and computes basic parameters
//  of the curvilinear helix domain

#ifndef _INPUTDATA_H_
#define _INPUTDATA_H_ 1

#include "DefsMTUsed.h"

#ifndef NSLICES
#define NSLICES 14
#endif

#ifndef NRADIANS
#define NRADIANS 15
#endif

/*
Angle = 60
Center =10
FiberContent = 40
Ef = 2000
Em = 0.25
Positions=40_r0102_C_40.dat
*/
class InputParams
{
public:
  int MeshNSlices, MeshNRadians, TotalNodesNum;
  double theta, a, Cf;
  double Ef, Em;
  double r, nu;
  double b, gam, k, t, AppliedExternalStrain;
public:
  void ReadInputDataFromFile(const char* InputFile);
  void SetParams();
  InputParams();

  // overload constructor
  InputParams(const char *InputFile)
  {
    // default parameter values
    theta=89.0; // 50-89
    r=1.0;      // const
    a=10;       // >r
    Ef=100;     // 200 -4000
    Em=0.1;     // 0.01-5
    nu=0.3;     // const
    Cf=50;      // range 35-85
    ReadInputDataFromFile(InputFile);
    SetParams();
  }
  ~InputParams() {};
};

void InputParams::SetParams()
{
  // Compute remaining
  b=a*tan(theta*M_PI/180);
  gam=sqrt(pow(a,2)+pow(b,2));
  k=a/pow(gam,2);
  t=b/pow(gam,2);
  AppliedExternalStrain=-0.01;
  MeshNSlices=NSLICES;   // Predefined mesh densities should not be modified
  MeshNRadians=NRADIANS;   // Predefined mesh densities should not be modified

  // Compute the total Nodes corresponding to the Mesh Specification
  TotalNodesNum=1;
  for (int i=1; i<=MeshNRadians; ++i)
  {
    TotalNodesNum=TotalNodesNum+i*MeshNSlices;
  }
}

void InputParams::ReadInputDataFromFile(const char* InputFile)
{
  // Define the input data file
  FILE *infile;
  char InpData[256];
  string DummyStr, ParamS, ValueS;
  int posit;

  // Check and put results in the output file
  infile=fopen(InputFile,"r");
  if (infile==NULL)
  {
    cout << "Input parameter file does not exist or could not be opened" <<  endl;
    exit(1);
  } else {
    while(fgets(InpData, 256, infile) != NULL)
    {
      // Get the number out of the string
      DummyStr=InpData;
      posit=DummyStr.find("=");
      ParamS=DummyStr.substr(0, posit);
      ValueS=DummyStr.substr(posit+1);

      ParamS.erase(std::remove(ParamS.begin(), ParamS.end(), ' '), ParamS.end());
      ValueS.erase(std::remove(ValueS.begin(), ValueS.end(), ' '), ValueS.end());

      if (ParamS.compare("Angle")==0)
      {
        theta=strtod(ValueS.c_str(), NULL);
        if ((theta < 50) || (theta > 89))
        {
          printf("Error in the input parameters: theta\n");
          exit(1);
        }
      }
      else if (ParamS.compare("Center")==0)
      {
        a=strtod(ValueS.c_str(), NULL);
        if (a < r)
        {
          printf("Error in the input parameters: r > a\n");
          exit(1);
        }
      }
      else if (ParamS.compare("FiberContent")==0)
      {
        Cf=strtod(ValueS.c_str(), NULL);
        if ((Cf < 35) || (Cf > 65))
        {
          printf("Error in the input parameters: Cf\n");
          exit(1);
        }
      }
      else if (ParamS.compare("Ef")==0)
      {
        Ef=strtod(ValueS.c_str(), NULL);
        if ((Ef < 200) || (Ef > 4000))
        {
          printf("Error in the input parameters: Ef\n");
          exit(1);
        }
      }
      else if (ParamS.compare("Em")==0)
      {
        Em=strtod(ValueS.c_str(), NULL);
        if ((Em < 0.01) || (Em > 5))
        {
          printf("Error in the input parameters: Em\n");
          exit(1);
        }
      }
    }
    fclose(infile);
  }
}

#endif
