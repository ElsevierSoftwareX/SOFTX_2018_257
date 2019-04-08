/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */
// Computes the mesh of the fascicle, namely the nodal coordinates to be used
// in the finite element discretization and their connectivity

#ifndef _MESHINGDATA_H_
#define _MESHINGDATA_H_ 1

#include "DefsMTUsed.h"
struct ElemIntegrCoord
{
    double Xic, Etac, Ac;
};

class MeshingData
{
public:
  int N_Nodes, N_dofs, N_Elements;
  Double2DVec coords, coords_Proj, Connectivity;
  vector<ElemIntegrCoord> IntegrNodeData;

public:
  MeshingData();
  // Define the necessary routines for the meshing
  void ComputeMesh(int N_Radians, int N_slices, int TotalNodesNum, double Hel_Radius);
  Double2DVec ComputeProjectMesh(Double2DVec coords, int TotalNumberNodes, double a, double r, double theta);
  int NumberOfElements(int N_Slices, int N_Radians);
  Double2DVec ComputeConnectiv(int N_Slices, int N_Radians, int MeshNElem);
  vector<ElemIntegrCoord> AreaCentersComp();

  // overload constructor
  MeshingData(InputParams InputParamData):N_Nodes(InputParamData.TotalNodesNum), N_dofs(3*InputParamData.TotalNodesNum),
    coords(InputParamData.TotalNodesNum, vector<double>(2,0)),coords_Proj(InputParamData.TotalNodesNum, vector<double>(2,0))
  {
    //Call different routines for the mesh and connectivity construction
    ComputeMesh(InputParamData.MeshNRadians, InputParamData.MeshNSlices, InputParamData.TotalNodesNum, InputParamData.r);
    N_Elements=NumberOfElements(InputParamData.MeshNSlices, InputParamData.MeshNRadians);
    Connectivity=ComputeConnectiv(InputParamData.MeshNSlices, InputParamData.MeshNRadians, N_Elements);
    IntegrNodeData=AreaCentersComp();
  }
  ~MeshingData() {};
};

void MeshingData::ComputeMesh(int N_Radians, int N_slices, int TotalNodesNum, double Hel_Radius)
{
  int counter;
  double phi, phi0, Lr;

  // Initialize and compute mesh coordinates
  counter=0;
  phi=2*M_PI/N_slices;
  coords[0][0]=0;
  coords[0][1]=0;
  Lr=Hel_Radius/N_Radians;
  for (int i=1; i<=N_Radians; ++i)
  {
    phi0=0;
    for (int j=1; j<=N_slices*i; ++j)
    {
      coords[counter+1][0]=i*Lr*cos(phi0/i);
      coords[counter+1][1]=i*Lr*sin(phi0/i);
      counter=counter+1;
      phi0=phi0+phi;
    }
  }
}

Double2DVec MeshingData::ComputeProjectMesh(Double2DVec coords, int TotalNumberNodes, double a, double r, double theta)
{
  Double2DVec Proj_Nodes(TotalNumberNodes, vector<double>(2));
  return Proj_Nodes;
}

int MeshingData::NumberOfElements(int N_Slices, int N_Radians)
{
  int NTrianglesPerSlice, TotalNumberElem;

  // iterate to compute the number of elements
  TotalNumberElem=0;
  for (int i=0; i<N_Radians; ++i)
  {
    // Compute the number of triangles per slide
    NTrianglesPerSlice=1+(i)*2;
    if (i==0)
    {
      for (int j=0; j<N_Slices; ++j)
      {
        TotalNumberElem=TotalNumberElem+1;
      }
    } else {
      for (int j=0; j<N_Slices; ++j)
      {
        for (int k=0; k<NTrianglesPerSlice; ++k)
        {
          TotalNumberElem=TotalNumberElem+1;
        }
      }
    }
  } // End of i

  return TotalNumberElem;
}

Double2DVec MeshingData::ComputeConnectiv(int N_Slices, int N_Radians, int MeshNElem)
{
  int NumberOfTriangle, in, out;
  int factor, innerstart, outerstart, NTrianglesPerSlice;
  Double2DVec MeshConnect(MeshNElem, vector<double>(3));

  // Paramater initialization and
  factor=0;
  NumberOfTriangle=0;
  for (int i=1; i<=N_Radians; ++i)
  {
    // Compute the number of triangles per slide
    NTrianglesPerSlice=1+(i-1)*2;
    if (i==1)
    {
      for (int j=1; j<=N_Slices; ++j)
      {
        MeshConnect[NumberOfTriangle][0]=i;
        MeshConnect[NumberOfTriangle][1]=i+j;
        MeshConnect[NumberOfTriangle][2]=i+j+1;
        if (j==(N_Slices))
        {
          MeshConnect[NumberOfTriangle][2]=(i-1)*N_Radians+2;
        }
        NumberOfTriangle=NumberOfTriangle+1;
      }
    } else {  // After the first radian
      innerstart=2+factor;
      factor=factor+N_Slices*(i-1);
      outerstart=2+factor;
      in=0;
      out=0;

      // Iterate over the slices
      for (int k=1; k<=N_Slices; ++k)
      {
        // Regular numbering per piece in each radian
        for (int piece=1; piece<=NTrianglesPerSlice; ++piece)
        {
          if (piece%2==1)
          {
            MeshConnect[NumberOfTriangle][0]=innerstart+in;
            MeshConnect[NumberOfTriangle][1]=outerstart+out;
            MeshConnect[NumberOfTriangle][2]=outerstart+out+1;
            out=out+1;
          } else {
            MeshConnect[NumberOfTriangle][0]=innerstart+in;
            MeshConnect[NumberOfTriangle][2]=innerstart+in+1;
            MeshConnect[NumberOfTriangle][1]=outerstart+out;
            in=in+1;
          }

          // Special case for the last slice
          if (k==(N_Slices))
          {
            if (piece==(NTrianglesPerSlice-1))
            {
              MeshConnect[NumberOfTriangle][0]=innerstart+in-1;
              MeshConnect[NumberOfTriangle][2]=innerstart+in-(i-1)*N_Slices;
              MeshConnect[NumberOfTriangle][1]=outerstart+out;
              in=in+1;
            }
            if (piece==(NTrianglesPerSlice))
            {
              MeshConnect[NumberOfTriangle][0]=innerstart+in-1-(i-1)*N_Slices;
              MeshConnect[NumberOfTriangle][1]=outerstart+out-1;
              MeshConnect[NumberOfTriangle][2]=outerstart+out-(i)*N_Slices;
              out=out+1;
            }
          } // End of special case correction
          NumberOfTriangle=NumberOfTriangle+1;
          if (NumberOfTriangle<20)
          {
          }
        }
      }
    }
  } // End of i...Nradians

  return MeshConnect;
}

vector<ElemIntegrCoord> MeshingData::AreaCentersComp()
{
  int Node_1, Node_2, Node_3;
  double x1, x2, x3, y1, y2, y3;
  vector<ElemIntegrCoord> Element_Data(Connectivity.size());

  // Compute element centers and Areas of each element
  for (size_t i=1; i<=Connectivity.size(); ++i)
  {
    // Locate Nodes
    Node_1=Connectivity[i-1][0]-1;
    Node_2=Connectivity[i-1][1]-1;
    Node_3=Connectivity[i-1][2]-1;

    // Extract coordinates
    x1=coords[Node_1][0];
    x2=coords[Node_2][0];
    x3=coords[Node_3][0];
    y1=coords[Node_1][1];
    y2=coords[Node_2][1];
    y3=coords[Node_3][1];
    Element_Data[i-1].Xic=(x1+x2+x3)/3;
    Element_Data[i-1].Etac=(y1+y2+y3)/3;
    Element_Data[i-1].Ac=(0.5)*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
  }

  return Element_Data;
}

#endif
