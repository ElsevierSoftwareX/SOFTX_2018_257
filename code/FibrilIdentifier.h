/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */
// Loads the fiber position for the fascicle geometry depending on the input
// data so that a certain fiber content is created. The content is computed
// In increments of 5%
#ifndef _FIBRILIDENTIFIER_H_
#define _FIBRILIDENTIFIER_H_ 1

double Cf_r_map[7][2]={ {35,0.0935},{40, 0.102},
                  {45, 0.1060},
                  {50, 0.1125},
                  {55, 0.11726},  
                  {60, 0.123},  
                  {65, 0.12747}}; //

double *LoadFibrils(char *filename)
{
  FILE *fp=fopen(filename, "r");
  if (fp==NULL)
  {
    printf("File %s not found. Exiting!\n", filename);
    exit(1);
  }

  double *FibrilCenterlines=(double *)malloc(40*2*sizeof(double));
  for (int i=0; i<40; i++)
    if (fscanf(fp, "%lf %lf", &FibrilCenterlines[i*2+0], &FibrilCenterlines[i*2+1]) == EOF)
    {
      printf("Error in reading data from %s. Exiting!", filename);
      exit(1);
    }

  fclose(fp);

  return FibrilCenterlines;
}

VectorXi FibrilIdentifier(MeshingData &MeshingDataC, double Cf)
{
  // The input number of fibers is 40
  int l_FibrilCenterLines=40;

  // round up Cf to multiple of 5
  double multiple = 5;
  double Cf_rounded = ((int)((Cf + multiple/2) / multiple)) * multiple;
  cout<< "The rounded content is: " << Cf_rounded << endl;
  double rf = 0;

  // get rf value from Cf
  for (int i=0; i<7; i++)
    if (Cf_rounded == Cf_r_map[i][0])
    {
      rf = Cf_r_map[i][1];
      break;
    }

  if (rf == 0)
  {
     printf("Fatal error in computing radius from Cf. Exiting...\n");
     exit(1); 
  }

  // get fiber positions filename from number of fibers (40) and rounded Cf value
  int Cf_int = (int) Cf_rounded;

  char filename[128];  
  sprintf(filename, "../data/%d_C_%d.dat", l_FibrilCenterLines, Cf_int);  //

  double *FibrilCenterlines;
  FibrilCenterlines = LoadFibrils(filename);

  int l_Xnc=MeshingDataC.Connectivity.size();
  double Xnc[l_Xnc];
  double Xbc[l_Xnc];

  for (int i=0; i<l_Xnc; i++)
    Xnc[i]=MeshingDataC.IntegrNodeData[i].Xic;

  for (int i=0; i<l_Xnc; i++)
    Xbc[i]=MeshingDataC.IntegrNodeData[i].Etac;

  double Xbounds[l_FibrilCenterLines*2], Ybounds[l_FibrilCenterLines*2];  // 40*2

  VectorXi Identities(l_Xnc);
  for (int i=0; i<l_Xnc; i++) Identities[i]=0;

  // Run over the centerlines of the fibrils
  for (int i=0; i<l_FibrilCenterLines; i++)
  {
    // Form margins for each of the fibrils
    Xbounds[i*2+0]=FibrilCenterlines[i*2+0]-rf;
    Xbounds[i*2+1]=FibrilCenterlines[i*2+0]+rf;

    Ybounds[i*2+0]=FibrilCenterlines[i*2+1]-rf;
    Ybounds[i*2+1]=FibrilCenterlines[i*2+1]+rf;

    // Run over the elements to identify
    for (int j=0; j<l_Xnc; j++)
    {
      double rpoint=sqrt(pow(Xnc[j]-FibrilCenterlines[i*2+0],2)+pow(Xbc[j]-FibrilCenterlines[i*2+1],2));
      if ((Xnc[j]>Xbounds[i*2+0]) && (Xnc[j]<Xbounds[i*2+1]) && (Xbc[j]>Ybounds[i*2+0]) && (Xbc[j]<Ybounds[i*2+1]) && (rpoint<rf))
        Identities[j]=1;
    }
  }

  free(FibrilCenterlines);
  return Identities;
}


#endif
