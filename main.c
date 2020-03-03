// This is the main function to do the post-processing on tagged particles
// © Shahram @ 2019

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "GlobalVars.h"
#include "proto.h"

#ifdef DoParallel
#include <mpi.h>
#endif

//#define DOLD_HDF5

int main(int argc, char *argv[])//char **argv)// *argv[])
{

sprintf(CurrentFile,"main.c");
int si=0,sf=0,step=0;

#ifdef DoParallel
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
MPI_Comm_size(MPI_COMM_WORLD, &NTask);
MPI_Barrier(MPI_COMM_WORLD);

 //MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &MPI_CommLocal);
 // MPI_Comm_rank(MPI_CommLocal, &ThisTask);
 // MPI_Comm_size(MPI_CommLocal, &NTask);

if(ThisTask==0)
{
#endif

printf("\n\n∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮☆ ∮\n");
printf("\ttagPP v1.0\n\t ©2019 Shahram ☻  @UA\n...\n");//%1.1f\n",Version);
printf("\ntagPP is for post-processing tag files!\n\n");

//printf("argc:%d\n",argc);
if(argc <2)
{
	printf("Prameters are missing\n");
	EndRun(27,CurrentFile);
}
fflush(stdout);

#ifdef DoParallel
}
#endif
strcpy(ParametersFile, argv[1]);
initialize();

#ifdef DoParallel
MPI_Barrier(MPI_COMM_WORLD);
// check # of cores and snaps ... snaps># of cores and snaps%#of cores!=0
si=GP.FirstSnap+ThisTask;
sf=GP.LastSnap;
step=NTask;
#else
si=GP.FirstSnap;
sf=GP.LastSnap;
step=1;
#endif

//printf("parameters: SageDir:%s, FirstSnap:%d\n",GP.SageDir,GP.FirstSnap);
//have to do on all snaps
//
for(TargetSnap=si; TargetSnap<=sf;TargetSnap+=step)
{
	//TargetSnap=GP.FirstSnap;
	ReadSage(TargetSnap);
	//#ifdef DoParallel
	//MPI_Barrier(MPI_COMM_WORLD);
	//#endif

	ReadTag(TargetSnap);
	//#ifdef DoParallel
	//fflush(stdout);
	//MPI_Barrier(MPI_COMM_WORLD);
	//#endif

	PaintStars(TargetSnap);
  //      #ifdef DoParallel
        //fflush(stdout);
        //MPI_Barrier(MPI_COMM_WORLD);
    //    #endif
	WriteTag(TargetSnap);
        //#ifdef DoParallel
        //fflush(stdout);
        //MPI_Barrier(MPI_COMM_WORLD);
        //#endif
	//free(AllStars);
	//free(AllStarsPre);
}

fflush(stdout);

//after for
//just write one tagged snap

//#ifdef DoParallel
//MPI_Barrier(MPI_COMM_WORLD);
//#endif
HashTags(GP.LastSnap);


//#ifdef DoParallel
//MPI_Barrier(MPI_COMM_WORLD);
//#endif
//WriteTaggedSnap(GP.LastSnap);

//free(P);
//free(StellarHaloAllSnaps);
//DeleteTable(table);

#ifdef DoParallel
MPI_Barrier(MPI_COMM_WORLD);
fflush(stdout);
if(ThisTask==0)
#endif
printf("end of the run!\n");

#ifdef DoParallel
MPI_Finalize();
#endif


return 0;
}



