// © Shahram Talei @ 2019
// Reading stored data in a sage file
//
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <ctype.h>

#ifdef DoParallel
#include <mpi.h>
#endif
//
#include "GlobalVars.h"
#include "proto.h"
/////// Configuration
//#define TestConfig

void ReadSage(int snap)
{
sprintf(CurrentFile,"ReadSage.c");
int i;
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("Extracting information from the sage file.\n");
//1
LoadSageFiles(snap);
printf("Sage file path is loaded for snap %d.\n",snap);
fflush(stdout);
//2
NumGalaxies = ReadSageHeader(SageFilesCount,SageFilesPath);
NumGalaxiesPre = ReadSageHeader(SageFilesCountPre,SageFilesPathPre);

printf("Recovered %d galaxies for snapshot %d. and ",NumGalaxies,snap);
printf("Recovered %d galaxies for snapshot %d.\n",NumGalaxiesPre,snap-1);


if((SageOutput = (struct SageGalaxies*)malloc(NumGalaxies * sizeof(struct SageGalaxies))) == NULL)
    {
        printf("Failed to allocate for target sage...");
        return;
    }
if((SageOutputPre = (struct SageGalaxies*)malloc(NumGalaxiesPre * sizeof(struct SageGalaxies))) == NULL)
    {
        printf("Failed to allocate for previous sage...");
        return;
    }


ReadSageModel(SageFilesCount, SageFilesPath, SageOutput);
ReadSageModel(SageFilesCountPre, SageFilesPathPre, SageOutputPre);
printf("Galaxy info loaded for snapshot:%d\n",snap);

//#ifdef DoParallel
//if(ThisTask==0)
//#endif
//if(snap==GP.FirstSnap)
//{
printf("Sample galaxy for snapshot:%d in processor %d\n",snap,ThisTask);
for(i=0;i<NumGalaxies;i++)
   PrintGalaxyInfo(SageOutput,i);
//}
//printf("Sample galaxy for snapshot:%d\n",snap);
//PrintGalaxyInfo(SageOutput,10);



//printf("Sample galaxy for snapshot:%d\n",snap-1);
//PrintGalaxyInfo(SageOutputPre,0);


/*
   mal_gal(totmal);
   printf("Reading all Sage files\n");
   read_sage();
   malallgal(totmal);
   read_trees(totmal);
   fill_gal(totmal);
   free(treepath);
*/
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
printf("Reading the sage file for snap %d,in process %d is done.\n",snap,ThisTask);
//#ifdef DoParallel
//}
//#endif

return;
}

/////// 1
void LoadSageFiles(int snap)
{
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
 //printf("Reading the sage file.\n");
    //int numfiles;
	//int NumSfiles; // change to global later!
    SageFilesCount = CountSageFiles(snap);
    SageFilesCountPre=CountSageFiles(snap-1);
    //NumSfiles = numfiles;
//   #ifdef DoParallel
//  if(ThisTask==0)
//  #endif
  //  printf("I found:\n %d file(s) for snapshot %d.\n %d file(s) for snapshot %d\n", SageFilesCount,snap,SageFilesCountPre,snap-1);
    SageFilesPath = (struct Path_Names*)malloc(SageFilesCount * sizeof(struct Path_Names));
    SageFilesPathPre = (struct Path_Names*)malloc(SageFilesCountPre * sizeof(struct Path_Names));

	  //  {
        //printf("Failed to allocate for sage...");
        //return;
 	    //}

    ReadSageFNames(snap,SageFilesPath);
    ReadSageFNames(snap-1,SageFilesPathPre);
//#ifdef DoParallel
//if(ThisTask==0){
//#endif
//#ifdef DoParallel
//if(ThisTask==0){
//#endif
   //printf("first path for snap %d is:%s\n",snap,SageFilesPath[0].paths);
   //printf("first path for snap %d is:%s\n",snap-1,SageFilesPathPre[0].paths);

//#ifdef DoParallel
//}
//#endif


    return;
}
    

void ReadSageFNames(int snap,struct Path_Names *SageFile)
{
	//printf("Reading sage file-name(s) for snap:%d\n",snap);
    int i=0;
    char line[10000];
    char file1[10000];
    char ABSpath[10000];//2
    char model[500];
    FILE *fd; 
    sprintf(file1, "%s/sagefile_%03d.txt", GP.SageDir, snap);     
    fd = fopen(file1, "r");
    while(EOF != fscanf(fd, "%s%*[^\n]", line))
    {
       strcpy(ABSpath,line);//2
       strcpy(model,strrchr(ABSpath,'/'));
       //strcpy(SageFilesPath[i].paths, line);//1
       if (model[0] == '/') 
         memmove(model, model+1, strlen(model));
       sprintf(SageFile[i].paths,"%s/%s",GP.SageDir,model);
       i++;
    }
   fclose(fd);
   return; 
}

int CountSageFiles(int snap)
{
    //int i=0;
    int n=0;
    char line[10000];
    char file1[10000];
    FILE *fd; 
    sprintf(file1, "%s/sagefile_%03d.txt",GP.SageDir, snap);     
    fd = fopen(file1, "r");
    while(fgets(line, sizeof(line), fd)!=NULL) 
    {
       n++;// 1;
    }
    fclose(fd);
    return n;
}
//////// 2

int ReadSageHeader(int FileCount,struct Path_Names *SageFile)
{

    int Ntrees;
    int NtotGal;
    int totmal; 
    int i;
    FILE *fd;
    char file1[1000];
    totmal = 0;
    for(i=0; i<FileCount; i++)
    {
         sprintf(file1, "%s", SageFile[i].paths);
         fd = fopen(file1, "rb");
         if(NULL == fd)
         {
            printf("Cannot open sage file");
            return(-1);
         }
         fread(&Ntrees, sizeof(int), 1, fd);
         fread(&NtotGal, sizeof(int), 1, fd);
         totmal = totmal + NtotGal;
         fclose(fd);
    }
  return totmal;
}
void ReadSageModel(int FileCount,struct Path_Names *SageFile,struct SageGalaxies *Output)
{  
   int Ntrees;
   int NtotGal;    
   int *galpertree;
   char file1[1000];
   int offset; 
   int i;
   FILE *fd;
   offset = 0;  
   for(i=0; i<FileCount; i++)
   {
      sprintf(file1, "%s", SageFile[i].paths);
      //printf("Loading galaxies for model = %s\n", file1);
      fd = fopen(file1, "rb");
      if(NULL == fd)
      {
         printf("Cannot open sage for reading");
         return;
      }                                                                      
      fread(&Ntrees, sizeof(int), 1, fd);
      fread(&NtotGal, sizeof(int), 1, fd);
      galpertree = (int*)malloc(Ntrees*sizeof(int));   
      fread(&galpertree[0], sizeof(int), Ntrees, fd);
      fread(&Output[offset], sizeof(struct SageGalaxies), NtotGal, fd);        
      fclose(fd);                                                                           
      offset = offset + NtotGal;                                                                            }
      free(galpertree);  
      return;
}   

/*
void malallgal(int totmal)
{   
    if((AllGal = (struct All_Gal_info*)malloc(totmal * sizeof(struct All_Gal_info))) == NULL)
    {
        printf("Failed to allocate for AllGal...");
        return;
    }
  return;
}

void fill_gal(int totmal)
{
    int i;
    int a;

    for(i=0; i<totmal; i++)
    {
       AllGal[i].Mvir = SageOutput[i].Mvir;

       for(a=0; a<3; a++)
       {
            AllGal[i].Pos[a] = SageOutput[i].Pos[a];
            AllGal[i].Vel[a] = SageOutput[i].Vel[a] * All.Time;
       }
   
       AllGal[i].Rvir = SageOutput[i].Rvir;
       
    }
  free(SageOutput);
  return;
}
*/
void PrintGalaxyInfo(struct SageGalaxies *Output,int index)
{
int i;
i=index;
//#ifdef DoParallel
//if(ThisTask==0)
//{
//#endif
printf("~~~~~~~~~~~~~~~~~~~~~~∮∮∮∮∮∮∮∮∮∮∮∮~~~~~~~~~~~~~~~~~~~~~~\n");
printf("Galaxy Information- galaxy:%d, Type:%d, FileNr:%d\n",i,Output[i].Type,Output[i].FileNr);
printf("GalIndex:%lld,HaloIndex:%d,TreeIndex:%d,CentralGal:%d,CentralMvir:%f\n",Output[i].GalaxyIndex,Output[i].HaloIndex,Output[i].TreeIndex,Output[i].CentralGal,Output[i].CentralMvir);

printf("Pos(x,y,z):(%f,%f,%f)\nVel(vx,vy,vz):(%f,%f,%f)\nSpin:(%f,%f,%f),Len:%d\n",Output[i].Pos[0],Output[i].Pos[1],Output[i].Pos[2],Output[i].Vel[0],Output[i].Vel[1],Output[i].Vel[2],Output[i].Spin[0],Output[i].Spin[1],Output[i].Spin[2],Output[i].Len);

printf("Mvir:%f,Rvir:%f,Vvir:%f,Vmax:%f,VelDisp:%f\n",Output[i].Mvir,Output[i].Rvir,Output[i].Vvir,Output[i].Vmax,Output[i].VelDisp);

printf("ColdGas:%f,StellarMass:%f,BulgeMass:%f,HotGas:%f\nEjectedMass:%f,BHMass:%f,ICS:%f\n",Output[i].ColdGas,Output[i].StellarMass,Output[i].BulgeMass,Output[i].HotGas,Output[i].EjectedMass,Output[i].BlackHoleMass,Output[i].ICS);

printf("\nMetals:\nColdGas:%f,StellarMass:%f,BulgeMass:%f,HotGas:%f\nEjectedMass:%f,ICS:%f\n\n",Output[i].MetalsColdGas,Output[i].MetalsStellarMass,Output[i].MetalsBulgeMass,Output[i].MetalsHotGas,Output[i].MetalsEjectedMass,Output[i].MetalsICS);


printf("Sfr:%f,SfrBulge:%f,SfrICS:%f\nRd:%f,Cooling:%f,Heating:%f\nLastMajorMerger:%f,OutflowRate:%f\n",Output[i].Sfr,Output[i].SfrBulge,Output[i].SfrICS,Output[i].DiskScaleRadius,Output[i].Cooling,Output[i].Heating,Output[i].LastMajorMerger,Output[i].OutflowRate);

printf("infallMvir:%f,infallVvir:%f,infallVmax:%f,r_heat:%f\n",Output[i].infallMvir,Output[i].infallVvir,Output[i].infallVmax,Output[i].r_heat);
//printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
//printf("∮∮∮∮∮∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮ ∮∮∮∮∮∮\n");
printf("~~~~~~~~~~~~~~~~~~~~~~∮∮∮∮∮∮∮∮∮∮∮∮~~~~~~~~~~~~~~~~~~~~~~\n");

//#ifdef DoParallel
//}
//#endif
return;

}
