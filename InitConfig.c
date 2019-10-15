// © Shahram @ 2019
// Initial setup and configuration
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
#include<mpi.h>
#endif

//
#include "GlobalVars.h"
#include "proto.h"

void initialize(void)
{
sprintf(CurrentFile,"InitConfig.c");
int outputfolder_status;
struct stat st={0};

#ifdef DoParallel
if(ThisTask==0){
#endif
printf("\n Initial configuration.\n");
ReadParameters(ParametersFile);
// create output folder
if (stat(GP.OutputDir, &st) == -1) 
	outputfolder_status=mkdir(GP.OutputDir,02755);
if(outputfolder_status<0)
	{
	printf("Can't create the output folder!\n");
	EndRun(32,CurrentFile);
	}
#ifdef DoParallel
}
#endif

SetGP(); //set default value for general parameters!
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("Finished the initial setup.\n");
return;
}
void ReadParameters(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

//#ifdef DoParallel
//if(ThisTask==0)
//#endif
printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.\n");

int id[MAXTAGS];
void *addr[MAXTAGS];
char tag[MAXTAGS][50];

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int errorFlag = 0; // set not used

//#ifdef DoParallel
//if(ThisTask==0){
//#endif
nt=0;

strcpy(tag[nt],"SageDir");
addr[nt]=GP.SageDir;
id[nt++]=STRING;

strcpy(tag[nt],"TagDir");
addr[nt]=GP.TagDir;
id[nt++]=STRING;

strcpy(tag[nt],"SnapDir");
addr[nt]=GP.SnapDir;
id[nt++]=STRING;


strcpy(tag[nt],"OutputDir");
addr[nt]=GP.OutputDir;
id[nt++]=STRING;

strcpy(tag[nt],"FirstSnap");
addr[nt]=&GP.FirstSnap;
id[nt++]=INT;

strcpy(tag[nt],"LastSnap");
addr[nt]=&GP.LastSnap;
id[nt++]=INT;

strcpy(tag[nt],"BufferSize");
addr[nt]=&GP.BufferSize;
id[nt++]=REAL;

strcpy(tag[nt],"TimeLimit");
addr[nt]=&GP.TimeLimit;
id[nt++]=REAL;

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	  #ifdef DoParallel
	   if(ThisTask==0)
	  #endif
	      printf("Obtaining parameters from file '%s':\n", fname);
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;
		  if(buf1[0] == '%')
		    continue;
		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }
		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case REAL:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  fprintf(stdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy((char *) addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  fprintf(stdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  fprintf(stdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
		    }
		}
	      fclose(fd);
	      fclose(fdout);
	      printf("\n");
	      i = strlen(GP.SageDir);
	      if(i > 0)
		if(GP.SageDir[i - 1] != '/')
		  strcat(GP.SageDir, "/");
	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", GP.SageDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
	    }
	printf("Finished reading parameters.\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.\n");
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	EndRun(136,CurrentFile);
	}

if(errorFlag==1) 
	printf("Error in reading parameter.\n");
//#ifdef DoParallel
//MPI_Bcast(&GP, sizeof(struct GlobalParameters), MPI_BYTE, 0, MPI_COMM_WORLD);
//}
//#endif

return;
}
void SetGP()
{
GP.TotNumStarsAllSnaps=0;
GP.TotNumTagsAllSnaps=0;
GP.CoolingOn=0;
GP.StarformationOn=0;

#ifdef DoParallel
MPI_Bcast(&GP, sizeof(struct GlobalParameters), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

return;
}
