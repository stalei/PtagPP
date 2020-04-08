// © Shahram Talei @ 2019
// // Reading stored data in a sage file
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
//
//
// //
#include "GlobalVars.h"
#include "proto.h"

#ifdef DoParallel
#include <mpi.h>
#endif


void PaintStars(int snap)
{
sprintf(CurrentFile,"PaintStars.c");
#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d Started painting snap %d! \n",ThisTask,snap);
#else
printf("I Started painting snap %d! \n",snap);
#endif
double Tf,Ti;
unsigned long int id;
int GalID=-1;
id=0; //particle id
fflush(stdout);
printf("test time:%g\n",AllStars[0].Time);
Tf=GetAge(FindTime(AllStars));
Ti=GetAge(FindTime(AllStarsPre));
printf("Time for snap:%d is %g & for snap:%d is %g.\nTime difference is: %g\n",snap,Tf,snap-1,Ti,Tf-Ti);
//

struct GalHashTable *gtable=GalEmptyTable(NumGalaxies);
if(gtable==NULL)
        {
        printf("can't allocate memory for hashtable!\n");
        EndRun(48,CurrentFile);
        }

//



int iter,iteration=4;
//for on all tagged particles
for(id=0;id<NumOfStars;id++)
{
iter=1;
do{
	GalID= LookupGalaxy(snap,id, SageOutput,NumGalaxies,2/iter);
        iter++;
  }while(GalID<0 && iter<iteration);
if(GalID>=0){
        GalInsertKey(gtable,GalID,&AllStars[id]);//insert tags in a galaxy hash table
	// write tag.galid here so I can use it later
	AllStars[id].GalNo=GalID;
        //CalculateStellarProperties(Ti,Tf, GalID,id);
        }
else
	printf("couldn't find associated galaxy for star %ld up to %gxRvir\n",id,(float)iter/2);
}
// now tags are in galaxy hash table, we can start painting them galaxy by galaxy
long int StCount=0;
double ECutoff;
for(id=0;id<NumOfStars;id++)
	{
	GalID=AllStars[id].GalNo;
	StCount=AllStars[id].Len;//CountStarsInGal(gtable,id);
	if(GP.f_mb<10)
		{
		AllStars[id].Len=(StCount*10.0)*(GP.f_mb/100);
		ECutoff=GalBndELimit(gtable,id,&AllStars,StCount,GP.f_mb);
		}
	else if(GP.f_mb==10)
		{
		//AllStars[id].Len=StCount;
		ECutoff=0;
		}
	else
		{
		printf("Invalid value for f_mb, we set default values (10 percent)\n");
                //AllStars[id].Len=StCount;
                ECutoff=0;
		}
	CalculateStellarProperties(Ti,Tf, GalID,id,ECutoff,StCount);
	}


#ifdef DoParallel
fflush(stdout);
if(ThisTask==0)
{
#endif
printf("finished calculation of stellar properties for all stars!\n");
PrintStar(0);
#ifdef DoParallel
}
#endif

return;
}
int LookupGalaxy(int snap, unsigned int id, struct SageGalaxies *Galaxy,int count,float ratio)
{
int i,j,galaxy=-1;
double rG,rS,rS2=0,rG2=0;
i=0;
for(j=0;j<3;j++)
{
	rS2+=AllStars[id].Pos[j]*AllStars[id].Pos[j];
}
rS=sqrt(rS2);
//printf("count:%d\n",count);
//printf("star %u is in %g\n",id,rS);
while(i<count && galaxy==-1)
   {
	rG=0;
	rG2=0;
	for(j=0;j<3;j++)
	{
           rG2+=Galaxy[i].Pos[j]*Galaxy[i].Pos[j];
	}
	rG=sqrt(rG2);
	//printf("galaxy %d is in %g\n",i,rG);
	if(abs(rG-rS)<(Galaxy[i].Rvir/ratio))
		galaxy=i;
	
//	if(AllStars[id].HaloIndex==SageOutput[i].GalaxyId)
//	{
//		galaxy=i;
//	}
	i++;
   }

return galaxy;
}
void CalculateStellarProperties(double ti,double tf, int galaxy, unsigned long int id,double BECut, long int count)
{
// I have to add more sophisticated calculations but let's start with simple method
// we have to devide this mass between all particles
// we also have to convert between units!
//if(AllStars[id].BindingEnergy < BECut)
//{
AllStars[id].StellarMass=1.0e9*(GetAge(tf)-GetAge(ti))*SageOutput[galaxy].Sfr/AllStars[id].Len;
//AllStars[id].GalNo=galaxy;//SageOutput[galaxy].
AllStars[id].TreeIndex=SageOutput[galaxy].TreeIndex;
AllStars[id].ZZ=SageOutput[galaxy].MetalsStellarMass/AllStars[id].Len;//lower than expected
AllStars[id].Mvir=SageOutput[galaxy].Mvir;
AllStars[id].Rvir=SageOutput[galaxy].Rvir;
AllStars[id].infallMvir=SageOutput[galaxy].infallMvir;
AllStars[id].Age=GetAge(AllStars[id].Time);//this makes sense
AllStars[id].LastMajorMerger=SageOutput[galaxy].LastMajorMerger;
//}
return;
}
double FindTime(struct tagged_particle *Stars)
{
return Stars[0].Time;
}
double GetAge(double a)
{
//high redshift approximation:
//double z=1/a-1;
//double H0=100*GP.HubbleParam;
//double t=2/(3*H0*sqrt(GP.Omega0+GP.OmegaLambda)*pow(1+z,3/2));

//or numerical integral
//get aMax from previous files -> GP.aMax
double age,ageGyr,t = 0;
int  n=1000,i;         // number of points in integrals
double WR=0,ai,adot;
double H0=GP.HubbleParam*100;
double Tyr = 977.8;    // coefficent for converting 1/H into Gyr
for(i=0;i<n;i++)
{
    ai = a*(i+0.5)/n;
    adot = sqrt((GP.Omega0/ai)+(WR/(ai*ai))+(GP.OmegaLambda*ai*ai));
    t+= 1./adot;
}

  age = a*t/n;
  ageGyr = (Tyr/H0)*age;

//t/=3.086e+19; //Mpc ->km
//t/=(3600*24*365.26); //s ->year
return ageGyr;
}
long int CountStarsInGal(struct GalHashTable *gtable,int id)
{
long int c=0;
struct GalLinkedList *glist;
glist=gtable->gtable[id];
while(glist->next)
{
        //printf("snap in IsTagged:%d\n",list->next->star->Snap);
        c++;
        glist=glist->next;
}

return c;
}

int subfind_compare_binding_energy(const void *a, const void *b)
{
  if(*((double *) a) > *((double *) b))
    return -1;

  if(*((double *) a) < *((double *) b))
    return +1;

  return 0;
}


double GalBndELimit(struct GalHashTable *gtable, int id, struct tagged_particle **Stars, long int count, double f_mb)
{
double BndELim=0;
long int c=0,NumLimit; // we tag 10 percent in sim but now we can use different fraction
NumLimit=(count*10.0)*(f_mb/100); //we tag 10 percent in CoSANG by default
// so there is a new limit for total number of tags for each galaxy
struct GalLinkedList *glist;
//sort binding energy & find energy cut-off
double *BndEnergy;
BndEnergy = (double *) malloc(count * sizeof(double));
glist=gtable->gtable[id];
while(glist->next)
{
        //printf("snap in IsTagged:%d\n",list->next->star->Snap);
	BndEnergy[c]=glist->next->star->BindingEnergy;
	c++;
        //
        glist=glist->next;
}
qsort(BndEnergy, count, sizeof(double), subfind_compare_binding_energy);
//now we can set the energy limit
BndELim=BndEnergy[NumLimit];
return BndELim;
}
