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
int GalID=-1,flag;
int *hid,*subhid;
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
	fflush(stdout);
        EndRun(48,CurrentFile);
        }

//



int iter,iteration=4;
//for on all tagged particles
long int *StCount;
StCount=(long int *) malloc(NumGalaxies * sizeof(long int));
double *ECutoff;
ECutoff= (double *) malloc(NumGalaxies * sizeof(double));
printf("got to ECutoff\n");
for(GalID=0;GalID<NumGalaxies;GalID++)
{
	ECutoff[GalID]=0;
	StCount[GalID]=0;
}
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
	StCount[GalID]++;
        //CalculateStellarProperties(Ti,Tf, GalID,id);
        }
else
	printf("couldn't find associated galaxy for star %ld up to %gxRvir\n",id,(float)iter/2);
}
// now tags are in galaxy hash table, we can start painting them galaxy by galaxy
//long int StCount=0;
//double *ECutoff;
//Ecutoff= (double *) malloc(NumGalaxies * sizeof(double));
for(GalID=0;GalID<NumGalaxies;GalID++)
	{
	//GalID=AllStars[id].GalNo;
	//StCount=AllStars[id].Len;//CountStarsInGal(gtable,id);
	if(GP.f_mb<10)
		{
		//StCount=(StCount*10.0)*(GP.f_mb/100);
		//printf("Before ECutoff\n");
		//fflush(stdout);
		ECutoff[GalID]=GalBndELimit(gtable,GalID,StCount[GalID],GP.f_mb);
                //printf("after ECutoff\n");
                //fflush(stdout);
		//StCount=(StCount*10.0)*(GP.f_mb/100);
		}
	else if(GP.f_mb==10)
		{
		//AllStars[id].Len=StCount;
		//printf("Before ECutoff\n");
                //fflush(stdout);
		ECutoff[GalID]=0;
		//printf("After ECutoff\n");
                //fflush(stdout);
		}
	else
		{
		printf("Invalid value for f_mb, we set default values (10 percent)\n");
                //AllStars[id].Len=StCount;
                ECutoff[GalID]=0;
		}
	//CalculateStellarProperties(Ti,Tf, GalID,id,ECutoff,StCount);
	}
//////////////////////////
if(GP.SubhaloSelectionOn !=0)
{
//int *hid,*subhid;
hid=malloc(sizeof(int));
*hid=-1;
subhid=malloc(sizeof(int));
*subhid=-1;
ReadTargetHalo(GP.SubhaloFile,snap,hid,subhid);
flag=(GP.SubhaloSelectionOn==1)?0:1;
}
////////////////////
for(id=0;id<NumOfStars;id++)
        {
	GalID=AllStars[id].GalNo;
	//printf("%ld after GalID\n",id);
	if(GP.f_mb<10)
		//StCount[GalID]=AllStars[id].Len;
		StCount[GalID]=(StCount[GalID]*10.0)*(GP.f_mb/100);
	//printf("%ld after StCount\n",id);
        //ECutoff[GalID]=GalBndELimit(gtable,GalID,StCount,GP.f_mb);
        if(GP.SubhaloSelectionOn !=0) 
	  CalculateStellarPropertiesSubSelection(Ti,Tf, GalID,id,ECutoff[GalID],StCount[GalID],*hid,*subhid,flag);
	else 
         CalculateStellarProperties(Ti,Tf, GalID,id,ECutoff[GalID],StCount[GalID]);
        }
/////////////////////////
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
long int Len;
Len=count;//AllStars[id].Len;
//if(G.f_mb<10)
if(AllStars[id].BindingEnergy <= BECut && Len>0)
{
printf("SMass:%g,SMassPre:%g\n",SageOutput[galaxy].StellarMass, SageOutputPre[galaxy].StellarMass);
//if(Len>0)
//	{
//printf("I got inside painting function.\n");
//fflush(stdout);
//AllStars[id].StellarMass=1.0e10*SageOutput[galaxy].StellarMass/Len;
//AllStars[id].StellarMass=1.0e10*((SageOutput[galaxy].StellarMass-SageOutputPre[galaxy].StellarMass)/Len); //1.0e9*(GetAge(tf)-GetAge(ti))*SageOutput[galaxy].Sfr/Len;
AllStars[id].StellarMass=1.0e10*((SageOutput[galaxy].Stars)/Len);
//AllStars[id].GalNo=galaxy;//SageOutput[galaxy].
AllStars[id].TreeIndex=SageOutput[galaxy].TreeIndex;
AllStars[id].ZZ=SageOutput[galaxy].MetalsStellarMass/SageOutput[galaxy].StellarMass;//Len;//lower than expected
AllStars[id].Mvir=SageOutput[galaxy].Mvir;
AllStars[id].Rvir=SageOutput[galaxy].Rvir;
AllStars[id].infallMvir=SageOutput[galaxy].infallMvir;
AllStars[id].Age=GetAge(AllStars[id].Time);//this makes sense
AllStars[id].LastMajorMerger=SageOutput[galaxy].LastMajorMerger;
AllStars[id].HaloIndex=SageOutput[galaxy].HaloIndex;
AllStars[id].SubhaloIndex=SageOutput[galaxy].FOFHaloIndex;
AllStars[id].GalIndex=SageOutput[galaxy].CentralGal;
}
/*else if(GP.f_mb==10 && Len>0)
//printf("BE:%g,BECut:%g\n",AllStars[id].BindingEnergy, BECut);
//if(Len>0)
        {
//printf("I got inside painting function.\n");
//fflush(stdout);
//AllStars[id].StellarMass=1.0e10*SageOutput[galaxy].StellarMass/Len;
AllStars[id].StellarMass=1.0e10*((SageOutput[galaxy].StellarMass-SageOutputPre[galaxy].StellarMass)/Len); //1.0e9*(GetAge(tf)-GetAge(ti))*SageOutput[galaxy].Sfr/Len;
//AllStars[id].GalNo=galaxy;//SageOutput[galaxy].
AllStars[id].TreeIndex=SageOutput[galaxy].TreeIndex;
AllStars[id].ZZ=SageOutput[galaxy].MetalsStellarMass/SageOutput[galaxy].StellarMass;//Len;//lower than expected
AllStars[id].Mvir=SageOutput[galaxy].Mvir;
AllStars[id].Rvir=SageOutput[galaxy].Rvir;
AllStars[id].infallMvir=SageOutput[galaxy].infallMvir;
AllStars[id].Age=GetAge(AllStars[id].Time);//this makes sense
AllStars[id].LastMajorMerger=SageOutput[galaxy].LastMajorMerger;
}
*/

else // if their energy is above the limit, unbind them
{
AllStars[id].BindingEnergy=0;
AllStars[id].StellarMass=1;
AllStars[id].ZZ=0;
}

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


double GalBndELimit(struct GalHashTable *gtable, int id, long int count, double f_mb)
{
double BndELim=0;
long int c=0,NumLimit; // we tag 10 percent in sim but now we can use different fraction
NumLimit=(count*10.0)*(f_mb/100); //we tag 10 percent in CoSANG by default
// so there is a new limit for total number of tags for each galaxy
struct GalLinkedList *glist2;
//sort binding energy & find energy cut-off
double *BndEnergy;
BndEnergy = (double *) malloc(count * sizeof(double));
glist2=gtable->gtable[id];
while(glist2->next)
{
        //printf("snap in IsTagged:%d\n",list->next->star->Snap);
	BndEnergy[c]=glist2->next->star->BindingEnergy;
	c++;
	printf("StCount:%ld,while count:%ld\n",count,c);
        //
        glist2=glist2->next;
}
qsort(BndEnergy, count, sizeof(double), subfind_compare_binding_energy);
//now we can set the energy limit
BndELim=BndEnergy[NumLimit];
return BndELim;
}

void CalculateStellarPropertiesSubSelection(double ti,double tf, int galaxy, unsigned long int id,double BECut, long int count,int hid,int subhid,int flag)
{
//if(flag) //exclude
//SageOutput[galaxy].HaloIndex
//if( (flag && SageOutput[galaxy].FOFHaloIndex !=hid && SageOutput[galaxy].HaloIndex !=subhid) || (flag==2 && SageOutput[galaxy].FOFHaloIndex ==hid && SageOutput[galaxy].HaloIndex ==subhid) )
long int Len;
Len=count;//AllStars[id].Len;
//if(G.f_mb<10)
//if( (flag && SageOutput[galaxy].FOFHaloIndex ==hid && SageOutput[galaxy].HaloIndex ==subhid) || (flag==2 && SageOutput[galaxy].FOFHaloIndex !=hid && SageOutput[galaxy].HaloIndex !=subhid) )
//hid=0;
//subhid=0;
//if( (flag && SageOutput[galaxy].FOFHaloIndex ==subhid && SageOutput[galaxy].HaloIndex ==hid) || (flag==2 && SageOutput[galaxy].FOFHaloIndex !=subhid && SageOutput[galaxy].HaloIndex !=hid) )
if(AllStars[id].BindingEnergy <= BECut && Len>0 && hid !=-1 && subhid !=-1)
{//2
if(SageOutput[galaxy].HaloIndex !=subhid){
//printf("SMass:%g,SMassPre:%g\n",SageOutput[galaxy].StellarMass, SageOutputPre[galaxy].StellarMass);
AllStars[id].StellarMass=1.0e10*((SageOutput[galaxy].Stars)/Len);
//AllStars[id].GalNo=galaxy;//SageOutput[galaxy].
AllStars[id].TreeIndex=SageOutput[galaxy].TreeIndex;
AllStars[id].ZZ=SageOutput[galaxy].MetalsStellarMass/SageOutput[galaxy].StellarMass;//Len;//lower than expected
AllStars[id].Mvir=SageOutput[galaxy].Mvir;
AllStars[id].Rvir=SageOutput[galaxy].Rvir;
AllStars[id].infallMvir=SageOutput[galaxy].infallMvir;
AllStars[id].Age=GetAge(AllStars[id].Time);//this makes sense
AllStars[id].LastMajorMerger=SageOutput[galaxy].LastMajorMerger;
AllStars[id].HaloIndex=SageOutput[galaxy].HaloIndex;
AllStars[id].SubhaloIndex=SageOutput[galaxy].FOFHaloIndex;
AllStars[id].GalIndex=SageOutput[galaxy].CentralGal;
}

else //flag // if their energy is above the limit, unbind them
{
AllStars[id].BindingEnergy=0;
AllStars[id].StellarMass=0;
AllStars[id].ZZ=0;
}//
}
else // if they don't/do belong to the target halo
{
AllStars[id].BindingEnergy=0;
AllStars[id].StellarMass=0;
AllStars[id].ZZ=0;
}
AllStars[id].TreeIndex=SageOutput[galaxy].TreeIndex;
AllStars[id].HaloIndex=SageOutput[galaxy].HaloIndex;
AllStars[id].SubhaloIndex=SageOutput[galaxy].FOFHaloIndex;
return;
}
void ReadTargetHalo(char *FileName, int snap, int *hid,int *subhid)
{
  FILE *f;
  int s,h,sub,status=-1;
  f=fopen(FileName,"rb");
  printf("opened the file\n");
  if(f!=NULL)
  {
    do
    {
    status=fread(&s,sizeof(s),1,f);
    printf("read:s=%d\n",s);
    fread(&h,sizeof(h),1,f);
    printf("read h=%d\n",h);
    fread(&sub,sizeof(sub),1,f);
    printf("read:sub=%d\n",sub);
    if(s==snap)
     {
      *hid=h;
      //printf("read:sub=%d\n",*subhid);
      *subhid=sub;
      return;
     }
     //else
     //{
      // *hid=-1;
       //*subhid=-1;
     //}
   }while(status==1);
  }

fclose(f);
}

