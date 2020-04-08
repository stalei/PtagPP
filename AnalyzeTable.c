// © Shahram Talei @ 2019
// // Hash table analysis
// // 
// // 
// //
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
//#include <stdbool.h>
#include "hdf5.h"
#include <dirent.h>

#include "GlobalVars.h"
#include "proto.h"

#ifdef DoParallel
#include <mpi.h>
#endif


void AnalyzeHashTable(int snap,struct HashTable *table)
{
sprintf(CurrentFile,"AnalyzeTable.c");
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("Let's do some analysis!.☻ \n");
// 17940
//printf("test isstar: %lld \n",table->table[0]->key);
//printf("test isstar: %d\n",IsStar(table->table[17940]));
//

#ifdef DoParallel
if(ThisTask==0)
#endif
printf("table address inside analysis:%p\n",table);

long long *IdList;
IdList=(long long *)malloc(GP.TotNumPart);

//long long UniqueStarCount;
//UniqueStarCount
GP.TotNumStars=CountUniqueStars(table,IdList);
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("I found %lld unique stars in hash table!\n",GP.TotNumStars);

//let's extract final stellar halo information
struct tagged_particle *FinalStellarHalo;

if((FinalStellarHalo=(struct tagged_particle *)malloc(GP.TotNumStars*sizeof(struct tagged_particle)))==NULL)
	printf("couldn't allocate memory for final stellar halo\n");

char SnapFile[500];
sprintf(SnapFile,"%s/snap_%03d",GP.SnapDir,snap);
ReadSnap(SnapFile);


ExtractFinalHalo(IdList,table,FinalStellarHalo);


WriteFinalTag(FinalStellarHalo,GP.TotNumStars);

return;
}
bool IsStar(struct LinkedList *list)
{
//if(list->next ==0) printf("link star is zero!\n");
printf("the key in isstar:%lld, snap:%d\n",list->next->key,list->next->star->Snap);
if(list->next ==NULL)
	return 0;
else
	return 1;
//if(tag->key !=0)
//	return 1;
//else
//	return 0;
//return (tag->key !=0);
}
int IsTagged(struct LinkedList *list)
{
int c=0;
while(list->next)
{
	printf("snap in IsTagged:%d\n",list->next->star->Snap);
        c++;
	list=list->next;
}
return c;
}

long long CountUniqueStars(struct HashTable *table,long long *IdList)
{
long long c,i;
c=0;
for(i=0;i<GP.TotNumPart;i++)
	{
	//printf("counting:%lld\n",i);
	if(table->table[i]->next !=0)//NULL)
	   {
	//if(IsStar(table->table[i]))
		IdList[c]=i;
		c++;
	   }
	}
return c;
}

void ExtractFinalHalo(long long *IdList,struct HashTable *table,struct tagged_particle *FinalStellarHalo)
{
long long i;
int j;
struct LinkedList *list;
for(i=0;i<GP.TotNumStars;i++)
{
	//printf("id list:%lld\n",IdList[i]);
	list=table->table[IdList[i]];
	//printf("list in extract:%p\n",list->next);
	FinalStellarHalo[i].StellarMass=0;
	FinalStellarHalo[i].ZZ=0;
	while(list->next)
	{
		//printf("i:%lld, mass:%g\n",i,list->next->star->StellarMass);
		FinalStellarHalo[i].StellarMass += list->next->star->StellarMass;
		//FinalStellarHalo[i].ZZ+=list->next->star->ZZ;
		FinalStellarHalo[i].ZZ+=(list->next->star->ZZ)*(list->next->star->StellarMass);
		FinalStellarHalo[i].Age=list->next->star->Age; //not += just take the laste age!
		//if(list->next->star->Snap==GP.LastSnap) may not be tagged in the last snap
		//{
		//this is where they tagged last time not their position at the last snapshot
		//for(j=0;j<3;j++)
		//	FinalStellarHalo[i].Pos[j]=list->next->star->Pos[j];
		//FinalStellarHalo[i].Pos[0]=list->next->star->Pos[0];
		//FinalStellarHalo[i].Pos[1]=list->next->star->Pos[1];
		//FinalStellarHalo[i].Pos[2]=list->next->star->Pos[2];
		//}
		//or their position in the last snapshot
		for(j=0;j<3;j++)
		{
			if(P[IdList[i]].Type==1)// && P[IdList[i]].Pos[j] !=0 )//just dm particles
			{
			FinalStellarHalo[i].Pos[j]=P[IdList[i]].Pos[j];
			FinalStellarHalo[i].Vel[j]=P[IdList[i]].Vel[j];
			}
			//FinalStellarHalo[i].Pos[j]=P[list->next->star->PID].Pos[j];
                	//FinalStellarHalo[i].Vel[j]=P[list->next->star->PID].Vel[j];
		}
		//
		list=list->next;
	}
	FinalStellarHalo[i].ZZ/=FinalStellarHalo[i].StellarMass;

}
return;
}
