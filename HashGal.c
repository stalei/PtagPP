// © Shahram Talei @ 2020
// // Build the hash table of tagged particles but with Galaxy ID as key, This is used for each snapshot analysis
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


struct GalHashTable *GalEmptyTable(size_t gsize)
{
size_t i;
struct GalHashTable *gtable= (struct GalHashTable*)malloc(sizeof(struct GalHashTable));
//if((table= (struct HashTable*)malloc(sizeof(struct HashTable)))==NULL)
//	printf("Couldn't allocate memory for hash table!\n");
gtable->gtable=(struct GalLinkedList**)malloc(gsize*sizeof(struct GalLinkedList*));
for(i=0;i<gsize;i++)
	{
	gtable->gtable[i]=GalNewLinkedList();
	}
gtable->gsize=gsize;
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
//printf("table inside empty hashtable:%p\n",table);
return gtable;
}
struct GalLinkedList *GalNewLinkedList()
{
struct GalLinkedList *sentinel;
if((sentinel=(struct GalLinkedList*)malloc(sizeof(struct GalLinkedList)))==NULL)
	printf("Couldn't allocate memory for the first linked list!\n");
sentinel->gkey=0;
sentinel->next=0;
sentinel->star=0;
return sentinel;
}
void GalInsertKey(struct GalHashTable *gtable,long long gkey, struct tagged_particle *tag)
{
long long index;
index=gkey;//tag->PID; in this case both are the same but nut in general
//printf("the key in insert:%lld\n",key);
if(!GalContainsElement(gtable->gtable[index],gkey))
	{
//	printf("the key in insert- contains:%lld\n",key);
	GalAddElement(gtable->gtable[index],index,tag);
	}
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
//printf("table address inside insert key:%p\n",table);
return;
}
bool GalContainsElement (struct GalLinkedList *glist,long long gkey)
{
return GalGetPreviousLink(glist,gkey)!=0;
}
void GalAddElement(struct GalLinkedList *glist,long long gkey, struct tagged_particle *tag)
{
struct GalLinkedList *glink;
if((glink= (struct GalLinkedList*)malloc(sizeof(struct GalLinkedList)))==NULL)
	printf("Couldn't allocate memory for linked list!\n");
glink->gkey= gkey;
//printf("the key in add:%lld\n",link->key);
glink->star=tag;
if(tag==NULL)
	printf("tag is null!\n");
//if(link->star==NULL)
//        printf("link star is null!\n");
//if(link->star->PID==17940)
//	printf("I got your star! snap:%d\n",link->star->Snap);
//printf("snap in hash is :%d\n",link->star->Snap);
//memcpy(link->star,tag,sizeof(*tag));
//link->star->Snap=tag->Snap;

glink->next=glist->next;
glist->next=glink;

//if(link->star->PID==17940 && table->table[17940]->star !=NULL)
  //      printf("I got your star! snap:%d, key:%lld,table_snap:%d\n",link->star->Snap,key,table->table[17940]->star->Snap);
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
//printf("table address inside add element:%p\n",table);

return;
}
struct GalLinkedList *GalGetPreviousLink(struct GalLinkedList *glist,long long gkey) //I removed static
{
while(glist->next)
	{
	//if(list->next->key==17940)
	//	printf("inside while snap:%d\n",list->next->star->Snap);
	if(glist->next->gkey==gkey)
		return glist;
	glist=glist->next;
	}
return 0;
}

void GalDeleteTable(struct GalHashTable *gtable)
{
size_t i;
for(i=0;i<gtable->gsize;i++)
{
	GalDeleteLinkedList(gtable->gtable[i]);
}
free(gtable->gtable);
free(gtable);
return;
}
void GalDeleteLinkedList(struct GalLinkedList *glist)
{
while(glist)
	{
	struct GalLinkedList *next=glist->next;
	free(glist);
	glist=next;
	}
}
