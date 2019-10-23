// © Shahram @ 2019
// // Build the hash table of tagged particles
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


void HashTags(int snap)
{
sprintf(CurrentFile,"HashTag.c");
long long int TotalCount=0;
#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d started hashing tagged particles in parallel mode for snap %d.☻ \n",ThisTask,snap);
#else
printf("I started hashing tagged particles in serial mode for snap %d.☻ \n",snap);
#endif
fflush(stdout);
//find file name, let's set it manually!
char SnapFile[500];
sprintf(SnapFile,"%s/snap_%03d",GP.SnapDir,snap);
ReadSnapHeader(SnapFile);

fflush(stdout);

//GP.TotNumStarsAllSnaps=0;
//struct HashTable *table=EmptyTable(GP.TotNumPart);
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("table address before load:%p\n",table);
TotalCount=LoadAllTags(GP.FirstSnap,snap);
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("table address after load:%p\n",table);

#ifdef DoParallel
if(ThisTask==0)
#endif
if(TotalCount != GP.TotNumTagsAllSnaps)
	printf("mismatch in the number of all tags after loading!count is: %lld and it was: %lld\n",TotalCount,GP.TotNumTagsAllSnaps);

//struct HashTable *table= (struct HashTable*)malloc(sizeof(struct HashTable));

//ConstructHashTable(table,TotalCount,StellarHaloAllSnaps);//this loads all tags, not an efficient way

#ifdef DoParallel
MPI_Barrier(MPI_COMM_WORLD);
#endif
fflush(stdout);

AnalyzeHashTable(table);//sends this hashtable to the analysis code!

//////////////
#ifdef DoParallel
if(ThisTask==0)
printf("Task %d finished hashing.☻ \n",ThisTask);
#else
printf("I finished hashing.☻ \n");
#endif

return;
}

long long ReadSnapHeader(char *fname)
{//0
FILE *fd=0;
int blksize1,blksize2;;
int i;//,bnr,type;
//enum iofields blocknr;
//char buf[500];
//int bytes_per_blockelement;
//int  offset=0,typelist[6];
//size_t blockmaxlen;
//int nread,pc,nstart;
//double bytes_tot = 0;
//size_t bytes;
#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}


if(!(fd=fopen(fname,"r")))
 {
	printf("can't open snapshot <%s>\n",fname);
	EndRun(47,CurrentFile);
 }

SKIP;
my_fread(&header, sizeof(header), 1, fd);
SKIP2;
if(blksize1 != 256 || blksize2 != 256)
{
	printf("incorrect header format\n");
	fflush(stdout);
	EndRun(63,CurrentFile);// Probable error is wrong size of fill[] in header file. Needs to be 256 bytes in total.
}
GP.TotNumPart=0;
GP.TotNumStars=0;
if(header.num_files <= 1)
	for(i = 0; i < 6; i++)
	  {
	    header.npartTotal[i] = header.npart[i];
	    GP.TotNumPart+=header.npart[i];
	    printf("# of type %d from the header:%u \n",i,header.npartTotal[i]);
	  }

  return GP.TotNumPart;
}

long long LoadAllTags(int snapi,int snapf) //return could be used just to check if we coverd all tags correctly otherwise we don't need the count as we know it from previous part.
{
//long long id,tagid;
int tf,i;
long long c=0;
printf("Total count of tags (including all duplicated particles:%lld\n",GP.TotNumTagsAllSnaps);
struct HashTable *table=EmptyTable(GP.TotNumPart);
if(table==NULL)
	{
	printf("can't allocate memory for hashtable!\n");
        EndRun(141,CurrentFile);
	}
if((StellarHaloAllSnaps=(struct tagged_particle *)malloc(GP.TotNumTagsAllSnaps*sizeof(struct tagged_particle)))==NULL)
	{	
	printf("can't allocate memory for all tags!\n");
	EndRun(146,CurrentFile);
	}
c=0;
struct tagged_particle *StellarHalo;
for(tf=snapi;tf<=snapf;tf++)
	{/*B*/
	hid_t file,dataset,TagDatatype,dataspace;
	size_t size;
	hsize_t dims_out[2];
	herr_t status;
	long long rows;
	char TagFile[500];
	char *DSName="FullTag";
	//struct tagged_particle *StellarHalo;
	int rank,status_n;
	sprintf(TagFile,"%s/tag_%03d.h5",GP.OutputDir,tf);
	/*strcpy(TagFile,"tag_098.h5");*/
	/*printf("openning tagfile: %s\n",TagFile);*/
	file = H5Fopen(TagFile, H5F_ACC_RDONLY, H5P_DEFAULT);
	/*printf("file:%d\n",file);*/
	fflush(stdout);
	dataset = H5Dopen(file,DSName,H5P_DEFAULT);
	//printf("dataset:%d\n",dataset);
	TagDatatype=H5Dget_type(dataset);
	size  = H5Tget_size(TagDatatype);
	if(size != sizeof(struct tagged_particle))
	{
	    printf("size mismatch, data size:%d, struct tag size:%lu d\n",(int)size,sizeof(struct tagged_particle));
	}
	dataspace = H5Dget_space(dataset);    /* dataspace handle*/ 
	rank = H5Sget_simple_extent_ndims(dataspace);
	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
	if (status_n<0) printf("Error in reading dimensions\n");
	rows = dims_out[0];
	if(rank<0)
	    printf("I coulldn't read the file:%s\n",TagFile);
	if((StellarHalo=(struct tagged_particle *)malloc(rows*size))==NULL)
		{
        	printf("can't allocate memory for stellar halo for snapshot %d!\n",tf);
        	EndRun(185,CurrentFile);
	        }

	status= H5Dread(dataset,TagDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, StellarHalo);
	if(status<0) printf("Error in reading data from the file:%s\n",TagFile);
	for(i=0;i<rows;i++)
		StellarHaloAllSnaps[c+i]=StellarHalo[i];
	//ConstructHashTable(rows,StellarHalo);
	//#ifdef DoParallel
	//if(ThisTask==0)
	//#endif
	//printf("I read the tag file, let's put it in the table!\n");
	//for(i=0;i<rows;i++)
	//	InsertKey(table,StellarHalo[i].PID,&StellarHalo[i]);
	//printf("Task %d is done with snap:%d!\n",ThisTask,tf);
	//#ifdef DoParallel
	//MPI_Barrier(MPI_COMM_WORLD);
	//#endif
	//free(StellarHalo);
	//H5Tclose(TagDatatype);
	//H5Dclose(dataset);
	//H5Sclose(dataspace);
	//H5Fclose(file);
	
	c+=rows;
	}/*B*/
//struct tagged_particle *t=*table->table[17940]->star;
//printf("key test:%lld, IsStar:%d, snap:%d \n",table->table[17940]->key,IsStar(table->table[17940]),table->table[17940]->star->Snap);
//struct HashTable *table=EmptyTable(GP.TotNumPart);
for(i=0;i<GP.TotNumTagsAllSnaps;i++)
	{
	printf("Tag add:%lld\n",StellarHaloAllSnaps[i].PID);
	InsertKey(table,StellarHaloAllSnaps[i].PID,&StellarHaloAllSnaps[i]);
	}
//ConstructHashTable(GP.TotNumTagsAllSnaps,StellarHaloAllSnaps);
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("table inside load:%p\n",table);

//printf("key test:%lld, IsStar:%d \n",table->table[17940]->key,IsStar(table->table[17940]));
if(table->table[17940]->star !=NULL)
   printf("key test:%lld, IsStar:%d, snap:%d \n",table->table[17940]->key,IsStar(table->table[17940]),table->table[17940]->star->Snap);
else
   printf("star is null!\n");
//printf("test tag count is:%d\n",IsTagged(table->table[17940]));

#ifdef DoParallel
if(ThisTask==0)
printf("Task %d finished loading all tagged files in hash table!\n",ThisTask);
#else
printf("I finished loading all tagged files in hash table!\n");
#endif

return c;
}

void ConstructHashTable(struct HashTable *table,long long TagCount,struct tagged_particle *Halo)
{
long long i;

//struct HashTable *table=EmptyTable(GP.TotNumPart);
for(i=0;i<TagCount;i++)
	InsertKey(table,i,&Halo[i]);

return;
}

struct HashTable *EmptyTable(size_t size)
{
size_t i;
//if((struct HashTable *table= (struct HashTable*)malloc(sizeof(struct HashTable)))==NULL)
if((table= (struct HashTable*)malloc(sizeof(struct HashTable)))==NULL)
	printf("Couldn't allocate memory for hash table!\n");
table->table=(struct LinkedList**)malloc(size*sizeof(struct LinkedList*));
for(i=0;i<size;i++)
	{
	table->table[i]=NewLinkedList();
	}
table->size=size;
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
//printf("table inside empty hashtable:%p\n",table);
return table;
}
struct LinkedList *NewLinkedList()
{
struct LinkedList *sentinel;
if((sentinel=(struct LinkedList*)malloc(sizeof(struct LinkedList)))==NULL)
	printf("Couldn't allocate memory for the first linked list!\n");
sentinel->key=0;
sentinel->next=0;
sentinel->star=0;
return sentinel;
}
void InsertKey(struct HashTable *table,long long key, struct tagged_particle *tag)
{
long long index;
index=key;//tag->PID; in this case both are the same but nut in general
printf("the key in insert:%lld\n",key);
if(!ContainsElement(table->table[index],key))
	{
	printf("the key in insert- contains:%lld\n",key);
	AddElement(table->table[index],index,tag);
	}
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
//printf("table address inside insert key:%p\n",table);
return;
}
bool ContainsElement (struct LinkedList *list,long long key)
{
return GetPreviousLink(list,key);
}
void AddElement(struct LinkedList *list,long long key, struct tagged_particle *tag)
{
struct LinkedList *link;
if((link= (struct LinkedList*)malloc(sizeof(struct LinkedList)))==NULL)
	printf("Couldn't allocate memory for linked list!\n");
link->key= key;
printf("the key in add:%lld\n",link->key);
link->star=tag;
if(tag==NULL)
	printf("tag is null!\n");
if(link->star==NULL)
        printf("link star is null!\n");
if(link->star->PID==17940)
	printf("I got your star! snap:%d\n",link->star->Snap);
//printf("snap in hash is :%d\n",link->star->Snap);
//memcpy(link->star,tag,sizeof(*tag));
//link->star->Snap=tag->Snap;

link->next=list->next;
list->next=link;

//if(link->star->PID==17940 && table->table[17940]->star !=NULL)
  //      printf("I got your star! snap:%d, key:%lld,table_snap:%d\n",link->star->Snap,key,table->table[17940]->star->Snap);
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
//printf("table address inside add element:%p\n",table);

return;
}
struct LinkedList *GetPreviousLink(struct LinkedList *list,long long key) //I removed static
{
while(list->next)
	{
	if(list->next->key==17940)
		printf("inside while snap:%d\n",list->next->star->Snap);
	if(list->next->key==key)
		return list;
	list=list->next;
	}
return 0;
}

void DeleteTable(struct HashTable *table)
{
size_t i;
for(i=0;i<table->size;i++)
{
	DeleteLinkedList(table->table[i]);
}
free(table->table);
free(table);
return;
}
void DeleteLinkedList(struct LinkedList *list)
{
while(list)
	{
	struct LinkedList *next=list->next;
	free(list);
	list=next;
	}
}
