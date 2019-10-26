// © Shahram @ 2019
// Reading stored data in a tag files
// 
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

#include "hdf5.h"
#include <dirent.h>
//
#include "GlobalVars.h"
#include "proto.h"

#ifdef DoParallel
#include <mpi.h>
#endif

void ReadTag(int snap)
{
sprintf(CurrentFile,"ReadTag.c");
long long counter,countersum;
counter=0;
countersum=0;
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("I started extracting information from the tag files.☻ \n");
TagFilesCount = CountTagFiles(snap);
TagFilesCountPre = CountTagFiles(snap-1);
if (TagFilesCount<1){
	printf("No tag file for snap %d!\n",snap);
	EndRun(30,CurrentFile);
	}
else
   printf("I found %d tag file(s) for snapshot %d.\n", TagFilesCount,snap);

if (TagFilesCountPre<1){
        printf("No tag file for snap %d!\n",snap-1);
        EndRun(30,CurrentFile);
        }
else
   printf("I found %d tag file(s) for snapshot %d.\n", TagFilesCountPre,snap-1);


TagFilesPath = (struct Path_Names*)malloc(TagFilesCount * sizeof(struct Path_Names));
TagFilesPathPre = (struct Path_Names*)malloc(TagFilesCountPre * sizeof(struct Path_Names));


ReadTagFNames(snap,TagFilesPath);// note that 0th file is the last!
ReadTagFNames(snap-1,TagFilesPathPre);// it is the same as previous, but anyway!
fflush(stdout);
//printf("sample tag file:%s\n",TagFilesPath[0].paths);
NumOfStars=0;
NumOfStarsPre=0;

counter=CountStars(TagFilesCount, TagFilesPath); //counter of stars for current snapshot current process
//counter+=NumOfStars;
#ifdef DoParallel
//MPI_Barrier(MPI_COMM_WORLD);
MPI_Allreduce(&counter, &countersum, 1, MPI_LONG_LONG_INT, MPI_SUM,MPI_COMM_WORLD);
#else
countersum=counter;
#endif
GP.TotNumTagsAllSnaps=countersum;
NumOfStars=counter;
printf("Total Number of Tagged Stars: all snapshots=%ld for snap %d=%ld\n",GP.TotNumTagsAllSnaps,snap,NumOfStars);

NumOfStarsPre=CountStars(TagFilesCountPre, TagFilesPathPre);


#ifdef DoParallel
if(ThisTask==0)
#endif
printf("Total Number of Tagged Stars:%ld for snap %d\n",NumOfStarsPre,snap-1);


if((AllStars=(struct tagged_particle *)malloc(NumOfStars*sizeof(struct tagged_particle)))==NULL)
    {
	printf("Failed to allocate memory for stars!!\n");
	EndRun(66,CurrentFile);
    }

if((AllStarsPre=(struct tagged_particle *)malloc(NumOfStarsPre*sizeof(struct tagged_particle)))==NULL)
    {
        printf("Failed to allocate memory for stars!!\n");
        EndRun(72,CurrentFile);
    }


ReadCombineTags(TagFilesCount,TagFilesPath,AllStars);

#ifdef DoParallel
if(ThisTask==0)
#endif
printf("All tags are loaded for snap:%d\n",snap);

ReadCombineTags(TagFilesCountPre,TagFilesPathPre,AllStarsPre);

#ifdef DoParallel
if(ThisTask==0){
#endif
printf("All tags are loaded for snap:%d\n",snap-1);

printf("sample star:%g\n",AllStars[500].Pos[0]);
//PrintStar(673);
PrintStar(0);
#ifdef DoParallel
}
#endif
//PrintStar(10000);
//SavePositions();//just a test
////////////////// Just remove this later
//
//free(AllStars);
///////////////////
fflush(stdout);

return;
}

int CountTagFiles(int snap)
{
//printf("Let me count tag files for snap %d:\n",snap);
char dir[500];
int c=0;
//strcpy(dir,GP.TagDir);
sprintf(dir,"%s/tag_%03d",GP.TagDir,snap);

printf("Counting files in %s\n",dir);
struct dirent *dp;
DIR *fd;

  if ((fd = opendir(dir)) == NULL) {
    fprintf(stderr, "listdir: can't open %s\n", dir);
    return -1;
  }
  while ((dp = readdir(fd)) != NULL) {
  if (!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, ".."))
    continue;    /* skip self and parent */
  //printf("%s/%s\n", dir, dp->d_name);
  c++;
  }
  closedir(fd);
return c;
}

void ReadTagFNames(int snap, struct Path_Names *TagFile)
{
printf("Loading tag file-name(s) for snap %d:\n",snap);
char dir[500];
int c=0;
//strcpy(dir,GP.TagDir);
sprintf(dir,"%s/tag_%03d",GP.TagDir,snap);
//printf("Counting files in %s\n",dir);
struct dirent *dp;
  DIR *fd;

  if ((fd = opendir(dir)) == NULL) {
    fprintf(stderr, "listdir: can't open %s\n", dir);
    return;
  }
  while ((dp = readdir(fd)) != NULL) {
  if (!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, ".."))
    continue;    /* skip self and parent */
  sprintf(TagFile[c].paths,"%s/%s", dir, dp->d_name);
  c++;
  }
  closedir(fd);
#ifdef DoParallel
if(ThisTask==0)
#endif
printf("I Finished reading tag files-name(s)\n");
return;
}

void InitializeTagDatatype()
{
hid_t TagDatatype;
typedef struct tagged_particle tagged_tmp;

TagDatatype=H5Tcreate(H5T_COMPOUND,sizeof(struct tagged_particle));//124);// sizeof(struct tagged_p));//H5Tcopy(H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "X", HOFFSET(tagged_tmp, Pos[0]), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Y", HOFFSET(tagged_tmp, Pos[1]), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Z", HOFFSET(tagged_tmp, Pos[2]), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Vx", HOFFSET(tagged_tmp, Vel[0]), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Vy", HOFFSET(tagged_tmp, Vel[1]), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Vz", HOFFSET(tagged_tmp, Vel[2]), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "ID", HOFFSET(tagged_tmp, ID), H5T_NATIVE_LLONG);
H5Tinsert(TagDatatype, "Snap", HOFFSET(tagged_tmp, Snap), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "PID", HOFFSET(tagged_tmp, PID), H5T_NATIVE_LLONG);
H5Tinsert(TagDatatype, "HaloIndex", HOFFSET(tagged_tmp, HaloIndex), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "SubhaloIndex", HOFFSET(tagged_tmp, SubhaloIndex), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "GalIndex", HOFFSET(tagged_tmp, GalIndex), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "AA", HOFFSET(tagged_tmp, AA), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "Age", HOFFSET(tagged_tmp, Age), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "ZZ", HOFFSET(tagged_tmp, ZZ), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "StellarMass", HOFFSET(tagged_tmp, StellarMass), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Time", HOFFSET(tagged_tmp, Time), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Len", HOFFSET(tagged_tmp, Len), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "MBID", HOFFSET(tagged_tmp, MBID), H5T_NATIVE_LLONG);
H5Tinsert(TagDatatype, "BindingEnergy", HOFFSET(tagged_tmp, BindingEnergy), H5T_NATIVE_DOUBLE);
H5Tclose(TagDatatype);
return;
}


long int CountStars(int count, struct Path_Names *FilePath)
{
long int c=0;
printf("Counting number of tagged particles.\n...\n");
//long int c=0;
//int i;
hid_t       file=0, dataset,TagDatatype;         /* handles */
hid_t        dataspace;   
//H5T_class_t class;                 /* data type class */
//H5T_order_t order;                 /* data order */
size_t      size;                  /*
				        * size of the data element	       
				        * 				        * stored in file
				        * 				        				        */
    hsize_t     dims_out[2];           /* dataset dimensions */      
    //herr_t      status;                             
    //int          **data_out; 
    long unsigned int          rows, cols;
    char        TagFile[500];//="/stuperm/stalei/run/hdftest/tag_098/tag_serial_030.h5";
    char        DSName[]="MiniTag";
    int          i, status_n, rank;
   for(i=0,c=0;i<count;i++)
{
    //sprintf(TagFile,"%s",TagFilesPath[0].paths);
    strcpy(TagFile,FilePath[i].paths);
    //printf("compare:%d\n", strncmp(TagFile,"/stuperm/stalei/run/hdftest/tag_098/tag_serial_030.h5",55));
    //printf("TagFile:-%s-\npath:-%s-\n",TagFile,TagFilesPath[i].paths);
     //InitializeTagDatatype();
    /*  Open the file and the dataset  */

       // printf("before openning the file: %s\n",TagFile);
        //fflush(stdout);

    file = H5Fopen(TagFile, H5F_ACC_RDONLY, H5P_DEFAULT);
	//printf("just opened the file: %s\n",TagFile);
	//fflush(stdout);
    //H5Fclose(file);
    dataset = H5Dopen(file, DSName,H5P_DEFAULT);
    TagDatatype=H5Dget_type(dataset);
    GTagDatatype=H5Dget_type(dataset);
    /* Get datatype and dataspace handles and then query
 *        dataset class, order, size, rank and dimensions  */
    //datatype  = H5Dget_type(dataset);     /* datatype handle */ 
    //class     = H5Tget_class(TagDatatype);
    //if (class == H5T_INTEGER) printf("Data set has INTEGER type \n");
    //order     = H5Tget_order(TagDatatype);
    //if (order == H5T_ORDER_LE) printf("Little endian order \n");
	//InitializeTagDatatype();
    size  = H5Tget_size(TagDatatype);
    if(size != sizeof(struct tagged_particle))
      {
      printf("size mismatch, data size:%d, struct tagsize:%lu d\n",(int)size,sizeof(struct tagged_particle));
    //H5Tclose(TagDatatype);
    //H5Dclose(dataset);
    //H5Sclose(dataspace);
    //H5Fclose(file);
     // EndRun(215,CurrentFile);
      }

    //printf(" Data size is %d \n", size);

    dataspace = H5Dget_space(dataset);    /* dataspace handle */
    //H5Dclose(dataset);
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if (status_n<0) printf("Error in reading dimensions\n");
    //printf("rank %d, dimensions %lu x %lu , each data size:%d \n", rank,
	//   (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]),size);
    rows = dims_out[0];
    cols = dims_out[1];
    printf("Found %lu star(s) in:<%s> (%d x %lu, size:%d) \n",rows,TagFile, rank,cols,(int)size);
    TagCols=cols;
    TagSize=size;
    if(rank>0)
	c+=rows;
    else
	printf("I coulldn't read the file:%s\n",TagFile);
   // printf("Num of records, Size:%d\n",dataset.Tables[0].Select("ID is not null").Length);
    H5Tclose(TagDatatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);
    //printf(".");
} // end of for i

printf("\nCounting is done!\n");

return c;
}

void ReadCombineTags(int count,struct Path_Names *TagPath,struct tagged_particle *Stars)
{

printf("Starting reading and combining tags!\n...\n");
int c;
hid_t       file, dataset,TagDatatype;         /* handles */
hid_t        dataspace;
//H5T_class_t class;                 /* data type class */
//H5T_order_t order;                 /* data order */
size_t      size;                  /* * size of the data element *** stored in file***/
hsize_t     dims_out[2];           /* dataset dimensions */
herr_t      status;
//int          **data_out;
long unsigned int          rows;//, cols;
char        TagFile[500];//="/stuperm/stalei/run/hdftest/tag_098/tag_serial_030.h5";
char        DSName[]="MiniTag";
int          i, j, status_n, rank;
c=0;
i=0;
struct tagged_particle *StellarHalo;
// Allocate memory
//MyStar=(struct tagged_particle **)malloc(NumOfStars*sizeof(struct tagged_particle*));
//AllStars=(struct tagged_particle *)malloc(count*sizeof(struct tagged_particle));
//MyStar[c]=(struct tagged_particle *)malloc(TagCols*TagSize);

//for (i=1; i < NumOfStars; i++) MyStar[i] = MyStar[0]+i*TagCols;

for(i=count-1;i>=0;i--)
   //for(i=0,c=0;i<TagFilesCount;i++)
{
    strcpy(TagFile,TagPath[i].paths);
    //printf("TagFile:-%s- Done!!!-\n",TagFile);
    /*  Open the file and the dataset  */
    file = H5Fopen(TagFile, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen(file, DSName,H5P_DEFAULT);
    TagDatatype=H5Dget_type(dataset);
    /* Get datatype and dataspace handles and then query
 **dataset class, order, size, rank and dimensions  */
    //datatype  = H5Dget_type(dataset);     /* datatype handle */ 
    //class     = H5Tget_class(TagDatatype);
    //if (class == H5T_INTEGER) printf("Data set has INTEGER type \n");
    //order     = H5Tget_order(TagDatatype);
    //if (order == H5T_ORDER_LE) printf("Little endian order \n");
    size  = H5Tget_size(TagDatatype);
    if(size != sizeof(struct tagged_particle))
      {
      printf("size mismatch, data size:%d, struct tag size:%lu d\n",(int)size,sizeof(struct tagged_particle));
      //EndRun(280,CurrentFile);
      }
    dataspace = H5Dget_space(dataset);    /* dataspace handle */
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if (status_n<0) printf("Error in reading dimensions\n");
    //printf("rank %d, dimensions %lu x %lu , each data size:%d \n", rank,
      //     (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]),size);

    rows = dims_out[0];
//    cols = dims_out[1];
  //  printf("rank %d, dimensions %lu x %lu , each data size:%d \n", rank,rows,cols,(int)size);


    if(rank>0)
        {}//c+=rows;
    else
        printf("I coulldn't read the file:%s\n",TagFile);
//    printf("Num of records, Size:%d\n",dataset.Tables[0].Select("ID is not null").Length);

/*************************** BEGIN  *******************************/

/* Allocate memory for new integer array[row][col]. First
 *    allocate the memory for the top-level array (rows).  Make
 *       sure you use the sizeof a *pointer* to your data type. */

    //data_out = (int**) malloc(rows*sizeof(int*));
   // MyStar=(struct tagged_particle **)malloc(NumOfStars*sizeof(struct tagged_particle*));
/* Allocate a contiguous chunk of memory for the array data values.  
 *    Use the sizeof the data type. */

    //data_out[0] = (int*)malloc( cols*rows*sizeof(int) );
    //MyStar[c]=(struct tagged_particle *)malloc(TagCols*TagSize);
/* Set the pointers in the top-level (row) array to the
 *    correct memory locations in the data value chunk. */

//    for (j=1; j < rows; j++) data_out[j] = data_out[0]+j*cols;
	
/************************* END *************************************/
StellarHalo=(struct tagged_particle *)malloc(rows*size);
//status = H5Dread(dataset, TagDatatype, H5S_ALL, H5S_ALL,H5P_DEFAULT, &data_out[0][0]);
status= H5Dread(dataset,TagDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, StellarHalo);
if(status<0) printf("Error in reading data from the file:%s\n",TagFile);
//status = H5Dread(dataset, TagDatatype, H5S_ALL, H5S_ALL,H5P_DEFAULT, MyStar[c][0]);
//printf("data:%g\n",StellarHalo[0].Pos[0]);
for(j=0;j<rows;j++)
    Stars[c+j]=StellarHalo[j];

H5Tclose(TagDatatype);
H5Dclose(dataset);
H5Sclose(dataspace);
H5Fclose(file);
free(StellarHalo);
c+=rows;
//#ifdef DoParallel
//if(ThisTask==0)
//#endif
printf("☆ ★ %d...✓ \n ",i);
}  //end of for i

//#ifdef DoParallel
//if(ThisTask==0)
//#endif
printf("\n\nRead and Combine is finished!\n");
return;
}
void PrintStar(long int index)
{
if(index>=0 && index<NumOfStars)
{
	int i=index;
	printf("\n\n☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆\n");
	printf("Information of the star %li\n",index);
	printf("Time:%g, Snap:%d\n",AllStars[i].Time,AllStars[i].Snap);
        printf("Pos:(%g,%g,%g)\n",AllStars[i].Pos[0],AllStars[i].Pos[1],AllStars[i].Pos[2]);
        printf("Vel:(%g,%g,%g)\n",AllStars[i].Vel[0],AllStars[i].Vel[1],AllStars[i].Vel[2]);
        printf("Halo:%d, Subhalo:%d,Gal:%d\n",AllStars[i].HaloIndex,AllStars[i].SubhaloIndex,AllStars[i].GalIndex);
        printf("AA:%g,Age:%g\n",AllStars[i].AA,AllStars[i].Age);
        //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        printf("☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆ ☆\n\n");
return;

}
else
{
	printf("Star index is out of range!");
	return;
}


}

void SavePositions()
{
int i;
FILE *fp;
char file[500];//="StarPos";
sprintf(file,"%s/StarPos",GP.OutputDir);
fp=fopen(file,"wa");
for(i=0;i<NumOfStars;i++)
{
fprintf(fp,"%f,%f,%f\n",AllStars[i].Pos[0],AllStars[i].Pos[1],AllStars[i].Pos[2]);
} 
fclose(fp);
return;
}
