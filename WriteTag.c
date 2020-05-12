// © Shahram @ 2019
// // Reading stored data in a tag files
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

#include "hdf5.h"
#include <dirent.h>

#include "GlobalVars.h"
#include "proto.h"

#ifdef DoParallel
#include<mpi.h>
#endif

void WriteTag(int snap)
{
sprintf(CurrentFile,"WriteTag.c");
#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d started writing the tag file for snapshot %d.☻ \n",ThisTask,snap);
#else
printf("I started writing the tag file for snapshot %d.☻ \n",snap);
#endif

// we have TagDatatype=H5Dget_type(dataset); from ReadTag.c
hid_t Tfile_id;
hid_t  Tdset_id;//, Tdatatype, Tgroup;
hid_t Tfilespace;//, Tmemspace;
hsize_t Tdimsf[2]; //1D array of struct, total size
hid_t Tplist_id;
herr_t Tstatus;
char TagFile[500];
char Tdataset_name[]="FullTag";//, group_name[100];

sprintf(TagFile,"%s/tag_%03d.h5",GP.OutputDir,snap);
Tplist_id=H5Pcreate(H5P_FILE_ACCESS);
Tfile_id=H5Fcreate(TagFile, H5F_ACC_TRUNC, H5P_DEFAULT,Tplist_id);

Tdimsf[0]=NumOfStars;
Tdimsf[1]=1;
//typedef struct tagged_particle tagged_tmp;

Tfilespace=H5Screate_simple(2,Tdimsf,NULL);
//Tdatatype=H5Tcreate(H5T_COMPOUND,sizeof(struct tagged_particle));
//sprintf(group_name,"/tagHalo_%d_Subhalo_%d",tagged_p.HaloIndex,tagged_p.SubhaloIndex);
//sprintf(Tdataset_name,"%s/gr%d", group_name, tagged_p.GalIndex);
//        Tgroup=H5Gcreate(Tfile_id,group_name,H5P_DEFAULT);
//
Tdset_id=H5Dcreate(Tfile_id,Tdataset_name,GTagDatatype, Tfilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
Tstatus=H5Dwrite(Tdset_id, GTagDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, AllStars);
if(Tstatus<0) printf("Failed to save full tag file for snap %d.\n",snap);
//H5Gclose(Tgroup);
H5Dclose(Tdset_id);
H5Tclose(GTagDatatype);
H5Sclose(Tfilespace);
H5Pclose(Tplist_id);
H5Fclose(Tfile_id);


#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d finished writing the tag file for snapshot %d.☻ \n",ThisTask,snap);
#else
printf("I finished writing the tag file for snapshot %d.☻ \n",snap);
#endif

return;
}

void WriteFinalTag(struct tagged_particle *halo,long long count)
{
#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d started writing final tag.☻ \n",ThisTask);
#else
printf("I started writing final tag.☻ \n");
#endif

// we have TagDatatype=H5Dget_type(dataset); from ReadTag.c
hid_t Tfile_id;
hid_t  Tdset_id, TagDatatype;//, Tgroup;
hid_t Tfilespace;//, Tmemspace;
hsize_t Tdimsf[2]; //1D array of struct, total size
hid_t Tplist_id;
herr_t Tstatus;
char TagFile[500];
char Tdataset_name[]="FinalTag";//, group_name[100];

sprintf(TagFile,"%s/StellarHalo.h5",GP.OutputDir);
Tplist_id=H5Pcreate(H5P_FILE_ACCESS);
Tfile_id=H5Fcreate(TagFile, H5F_ACC_TRUNC, H5P_DEFAULT,Tplist_id);

Tdimsf[0]=count;
Tdimsf[1]=1;
//typedef struct tagged_particle tagged_tmp;

Tfilespace=H5Screate_simple(2,Tdimsf,NULL);
//Tdatatype=H5Tcreate(H5T_COMPOUND,sizeof(struct tagged_particle));
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
H5Tinsert(TagDatatype, "GalNo", HOFFSET(tagged_tmp, GalNo), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "TreeIndex", HOFFSET(tagged_tmp, TreeIndex), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "AA", HOFFSET(tagged_tmp, AA), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "Age", HOFFSET(tagged_tmp, Age), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "ZZ", HOFFSET(tagged_tmp, ZZ), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "StellarMass", HOFFSET(tagged_tmp, StellarMass), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Time", HOFFSET(tagged_tmp, Time), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Len", HOFFSET(tagged_tmp, Len), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "MBID", HOFFSET(tagged_tmp, MBID), H5T_NATIVE_LLONG);
H5Tinsert(TagDatatype, "BindingEnergy", HOFFSET(tagged_tmp, BindingEnergy), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Mvir", HOFFSET(tagged_tmp, Mvir), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Rvir", HOFFSET(tagged_tmp, Rvir), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "infallMvir", HOFFSET(tagged_tmp, infallMvir), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "LastMajorMerger", HOFFSET(tagged_tmp, LastMajorMerger), H5T_NATIVE_FLOAT);

//sprintf(group_name,"/tagHalo_%d_Subhalo_%d",tagged_p.HaloIndex,tagged_p.SubhaloIndex);
//sprintf(Tdataset_name,"%s/gr%d", group_name, tagged_p.GalIndex);
//        Tgroup=H5Gcreate(Tfile_id,group_name,H5P_DEFAULT);
//
Tdset_id=H5Dcreate(Tfile_id,Tdataset_name,TagDatatype, Tfilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
Tstatus=H5Dwrite(Tdset_id, TagDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, halo);
if(Tstatus<0) printf("Failed to save final full tag.\n");
//H5Gclose(Tgroup);
H5Dclose(Tdset_id);
H5Tclose(TagDatatype);
H5Sclose(Tfilespace);
H5Pclose(Tplist_id);
H5Fclose(Tfile_id);


#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d finished writing final tag.☻ \n",ThisTask);
#else
printf("I finished writing final tag.☻ \n");
#endif

return;
}

void WriteFirstTagged(struct tagged_particle *halo,long long count)
{
#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d started writing first tagged.☻ \n",ThisTask);
#else
printf("I started writing first tagged.☻ \n");
#endif

// we have TagDatatype=H5Dget_type(dataset); from ReadTag.c
hid_t Tfile_id;
hid_t  Tdset_id, TagDatatype;//, Tgroup;
hid_t Tfilespace;//, Tmemspace;
hsize_t Tdimsf[2]; //1D array of struct, total size
hid_t Tplist_id;
herr_t Tstatus;
char TagFile[500];
char Tdataset_name[]="FirstTagged";//, group_name[100];

sprintf(TagFile,"%s/FirstTagged.h5",GP.OutputDir);
Tplist_id=H5Pcreate(H5P_FILE_ACCESS);
Tfile_id=H5Fcreate(TagFile, H5F_ACC_TRUNC, H5P_DEFAULT,Tplist_id);

Tdimsf[0]=count;
Tdimsf[1]=1;
//typedef struct tagged_particle tagged_tmp;

Tfilespace=H5Screate_simple(2,Tdimsf,NULL);
//Tdatatype=H5Tcreate(H5T_COMPOUND,sizeof(struct tagged_particle));
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
H5Tinsert(TagDatatype, "GalNo", HOFFSET(tagged_tmp, GalNo), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "TreeIndex", HOFFSET(tagged_tmp, TreeIndex), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "AA", HOFFSET(tagged_tmp, AA), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "Age", HOFFSET(tagged_tmp, Age), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "ZZ", HOFFSET(tagged_tmp, ZZ), H5T_NATIVE_FLOAT);
H5Tinsert(TagDatatype, "StellarMass", HOFFSET(tagged_tmp, StellarMass), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Time", HOFFSET(tagged_tmp, Time), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Len", HOFFSET(tagged_tmp, Len), H5T_NATIVE_INT);
H5Tinsert(TagDatatype, "MBID", HOFFSET(tagged_tmp, MBID), H5T_NATIVE_LLONG);
H5Tinsert(TagDatatype, "BindingEnergy", HOFFSET(tagged_tmp, BindingEnergy), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Mvir", HOFFSET(tagged_tmp, Mvir), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "Rvir", HOFFSET(tagged_tmp, Rvir), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "infallMvir", HOFFSET(tagged_tmp, infallMvir), H5T_NATIVE_DOUBLE);
H5Tinsert(TagDatatype, "LastMajorMerger", HOFFSET(tagged_tmp, LastMajorMerger), H5T_NATIVE_FLOAT);

//sprintf(group_name,"/tagHalo_%d_Subhalo_%d",tagged_p.HaloIndex,tagged_p.SubhaloIndex);
//sprintf(Tdataset_name,"%s/gr%d", group_name, tagged_p.GalIndex);
//
Tdset_id=H5Dcreate(Tfile_id,Tdataset_name,TagDatatype, Tfilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
Tstatus=H5Dwrite(Tdset_id, TagDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, halo);
if(Tstatus<0) printf("Failed to save final full tag.\n");
//H5Gclose(Tgroup);
H5Dclose(Tdset_id);
H5Tclose(TagDatatype);
H5Sclose(Tfilespace);
H5Pclose(Tplist_id);
H5Fclose(Tfile_id);


#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d finished writing first tagged.☻ \n",ThisTask);
#else
printf("I finished writing first tagged.☻ \n");
#endif

return;
}

