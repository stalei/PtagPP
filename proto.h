// © Shahram Talei @ 2019
#ifndef ALLVARS_H
#include "GlobalVars.h"
#endif
#include <stdbool.h>

void EndRun(int,char[]);
void initialize(void);
void ReadParameters(char *fname);
void SetGP(void);
void ReadSage(int);
void LoadSageFiles(int);
int CountSageFiles(int);
void ReadSageFNames(int, struct Path_Names *);
int ReadSageHeader(int, struct Path_Names *);
void ReadSageModel(int, struct Path_Names *, struct SageGalaxies *);
void PrintGalaxyInfo(struct SageGalaxies *,int);
void ReadTag(int);
int CountTagFiles(int);
void ReadTagFNames(int,struct Path_Names *);

long int CountStars(int, struct Path_Names *);

void InitializeTagDatatype(void);
void ReadCombineTags(int, struct Path_Names *,struct tagged_particle *);
void PrintStar(long int);
void SavePositions(void);

void PaintStars(int);
int LookupGalaxy(int, unsigned int, struct SageGalaxies *,int,float);
void CalculateStellarProperties(double,double,int,unsigned long int);
//double FindTime(int);
double FindTime(struct tagged_particle *);
void WriteTag(int);

void WriteTaggedSnap(int);
void ReadSnap(char *);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
int blockpresent(enum iofields blocknr);
void get_dataset_name(enum iofields blocknr, char *buf);
int get_bytes_per_blockelement(enum iofields blocknr, int mode);
int get_particles_in_block(enum iofields blocknr, int *typelist);
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);

#ifdef DoParallel
long long CountTotStarsPar(int,int);
#else
long long CountTotStars(int,int);
#endif
void SaveInBinary(void);
double GetAge(double);
void HashTags(int);
long long ReadSnapHeader(char *);
long long CountAllTags(int,int);
long long LoadAllTags(int,int,struct HashTable *);
void ConstructHashTable(struct HashTable *,long long,struct tagged_particle *);
struct HashTable *EmptyTable(size_t);
struct LinkedList *NewLinkedList(void);
void InsertKey(struct HashTable *,long long,struct tagged_particle *);
bool ContainsElement (struct LinkedList *,long long);
void AddElement(struct LinkedList *,long long,struct tagged_particle *);
struct LinkedList *GetPreviousLink(struct LinkedList *,long long); // I removed static
void DeleteTable(struct HashTable *);
void DeleteLinkedList(struct LinkedList *);
//struct HashTable *EmptyTable(size_t);
//struct LinkedList *NewLinkedList();



void AnalyzeHashTable(int ,struct HashTable *);
bool IsStar(struct LinkedList *);
int IsTagged(struct LinkedList *);
long long CountUniqueStars(struct HashTable *,long long *);
void ExtractFinalHalo(long long *,struct HashTable *,struct tagged_particle *);
void WriteFinalTag(struct tagged_particle *, long long);

// galaxy tag
struct GalHashTable *GalEmptyTable(size_t);
struct GalLinkedList *GalNewLinkedList(void);
void GalInsertKey(struct GalHashTable *,long long,struct tagged_particle *);
bool GalContainsElement (struct GalLinkedList *,long long);
void GalAddElement(struct GalLinkedList *,long long,struct tagged_particle *);
struct GalLinkedList *GalGetPreviousLink(struct GalLinkedList *,long long); // I removed static
void GalDeleteTable(struct GalHashTable *);
void GalDeleteLinkedList(struct GalLinkedList *);
