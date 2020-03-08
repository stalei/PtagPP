// © Shahram @ 2019
/*! \file GlobalVars.c
 *  \brief declares global variables.
 * 
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 * 
 *     - Erase all #define statements
 *     - add #include "allvars.h"
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */
#include "GlobalVars.h"

char CurrentFile[100];
char ParametersFile[100];
struct GlobalParameters GP;

void *CommBuffer;		/*!< points to communication buffer, which is used at a few places */

int TargetSnap;

struct Path_Names *SageFilesPath,*SageFilesPathPre,*TagFilesPath,*TagFilesPathPre;
int SageFilesCount,SageFilesCountPre;
int NumGalaxies;
int NumGalaxiesPre;
struct SageGalaxies *SageOutput,*SageOutputPre;
struct tagged_particle tagged_p,*AllStars,*AllStarsPre,*StellarHaloAllSnaps;
int TagFilesCount,TagFilesCountPre;
long int TaggedParticlesCount;
long int NumOfStars, NumOfStarsPre;
hid_t GTagDatatype;
int TagCols;
int TagSize;

struct io_header header;	/*!< holds header for snapshot files */

/*! This structure holds all the information that is
 *  * stored for each particle of the simulation.
 *   */
struct particle_data *P,	/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */

/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 *  * variables.
 *   */
struct sph_particle_data *SphP,	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;	

#ifdef DoParallel
int ThisTask;			/*!< the number of the local processor  */
int NTask;			/*!< number of processors */
MPI_Comm MPI_CommLocal;
#endif

struct HashTable *table;
//struct HashTable *EmptyTable(size_t);
struct LinkedList;
//struct LinkedList *NewLinkedList();
struct GalHashTable *gtable;
struct GalLinkedList;

