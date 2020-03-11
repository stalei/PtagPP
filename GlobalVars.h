// © Shahram Talei @ 2019
/*! \file GlobalVars.h
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
#ifndef ALLVARS_H
#define ALLVARS_H

#include<stdio.h>
#include "hdf5.h"


#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif


#ifndef DOUBLEPRECISION     /* default is single-precision */
typedef float  MyFloat;
typedef float  MyDouble;
#else
#if (DOUBLEPRECISION+0) == 2 
typedef float   MyFloat;
typedef double  MyDouble;
#else                        /* everything double-precision */
typedef double  MyFloat;
typedef double  MyDouble;
#endif
#endif


#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif

/*Determines the maximum size of arrays related to the number of CR populations */
#ifndef NUMCRPOP   /*!< Number of CR populations pressent in parameter file */
#define NUMCRPOP 1
#endif

#ifndef GDE_TYPES 
#define GDE_TYPES 2
#endif

#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif


#ifdef   ENLARGE_DYNAMIC_RANGE_IN_TIME
typedef  long long integertime;
#define  TIMEBINS        60
#define  TIMEBASE        (((long long)1)<<TIMEBINS)     /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                                         *   where TIMESPAN needs to be a power of 2.
                                                         *                                                            */
#else
typedef  int integertime;
#define  TIMEBINS        29
#define  TIMEBASE        (1<<TIMEBINS)  /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                         *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
                                         *                                            *   to 2^29
                                         *                                                                                     */
#endif

#ifdef FLTROUNDOFFREDUCTION
#define FLT(x) ((MyFloat)(x))
#ifdef SOFTDOUBLEDOUBLE      /* this requires a C++ compilation */
#include "dd.h"
typedef dd MyLongDouble;
#else
typedef long double MyLongDouble;
#endif
#else  /* not enabled */
#define FLT(x) (x)
typedef MyFloat MyLongDouble;
#endif  /* end FLTROUNDOFFREDUCTION */

#ifndef  GRAVCOSTLEVELS
#define  GRAVCOSTLEVELS      6
#endif

/* some flags for the field "flag_ic_info" in the file header */
#define FLAG_ZELDOVICH_ICS     1
#define FLAG_SECOND_ORDER_ICS  2
#define FLAG_EVOLVED_ZELDOVICH 3
#define FLAG_EVOLVED_2LPT      4
#define FLAG_NORMALICS_2LPT    5


// compiler specific data alignment hints
// // XLC compiler
#if defined(__xlC__)
#define ALIGN(n) __attribute__((__aligned__(n)))
// // GNU compiler 
#elif defined(__GNUC__)
#define ALIGN(n) __attribute__((__aligned__(n)))
// // Intel Compiler
#elif defined(__INTEL_COMPILER)
// // GNU Intel Compiler
#define ALIGN(n) __declspec(align(n))
// // Unknown Compiler
#else
#define ALIGN(n) 
#endif

#if defined (BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING) || defined(LT_STELLAREVOLUTION)
#define PPP P
#else
#define PPP SphP
#endif

extern char CurrentFile[100];
extern char ParametersFile[100];


extern int TargetSnap;

extern void *CommBuffer;	/*!< points to communication buffer, which is used at a few places */

extern struct GlobalParameters
{
int FirstSnap; 
int LastSnap;
long int TimeLimit;
char SageDir[500];
char TagDir[500];
char SnapDir[500];
char OutputDir[500];
double f_mb;
int BufferSize;
double Omega0;
double OmegaLambda;
double OmegaBaryon;
double HubbleParam;
//int SnapFormat;
////
int CoolingOn;
int StarformationOn;
long long TotNumPart;
long long TotNumStars;
long long TotNumStarsAllSnaps;
long long TotNumTagsAllSnaps;
  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
 *    *  value, * allowing the size of the snapshot files to be reduced
 *       */
  double MassTable[6];

} 
GP;

extern struct Path_Names
{
    char paths[100];
}*SageFilesPath,*SageFilesPathPre,*TagFilesPath,*TagFilesPathPre;

extern int SageFilesCount,SageFilesCountPre;
extern int NumGalaxies;
extern int NumGalaxiesPre;

// This structune holds all of the information from the Sage files
extern struct SageGalaxies
{
  int   Type;
  int   FileNr;
  long long   GalaxyIndex;
  int   HaloIndex;
  int   FOFHaloIndex;
  int   TreeIndex;
  
  int   SnapNum;
  int   CentralGal;
  float CentralMvir;

   //properties of subhalo at the last time this galaxy was a central galaaxy 
  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;   
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float VelDisp;

   //baryonic reservoirs 
  float ColdGas;
  float StellarMass;
  float BulgeMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float ICS;

   //metals
  float MetalsColdGas;
  float MetalsStellarMass;
  float MetalsBulgeMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsICS;
 
   //misc 
  float Sfr;
  float SfrBulge;
  float SfrICS;
  float DiskScaleRadius;
  float Cooling;
  float Heating;
  float LastMajorMerger;
  float OutflowRate;
  
  float infallMvir;  //infall properties
  float infallVvir;
  float infallVmax;
  float r_heat;

}*SageOutput,*SageOutputPre;

extern struct tagged_particle
{
double Pos[3]; // float or ALIGN(32) MyDouble Pos[3];
double Vel[3]; // same as Pos
long long ID; //MyIDType ID;
int Snap;
long long PID; //unsigned int PID;
int HaloIndex;
int SubhaloIndex;
int GalIndex;
int GalNo;
int TreeIndex;
float AA;
float Age;
float ZZ;
double StellarMass;
double Time;
int Len;
//double Pos[3];
////double Vel[3];
long long MBID; //most bound ID  long long  MostBoundID;
double BindingEnergy; //Myfloat
double Mvir;
double Rvir;
double infallMvir;
float LastMajorMerger;
}tagged_p,*AllStars,*AllStarsPre,*StellarHaloAllSnaps;//,*tagged_list;
//


extern int TagFilesCount, TagFilesCountPre;
extern long int TaggedParticlesCount;
extern long int NumOfStars, NumOfStarsPre;
extern hid_t GTagDatatype;
extern int TagCols;
extern int TagSize;




/*! Header for the standard file format.
 *  */
extern struct io_header
{
  int npart[6];			/*!< number of particles of each type in this file */
  double mass[6];		/*!< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;			/*!< time of snapshot file */
  double redshift;		/*!< redshift of snapshot file */
  int flag_sfr;			/*!< flags whether the simulation was including star formation */
  int flag_feedback;		/*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
				   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;		/*!< flags whether cooling was included  */
  int num_files;		/*!< number of files in multi-file snapshot */
  double BoxSize;		/*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;		/*!< matter density in units of critical density */
  double OmegaLambda;		/*!< cosmological constant parameter */
  double HubbleParam;		/*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;		/*!< flags whether the file contains formation times of star particles */
  int flag_metals;		/*!< flags whether the file contains metallicity values for gas and star
				   particles */
  unsigned int npartTotalHighWord[6];	/*!< High word of the total number of particles of each type */
  int flag_entropy_instead_u;	/*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;	/*!< flags that snapshot contains double-precision instead of single precision */

  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                     or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                     For snapshots files, the value informs whether the simulation was evolved from
                                     Zeldoch or 2lpt ICs. Encoding is as follows:
                                        FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                        FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                        FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                        FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                        FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */

  char fill[18];		/*!< fills to 256 Bytes */

  char names[15][2];
}
header;	

enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_SECONDORDERMASS,
  IO_U,
  IO_RHO,
  IO_NE,
  IO_NH,
  IO_HSML,
  IO_VALPHA,
  IO_SFR,
  IO_AGE,
  IO_HSMS,
  IO_ACRS,
  IO_Z,
  IO_BHMASS,
  IO_BHMDOT,
  IO_BHPROGS,
  IO_BHMBUB,
  IO_BHMINI,
  IO_BHMRAD,
  IO_ACRB,
  IO_POT,
  IO_ACCEL,
  IO_CR_C0,
  IO_CR_Q0,
  IO_CR_P0,
  IO_CR_E0,
  IO_CR_n0,
  IO_CR_ThermalizationTime,
  IO_CR_DissipationTime,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_HeIII,
  IO_H2I,
  IO_H2II,
  IO_HM,
  IO_HD,
  IO_DI,
  IO_DII,
  IO_HeHII,
  IO_DTENTR,
  IO_STRESSDIAG,
  IO_STRESSOFFDIAG,
  IO_STRESSBULK,
  IO_SHEARCOEFF,
  IO_TSTP,
  IO_BFLD,
  IO_BSMTH,
  IO_DBDT,
  IO_DIVB,
  IO_ABVC,
  IO_AMDC,
  IO_VTURB,
  IO_LTURB,
  IO_ALFA2_DYN,
  IO_ETA2_DYN,
  IO_PHI,
  IO_XPHI,
  IO_GRADPHI,
  IO_ROTB,
  IO_SROTB,
  IO_COOLRATE,
  IO_CONDRATE,
  IO_DENN,
  IO_EGYPROM,
  IO_EGYCOLD,
  IO_MACH,
  IO_MACH1,
  IO_MACH2,
  IO_MACH3,
  IO_RUP,
  IO_RDOWN,
  IO_PUP,
  IO_PDOWN,
  IO_SHSPEED,
  IO_SHCOMPRESS,
  IO_SHNORMAL,
  IO_VUP,
  IO_VDOWN,
  IO_DTENERGY,
  IO_PRESHOCK_CSND,
  IO_PRESHOCK_DENSITY,
  IO_PRESHOCK_ENERGY,
  IO_PRESHOCK_XCR,
  IO_DENSITY_JUMP,
  IO_ENERGY_JUMP,
  IO_CRINJECT,
  IO_TIDALTENSORPS,
  IO_DISTORTIONTENSORPS,
  IO_EULERA,
  IO_EULERB,
  IO_FLOW_DETERMINANT,
  IO_PHASE_SPACE_DETERMINANT,
  IO_ANNIHILATION_RADIATION,
  IO_STREAM_DENSITY,
  IO_EOSTEMP,
  IO_EOSXNUC,
  IO_PRESSURE,
  IO_RADGAMMA,
  IO_RAD_ACCEL,
  IO_EDDINGTON_TENSOR,
  IO_LAST_CAUSTIC,
  IO_SHEET_ORIENTATION,
  IO_INIT_DENSITY,
  IO_CAUSTIC_COUNTER,
  IO_DMHSML,                    /* for 'SUBFIND_RESHUFFLE_CATALOGUE' option */
  IO_DMDENSITY,
  IO_DMVELDISP,
  IO_DMHSML_V,                 /* for 'SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI' option */
  IO_DMDENSITY_V,
  IO_VRMS,
  IO_VBULK,
  IO_VRAD,
  IO_VTAN,
  IO_TRUENGB,
  IO_VDIV,
  IO_VROT,
  IO_VORT,
  IO_DPP,
  IO_BPCR_pNORM,
  IO_BPCR_eNORM,
  IO_BPCR_pSLOPE,
  IO_BPCR_eSLOPE,
  IO_BPCR_pE,
  IO_BPCR_pN,
  IO_BPCR_ePRESSURE,
  IO_BPCR_pPRESSURE,

  IO_iMass,
  IO_Zs,
  IO_ZAGE,
  IO_ZAGE_LLV,
  IO_CLDX,
  IO_HTEMP,
  IO_TEMP,
  IO_CONTRIB,
  IO_ZSMOOTH,
  IO_allZSMOOTH,
  IO_CHEM,
  IO_DELAYTIME,

  IO_ABUNDANCE,
  IO_ABND_GRAD,
  IO_ABND_SPH,
  IO_DIFFUSING_CB,
  IO_DIFFUSING_CA,
  IO_DIFFUSING_CD,
  IO_DIFF_TSTEP,
  
  IO_PSUM,
  IO_SIDMNUMNGB, 
  IO_NUMTOTALSCATTER,
  IO_SIDMHSML,
  IO_SIDMDENSITY,
  IO_SIDMVELDISP,

  IO_AGS_SOFT,
  IO_AGS_DENS,
  IO_AGS_ZETA,
  IO_AGS_OMEGA,
  IO_AGS_CORR,
  IO_AGS_NGBS,

  IO_VSTURB_DISS,
  IO_VSTURB_DRIVE,
  
  IO_MG_PHI,
  IO_MG_ACCEL,
  
  IO_LASTENTRY			/* This should be kept - it signals the end of the list */
};


/*! This structure holds all the information that is
 *  * stored for each particle of the simulation.
 *   */
extern ALIGN(32) struct particle_data
{
  short int Type;		/*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  short int TimeBin;
  MyIDType ID;

  integertime Ti_begstep;		/*!< marks start of current timestep of particle on integer timeline */
  integertime Ti_current;		/*!< current time of the particle */

  ALIGN(32) MyDouble Pos[3];   /*!< particle position at its current time */
  MyDouble Mass;     /*!< particle mass */
  MyDouble InitMass;
  
  MyDouble Vel[3];   /*!< particle velocity at its current time */
  MyDouble dp[3];

  union
  {
    MyFloat       GravAccel[3];		/*!< particle acceleration due to gravity */
    MyLongDouble dGravAccel[3];
  } g;
#ifdef PMGRID
  MyFloat GravPM[3];		/*!< particle acceleration due to long-range PM gravity force */
#endif
#ifdef FORCETEST
  MyFloat GravAccelDirect[3];	/*!< particle acceleration calculated by direct summation */
#endif
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
  union
  {
    MyFloat       Potential;		/*!< gravitational potential */
    MyLongDouble dPotential;
  } p;
#endif

#ifdef FS_TURB_ESTIM
  MyFloat StrFnc[FS_BINS];
  long int StrFnc_count[FS_BINS];
#endif
#ifdef MODGRAV
  MyDouble MassOrig;  /* the mass that particles would have without mass correction */ 
  MyDouble mg_phi;    /* value of the scalar field fluctuation */
  MyDouble mg_grad_phi[3];    /* value of the scalar field fluctuation */
  MyDouble mg_beta;   /* local scalar coupling */
  MyDouble ModGravAccel[3];  /* acceleration due to scalar field fluctuations */
  int mg_treelevel;   /* level of parent AMR node */
#endif

#ifdef DISTORTIONTENSORPS
  MyBigFloat distortion_tensorps[6][6];               /*!< phase space distortion tensor */
  MyBigFloat last_determinant;                        /*!< last real space distortion tensor determinant */
  MyBigFloat stream_density;                          /*!< physical stream density that is going to be integrated */
  double tidal_tensorps[3][3];                        /*!< tidal tensor (=second derivatives of grav. potential) */
  float caustic_counter;                              /*!< caustic counter */
#ifndef GDE_LEAN
  MyBigFloat annihilation;                            /*!< integrated annihilation rate */
  MyBigFloat analytic_annihilation;                   /*!< analytically integrated annihilation rate */
  MyBigFloat rho_normed_cutoff_current;               /*!< current and last normed_cutoff density in rho_max/rho_init * sqrt(sigma) */
  MyBigFloat rho_normed_cutoff_last;
  MyBigFloat s_1_current, s_2_current, s_3_current;   /*! < current and last stretching factor */
  MyBigFloat s_1_last, s_2_last, s_3_last;
  MyBigFloat second_deriv_current;                    /*! < current and last second derivative */
  MyBigFloat second_deriv_last;
  double V_matrix[3][3];                              /*!< initial orientation of CDM sheet the particle is embedded in */
  float init_density;                                 /*!< initial stream density */
  float analytic_caustics;                            /*!< number of caustics that were integrated analytically */
  float a0;
#endif

#ifdef OUTPUT_LAST_CAUSTIC
  MyFloat lc_Time;                                  /*!< time of caustic passage */
  MyFloat lc_Pos[3];                                /*!< position of caustic */
  MyFloat lc_Vel[3];                                /*!< particle velocity when passing through caustic */
  MyFloat lc_rho_normed_cutoff;                     /*!< normed_cutoff density at caustic */
  MyFloat lc_Dir_x[3];                              /*!< principal axis frame of smear out */
  MyFloat lc_Dir_y[3];
  MyFloat lc_Dir_z[3];
  MyFloat lc_smear_x;                               /*!< smear out length */
  MyFloat lc_smear_y;
  MyFloat lc_smear_z;
#endif
#ifdef PMGRID
  double tidal_tensorpsPM[3][3];	            /*!< for TreePM simulations, long range tidal field */
#endif
#endif

  MyFloat OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening
                                          criterion */
#if defined(EVALPOTENTIAL) && defined(PMGRID)
  MyFloat PM_Potential;
#endif

#ifdef STELLARAGE
  MyFloat StellarAge;		/*!< formation time of star particle */
#endif
#ifdef METALS
  MyFloat Metallicity;		/*!< metallicity of gas or star particle */
#endif				/* closes METALS */

#if defined (BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING) || defined(LT_STELLAREVOLUTION) || defined(LT_MV_CHEMICALDIFFUSION) || defined(UM_CHEMISTRY)
  MyFloat Hsml;

  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#if defined(RADTRANSFER) || defined(SNIA_HEATING)
  MyFloat DensAroundStar;
#endif
#endif


#if defined(BLACK_HOLES)
  MyIDType SwallowID;
#if !defined(DETACH_BLACK_HOLES)
#ifdef BH_COUNTPROGS
  int BH_CountProgs;
#endif
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
  int     BH_TimeBinGasNeighbor;
#ifdef BH_BUBBLES
  MyFloat BH_Mass_bubbles;
  MyFloat BH_Mass_ini;
#ifdef UNIFIED_FEEDBACK
  MyFloat BH_Mass_radio;
#endif
#endif
  union
  {
    MyFloat BH_Density;
    MyLongDouble dBH_Density;
  } b1;
  union
  {
    MyFloat BH_Entropy;
    MyLongDouble dBH_Entropy;
  } b2;
  union
  {
    MyFloat BH_SurroundingGasVel[3];
    MyLongDouble dBH_SurroundingGasVel[3];
  } b3;
  union
  {
    MyFloat BH_accreted_Mass;
    MyLongDouble dBH_accreted_Mass;
  } b4;
  union
  {
    MyFloat BH_accreted_BHMass;
    MyLongDouble dBH_accreted_BHMass;
  } b5;
  union
  {
    MyFloat BH_accreted_momentum[3];
    MyLongDouble dBH_accreted_momentum[3];
  } b6;
#ifdef BH_BUBBLES
  union
  {
    MyFloat BH_accreted_BHMass_bubbles;
    MyLongDouble dBH_accreted_BHMass_bubbles;
  } b7;
#ifdef UNIFIED_FEEDBACK
  union
  {
    MyFloat BH_accreted_BHMass_radio;
    MyLongDouble dBH_accreted_BHMass_radio;
  } b8;
#endif
#endif
#ifdef KD_FRICTION
  MyFloat BH_SurroundingVel[3];
  MyFloat BH_SurroundingDensity;
#endif
#ifdef KD_FRICTION_DYNAMIC
  MyFloat BH_sigma;
  MyFloat BH_bmax;
#endif
#ifdef KD_TAKE_CENTER_OF_MASS_FOR_BH_MERGER
  MyFloat BH_SwallowPos[3];
#endif
#ifdef LT_DF_BH_BHAR_SWITCH
  MyFloat BlackHoleFeedbackFactor;
#endif
#ifdef LT_BH_GUESSHSML
  MyFloat mean_hsml;
  MyFloat mean_rho;
#endif
#ifdef REPOSITION_ON_POTMIN
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
#ifdef BH_KINETICFEEDBACK
  MyFloat ActiveTime;
  MyFloat ActiveEnergy;
#endif
#endif  /* if !defined(DETACH_BLACK_HOLES) */
#endif  /* if deffined(BLACK_HOLES) */ 

#ifdef SIDM
  int ShouldScatterInStep;
  int DidScatterInStep;
  MyIDType ScatterID;
  MyDouble RandX; 
  MyDouble sidm_PSum;
  MyDouble sidm_NumTotalScatter;
  MyDouble sidm_Hsml;
  MyDouble sidm_Density;
  MyDouble sidm_VelDisp;
  MyDouble sidm_NumNgb;
#endif

  
#if defined(DISKPOT) && !defined(SUBFIND) && !defined(ORDER_SNAPSHOTS_BY_ID)
	int GrNr;
#endif
  
#if defined(SUBFIND) 
  int GrNr;
  int SubNr;
  int DM_NumNgb;
  unsigned short targettask, origintask2;
  int origintask, submark, origindex;
  MyFloat DM_Hsml;
  union
  {
    MyFloat DM_Density;
    MyFloat DM_Potential;
  } u;
  union
  {
    MyFloat DM_VelDisp;
    MyFloat DM_BindingEnergy;
  } v;
#ifdef DENSITY_SPLIT_BY_TYPE
  union
  {
    MyFloat int_energy;
    MyFloat density_sum;
  } w;
#endif

#ifdef SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI
  MyFloat DM_Hsml_V;
  MyFloat DM_Density_V;
#endif

#ifdef SAVE_HSML_IN_IC_ORDER
  MyIDType ID_ic_order;
#endif
#ifdef SUBFIND_ALTERNATIVE_COLLECTIVE
  peanokey Key;
#endif
#endif

#if defined(ORDER_SNAPSHOTS_BY_ID) && !defined(SUBFIND)
  int     GrNr;
  int     SubNr;
#endif


#ifdef CS_MODEL
  MyFloat Zm[12];
  MyFloat ZmReservoir[12];
#ifdef CS_FEEDBACK
  MyFloat EnergySN;
  MyFloat EnergySNCold;
#endif
#endif

  float GravCost[GRAVCOSTLEVELS];   /*!< weight factor used for balancing the work-load */

#ifdef WAKEUP
  int dt_step;
#endif

#if defined(LT_STELLAREVOLUTION) || defined(DETACH_BLACK_HOLES)
  union
  {
    unsigned int BHID;
    unsigned int MetID;
  } pt;
#endif

#ifdef SCF_HYBRID
  MyDouble GravAccelSum[3];
  MyFloat MassBackup;
#endif

#ifdef MOL_CLOUDS
  MyFloat MOL_CLOUDS_TimeBorn;
  MyFloat MOL_CLOUDS_LifeTime;

  unsigned int MOL_CLOUDS_index;
#endif

#if defined(JD_VTURB) || defined(HIGH_ORDER_INDUCTION) || defined(KD_RESTRICT_NEIGHBOURS)
  int TrueNGB;			/*!< Number of neighbours inside hsml */
#endif

#ifdef ADAPTGRAVSOFT
  MyFloat AGS_Density;		/* !< (mass/number) density of particle */
#ifdef AGS_OUTPUTNGBS
  int AGS_REALNumNgb;
#endif
  MyFloat AGS_NumNgb;
  MyFloat AGS_Hsml;
  MyFloat AGS_zeta;           /*!< factor in the correction term */
  MyFloat AGS_omega;          /*!< factor in the correction term */
#ifdef AGS_OUTPUTCORR
  MyFloat AGS_corr;          /*!< factor in the correction term */
#endif
#endif

#ifdef PARTICLE_TAGGING
 double StellarAge;
 double Metallicity;
 double StellarMass;
#endif

}
 *P,				/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */



/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 *  * variables.
 *   */
extern struct sph_particle_data
{
  MyDouble Entropy;		/*!< entropy (actually entropic function) of particle */
  MyFloat  EntropyPred;         /*!< predicted value of the entropy at the current time */
  MyFloat  Pressure;		/*!< current pressure */
  MyFloat  VelPred[3];		/*!< predicted SPH particle velocity at the current time */
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat MinViscousDt;
#else
  MyFloat MaxSignalVel;           /*!< maximum signal velocity */
#endif

#ifdef VORONOI
  MyFloat MaxDelaunayRadius;
  MyFloat Volume;
  MyFloat Center[3];
#ifdef VORONOI_SHAPESCHEME
  MyFloat W;
#endif
#endif

  union
  {
    MyFloat       Density;		/*!< current baryonic mass density of particle */
    MyLongDouble dDensity;
  } d;
  union
  {
    MyFloat       DtEntropy;		/*!< rate of change of entropy */
    MyLongDouble dDtEntropy;
  } e;
  union
  {
    MyFloat       HydroAccel[3];	/*!< acceleration due to hydrodynamical force */
    MyLongDouble dHydroAccel[3];
  } a;
  union
  {
    MyFloat       DhsmlDensityFactor;	/*!< correction factor needed in entropy formulation of SPH */
    MyLongDouble dDhsmlDensityFactor;
  } h;
  union
  {
    MyFloat       DivVel;		/*!< local velocity divergence */
    MyLongDouble dDivVel;
  } v;
#ifndef NAVIERSTOKES
  union
  {
    MyFloat CurlVel;     	        /*!< local velocity curl */
    MyFloat       Rot[3];		/*!< local velocity curl */
    MyLongDouble dRot[3];
  } r;
#else
  union
  {
    MyFloat DV[3][3];
    struct
    {
      MyFloat DivVel;
      MyFloat CurlVel;
      MyFloat StressDiag[3];
      MyFloat StressOffDiag[3];
#ifdef NAVIERSTOKES_BULK
      MyFloat StressBulk;
#endif
    } s;
  } u;
#endif

#if defined(AB_TURB) || defined(VS_TURB) || defined(OUTPUT_VORTICITY) 
   MyFloat Vorticity[3];
   MyFloat SmoothedVel[3];
#endif

#if !(defined(BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING) || defined(LT_STELLAREVOLUTION)) || defined(UM_CHEMISTRY)
  MyFloat Hsml;			/*!< current smoothing length */
  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
  union
  {
    MyFloat       Injected_BH_Energy;
    MyLongDouble dInjected_BH_Energy;
  } i;
#endif

#ifdef COOLING
#if !defined(UM_CHEMISTRY)
  MyFloat Ne;  /*!< electron fraction, expressed as local electron number
		    density normalized to the hydrogen number density. Gives
		    indirectly ionization state and mean molecular weight. */
#endif
#endif
#ifdef SFR
  MyFloat Sfr;
#endif
#ifdef WINDS
  MyFloat DelayTime;		/*!< remaining maximum decoupling time of wind particle */
#ifdef VARIABLE_WINDS
  MyFloat HostHaloMass;
#endif
#endif

#ifdef JD_VTURB
  MyFloat Vrms;		    /*!< RMS velocity inside kernel around Vbulk */
  MyFloat Vbulk[3];	    /*!< Mean velocity inside kernel */
#ifdef JD_DECOMPOSE_VTURB 
  MyFloat Vtan;         /*!< RMS of tangential component of velocity inside kernel around Vbulk */
  MyFloat Vrad;         /*!< RMS of radial component of velocity inside kernel around Vbulk */
#endif
#ifdef JD_DPP
  MyFloat Dpp;			/*!< Reacceleration Coefficient as (Cassano+ '04) */
#endif
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  MyDouble DuDt_diss;
  MyDouble DuDt_drive;
  MyDouble EgyDiss;
  MyDouble EgyDrive;
  MyDouble TurbAccel[3];
#endif
#if defined(TURB_DRIVING)
  MyDouble TurbAccel[3];
#endif

#ifdef AB_SHOCK
  MyFloat Shock_N1[3];           /*!< Parallel to shock >*/
  MyFloat Shock_Up_V1[3];        /*!< Upstream parallel velocity >*/
  MyFloat Shock_Down_V1[3];      /*!< Downstream parallel velocity >*/
  MyFloat Shock_Up_Signal;       /*!< Upstream signal velocity >*/
  MyFloat Shock_Down_Signal;     /*!< Downstream signal velocity >*/
  MyFloat Shock_Mach;            /*!< Machnumber >*/
  MyFloat Shock_Mach1;            /*!< Machnumber >*/
  MyFloat Shock_Mach2;            /*!< Machnumber >*/
  MyFloat Shock_Mach3;            /*!< Machnumber >*/
  MyFloat Shock_Speed;           /*!< Shock speed >*/
  MyFloat Shock_Compress;        /*!< Compression ratio >*/
  MyFloat Shock_Up_Weight;       /*!< Upstream weight >*/
  MyFloat Shock_Down_Weight;     /*!< Downstream weight >*/
  MyFloat Shock_Up_Rho;          /*!< Upstream density >*/
  MyFloat Shock_Down_Rho;        /*!< Downstream density >*/
  MyFloat Shock_Up_Pressure;     /*!< Upstream pressure >*/
  MyFloat Shock_Down_Pressure;   /*!< Downstream pressure >*/
#ifdef AB_SHOCK_VELDIV
  MyFloat Shock_N2[3];           /*!< Perpendicular to shock >*/
  MyFloat Shock_N3[3];           /*!< Perpendicular to shock >*/
  MyFloat Shock_Up_V2[3];        /*!< Upstream perpendicular velocity >*/
  MyFloat Shock_Down_V2[3];      /*!< Downstream perpendicular velocity >*/
  MyFloat Shock_Up_V3[3];        /*!< Upstream perpendicular velocity >*/
  MyFloat Shock_Down_V3[3];      /*!< Downstream perpendicular velocity >*/
  MyFloat Shock_Up_Weight23[2];  /*!< Upstream perpendicular weights >*/
  MyFloat Shock_Down_Weight23[2];/*!< Downstream perpendicular weights >*/
#endif
#endif /* AB_SHOCK */

#ifdef VISCOSITY_SUPPRESSION
  MyFloat NV_R;
  MyFloat NV_DivVel;
  MyFloat NV_dt_DivVel;
  MyFloat NV_A[3][3];
  MyFloat NV_D[3][3];
  MyFloat NV_T[3][3];
  MyFloat NV_trSSt;
  MyFloat alphaloc, alpha;
#endif

#ifdef FS_TURB_ESTIM
  MyFloat Vturb;
  MyFloat Lturb;
#endif
#ifdef MAGNETIC
  union
  {
    MyFloat B[3];
#if defined(EULERPOTENTIALS)
    MyFloat dEulerA[3];
#endif
  } b1;

  union
  {  
    MyFloat BPred[3];
#if defined(EULERPOTENTIALS)
    MyFloat dEulerB[3];
#endif
  } b2;

#ifdef MAGNETIC_SN_SEEDING
  MyFloat MagSeed[3];
#endif

#ifdef HIGH_ORDER_INDUCTION
  MyFloat Chi[3][3];
  MyFloat Xix[3], Xiy[3], Xiz[3];
#endif
#ifdef FS_ALFA2_DYN
  MyFloat alfa2;
#endif
#ifdef FS_ETA2_DYN
  MyFloat eta2;
#endif
#ifdef DIVBFORCE3
  MyFloat magacc[3];
  MyFloat magcorr[3];
#endif
#ifdef EULERPOTENTIALS
  MyFloat EulerA,EulerB;
#ifdef EULER_DISSIPATION
  MyFloat DtEulerA,DtEulerB;
#endif
#endif
#if !defined(EULERPOTENTIALS)
  MyFloat DtB[3];
#endif
#if defined(TRACEDIVB) || defined(TIME_DEP_MAGN_DISP) || defined(DIVBCLEANING_DEDNER)
  MyFloat divB;
#endif
#if defined(BSMOOTH) || defined(BFROMROTA) || defined(BSMOOTH_TIME)
  MyFloat BSmooth[3];
#endif
#ifdef TIME_DEP_MAGN_DISP
  MyFloat Balpha, DtBalpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  MyFloat Phi, PhiPred, DtPhi;
  MyFloat GradPhi[3];
#ifdef SMOOTH_PHI
  MyFloat SmoothPhi;
#endif
#endif
#if defined(SMOOTH_DIVB) || defined(DIVBCLEANING_DEDNER)
  MyFloat SmoothDivB;
#endif

#if defined(ROT_IN_MAG_DIS) || defined(OUTPUT_ROTB) || defined(FS_ETA2_DYN)
  MyFloat RotB[3];
#ifdef SMOOTH_ROTB
  MyFloat SmoothedRotB[3];
#endif
#endif

#endif /* MAGNETIC */

#ifdef VSMOOTH
  MyFloat VSmooth[3];
  int SmoothNgb;
#endif

#if (defined(SMOOTH_DIVB) || defined(DIVBCLEANING_DEDNER) || defined(BSMOOTH) || defined(SMOOTH_ROTB) || ((defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)) || defined(LT_SMOOTH_Z) || defined(LT_SMOOTH_XCLD) || defined(LT_TRACK_WINDS))) || defined(VSMOOTH)
  MyFloat DensityNorm;
#endif

#ifdef TIME_DEP_ART_VISC
  MyFloat alpha, Dtalpha;
#endif
#ifdef NS_TIMESTEP
  MyFloat ViscEntropyChange;
#endif
#ifdef CONDUCTION_SATURATION
  MyFloat GradEntr[3];
#endif

#ifdef MHM
  MyFloat FeedbackEnergy;
#endif

#ifdef COSMIC_RAYS
  MyFloat CR_C0[NUMCRPOP];			/*!< Cosmic ray amplitude adiabatic invariable */
  MyFloat CR_q0[NUMCRPOP];			/*!< Cosmic ray cutoff adiabatic invariable */
  MyFloat CR_E0[NUMCRPOP];			/*!< Specific Energy at Rho0 */
  MyFloat CR_n0[NUMCRPOP];			/*!< baryon fraction in cosmic rays */

  MyFloat CR_DeltaE[NUMCRPOP];		/*!< Specific Energy growth during timestep */
  MyFloat CR_DeltaN[NUMCRPOP];		/*!< baryon fraction growth during timestep */
#ifdef MACHNUM
  MyFloat CR_Gamma0[NUMCRPOP];
#endif

#ifdef CR_OUTPUT_INJECTION
  MyFloat CR_Specific_SupernovaHeatingRate;
#endif
#endif				/* COSMIC_RAYS */

#ifdef MACHNUM
  MyFloat Shock_MachNumber;	/*!< Mach number */
  MyFloat Shock_DecayTime;	/*!< Shock decay time */
#ifdef COSMIC_RAYS
  MyFloat Shock_DensityJump;	/*!< Density jump at the shock */
  MyFloat Shock_EnergyJump;	/*!< Energy jump at the shock */
  MyFloat PreShock_PhysicalDensity;	/*!< Specific energy in the preshock regime */
  MyFloat PreShock_PhysicalEnergy;	/*!< Density in the preshock regime */
  MyFloat PreShock_XCR;		/*!< XCR = PCR / Pth in the preshock regime */
#endif
#ifdef MACHSTATISTIC
  MyFloat Shock_DtEnergy;		/*!< Change of thermal specific energy at Shocks */
#endif
#ifdef OUTPUT_PRESHOCK_CSND
  MyFloat PreShock_PhysicalSoundSpeed;	/*!< Sound speed in the preshock regime */
  MyFloat PreShock_PhysicalDensity;	/*!< Specific energy in the preshock regime */
#endif
#endif				/* Mach number estimate */


#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
  MyFloat elec;
  MyFloat HI;
  MyFloat HII;

  MyFloat HeI;
  MyFloat HeII;
  MyFloat HeIII;

  MyFloat H2I;
  MyFloat H2II;

  MyFloat HM;

  MyFloat Gamma;
  MyFloat t_elec, t_cool;

#ifdef UM_CHEMISTRY
  MyFloat Um_MeanMolecularWeight;
  MyFloat HD;
  MyFloat DI;
  MyFloat DII;
  MyFloat HeHII;
#endif

#endif

#if defined(RADTRANSFER)
  MyFloat ET[6];                /* eddington tensor - symmetric -> only 6 elements needed */
  MyFloat Je[N_BINS];           /* emmisivity */

#ifndef UM_CHEMISTRY
  MyFloat HI;                  /* HI fraction */
  MyFloat HII;                 /* HII fraction */
  MyFloat elec;               /* electron fraction */
#ifdef RT_INCLUDE_HE
  MyFloat HeI;                 /* HeI fraction */
  MyFloat HeII;                 /* HeII fraction */
  MyFloat HeIII;                 /* HeIII fraction */
#endif
#endif

  MyFloat n_gamma[N_BINS];

#ifdef RADTRANSFER_FLUXLIMITER
  MyFloat Grad_ngamma[3][N_BINS];
#endif

#ifdef RT_RAD_PRESSURE
  MyFloat n[3];
  MyFloat RadAccel[3];
#endif

#ifdef SFR
  MyDouble DensitySfr;
  MyDouble HsmlSfr;
  MyDouble DhsmlDensityFactorSfr;
  MyDouble NgbSfr;
#endif

#endif //radtransfer 

#if defined CS_MODEL
  MyFloat DensityOld;
#ifdef CS_FEEDBACK
  union
  {
    MyFloat       DensityAvg;		/*!< current baryonic mass density of particle */
    MyLongDouble dDensityAvg;
  } da;
  union
  {
    MyFloat       EntropyAvg;		/*!< current baryonic mass density of particle */
    MyLongDouble dEntropyAvg;
  } ea;
  MyFloat HotHsml;
  int     HotNgbNum;
  MyFloat DensPromotion;
  MyFloat TempPromotion;
#endif
#endif

#ifdef EOS_DEGENERATE
  MyFloat u;                            /* internal energy density */
  MyFloat temp;                         /* temperature */
  MyFloat dpdr;							/* derivative of pressure with respect to density at constant entropy */
  MyFloat xnuc[EOS_NSPECIES];           /* nuclear mass fractions */
  MyFloat dxnuc[EOS_NSPECIES];          /* change of nuclear mass fractions */
  MyFloat xnucPred[EOS_NSPECIES];
#endif

#ifdef WAKEUP
  short int wakeup;             /*!< flag to wake up particle */
#endif

#ifdef BP_REAL_CRs
  MyFloat CRpNorm[BP_REAL_CRs];         /*!< normalization of CR protons spectrum */
  MyFloat CRpSlope[BP_REAL_CRs];        /*!< slope of CR protons spectrum */
  MyFloat CRpCut;                       /*!< cutoff of CR protons spectrum  */
  MyFloat CRpN[BP_REAL_CRs];            /*!< number of CR p */
  MyFloat CRpE[BP_REAL_CRs];            /*!< energy of CR p */
  MyFloat CRpPressure;                  /*!< pressure of CR p */
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
  MyFloat DtCRpE[BP_REAL_CRs];		/*!< time derivative of CR p energy */
  MyFloat DtCRpN[BP_REAL_CRs];		/*!< time derivative of CR p number */
#endif
  MyFloat CReNorm[BP_REAL_CRs];         /*!< normalization of CR electrons spectrum */
  MyFloat CReSlope[BP_REAL_CRs];        /*!< slope of CR electrons spectrum */
  MyFloat CReCut;                       /*!< cutoff of CR electrons spectrum  */
  MyFloat CReN[BP_REAL_CRs];            /*!< number of CR e */
  MyFloat CReE[BP_REAL_CRs];            /*!< energy of CR e */
  MyFloat CRePressure;                  /*!< pressure of CR e */
  MyFloat DensityOld;
#endif

#if (defined(SFR) && defined(MAGNETIC)) || defined(LT_STELLAREVOLUTION)
  float XColdCloud;
#endif
#if defined(LT_STELLAREVOLUTION)
  float Temperature;
#endif

#ifdef LT_STELLAREVOLUTION
  MyFloat  MassRes;
  MyDouble EgyRes;                      /*!< external (Sn) energy resorvoir */
  float    Metals[LT_NMetP];            /*!< H is not stored here */
#ifndef LT_LOCAL_IRA
  double   mstar;
#endif
#ifdef LT_TRACK_CONTRIBUTES
  Contrib  contrib;
#endif
#ifdef LT_ZAGE
  MyFloat  ZAge, ZAgeW;
#endif
#ifdef LT_ZAGE_LLV
  MyFloat  ZAge_llv, ZAgeW_llv;
#endif
#ifdef LT_TRACK_WINDS
  MyFloat  AvgHsml;
#endif

#endif

#ifdef LT_MV_CHEMICALDIFFUSION
  
  /* int DiffusionIntSteps; */
  /* double DiffusionIntTime; */
  
  MyFloat A, D, B[LT_NMetP];               /*!< The two coefficients in the Greif et al. approximation */
  MyFloat GradAbnd[LT_NMetP][3];
  MyFloat Abundance[LT_NMetP];             /*!< The quantity to be diffused */
  MyFloat deltaMetals[LT_NMetP];           /*!< needed to renormalize  */
  MyFloat AvgAbnd[LT_NMetP];               /*!< only for testing purposes  */

  MyDouble prevMass;
  union
  {
    MyFloat       vel_rms;
    MyDouble      dvel_rms;
  } v_rms;           /*!< rms velocity dispersion within its smoothing length */
  int NeighboringTo, TaskNeighTo, NeighboredBy, TaskNeighBy;
  int ChemFof;
  /*  int vrms_NGB;*/
#endif

  
#if defined(LT_SMOOTH_Z)
#if defined(LT_SMOOTH_ALLMETALS)
  MyFloat Zsmooth[LT_NMetP];
#else
  MyFloat Zsmooth;
  MyFloat Zsmooth_a;
  MyFloat Zsmooth_b;
#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
  float SmoothDens;
  float SmoothDens_b;
  int SmoothNgb;
#endif
#if defined(LT_SMOOTH_NGB)
  float SmoothHsml;
#endif
#endif                                 /* closes LT_METAL_COOLING_WAL  */
#endif                                 /* closes LT_SMOOTH_Z */

#ifdef LT_SMOOTH_XCLD                  /* smooth the cloud fraction */
  float XCLDsmooth;
#endif

#ifdef CHEMCOOL
  double TracAbund[TRAC_NUM];
#endif

#if defined(BLACK_HOLES) && defined(LT_BH_ACCRETE_SLICES)
  int NSlicesSwallowed;
#endif
}
  *SphP,				/*!< holds SPH particle data on local processor */
  *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */



#ifdef DoParallel
extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;		/*!< number of processors */
extern MPI_Comm MPI_CommLocal;
#endif

// for hashing tags
struct HashTable {
    struct LinkedList **table;
    size_t size;
};//*table;
struct HashTable *EmptyTable(size_t);

struct LinkedList {
    long long key;
    struct tagged_particle *star;
    struct LinkedList *next;
};
struct LinkedList *NewLinkedList(void);
//
// for hashing galaxies
struct GalHashTable {
    struct GalLinkedList **gtable;
    size_t gsize;
};//*table;
struct GalHashTable *GalEmptyTable(size_t);

struct GalLinkedList {
    long long gkey;
    struct tagged_particle *star;
    //long long ID; //particle ID
    //double BE; //binding energy
    //long long MostBndID;
    //long int rank;
    struct GalLinkedList *next;
};
struct GalLinkedList *GalNewLinkedList(void);




#endif // Don't put anything after this!!!
///////////////////////////////////////////////////
//////////////////////////////////////////////////
/////////////////////////////////////////////////
