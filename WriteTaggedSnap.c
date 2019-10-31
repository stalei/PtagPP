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
#include <mpi.h>
#endif


void WriteTaggedSnap(int snap)
{
sprintf(CurrentFile,"WriteTaggedSnap.c");
/*
long long int StarCount=0,TotalCount=0;
#ifdef DoParallel
//if(ThisTask==0)
printf("Task %d started writing tagged snap in parallel mode %d.☻ \n",ThisTask,snap);
#else
printf("I started writing tagged snap in serial mode %d.☻ \n",snap);
#endif
fflush(stdout);
//find file name, let's set it manually!
char SnapFile[500];
sprintf(SnapFile,"%s/snap_%03d",GP.SnapDir,snap);
ReadSnap(SnapFile);

fflush(stdout);

//GP.TotNumStarsAllSnaps=0;

#ifdef DoParallel
StarCount=CountTotStarsPar(GP.FirstSnap,snap); //counts and changes the type from 1 (halo) -> 4 (star)
#else
StarCount=CountTotStars(GP.FirstSnap,snap); //counts and changes the type from 1 (halo) -> 4 (star)
#endif
//int *pppp;
//pppp=0x00002b025d2d4010;
//printf("pppp:%d", *pppp);
MPI_Allreduce(&StarCount, &TotalCount, 1, MPI_LONG_LONG_INT, MPI_SUM,MPI_COMM_WORLD);
GP.TotNumStarsAllSnaps=TotalCount;

#ifdef DoParallel
if(ThisTask==0)
#endif
printf("total number of stars:%lld\n",GP.TotNumStarsAllSnaps);

int i;
#ifdef DoParallel
if(ThisTask==0)
#endif
for(i=0;i<30000;i++)
	if(P[i].Type !=1)
		printf("id:%d-type:%d-mass:%g- age:%g\n",i,P[i].Type,P[i].Mass,P[i].StellarAge,P[i].Metallicity);


#ifdef DoParallel
if(ThisTask==0){
#endif
printf("Writing binary files!\n");
SaveInBinary();

#ifdef DoParallel
}
//MPI_Finalize();
#endif
//////////////
#ifdef DoParallel
if(ThisTask==0)
#endif
*/
printf("I finished writing tagged snap %d.☻ \n",snap);
return;
}

void ReadSnap(char *fname)
{//0
FILE *fd=0;
int blksize1,blksize2;;
int i,bnr,type;
enum iofields blocknr;
char buf[500];
int bytes_per_blockelement;
int n_in_file,npart, offset=0,typelist[6];
size_t blockmaxlen;
int nread,pc,nstart;
double bytes_tot = 0;
size_t bytes;
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
//GP.TotNumPart=0;
//GP.TotNumStars=0;
if(header.num_files <= 1)
	for(i = 0; i < 6; i++)
	  {
	    header.npartTotal[i] = header.npart[i];
	    //GP.TotNumPart+=header.npart[i];
	    printf("from header:%u \n",header.npartTotal[i]);
	  }

//////////////////
//Allocate memory for particles
  if(GP.TotNumPart > 0)
    {
      if(!(P = (struct particle_data *) malloc(bytes = GP.TotNumPart * sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  EndRun(89,CurrentFile);
	}
      bytes_tot += bytes;
	printf("\nAllocated %g MByte for particle storage.\n\n", bytes_tot / (1024.0 * 1024.0));
    }
      if(!(CommBuffer = malloc( bytes = GP.BufferSize * 1024 * 1024)))
	{
	  printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  EndRun(101,CurrentFile);
	}
//for(i = 0, n_in_file = 0; i < 6; i++)
//	n_in_file += header.npart[i];
/////////////////


  for(bnr = 0; bnr < 1000; bnr++)
    {//A
	if(bnr>=3) break; /// just for test to stop at IDS.
      blocknr = (enum iofields) bnr;
      if(blocknr == IO_LASTENTRY)
	break;
	if(blockpresent(blocknr))
	{//B
	//printf("%s is saved in the file. \n",buf);
		//continue; //ignore all other blocks
	//printf("blocknr3:%d \n",blocknr);
	  if(blocknr == IO_HSMS)
	    continue;
	  get_dataset_name(blocknr, buf);
	  printf("reading block %d (%s)...\n", bnr, buf);
	  fflush(stdout);
//	printf("before bytes p element.\n");
	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 1);
	//printf("bytes per element:%d\n", bytes_per_blockelement);
	  blockmaxlen = (size_t) ((GP.BufferSize * 1024 * 1024) / bytes_per_blockelement);
	  npart = get_particles_in_block(blocknr, &typelist[0]);
	 //printf("npart:%d.\n",npart);
	  if(npart > 0)
	    {//C
	      if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
		SKIP;
		for(type = 0, offset = 0, nread = 0; type < 6; type++)
		{//E
		  n_in_file = header.npart[type];
		//printf("n_in_file:%d\n",n_in_file);
		if(typelist[type]==0)
		{
		 //n_in_file++;
		 offset +=n_in_file;
		}
		else
		{ // F 4.5
		 //n_in_file++;  
		do
		{// H 5
			pc=n_in_file;
			if(pc > (int) blockmaxlen)
				pc = blockmaxlen;
			printf("pc:%d\n",pc);
			if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
					{
					  printf("my_fread:%d\n",my_fread(CommBuffer, bytes_per_blockelement, pc, fd));
					  nread += pc;
					}
				      else
					{
					  nread += pc;
					}
			nstart=0;//no gas!
			empty_read_buffer(blocknr, nstart + offset, pc, type);
			offset += pc;
			n_in_file-=pc;
		printf("read file is done, nread:%d, offset:%d\n",nread,offset);
		}//H
		while(n_in_file > 0);
		}//F 4.5
		}// E 4
		if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
			SKIP2;
		if(blksize1 != blksize2)
		{
		printf("incorrect block-sizes detected!\n");
		printf("blocknr=%d  blksize1=%d  blksize2=%d\n", bnr,blksize1, blksize2);
		if(blocknr == IO_ID)
		{
		printf("Possible mismatch of 32bit and 64bit ID's in IC file and GADGET compilation !\n");
		}
		fflush(stdout);
		//EndRun(163,CurrentFile);
		}
	    }//C 3 
	}// B 2
    }//A 1
/////////////////////
for(type = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];
    }
fclose(fd);
}


size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if(size * nmemb == 0)
    return 0;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      if(feof(stream))
	printf("I/O error (fread) has occured: end of file\n");
      else
	printf("I/O error (fread) has occured:\n");
      fflush(stdout);
      EndRun(77,CurrentFile);
    }
  return nread;
}

/*! This function tells whether a block in the output file is present or not.
 *  */
int blockpresent(enum iofields blocknr)
{
  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
      return 1;			/* always present */

    case IO_NE:
    case IO_NH:
      if(GP.CoolingOn == 0)
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY) || defined(RADTRANSFER)
	return 1;
#else
	return 0;
#endif
      else
	return 1;
      break;

    case IO_RADGAMMA:
#ifdef RADTRANSFER
      return N_BINS;
#else
      return 0;
#endif
      break;

    case IO_RAD_ACCEL:
#ifdef RT_OUTPUT_RAD_ACCEL
      return 3;
#else
      return 0;
#endif

    case IO_HSMS:
#if defined(SUBFIND)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ACRS:
#ifdef LT_STELLAREVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_SFR:
    case IO_AGE:
    case IO_Z:
      if(GP.StarformationOn == 0)
	return 0;
      else
	{
#ifdef SFR
	  if(blocknr == IO_SFR)
	    return 1;
#endif
#ifdef STELLARAGE
	  if(blocknr == IO_AGE)
	    return 1;
#endif
#ifdef METALS
	  if(blocknr == IO_Z)
	    return 1;
#endif
	}
      return 0;
      break;

    case IO_DELAYTIME:
#ifdef WINDS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VALPHA:
#ifdef VORONOI_TIME_DEP_ART_VISC
      return 1;
#else
      return 0;
#endif
      break;

    case IO_EGYPROM:
    case IO_EGYCOLD:
#ifdef CS_FEEDBACK
      return 1;
#else
      return 0;
#endif
      break;

    case IO_HeI:
    case IO_HeII:
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY) ||defined(RADTRANSFER)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_HII:
    case IO_HeIII:
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
      return 1;
#else
      return 0;
#endif
      break;



    case IO_H2I:
    case IO_H2II:
    case IO_HM:
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_HD:
    case IO_DI:
    case IO_DII:
    case IO_HeHII:
#if defined (UM_CHEMISTRY)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_POT:
#if defined(OUTPUTPOTENTIAL) || defined(SUBFIND_RESHUFFLE_AND_POTENTIAL)
      return 1;
#else
      return 0;
#endif


    case IO_VSTURB_DISS:
    case IO_VSTURB_DRIVE:
#if defined (VS_TURB) || defined (AB_TURB)
      return 1;
#else
      return 0;
#endif
      break;


    case IO_ACCEL:
#ifdef OUTPUTACCELERATION
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DTENTR:
#ifdef OUTPUTCHANGEOFENTROPY
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSDIAG:
#ifdef OUTPUTSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSOFFDIAG:
#ifdef OUTPUTSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSBULK:
#ifdef OUTPUTBULKSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_SHEARCOEFF:
#ifdef OUTPUTSHEARCOEFF
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TSTP:
#ifdef OUTPUTTIMESTEP
      return 1;
#else
      return 0;
#endif


    case IO_BFLD:
#ifdef MAGNETIC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_BSMTH:
#ifdef OUTPUTBSMOOTH
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DBDT:
#ifdef DBOUTPUT
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VRMS:
#if defined(JD_VTURB) && !defined(FS_TURB_ESTIM)
    case IO_VBULK:
#endif
    case IO_TRUENGB:
#if defined(JD_VTURB)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VRAD:
    case IO_VTAN:
#if defined(JD_DECOMPOSE_VTURB)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VDIV:
    case IO_VROT:
#if defined(OUTPUT_DIV_CURL)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VORT:
#ifdef OUTPUT_VORTICITY
      return 1;
#else
      return 0;
#endif
      break;

    case IO_DPP:
#ifdef JD_DPP
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DIVB:
#ifdef TRACEDIVB
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ABVC:
#if defined(TIME_DEP_ART_VISC) || defined(VISCOSITY_SUPPRESSION)
      return 1;
#else
      return 0;
#endif
      break;


    case IO_AMDC:
#ifdef TIME_DEP_MAGN_DISP
      return 1;
#else
      return 0;
#endif
      break;

    case IO_LTURB:
#ifdef FS_TURB_ESTIM
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VTURB:
#ifdef FS_TURB_ESTIM
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ALFA2_DYN:
#ifdef FS_ALFA2_DYN
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ETA2_DYN:
#ifdef FS_ETA2_DYN
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PHI:
#ifdef OUTPUTDEDNER
      return 1;
#else
      return 0;
#endif
      break;

    case IO_XPHI:
#if (defined(SFR) && (defined(OUTPUT_XPHI) || defined(OUTPUTDEDNER)))
      return 1;
#else
      return 0;
#endif
      break;

    case IO_GRADPHI:
#ifdef OUTPUTDEDNER
      return 1;
#else
      return 0;
#endif
      break;


    case IO_ROTB:
#ifdef OUTPUT_ROTB
      return 1;
#else
      return 0;
#endif
      break;


    case IO_SROTB:
#ifdef OUTPUT_SROTB
      return 1;
#else
      return 0;
#endif
      break;

    case IO_EULERA:
#ifdef EULERPOTENTIALS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_EULERB:
#ifdef EULERPOTENTIALS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_COOLRATE:
#ifdef OUTPUTCOOLRATE
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CONDRATE:
#ifdef OUTPUTCONDRATE
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DENN:
#ifdef OUTPUTDENSNORM
      return 1;
#else
      return 0;
#endif
      break;


    case IO_ACRB:
    case IO_BHMASS:
    case IO_BHMDOT:
#ifdef BLACK_HOLES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BHPROGS:
#ifdef BH_COUNTPROGS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BHMBUB:
    case IO_BHMINI:
#ifdef BH_BUBBLES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MACH:
#if defined(MACHNUM) || defined(AB_SHOCK)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_SHSPEED:
    case IO_SHCOMPRESS:
    case IO_SHNORMAL:
#ifdef AB_SHOCK
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MACH1:
    case IO_MACH2:
    case IO_MACH3:
    case IO_RUP:
    case IO_RDOWN:
    case IO_PUP:
    case IO_PDOWN:
    case IO_VUP:
    case IO_VDOWN:
#if defined(AB_SHOCK_EXTRAOUTPUT)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_DTENERGY:
#ifdef MACHSTATISTIC
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PRESHOCK_CSND:
#ifdef OUTPUT_PRESHOCK_CSND
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CR_C0:
    case IO_CR_Q0:
#ifdef COSMIC_RAYS
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
#ifdef CR_OUTPUT_THERMO_VARIABLES
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
#ifdef CR_OUTPUT_TIMESCALES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PRESHOCK_DENSITY:
#if defined(CR_OUTPUT_JUMP_CONDITIONS) || defined(OUTPUT_PRESHOCK_CSND)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BPCR_pNORM:
    case IO_BPCR_eNORM:
    case IO_BPCR_pSLOPE:
    case IO_BPCR_eSLOPE:
    case IO_BPCR_pE:
    case IO_BPCR_pN:
    case IO_BPCR_pPRESSURE:
    case IO_BPCR_ePRESSURE:
#ifdef BP_REAL_CRs
      return 1;
#else
      return 0;
#endif
      break;

    case IO_CRINJECT:
#ifdef CR_OUTPUT_INJECTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TIDALTENSORPS:
#ifdef OUTPUT_TIDALTENSORPS
      return 1;
#else
      return 0;
#endif
    case IO_DISTORTIONTENSORPS:
#ifdef OUTPUT_DISTORTIONTENSORPS
      return 1;
#else
      return 0;
#endif

    case IO_CAUSTIC_COUNTER:
#ifdef DISTORTIONTENSORPS
      return 1;
#else
      return 0;
#endif

    case IO_FLOW_DETERMINANT:
#if defined(DISTORTIONTENSORPS) && !defined(GDE_LEAN)
      return 1;
#else
      return 0;
#endif

    case IO_STREAM_DENSITY:
#ifdef DISTORTIONTENSORPS
      return 1;
#else
      return 0;
#endif

    case IO_PHASE_SPACE_DETERMINANT:
#ifdef DISTORTIONTENSORPS
      return 1;
#else
      return 0;
#endif

    case IO_ANNIHILATION_RADIATION:
#if defined(DISTORTIONTENSORPS) && !defined(GDE_LEAN)
      return 1;
#else
      return 0;
#endif

    case IO_LAST_CAUSTIC:
#ifdef OUTPUT_LAST_CAUSTIC
      return 1;
#else
      return 0;
#endif

    case IO_SHEET_ORIENTATION:
#if defined(DISTORTIONTENSORPS) && (!defined(GDE_LEAN) || defined(GDE_READIC))
      return 1;
#else
      return 0;
#endif

    case IO_INIT_DENSITY:
#if defined(DISTORTIONTENSORPS) && (!defined(GDE_LEAN) || defined(GDE_READIC))
      return 1;
#else
      return 0;
#endif


    case IO_SECONDORDERMASS:
      if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
	return 1;
      else
	return 0;

    case IO_EOSTEMP:
    case IO_EOSXNUC:
    case IO_PRESSURE:
#ifdef EOS_DEGENERATE
      return 1;
#else
      return 0;
#endif

    case IO_EDDINGTON_TENSOR:
#if defined(RADTRANSFER) && defined(RT_OUTPUT_ET)
      return 1;
#else
      return 0;
#endif

    case IO_DMHSML:
    case IO_DMDENSITY:
    case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
#ifdef SIDM
      return 1;
#else
      return 0;
#endif
      break;





    case IO_DMHSML_V:
    case IO_DMDENSITY_V:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI) && defined(SUBFIND)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_Zs:
    case IO_iMass:
#ifdef LT_STELLAREVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_CLDX:
#ifdef LT_STELLAREVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_HTEMP:
#ifdef LT_STELLAREVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TEMP:
#ifdef LT_METAL_COOLING_WAL
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ZAGE:
#if defined(LT_ZAGE)
      return 1;
#else
      return 0;
#endif
      break;
    case IO_ZAGE_LLV:
#if defined(LT_ZAGE_LLV)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_CONTRIB:
#ifdef LT_TRACK_CONTRIBUTES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ZSMOOTH:
#if defined(LT_SMOOTH_Z) && !defined(LT_SMOOTH_ALLMETALS)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_allZSMOOTH:
#ifdef LT_SMOOTH_ALLMETALS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ABUNDANCE:
    case IO_ABND_GRAD:
    case IO_ABND_SPH:
    case IO_DIFFUSING_CB:
    case IO_DIFFUSING_CA:
    case IO_DIFFUSING_CD:
    case IO_DIFF_TSTEP:
#ifdef LT_MV_CHEMICALDIFFUSION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_CHEM:
#ifdef CHEMCOOL
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_SOFT:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTGRAVSOFT)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_DENS:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTGRAVNUMDENS)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_ZETA:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTZETA)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_OMEGA:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTOMEGA)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_CORR:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_NGBS:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTNGBS)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MG_PHI:
#ifdef MODGRAV
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MG_ACCEL:
#if defined(OUTPUTACCELERATION) && defined(MODGRAV)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_LASTENTRY:
      return 0;			/* will not occur */
    }


  return 0;			/* default: not present */
}

void get_dataset_name(enum iofields blocknr, char *buf)
{
  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;
    case IO_NE:
      strcpy(buf, "ElectronAbundance");
      break;
    case IO_NH:
      strcpy(buf, "NeutralHydrogenAbundance");
      break;
    case IO_RADGAMMA:
      strcpy(buf, "photon number density");
      break;
    case IO_RAD_ACCEL:
      strcpy(buf, "rad acceleration");
      break;
    case IO_HII:
      strcpy(buf, "HII");
      break;
    case IO_HeI:
      strcpy(buf, "HeI");
      break;
    case IO_HeII:
      strcpy(buf, "HeII");
      break;
    case IO_HeIII:
      strcpy(buf, "HeIII");
      break;
    case IO_H2I:
      strcpy(buf, "H2I");
      break;
    case IO_H2II:
      strcpy(buf, "H2II");
      break;
    case IO_HM:
      strcpy(buf, "HM");
      break;
    case IO_HD:
      strcpy(buf, "HD  ");
      break;
    case IO_DI:
      strcpy(buf, "DI  ");
      break;
    case IO_DII:
      strcpy(buf, "DII ");
      break;
    case IO_HeHII:
      strcpy(buf, "HeHp");
      break;
    case IO_DELAYTIME:
      strcpy(buf, "DelayTime");
      break;
    case IO_HSML:
      strcpy(buf, "SmoothingLength");
      break;
    case IO_VALPHA:
      strcpy(buf, "ArtificialViscosityV");
      break;
    case IO_SFR:
      strcpy(buf, "StarFormationRate");
      break;
    case IO_AGE:
      strcpy(buf, "StellarFormationTime");
      break;
    case IO_HSMS:
      strcpy(buf, "StellarSmoothingLength");
      break;
    case IO_ACRS:
      strcpy(buf, "StellarSpreadingLength");
      break;
    case IO_Z:
      strcpy(buf, "Metallicity");
      break;
    case IO_POT:
      strcpy(buf, "Potential");
      break;
    case IO_ACCEL:
      strcpy(buf, "Acceleration");
      break;
    case IO_DTENTR:
      strcpy(buf, "RateOfChangeOfEntropy");
      break;
    case IO_STRESSDIAG:
      strcpy(buf, "DiagonalStressTensor");
      break;
    case IO_STRESSOFFDIAG:
      strcpy(buf, "OffDiagonalStressTensor");
      break;
    case IO_STRESSBULK:
      strcpy(buf, "BulkStressTensor");
      break;
    case IO_SHEARCOEFF:
      strcpy(buf, "ShearCoefficient");
      break;
    case IO_TSTP:
      strcpy(buf, "TimeStep");
      break;
    case IO_BFLD:
      strcpy(buf, "MagneticField");
      break;
    case IO_BSMTH:
      strcpy(buf, "SmoothedMagneticField");
      break;
    case IO_DBDT:
      strcpy(buf, "RateOfChangeOfMagneticField");
      break;
    case IO_VRMS:
      strcpy(buf, "RMSVelocity");
      break;
    case IO_VBULK:
      strcpy(buf, "BulkVelocity");
      break;
    case IO_VRAD:
      strcpy(buf, "RMSRadialVelocity");
      break;
    case IO_VTAN:
      strcpy(buf, "RMSTangentialVelocity");
      break;
    case IO_TRUENGB:
      strcpy(buf, "TrueNumberOfNeighbours");
      break;
    case IO_DPP:
      strcpy(buf, "MagnetosReaccCoefficient");
      break;
    case IO_VDIV:
      strcpy(buf, "VelocityDivergence");
      break;
    case IO_VROT:
      strcpy(buf, "VelocityCurl");
      break;
    case IO_VORT:
      strcpy(buf, "Vorticity");
      break;
    case IO_DIVB:
      strcpy(buf, "DivergenceOfMagneticField");
      break;
    case IO_ABVC:
      strcpy(buf, "ArtificialViscosity");
      break;
    case IO_AMDC:
      strcpy(buf, "ArtificialMagneticDissipatio");
      break;
    case IO_VTURB:
      strcpy(buf, "TurbVelociy");
      break;
    case IO_LTURB:
      strcpy(buf, "TurbLengh");
      break;
    case IO_ALFA2_DYN:
      strcpy(buf, "AlfaDynamo");
      break;
    case IO_ETA2_DYN:
      strcpy(buf, "EtaDynamo");
      break;
    case IO_PHI:
      strcpy(buf, "DivBcleaningFunctionPhi");
      break;
    case IO_XPHI:
      strcpy(buf, "ColdGasFraction_PHI");
      break;
    case IO_GRADPHI:
      strcpy(buf, "DivBcleaningFunctionGadPhi");
      break;
    case IO_ROTB:
      strcpy(buf, "RotationB");
      break;
    case IO_SROTB:
      strcpy(buf, "SmoothedRotationB");
      break;
    case IO_EULERA:
      strcpy(buf, "EulerPotentialA");
      break;
    case IO_EULERB:
      strcpy(buf, "EulerPotentialB");
      break;
    case IO_COOLRATE:
      strcpy(buf, "CoolingRate");
      break;
    case IO_CONDRATE:
      strcpy(buf, "ConductionRate");
      break;
    case IO_DENN:
      strcpy(buf, "Denn");
      break;
    case IO_EGYPROM:
      strcpy(buf, "EnergyReservoirForFeeback");
      break;
    case IO_EGYCOLD:
      strcpy(buf, "EnergyReservoirForColdPhase");
      break;
    case IO_CR_C0:
      strcpy(buf, "CR_C0");
      break;
    case IO_CR_Q0:
      strcpy(buf, "CR_q0");
      break;
    case IO_CR_P0:
      strcpy(buf, "CR_P0");
      break;
    case IO_CR_E0:
      strcpy(buf, "CR_E0");
      break;
    case IO_CR_n0:
      strcpy(buf, "CR_n0");
      break;
    case IO_CR_ThermalizationTime:
      strcpy(buf, "CR_ThermalizationTime");
      break;
    case IO_CR_DissipationTime:
      strcpy(buf, "CR_DissipationTime");
      break;
    case IO_BHMASS:
      strcpy(buf, "BH_Mass");
      break;
    case IO_ACRB:
      strcpy(buf, "BH_AccreationLength");
      break;
    case IO_BHMDOT:
      strcpy(buf, "BH_Mdot");
      break;
    case IO_BHPROGS:
      strcpy(buf, "BH_NProgs");
      break;
    case IO_BHMBUB:
      strcpy(buf, "BH_Mass_bubbles");
      break;
    case IO_BHMINI:
      strcpy(buf, "BH_Mass_ini");
      break;
    case IO_BHMRAD:
      strcpy(buf, "BH_Mass_radio");
      break;
    case IO_MACH:
      strcpy(buf, "MachNumber");
      break;
    case IO_MACH1:
      strcpy(buf, "MachNumber Vel");
      break;
    case IO_MACH2:
      strcpy(buf, "MachNumber Press");
      break;
    case IO_MACH3:
      strcpy(buf, "MachNumber Rho");
      break;
    case IO_RUP:
      strcpy(buf, "Rho Upstream");
      break;
    case IO_RDOWN:
      strcpy(buf, "Rho Downstream");
      break;
    case IO_PUP:
      strcpy(buf, "Pressure Upstream");
      break;
    case IO_PDOWN:
      strcpy(buf, "Pressure Downstream");
      break;
    case IO_SHSPEED:
      strcpy(buf, "ShockSpeed");
      break;
    case IO_SHCOMPRESS:
      strcpy(buf, "ShockCompressionRatio");
      break;
    case IO_SHNORMAL:
      strcpy(buf, "ShockNormal");
      break;
    case IO_VUP:
      strcpy(buf, "Velocity Upstream");
      break;
    case IO_VDOWN:
      strcpy(buf, "Velocity Downstream");
      break;
    case IO_DTENERGY:
      strcpy(buf, "DtEnergy");
      break;
    case IO_PRESHOCK_CSND:
      strcpy(buf, "Preshock_SoundSpeed");
      break;
    case IO_PRESHOCK_DENSITY:
      strcpy(buf, "Preshock_Density");
      break;
    case IO_PRESHOCK_ENERGY:
      strcpy(buf, "Preshock_Energy");
      break;
    case IO_PRESHOCK_XCR:
      strcpy(buf, "Preshock_XCR");
      break;
    case IO_DENSITY_JUMP:
      strcpy(buf, "DensityJump");
      break;
    case IO_ENERGY_JUMP:
      strcpy(buf, "EnergyJump");
      break;
    case IO_CRINJECT:
      strcpy(buf, "CR_DtE");
      break;
    case IO_TIDALTENSORPS:
      strcpy(buf, "TidalTensorPS");
      break;
    case IO_DISTORTIONTENSORPS:
      strcpy(buf, "DistortionTensorPS");
      break;
    case IO_CAUSTIC_COUNTER:
      strcpy(buf, "CausticCounter");
      break;
    case IO_FLOW_DETERMINANT:
      strcpy(buf, "FlowDeterminant");
      break;
    case IO_STREAM_DENSITY:
      strcpy(buf, "StreamDensity");
      break;
    case IO_SECONDORDERMASS:
      strcpy(buf, "2lpt-mass");
      break;
    case IO_PHASE_SPACE_DETERMINANT:
      strcpy(buf, "PhaseSpaceDensity");
      break;
    case IO_ANNIHILATION_RADIATION:
      strcpy(buf, "AnnihilationRadiation");
      break;
    case IO_LAST_CAUSTIC:
      strcpy(buf, "LastCaustic");
      break;
    case IO_SHEET_ORIENTATION:
      strcpy(buf, "SheetOrientation");
      break;
    case IO_INIT_DENSITY:
      strcpy(buf, "InitDensity");
      break;
    case IO_EOSTEMP:
      strcpy(buf, "Temperature");
      break;
    case IO_EOSXNUC:
      strcpy(buf, "Nuclear mass fractions");
      break;
    case IO_PRESSURE:
      strcpy(buf, "Pressure");
      break;
    case IO_EDDINGTON_TENSOR:
      strcpy(buf, "EddingtonTensor");
      break;
    case IO_PSUM:
      strcpy(buf, "PSum");
      break;
    case IO_SIDMNUMNGB:
      strcpy(buf, "DMNumNgb");
      break;
    case IO_NUMTOTALSCATTER:
      strcpy(buf, "NumTotalScatter");
      break;
    case IO_SIDMHSML:
      strcpy(buf, "SIDMHsml");
      break;
    case IO_SIDMDENSITY:
      strcpy(buf, "SIDMRho");
      break;
    case IO_SIDMVELDISP:
      strcpy(buf, "SVelDisp");
      break;
    case IO_DMHSML:
      strcpy(buf, "DM Hsml");
      break;
    case IO_DMDENSITY:
      strcpy(buf, "DM Density");
      break;
    case IO_DMVELDISP:
      strcpy(buf, "DM Velocity Dispersion");
      break;
    case IO_DMHSML_V:
      strcpy(buf, "DM Hsml Voronoi");
      break;
    case IO_DMDENSITY_V:
      strcpy(buf, "DM Density Voronoi");
      break;
    case IO_Zs:
      strcpy(buf, "Mass of Metals");
      break;
    case IO_ZAGE:
      strcpy(buf, "Metallicity-averaged time");
      break;
    case IO_ZAGE_LLV:
      strcpy(buf, "long-living-Metallicity-averaged time");
      break;
    case IO_iMass:
      strcpy(buf, "SSPInitialMass");
      break;
    case IO_CLDX:
      strcpy(buf, "CloudFraction");
      break;
    case IO_HTEMP:
      strcpy(buf, "HotPhaseTemperature");
      break;
    case IO_TEMP:
      strcpy(buf, "Temperature");
      break;
    case IO_CONTRIB:
      strcpy(buf, "TrackContributes");
      break;
    case IO_ZSMOOTH:
      strcpy(buf, "smoothed metallicity");
      break;
    case IO_allZSMOOTH:
      strcpy(buf, "smoothed metallicities (all elements)");
      break;
    case IO_ABUNDANCE:
      strcpy(buf, "abundance by diffusion");
      break;
    case IO_ABND_GRAD:
      strcpy(buf, "abundance gradient");
      break;
    case IO_ABND_SPH:
      strcpy(buf, "SPH-smoothed abundance");
      break;
    case IO_DIFFUSING_CB:
      strcpy(buf, "Diffusion B-coefficient");
      break;
    case IO_DIFFUSING_CA:
      strcpy(buf, "Diffusion A-coeffient");
      break;
    case IO_DIFFUSING_CD:
      strcpy(buf, "Diffusion D-coefficient");
      break;
    case IO_DIFF_TSTEP:
      strcpy(buf, "Diffusion Time-Step");
      break;
    case IO_CHEM:
      strcpy(buf, "ChemicalAbundances");
      break;
    case IO_BPCR_pNORM:
      strcpy(buf, "BPCRpNormalization");
      break;
    case IO_BPCR_eNORM:
      strcpy(buf, "BPCReNormalization");
      break;
    case IO_BPCR_pSLOPE:
      strcpy(buf, "BPCRpSlope");
      break;
    case IO_BPCR_eSLOPE:
      strcpy(buf, "BPCReSlope");
      break;
    case IO_BPCR_pE:
      strcpy(buf, "BPCRpEnergy");
      break;
    case IO_BPCR_pN:
      strcpy(buf, "BPCRpNumber");
      break;
    case IO_BPCR_pPRESSURE:
      strcpy(buf, "BPCRpPressure");
      break;
    case IO_BPCR_ePRESSURE:
      strcpy(buf, "BPCRePressure");
      break;
    case IO_AGS_SOFT:
      strcpy(buf, "AGS-Softening");
      break;
    case IO_AGS_DENS:
      strcpy(buf, "AGS-Density");
      break;
    case IO_AGS_ZETA:
      strcpy(buf, "AGS-Zeta");
      break;
    case IO_AGS_OMEGA:
      strcpy(buf, "AGS-Omega");
      break;
    case IO_AGS_CORR:
      strcpy(buf, "AGS-Correction");
      break;
    case IO_AGS_NGBS:
      strcpy(buf, "AGS-Neighbours");
      break;
    case IO_VSTURB_DISS:
      strcpy(buf, "TurbulenceDissipation");
      break;
    case IO_VSTURB_DRIVE:
      strcpy(buf, "TurbulenceDriving");
      break;
    case IO_MG_PHI:
      strcpy(buf, "ModifiedGravityPhi");
      break;
    case IO_MG_ACCEL:
      strcpy(buf, "ModifiedGravityAcceleration");
      break;

    case IO_LASTENTRY:
	printf("Reached the last entry!\n");
      //EndRun(1581,CurrentFile);
      break;
    }
}


/*
utput file.
 */
int get_bytes_per_blockelement(enum iofields blocknr, int mode)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_VBULK:
    case IO_BFLD:
    case IO_GRADPHI:
    case IO_BSMTH:
    case IO_DBDT:
    case IO_ROTB:
    case IO_SROTB:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_RAD_ACCEL:
    case IO_VORT:
    case IO_SHNORMAL:
    case IO_VUP:
    case IO_VDOWN:
    case IO_MG_ACCEL:
      if(mode)
	bytes_per_blockelement = 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
      break;

    case IO_ID:
      bytes_per_blockelement = sizeof(MyIDType);
      break;

    case IO_BHPROGS:
    case IO_TRUENGB:
    case IO_AGS_NGBS:
      bytes_per_blockelement = sizeof(int);
      break;

    case IO_MASS:
    case IO_SECONDORDERMASS:
    case IO_U:
    case IO_RHO:
    case IO_NE:
    case IO_NH:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
    case IO_HD:
    case IO_DI:
    case IO_DII:
    case IO_HeHII:
    case IO_HSML:
    case IO_VALPHA:
    case IO_SFR:
    case IO_AGE:
    case IO_DELAYTIME:
    case IO_HSMS:
    case IO_ACRS:
    case IO_POT:
    case IO_DTENTR:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_TSTP:
    case IO_DIVB:
    case IO_VRMS:
    case IO_VRAD:
    case IO_VTAN:
    case IO_VDIV:
    case IO_VROT:
    case IO_DPP:
    case IO_ABVC:
    case IO_AMDC:
    case IO_VTURB:
    case IO_LTURB:
    case IO_ALFA2_DYN:
    case IO_ETA2_DYN:
    case IO_PHI:
    case IO_XPHI:
    case IO_EULERA:
    case IO_EULERB:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_EGYPROM:
    case IO_EGYCOLD:
    case IO_BHMASS:
    case IO_ACRB:
    case IO_BHMDOT:
    case IO_BHMBUB:
    case IO_BHMINI:
    case IO_BHMRAD:
    case IO_MACH:
    case IO_MACH1:
    case IO_MACH2:
    case IO_MACH3:
    case IO_RUP:
    case IO_RDOWN:
    case IO_PUP:
    case IO_PDOWN:
    case IO_SHSPEED:
    case IO_SHCOMPRESS:
    case IO_DTENERGY:
    case IO_PRESHOCK_CSND:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_CAUSTIC_COUNTER:
    case IO_FLOW_DETERMINANT:
    case IO_STREAM_DENSITY:
    case IO_PHASE_SPACE_DETERMINANT:
    case IO_EOSTEMP:
    case IO_PRESSURE:
    case IO_INIT_DENSITY:
    case IO_iMass:
    case IO_CLDX:
    case IO_HTEMP:
    case IO_TEMP:
    case IO_ZSMOOTH:
    case IO_BPCR_pPRESSURE:
    case IO_BPCR_ePRESSURE:
    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
    case IO_AGS_SOFT:
    case IO_AGS_DENS:
    case IO_AGS_ZETA:
    case IO_AGS_OMEGA:
    case IO_AGS_CORR:
    case IO_VSTURB_DISS:
    case IO_VSTURB_DRIVE:
    case IO_MG_PHI:
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      if(mode)
	bytes_per_blockelement = NUMCRPOP * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = NUMCRPOP * sizeof(MyOutputFloat);
      break;

    case IO_BPCR_pNORM:
    case IO_BPCR_eNORM:
    case IO_BPCR_pSLOPE:
    case IO_BPCR_eSLOPE:
    case IO_BPCR_pE:
    case IO_BPCR_pN:
#ifndef BP_REAL_CRs
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;
#else
      if(mode)
	bytes_per_blockelement = BP_REAL_CRs * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = BP_REAL_CRs * sizeof(MyOutputFloat);
      break;
#endif

    case IO_DMHSML:
    case IO_DMDENSITY:
    case IO_DMVELDISP:
    case IO_DMHSML_V:
    case IO_DMDENSITY_V:
      bytes_per_blockelement = sizeof(float);
      break;

    case IO_RADGAMMA:
#ifdef RADTRANSFER
      if(mode)
	bytes_per_blockelement = N_BINS * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = N_BINS * sizeof(MyOutputFloat);
#endif
      break;

    case IO_EDDINGTON_TENSOR:
      if(mode)
	bytes_per_blockelement = 6 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 6 * sizeof(MyOutputFloat);


    case IO_Z:
#ifndef CS_MODEL
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;
#else
      if(mode)
	bytes_per_blockelement = 12 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 12 * sizeof(MyOutputFloat);
      break;
#endif

    case IO_TIDALTENSORPS:
      if(mode)
	bytes_per_blockelement = 9 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
      break;

    case IO_DISTORTIONTENSORPS:
      if(mode)
	bytes_per_blockelement = 36 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 36 * sizeof(MyOutputFloat);
      break;

    case IO_ANNIHILATION_RADIATION:
      if(mode)
	bytes_per_blockelement = 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
      break;

    case IO_LAST_CAUSTIC:
      if(mode)
	bytes_per_blockelement = 20 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 20 * sizeof(MyOutputFloat);
      break;

    case IO_SHEET_ORIENTATION:
      if(mode)
	bytes_per_blockelement = 9 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
      break;

    case IO_EOSXNUC:
#ifdef EOS_DEGENERATE
      if(mode)
	bytes_per_blockelement = EOS_NSPECIES * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = EOS_NSPECIES * sizeof(MyOutputFloat);
      break;
#else
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;
#endif

    case IO_Zs:
#ifdef LT_STELLAREVOLUTION
      if(mode)
	bytes_per_blockelement = LT_NMetP * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = LT_NMetP * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_allZSMOOTH:
#ifdef LT_SMOOTH_ALLMETALS
      if(mode)
	bytes_per_blockelement = LT_NMetP * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = LT_NMetP * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_ZAGE:
    case IO_ZAGE_LLV:
#if defined(LT_ZAGE) || defined(LT_ZAGE_LLV)
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_CONTRIB:
#ifdef LT_TRACK_CONTRIBUTES
      bytes_per_blockelement = sizeof(Contrib);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_DIFFUSING_CB:
    case IO_ABUNDANCE:
    case IO_ABND_SPH:
#ifdef LT_MV_CHEMICALDIFFUSION
      if(mode)
	bytes_per_blockelement = LT_NMetP * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = LT_NMetP * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_DIFFUSING_CA:
    case IO_DIFFUSING_CD:
    case IO_DIFF_TSTEP:
#ifdef LT_MV_CHEMICALDIFFUSION
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_ABND_GRAD:
#ifdef LT_MV_CHEMICALDIFFUSION
      if(mode)
	bytes_per_blockelement = LT_NMetP * 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = LT_NMetP * 3 * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;



    case IO_CHEM:
#ifdef CHEMCOOL
      if(mode)
	bytes_per_blockelement = TRAC_NUM * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = TRAC_NUM * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;


    case IO_LASTENTRY:
      EndRun(1921,CurrentFile);
      break;
    }

  return bytes_per_blockelement;
}


/*! This function determines how many particles there are in a given block,
 *  *  based on the information in the header-structure.  It also flags particle
 *   *  types that are present in the block in the typelist array.
 *    */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, nsel, ntot_withmasses, ngas, ngasAlpha, nstars, nngb;

  nall = 0;
  nsel = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(GP.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  ngasAlpha = NUMCRPOP * header.npart[0];
  nstars = header.npart[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:
    case IO_SECONDORDERMASS:
    case IO_AGS_SOFT:
    case IO_AGS_DENS:
    case IO_AGS_ZETA:
    case IO_AGS_OMEGA:
    case IO_AGS_CORR:
    case IO_AGS_NGBS:
    case IO_MG_PHI:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(GP.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_RAD_ACCEL:
    case IO_RADGAMMA:
    case IO_EDDINGTON_TENSOR:
    case IO_U:
    case IO_RHO:
    case IO_NE:
    case IO_NH:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
    case IO_HD:
    case IO_DI:
    case IO_DII:
    case IO_HeHII:
    case IO_HSML:
    case IO_DELAYTIME:
    case IO_VALPHA:
    case IO_SFR:
    case IO_DTENTR:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_BSMTH:
    case IO_BFLD:
    case IO_DBDT:
    case IO_VRMS:
    case IO_VBULK:
    case IO_VTAN:
    case IO_VRAD:
    case IO_VDIV:
    case IO_VROT:
    case IO_VORT:
    case IO_DPP:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_VTURB:
    case IO_LTURB:
    case IO_ALFA2_DYN:
    case IO_ETA2_DYN:
    case IO_PHI:
    case IO_XPHI:
    case IO_GRADPHI:
    case IO_ROTB:
    case IO_SROTB:
    case IO_EULERA:
    case IO_EULERB:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_MACH:
    case IO_MACH1:
    case IO_MACH2:
    case IO_MACH3:
    case IO_RUP:
    case IO_RDOWN:
    case IO_PUP:
    case IO_PDOWN:
    case IO_SHSPEED:
    case IO_SHCOMPRESS:
    case IO_SHNORMAL:
    case IO_VUP:
    case IO_VDOWN:
    case IO_DTENERGY:
    case IO_PRESHOCK_CSND:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_EOSTEMP:
    case IO_EOSXNUC:
    case IO_PRESSURE:
    case IO_CHEM:
    case IO_BPCR_pNORM:
    case IO_BPCR_eNORM:
    case IO_BPCR_pSLOPE:
    case IO_BPCR_eSLOPE:
    case IO_BPCR_pE:
    case IO_BPCR_pN:
    case IO_BPCR_pPRESSURE:
    case IO_BPCR_ePRESSURE:
    case IO_VSTURB_DISS:
    case IO_VSTURB_DRIVE:
    case IO_ABUNDANCE:
    case IO_DIFFUSING_CD:
    case IO_DIFFUSING_CA:
    case IO_DIFFUSING_CB:
    case IO_ABND_GRAD:
    case IO_ABND_SPH:
    case IO_DIFF_TSTEP:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngas;
      break;

    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngasAlpha;
      break;

    case IO_AGE:
      for(i = 0; i < 6; i++)
#ifdef BLACK_HOLES
	if(i != 4 && i != 5)
	  typelist[i] = 0;
      return nstars + header.npart[5];
#else
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
#endif
      break;

    case IO_TRUENGB:
      nngb = ngas;
      for(i = 1; i < 4; i++)
	typelist[i] = 0;
#ifdef LT_STELLAREVOLUTION
      nngb += header.npart[4];
#else
      typelist[4] = 0;
#endif
#ifdef BLACK_HOLES
      nngb += header.npart[5];
#else
      typelist[5] = 0;
#endif
      return nngb;
      break;

    case IO_HSMS:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
      break;

    case IO_ACRS:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
      break;

    case IO_Z:
    case IO_EGYPROM:
    case IO_EGYCOLD:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;

    case IO_BHMASS:
    case IO_ACRB:
    case IO_BHMDOT:
    case IO_BHMBUB:
    case IO_BHMINI:
    case IO_BHMRAD:
    case IO_BHPROGS:
      for(i = 0; i < 6; i++)
	if(i != 5)
	  typelist[i] = 0;
      return header.npart[5];
      break;

    case IO_TIDALTENSORPS:
    case IO_DISTORTIONTENSORPS:
    case IO_CAUSTIC_COUNTER:
    case IO_FLOW_DETERMINANT:
    case IO_STREAM_DENSITY:
    case IO_PHASE_SPACE_DETERMINANT:
    case IO_ANNIHILATION_RADIATION:
    case IO_LAST_CAUSTIC:
    case IO_SHEET_ORIENTATION:
    case IO_INIT_DENSITY:
      for(i = 0; i < 6; i++)
	if(((1 << i) & (GDE_TYPES)))
	  nsel += header.npart[i];
	else
	  typelist[i] = 0;
      return nsel;
      break;


    case IO_DMHSML:
    case IO_DMDENSITY:
    case IO_DMVELDISP:
    case IO_DMHSML_V:
    case IO_DMDENSITY_V:
      for(i = 0; i < 6; i++)
	if(((1 << i) & (FOF_PRIMARY_LINK_TYPES)))
	  nsel += header.npart[i];
	else
	  typelist[i] = 0;
      return nsel;
      break;

    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
      for(i = 0; i < 6; i++)
#ifdef SIDM
	if(((1 << i) & (SIDM)))
	  nsel += header.npart[i];
	else
	  typelist[i] = 0;
#else
	typelist[i] = 0;
#endif
      return nsel;
      break;


    case IO_Zs:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;
    case IO_ZAGE:
    case IO_ZAGE_LLV:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;
    case IO_iMass:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
      break;
    case IO_CLDX:
    case IO_HTEMP:
      for(i = 0; i < 6; i++)
	if(i != 0)
	  typelist[i] = 0;
      return ngas;
      break;
    case IO_TEMP:
      for(i = 0; i < 6; i++)
	if(i != 0)
	  typelist[i] = 0;
      return ngas;
      break;
    case IO_ZSMOOTH:
      for(i = 0; i < 6; i++)
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
	if(i != 0)
#else
	if(i != 0 && i != 4)
#endif
	  typelist[i] = 0;
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
      return ngas;
#else
      return ngas + nstars;
#endif
      break;
    case IO_allZSMOOTH:
      for(i = 0; i < 6; i++)
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
	if(i != 0)
#else
	if(i != 0 && i != 4)
#endif
	  typelist[i] = 0;
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
      return ngas;
#else
      return ngas + nstars;
#endif
      break;
    case IO_CONTRIB:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;

    case IO_MG_ACCEL:
      for(i = 2; i < 6; i++)
	typelist[i] = 0;
      return ngas + header.npart[1];
      break;

    case IO_LASTENTRY:
      EndRun(2292,CurrentFile);
      break;
    }

  EndRun(2296,CurrentFile);
  return 0;
}




//! This function reads out the buffer that was filled with particle data.
//
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
  int n, k;
  MyInputFloat *fp;
  MyIDType *ip;
  float *fp_single;

#ifdef AUTO_SWAP_ENDIAN_READIC
  int vt, vpb;
  char *cp;
#endif

  fp = (MyInputFloat *) CommBuffer;
  fp_single = (float *) CommBuffer;
  ip = (MyIDType *) CommBuffer;
#ifdef LT_TRACK_CONTRIBUTES
  contrib = (Contrib *) CommBuffer;
#endif

#ifdef AUTO_SWAP_ENDIAN_READIC
  if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V
     && blocknr != IO_DMDENSITY_V)
    {
      cp = (char *) CommBuffer;
      vt = get_datatype_in_block(blocknr);
      vpb = get_values_per_blockelement(blocknr);
      if(vt == 2)
	swap_Nbyte(cp, pc * vpb, 8);
      else
	{
#ifdef INPUT_IN_DOUBLEPRECISION
	  if(vt == 1)
	    swap_Nbyte(cp, pc * vpb, 8);
	  else
#endif
	    swap_Nbyte(cp, pc * vpb, 4);
	}
    }
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#endif

  switch (blocknr)
    {
    case IO_POS:		// positions 
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  P[offset + n].Pos[k] = *fp++;

      for(n = 0; n < pc; n++)
	P[offset + n].Type = type;	// initialize type here as well 
      break;

    case IO_VEL:		// velocities 
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
#ifdef RESCALEVINI
	  // scaling v to use same IC's for different cosmologies 
	  if(RestartFlag == 0)
	    P[offset + n].Vel[k] = (*fp++) * All.VelIniScale;
	  else
	    P[offset + n].Vel[k] = *fp++;
#else
	  P[offset + n].Vel[k] = *fp++;
#endif
      break;

    case IO_ID:		// particle ID 
      for(n = 0; n < pc; n++)
	P[offset + n].ID = *ip++;
      break;

    case IO_MASS:		// particle mass 
      for(n = 0; n < pc; n++)
	P[offset + n].Mass = *fp++;
      break;


    case IO_SHEET_ORIENTATION:	// initial particle sheet orientation 
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
      for(n = 0; n < pc; n++)
	{
#ifndef GDE_LEAN
	  P[offset + n].V_matrix[0][0] = *fp++;
	  P[offset + n].V_matrix[0][1] = *fp++;
	  P[offset + n].V_matrix[0][2] = *fp++;
	  P[offset + n].V_matrix[1][0] = *fp++;
	  P[offset + n].V_matrix[1][1] = *fp++;
	  P[offset + n].V_matrix[1][2] = *fp++;
	  P[offset + n].V_matrix[2][0] = *fp++;
	  P[offset + n].V_matrix[2][1] = *fp++;
	  P[offset + n].V_matrix[2][2] = *fp++;
#else
	  *fp += 8;
#endif
	}
#endif
      break;

    case IO_INIT_DENSITY:	// initial stream density 
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
      for(n = 0; n < pc; n++)
	GDE_INITDENSITY(offset + n) = *fp++;
      break;
#endif

    case IO_CAUSTIC_COUNTER:	// initial caustic counter 
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
      for(n = 0; n < pc; n++)
	P[offset + n].caustic_counter = *fp++;
      break;
#endif

    case IO_SECONDORDERMASS:
      for(n = 0; n < pc; n++)
	{
	  P[offset + n].OldAcc = P[offset + n].Mass;	// use this to temporarily store the masses in the 2plt IC case 
	  P[offset + n].Mass = *fp++;
	}
      break;

    case IO_U:			// temperature
      for(n = 0; n < pc; n++)
	SphP[offset + n].Entropy = *fp++;
      break;

    case IO_RHO:		// density 
      for(n = 0; n < pc; n++)
	SphP[offset + n].d.Density = *fp++;
      break;

    case IO_NE:		// electron abundance
#if defined(COOLING) || defined(CHEMISTRY) || defined(UM_CHEMISTRY) || defined(RADTRANSFER)
      for(n = 0; n < pc; n++)
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY) || defined(RADTRANSFER)
	SphP[offset + n].elec = *fp++;
#else
	SphP[offset + n].Ne = *fp++;
#endif
#endif
      break;

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY) || defined(RADTRANSFER)
    case IO_NH:		// neutral hydrogen abundance 
      for(n = 0; n < pc; n++)
	SphP[offset + n].HI = *fp++;
      break;

    case IO_HII:		// ionized hydrogen abundance
      for(n = 0; n < pc; n++)
	SphP[offset + n].HII = *fp++;
      break;

    case IO_HeI:		// neutral Helium 
      for(n = 0; n < pc; n++)
	SphP[offset + n].HeI = *fp++;
      break;

    case IO_HeII:		// ionized Heluum 
      for(n = 0; n < pc; n++)
	SphP[offset + n].HeII = *fp++;

    case IO_HeIII:		// double ionised Helium 
      for(n = 0; n < pc; n++)
	SphP[offset + n].HeIII = *fp++;
      break;
#endif

    case IO_H2I:		// H2 molecule 
      for(n = 0; n < pc; n++)
	SphP[offset + n].H2I = *fp++;
      break;

    case IO_H2II:		// ionised H2 molecule
      for(n = 0; n < pc; n++)
	SphP[offset + n].H2II = *fp++;

    case IO_HM:		// H minus
      for(n = 0; n < pc; n++)
	SphP[offset + n].HM = *fp++;
      break;

    case IO_HeHII:		// HeH+
#if defined (UM_CHEMISTRY)
      for(n = 0; n < pc; n++)
	SphP[offset + n].HeHII = *fp++;
#endif
      break;

    case IO_HD:		// HD
#if defined (UM_CHEMISTRY) &&  defined (UM_HD_COOLING)
      for(n = 0; n < pc; n++)
	SphP[offset + n].HD = *fp++;
#endif
      break;

    case IO_DI:		// D //
#if defined (UM_CHEMISTRY) &&  defined (UM_HD_COOLING)
      for(n = 0; n < pc; n++)
	SphP[offset + n].DI = *fp++;
#endif
      break;

    case IO_DII:		// D plus //
#if defined (UM_CHEMISTRY) &&  defined (UM_HD_COOLING)
      for(n = 0; n < pc; n++)
	SphP[offset + n].DII = *fp++;
#endif
      break;

#else
    case IO_NH:		// neutral hydrogen abundance 
    case IO_HII:		// ionized hydrogen abundance 
    case IO_HeI:		// neutral Helium 
    case IO_HeII:		 //ionized Heluum 
    case IO_HeIII:		// double ionised Helium 
    case IO_H2I:		// H2 molecule 
    case IO_H2II:		// ionised H2 molecule
    case IO_HM:		// H minus 
    case IO_HeHII:		// HeH+
    case IO_HD:		// HD 
    case IO_DI:		// D 
    case IO_DII:		// D plus  
      break;
#endif

    case IO_HSML:		// SPH smoothing length 
      for(n = 0; n < pc; n++)
	PPP[offset + n].Hsml = *fp++;
      break;

    case IO_DELAYTIME:
#ifdef WINDS
      for(n = 0; n < pc; n++)
	SphP[offset + n].DelayTime = *fp++;
#endif
      break;

    case IO_AGE:		// Age of stars 
#ifdef STELLARAGE
      for(n = 0; n < pc; n++)
	P[offset + n].StellarAge = *fp++;
#endif
      break;

    case IO_Z:			// Gas and star metallicity
#ifdef METALS
#ifndef CS_MODEL
      for(n = 0; n < pc; n++)
	P[offset + n].Metallicity = *fp++;
#else
      for(n = 0; n < pc; n++)
	for(k = 0; k < 12; k++)
	  P[offset + n].Zm[k] = *fp++;
#endif
#endif
      break;

    case IO_EGYPROM:		// SN Energy Reservoir
#ifdef CS_FEEDBACK
      for(n = 0; n < pc; n++)
	P[offset + n].EnergySN = *fp++;
#endif
      break;

    case IO_EGYCOLD:		// Cold  SN Energy Reservoir
#ifdef CS_FEEDBACK
      for(n = 0; n < pc; n++)
	P[offset + n].EnergySNCold = *fp++;
#endif
      break;

    case IO_VRMS:		// Turbulence on kernel scale
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset + n].Vrms = *fp++;
#endif
      break;
    case IO_VBULK:
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  SphP[offset + n].Vbulk[k] = *fp++;
#endif
      break;
    case IO_VTAN:
#ifdef JD_DECOMPOSE_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset + n].Vtan = *fp++;
#endif
      break;
    case IO_VRAD:
#ifdef JD_DECOMPOSE_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset + n].Vrad = *fp++;
#endif
      break;
    case IO_VDIV:
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset + n].v.DivVel = *fp++;
#endif
      break;
    case IO_VROT:
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset + n].r.CurlVel = *fp++;
#endif

#ifdef ADJ_BOX_POWERSPEC
      if(RestartFlag == 6)
	{
	  for(n = 0; n < pc; n++)
	    SphP[offset + n].r.CurlVel = *fp++;
	}
#endif
      break;
    case IO_VORT:
#ifdef ADJ_BOX_POWERSPEC
      if(RestartFlag == 6)
	{
	  for(n = 0; n < pc; n++)
	    for(k = 0; k < 3; k++)
	      SphP[offset + n].Vorticity[k] = *fp++;
	}
#endif
      break;
    case IO_TRUENGB:
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	P[offset + n].TrueNGB = *fp++;
#endif
      break;
    case IO_DPP:
#ifdef JD_DPP
      for(n = 0; n < pc; n++)
	SphP[offset + n].Dpp = *fp++;
#endif
      break;

    case IO_BFLD:		// Magnetic field
#ifdef MAGNETIC
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  SphP[offset + n].b2.BPred[k] = *fp++;
#ifdef TRACEDIVB
      SphP[offset + n].divB = 0;
#endif
#ifdef MAGNETICZERO
      for(k = 0; k < 3; k++)
	SphP[offset + n].b2.BPred[k] = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
      SphP[offset + n].Phi = 0;
      SphP[offset + n].PhiPred = 0;
#endif
#endif
      break;

    case IO_CR_C0:		// Adiabatic invariant for cosmic rays
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; n++)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  SphP[offset + n].CR_C0[CRpop] = *fp++;
#endif
      break;

    case IO_CR_Q0:		// Adiabatic invariant for cosmic rays
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; n++)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  SphP[offset + n].CR_q0[CRpop] = *fp++;
#endif
      break;

    case IO_CR_P0:
      break;

    case IO_CR_E0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; n++)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  SphP[offset + n].CR_E0[CRpop] = *fp++;
#endif
      break;

    case IO_CR_n0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; n++)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  SphP[offset + n].CR_n0[CRpop] = *fp++;
#endif
      break;

    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      break;

    case IO_BHMASS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
#if !defined(DETACH_BLACK_HOLES)
	P[offset + n].BH_Mass = *fp++;
#else
	BHP[N_BH_idx + n].BH_Mass = *fp++;
      N_BH_idx += pc;
#endif
#endif
      break;

    case IO_BHMDOT:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
#if !defined(DETACH_BLACK_HOLES)
	P[offset + n].BH_Mdot = *fp++;
#else
	BHP[N_BH_idx + n].BH_Mdot = *fp++;
      N_BH_idx += pc;
#endif
#endif
      break;

    case IO_BHPROGS:
#ifdef BH_COUNTPROGS
      for(n = 0; n < pc; n++)
#if !defined(DETACH_BLACK_HOLES)
	P[offset + n].BH_CountProgs = *fp++;
#else
	BHP[N_BH_idx + n].BH_CountProgs = *fp++;
      N_BH_idx += pc;
#endif
#endif
      break;

    case IO_BHMBUB:
#ifdef BH_BUBBLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mass_bubbles = *fp++;
#endif
      break;

    case IO_BHMINI:
#ifdef BH_BUBBLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mass_ini = *fp++;
#endif
      break;

    case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mass_radio = *fp++;
#endif
      break;

    case IO_EOSXNUC:
#ifdef EOS_DEGENERATE
      for(n = 0; n < pc; n++)
	for(k = 0; k < EOS_NSPECIES; k++)
	  SphP[offset + n].xnuc[k] = *fp++;
#endif
      break;

    case IO_Zs:
#ifdef LT_STELLAREVOLUTION
      if(type == 4)
	{
	  for(n = 0; n < pc; n++, fp_single += LT_NMetP)
	    memcpy(MetP[N_star_idx + n].Metals, fp_single, LT_NMetP * sizeof(float));
	  N_star_idx += pc;
	}
      else if(type == 0)
	for(n = 0; n < pc; n++, fp_single += LT_NMetP)
	  memcpy(SphP[offset + n].Metals, fp_single, LT_NMetP * sizeof(float));
#endif
      break;

    case IO_ZAGE:
#ifdef LT_ZAGE
      if(type == 4)
	{
	  for(n = 0; n < pc; n++)
	    MetP[N_star_idx + n].ZAge = *fp++;
	  N_star_idx += pc;
	}
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    // note this is not the weight that was used when the snapshot has been written
	    SphP[offset + n].ZAgeW = get_metalmass(SphP[offset + n].Metals);
#ifndef LT_LOGZAGE
	    SphP[offset + n].ZAge = *fp++ * SphP[offset + n].ZAgeW;
#else
	    if(SphP[offset + n].ZAgeW > 0)
	      SphP[offset + n].ZAge = log10(*fp++ * SphP[offset + n].ZAgeW);
	    else
	      SphP[offset + n].ZAge = 0;
#endif
	  }
#endif
      break;

    case IO_ZAGE_LLV:
#ifdef LT_ZAGE_LLV
      if(type == 4)
	{
	  for(n = 0; n < pc; n++)
	    MetP[N_star_idx + n].ZAge_llv = *fp++;
	  N_star_idx += pc;
	}
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    // note this is not the weight that was used when the snapshot has been written
	    SphP[offset + n].ZAgeW_llv = SphP[offset + n].Metals[Iron];
#ifndef LT_LOGZAGE
	    SphP[offset + n].ZAge_llv = *fp++ * SphP[offset + n].ZAgeW_llv;
#else
	    if(SphP[offset + n].ZAgeW_llv > 0)
	      SphP[offset + n].ZAge_llv = log10(*fp++ * SphP[offset + n].ZAgeW_llv);
	    else
	      SphP[offset + n].ZAge_llv = 0;
#endif
	  }
#endif
      break;

    case IO_iMass:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; n++)
	MetP[N_star_idx + n].iMass = *fp++;
      N_star_idx += pc;
#endif
      break;

    case IO_CONTRIB:
#if defined(LT_STELLAREVOLUTION) && defined(LT_TRACK_CONTRIBUTES)
      if(type == 4)
	{
	  for(n = 0; n < pc; n++)
	    MetP[N_star_idx + n].contrib = *contrib++;
	  N_star_idx += pc;
	}
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].contrib = *contrib++;
#endif
      break;

    case IO_ABUNDANCE:
#ifdef LT_MV_CHEMICALDIFFUSION
      for(n = 0; n < pc; n++, fp += LT_NMetP)
	memcpy(SphP[offset + n].Abundance, fp, LT_NMetP * sizeof(float));
#endif
      break;

    case IO_RADGAMMA:
#ifdef RADTRANSFER
      if(RestartFlag != 2)
	{
	  for(n = 0; n < pc; n++)
	    for(k = 0; k < N_BINS; k++)
	      SphP[offset + n].n_gamma[k] = *fp++;
	}
#endif
      break;

    case IO_DMHSML:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].DM_Hsml = *fp_single++;
#endif
      break;

    case IO_DMDENSITY:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].u.DM_Density = *fp_single++;
#endif
      break;

    case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].v.DM_VelDisp = *fp_single++;
#endif
      break;

    case IO_DMHSML_V:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].DM_Hsml_V = *fp_single++;
#endif
      break;

    case IO_DMDENSITY_V:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].DM_Density_V = *fp_single++;
#endif
      break;

    case IO_EULERA:
#ifdef READ_EULER
      for(n = 0; n < pc; n++)
	SphP[offset + n].EulerA = *fp++;
#endif
      break;

    case IO_EULERB:
#ifdef READ_EULER
      for(n = 0; n < pc; n++)
	SphP[offset + n].EulerB = *fp++;
#endif
      break;

    case IO_ALFA2_DYN:
#if defined(FS_ALFA2_DYN) && !defined(FS_ALFA2_TURB)
      for(n = 0; n < pc; n++)
	SphP[offset + n].alfa2 = *fp++;
#endif
      break;

    case IO_ETA2_DYN:
#if defined(FS_ETA2_DYN) && !defined(FS_ETA2_TURB)
      for(n = 0; n < pc; n++)
	SphP[offset + n].eta2 = *fp++;
#endif
      break;


    case IO_CHEM:		// Chemical abundances
#ifdef CHEMCOOL
      for(n = 0; n < pc; n++)
	for(k = 0; k < TRAC_NUM; k++)
	  SphP[offset + n].TracAbund[k] = *fp++;
#endif
      break;

    case IO_CLDX:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; n++)
	SphP[offset + n].XColdCloud = *fp++;
#endif
      break;

    case IO_HTEMP:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; n++)
	SphP[offset + n].Temperature = *fp++;
#endif
      break;

    case IO_TEMP:
#ifdef LT_METAL_COOLING_WAL
      for(n = 0; n < pc; n++)
	SphP[offset + n].Temperature = *fp++;
#endif
      break;

      // the other input fields (if present) are not needed to define the 
 //          initial conditions of the code 

    case IO_SFR:
    case IO_ZSMOOTH:
    case IO_allZSMOOTH:
    case IO_POT:
    case IO_ACCEL:
    case IO_DTENTR:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_TSTP:
    case IO_DBDT:
    case IO_DIVB:
    case IO_ABVC:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_BSMTH:
    case IO_DENN:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_AMDC:
    case IO_PHI:
    case IO_XPHI:
    case IO_GRADPHI:
    case IO_TIDALTENSORPS:
    case IO_ROTB:
    case IO_SROTB:
    case IO_FLOW_DETERMINANT:
    case IO_STREAM_DENSITY:
    case IO_PHASE_SPACE_DETERMINANT:
    case IO_ANNIHILATION_RADIATION:
    case IO_EOSTEMP:
    case IO_PRESSURE:
    case IO_PRESHOCK_CSND:
    case IO_EDDINGTON_TENSOR:
    case IO_LAST_CAUSTIC:
    case IO_VALPHA:
    case IO_HSMS:
    case IO_ACRS:
    case IO_ACRB:
    case IO_BPCR_pNORM:
    case IO_BPCR_eNORM:
    case IO_BPCR_pSLOPE:
    case IO_BPCR_eSLOPE:
    case IO_BPCR_pE:
    case IO_BPCR_pN:
    case IO_BPCR_ePRESSURE:
    case IO_BPCR_pPRESSURE:
    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
    case IO_AGS_SOFT:
    case IO_AGS_DENS:
    case IO_AGS_ZETA:
    case IO_AGS_OMEGA:
    case IO_AGS_CORR:
    case IO_AGS_NGBS:
    case IO_VSTURB_DISS:
    case IO_VSTURB_DRIVE:
    case IO_MG_PHI:
    case IO_MG_ACCEL:
      break;

    case IO_LASTENTRY:
      EndRun(3067,CurrentFile);
      break;
    }
}

#ifdef DoParallel
long long CountTotStarsPar(int snapi,int snapf)
{
long long id,tagid;
long long GlobalSum=0;
int tf,i;
long long c=0;

printf("Task %d started counting!\n",ThisTask);
//if(ThisTask==0)
printf("Total count of tags (including all duplicated particles):%lld\n",GP.TotNumTagsAllSnaps);

if((StellarHaloAllSnaps=(struct tagged_particle *)malloc(GP.TotNumTagsAllSnaps*sizeof(struct tagged_particle)))==NULL)
        {
	//if(ThisTask==0)
        printf("can't allocate memory for all tags!\n");
        EndRun(3130,CurrentFile);
        }
//if(ThisTask==0)
//{
c=0;
for(tf=snapi;tf<snapf;tf++)
        {/*B*/
        hid_t file,dataset,TagDatatype,dataspace;
        size_t size;
        hsize_t dims_out[2];
        herr_t status;
        long unsigned int rows;
        char TagFile[500];
        char *DSName="FullTag";
        struct tagged_particle *StellarHalo;
        int rank,status_n;
        sprintf(TagFile,"%s/tag_%03d.h5",GP.OutputDir,tf);
        /*strcpy(TagFile,"tag_098.h5");*/
        printf("openning tagfile: %s\n",TagFile);
        file = H5Fopen(TagFile, H5F_ACC_RDONLY, H5P_DEFAULT);
        //printf("file:%d\n",file);
        fflush(stdout);
        dataset = H5Dopen(file,DSName,H5P_DEFAULT);
        //printf("dataset:%d\n",dataset);
        TagDatatype=H5Dget_type(dataset);
	//printf("got the type!\n");
        size  = H5Tget_size(TagDatatype);
        //printf("got the size!\n");
        if(size != sizeof(struct tagged_particle))
        {
	//if(ThisTask==0)
            printf("size mismatch, data size:%d, struct tag size:%lu d\n",(int)size,sizeof(struct tagged_particle));
        }
        dataspace = H5Dget_space(dataset);    /* dataspace handle*/
        //printf("got the space!\n");
        rank = H5Sget_simple_extent_ndims(dataspace);
        //printf("got the rank!\n");
        status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
        //printf("got the dims!\n");
        if (status_n<0) 
	  // if(ThisTask==0)
		printf("Error in reading dimensions\n");
        rows = dims_out[0];
        if(rank<0)
	//#ifdef DoParallel
	//if(ThisTask==0)
	//#endif
            printf("I coulldn't read the file:%s\n",TagFile);

        StellarHalo=(struct tagged_particle *)malloc(rows*size);
        status= H5Dread(dataset,TagDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, StellarHalo);
        if(status<0) 
        //if(ThisTask==0)
		printf("Error in reading data from the file:%s\n",TagFile);
        for(i=0;i<rows;i++)
                StellarHaloAllSnaps[c+i]=StellarHalo[i];
	//printf("made it to the end at task: %d!\n",ThisTask);
	free(StellarHalo);
	H5Sclose(dataspace);
        H5Tclose(TagDatatype);
        H5Dclose(dataset);
        H5Fclose(file);
	//printf("All are closed now!");
        c+=rows;
        }/*B*/
//}/* ThisTask==0   */
printf("Task %d finished loading all tagged files on memory!\n",ThisTask);
c=0;
//MPI_Barrier(MPI_COMM_WORLD);
//if(ThisTask==0)
	printf("TotNumPart:%lld, TotNumTagsAll:%lld\n",GP.TotNumPart,GP.TotNumTagsAllSnaps);

//if(ThisTask !=0)
//{// !=T0
for(id=ThisTask;id<GP.TotNumPart;id+=NTask)
        {//C
        for(tagid=0;tagid<GP.TotNumTagsAllSnaps;tagid++)
                {//D
                if(P[id].ID==StellarHaloAllSnaps[tagid].PID)
                        {//E
                        P[id].Type=4;
			P[id].StellarAge+=StellarHaloAllSnaps[tagid].Age;
			P[id].Metallicity+=StellarHaloAllSnaps[tagid].ZZ;
			P[id].StellarMass+=StellarHaloAllSnaps[tagid].StellarMass;
                        c++;
                        }//E
                }//D
        if(id%10000==ThisTask) 
	//if(ThisTask==0)
	     printf("★  %lld x 10000 particles are done with %lld stars. Task:%d ✓\n",id/10000,c,ThisTask);
        }//C

//}// !=T0
//long long GlobalSum=0;
fflush(stdout);
//MPI_Barrier(MPI_COMM_WORLD);
printf("before Allreduce!\n");
fflush(stdout);
//MPI_Allreduce(&c, &GlobalSum, 1, MPI_LONG_LONG_INT, MPI_SUM,MPI_COMM_WORLD);
//}//A
GlobalSum=c;
//fflush(stdout);
printf("before return in task %d!\n",ThisTask);
return GlobalSum;
}

#else

long long CountTotStars(int snapi,int snapf)
{
long long id,tagid;
int tf,i;
long long c=0;
printf("Total count of tags (including all duplicated particles:%lld\n",GP.TotNumTagsAllSnaps);

if((StellarHaloAllSnaps=(struct tagged_particle *)malloc(GP.TotNumTagsAllSnaps*sizeof(struct tagged_particle)))==NULL)
	{	
	printf("can't allocate memory for all tags!\n");
	EndRun(3130,CurrentFile);
	}
c=0;
for(tf=snapi;tf<snapf;tf++)
	{/*B*/
	hid_t file,dataset,TagDatatype,dataspace;
	size_t size;
	hsize_t dims_out[2];
	herr_t status;
	long unsigned int rows;
	char TagFile[500];
	char *DSName="FullTag";
	struct tagged_particle *StellarHalo;
	int rank,status_n;
	sprintf(TagFile,"%s/tag_%03d.h5",GP.OutputDir,tf);
	/*strcpy(TagFile,"tag_098.h5");*/
	/*printf("openning tagfile: %s\n",TagFile);*/
	file = H5Fopen(TagFile, H5F_ACC_RDONLY, H5P_DEFAULT);
	/*printf("file:%d\n",file);*/
	fflush(stdout);
	dataset = H5Dopen(file,DSName,H5P_DEFAULT);
	/*printf("dataset:%d\n",dataset);*/
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
	StellarHalo=(struct tagged_particle *)malloc(rows*size);
	status= H5Dread(dataset,TagDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, StellarHalo);
	if(status<0) printf("Error in reading data from the file:%s\n",TagFile);
	for(i=0;i<rows;i++)
		StellarHaloAllSnaps[c+i]=StellarHalo[i];
	
	free(StellarHalo);
	H5Tclose(TagDatatype);
	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Fclose(file);
	
	c+=rows;
	}/*B*/
printf("I finished loading all tagged files on memory!\n");
c=0;
for(id=0;id<GP.TotNumPart;id++)
	{/*C*/ 
	for(tagid=0;tagid<GP.TotNumTagsAllSnaps;tagid++)
		{/*D*/
		if(P[id].ID==StellarHaloAllSnaps[tagid].PID)
			{/*E*/
			P[id].Type=4;
			c++;
			}/*E*/
		}/*D*/
	if(id%10000==0) printf("%lld x 10000 particles are done with %lld stars\n",id/10000,c);
	}/*C*/

/*}A*/
return c;
}

#endif

void SaveInBinary()
{
int i;
char StellarMassFile[500],MetallicityFile[500];
sprintf(StellarMassFile,"%s/StellarMass",GP.OutputDir);
sprintf(MetallicityFile,"%s/Metallicity",GP.OutputDir);
FILE *SM_File, *M_File;
typedef struct StellarMassRecord{
MyDouble Pos[3];
MyDouble SM;
};//SMRecord;
//typedef struct StellarMassRecord SMRecord;
 if(!(SM_File=fopen(StellarMassFile,"ab")))
        printf("Can't open tagged particles file at %s\n",StellarMassFile);

 if(!(M_File=fopen(MetallicityFile,"ab")))
        printf("Can't open tagged particles file at %s\n",MetallicityFile);

for(i=0;i<GP.TotNumPart;i++)
	if(P[i].Type==4)
	{
//		SMRecord.Pos[0]=P[i].Pos[0];
//		SMRecord.Pos[1]=P[i].Pos[1];
//		SMRecord.Pos[2]=P[i].Pos[2];
//		SMRecord.SM=P[i].StellarMass;
//		fwrite(SMRecord,1,sizeof(SMRecord),SM_File);
	}


fclose(SM_File);
fclose(M_File);
return;
}
