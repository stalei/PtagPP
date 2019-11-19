/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


void read_trees(int totmal)
{
    int i;
    int fnr;
    int haloindex;
    int treeindex;
    long long mID;
    int subfilnum;
    int subindex;
    int MIdx;


    Treefiles();
    
    for(i=0; i<totmal; i++)
    { 
       fnr = SageOutput[i].FileNr;
       treeindex = SageOutput[i].TreeIndex;
       haloindex = SageOutput[i].HaloIndex;
       AllGal[i].DiskScaleRadius = SageOutput[i].DiskScaleRadius;
       AllGal[i].BulgeMass = SageOutput[i].BulgeMass;
       AllGal[i].BlackHoleMass = SageOutput[i].BlackHoleMass;
       AllGal[i].StellarMass = SageOutput[i].StellarMass;
       AllGal[i].ColdGas = SageOutput[i].ColdGas;

       malhalo(fnr);
       read_haloheader(fnr);

       read_halotree(fnr);

       mID = Find_mostbound(haloindex, treeindex);
       AllGal[i].mostbound = mID;

       subfilnum = Find_subfnr(haloindex, treeindex);
       AllGal[i].subfnum = subfilnum;

       subindex = Find_subhidx(haloindex, treeindex);
       AllGal[i].subhaloidx = subindex;

      MIdx = Find_MIdx(haloindex, treeindex);
 
      AllGal[i].MhaloIdx = MIdx;
   
       free(halo_tree);
       free(hptree);

       halo_tree = NULL;
       hptree = NULL;
     }
    
    return;
}

void malhalo(int fnr)
{
    int haloNtrees;
    int Ntothalos;
    char file1[1000];
    FILE *ht;

    sprintf(file1, "%s", treepath[fnr].tpaths);
        
    //Open Tree file
         
    ht = fopen(file1, "rb");
        
    if(NULL == ht)
    {
        perror("Cannot open tree file");
        return;
    }
    

    //Read in the header
    fread(&haloNtrees, sizeof(int), 1, ht);
    fread(&Ntothalos, sizeof(int), 1, ht);

    //Part of the header allocate the needed array size
    if((hptree = (struct Halos_per_tree*)malloc(haloNtrees*sizeof(struct Halos_per_tree))) == NULL)
    {
         perror("Failed to allocate for HaloTree header....");
         return;
    }

    if((halo_tree = (struct Halo_tree_read*)malloc(Ntothalos*sizeof(struct Halo_tree_read))) == NULL)
    {
         perror("Failed to allocate for HaloTree....");
         return;
    }
 
    fclose(ht);
 
    return;
}


void read_haloheader(int fnr)
{
    char file1[1000];
    FILE *ht;

    sprintf(file1, "%s", treepath[fnr].tpaths);
        
    //Open Tree file
         
    ht = fopen(file1, "rb");
        
    if(NULL == ht)
    {
        perror("Cannot open tree file");
        return;
    }

    halohead = (struct Halo_tree_header*)malloc(sizeof(struct Halo_tree_header));

    fread(&halohead[0], sizeof(struct Halo_tree_header), 1, ht);

    //Read in the header
    fread(&hptree[0], sizeof(struct Halos_per_tree),halohead[0].haloNtrees,ht);
    
    
    fclose(ht);

    free(halohead);
    halohead = NULL;
    return;

}


void read_halotree(int fnr)
{

    int haloNtrees;
    int Ntothalos;
    int *tree;
    char file1[1000];
    FILE *ht;

   
    sprintf(file1, "%s", treepath[fnr].tpaths);
        
    //Open Tree file
         
    ht = fopen(file1, "rb");
        
    if(NULL == ht)
    {
        perror("Cannot open tree file");
        return;
    }

    //Read in the header
    fread(&haloNtrees, sizeof(int), 1, ht);
    fread(&Ntothalos, sizeof(int), 1, ht);
    

    //Part of the header allocate the needed array size
    tree = (int*)malloc(haloNtrees*sizeof(int));

   

    //Read in the rest of the header into an array
     
    fread(&tree[0], sizeof(int),haloNtrees,ht);

    fread(&halo_tree[0], sizeof(struct Halo_tree_read), Ntothalos, ht);

   
    fclose(ht);
    
    free(tree);
    tree = NULL;
    return;

}

long long Find_mostbound(int haloindex, int treeindex)
{
    int numthalos;
    int i;
    int addhalos;
    int Shalo;
    long long mboundID;
      
    addhalos = 0;
    
    numthalos = hptree[treeindex].halospertree;
    
    for(i=0; i<treeindex+1; i++)
    {
        addhalos = hptree[i].halospertree + addhalos;
        
    }

    Shalo = addhalos - numthalos + haloindex;

    mboundID = halo_tree[Shalo].MostBoundID;


    return mboundID;

}

int Find_subhidx(int haloindex, int treeindex)
{
    int numthalos;
    int i;
    int addhalos;
    int Shalo;
    int subhaloindex;
      
    addhalos = 0;
    
    numthalos = hptree[treeindex].halospertree;
    
    for(i=0; i<treeindex+1; i++)
    {
        addhalos = hptree[i].halospertree + addhalos;
        
    }

    Shalo = addhalos - numthalos + haloindex;

    subhaloindex = halo_tree[Shalo].SubhaloIndex;

    //printf("Halo Read:  haloidx = %d, treeidx = %d, addhalo = %d, numthalos = %d, Shalo = %d, subhaloidx = %d\n", haloindex, treeindex, addhalos, numthalos, Shalo, subhaloindex);

    return subhaloindex;

}


int Find_subfnr(int haloindex, int treeindex)
{
    int numthalos;
    int i;
    int addhalos;
    int Shalo;
    int subf_fnr;
      
    addhalos = 0;
    
    numthalos = hptree[treeindex].halospertree;
    
    for(i=0; i<treeindex+1; i++)
    {
        addhalos = hptree[i].halospertree + addhalos;
        
    }

    Shalo = addhalos - numthalos + haloindex;  

    subf_fnr = halo_tree[Shalo].FileNr;

    return subf_fnr;

}

int Find_MIdx(int haloindex, int treeindex)
{
    int numthalos;
    int i;
    int addhalos;
    int Shalo;
    int MIdx;

    addhalos = 0;

    numthalos = hptree[treeindex].halospertree;

    for(i=0; i<treeindex+1; i++)
    {
        addhalos = hptree[i].halospertree + addhalos;

    }

    Shalo = addhalos - numthalos + haloindex;

    MIdx = halo_tree[Shalo].SubhaloIndex;

    return MIdx;
}
