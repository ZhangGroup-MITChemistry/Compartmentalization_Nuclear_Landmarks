/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Arben Jusufi, Axel Kohlmeyer (Temple U.)
------------------------------------------------------------------------- */

// Kartik Version Comment 
// 1. Edited this file to zero out the specific ab potential completely
// 2. Removes the redudancy of 7 atom types and just uses 4 atom types
// 3. Implements mie potential for lamina-chromatin

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_tanhlr_cut_domainab.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairTanhlrCutDomainab::PairTanhlrCutDomainab(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairTanhlrCutDomainab::~PairTanhlrCutDomainab()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(htanh);
    memory->destroy(sigmah);
    memory->destroy(rmh);
    memory->destroy(ptanh);
    memory->destroy(offset);
    memory->destroy(ideal_potential);
  }
}

/* ---------------------------------------------------------------------- */

void PairTanhlrCutDomainab::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,imol,jmol;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,rexp,utanh,factor_lj,tanr,alpha,beta,gamma;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double r2inv,r6inv,sig12nad,sig6nad,sig12spec,sig6spec;  //Used in LJ potential interaction if turned on
  double gamA,gamR,Cmie,sigRlad,sigAlad,rgamA,rgamR; //Used in the mie potential for lamina-chromatin
  //Change the gamA amd gamR in the offset and single() when you want to change the parameters
  double sigwca,cutwca,cutsqwca,sig12wca,sig6wca,offsetwca; //wca repulsion between 0.0 chromatin and lamina

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  //Kartik edit
  //assumed for chr-lamina interaction that sigma of pairs like 1-9,2-9,3-9... (i.e. types within chr do not matter) is the same and precomputes certain quantities.
  //If LJ interaction not used the following terms get computed but not applied.
  //Chr-Lamina
  //Type 1 is Chr Type 9 is Lamina. Mie Potential
  gamR     = 12.0;
  gamA     = 6.0;
  
  Cmie     = ((gamR/(gamR-gamA))*pow((gamR/gamA),(gamA/(gamR-gamA))));
  sigRlad  = pow(sigmah[1][6],gamR);
  sigAlad  = pow(sigmah[1][6],gamA);
  //Adding in the wca parameters 
  sigwca   = sigmah[1][6];
  cutwca   = pow(2.0,1.0/6.0)*sigwca;
  cutsqwca = pow(cutwca,2.0);
  sig12wca = pow(sigwca,12.0);
  sig6wca  = pow(sigwca,6.0);
  offsetwca = (pow(sigwca/cutwca,12.0)-pow(sigwca/cutwca,6.0));
  //Chr-nucelolus via NAD 
  sig12nad = pow(sigmah[1][5],12.0);
  sig6nad  = pow(sigmah[1][5],6.0); 
  //Quantities for chr-speckle interactions. Because sigma of chr-lamina, chr-nucleolus and chr-speckle interaction maybe different.
  //Chr-Speckle
  //Type 1 is Chr Type 10 is Speckles.
  sig12spec = pow(sigmah[1][7],12.0);
  sig6spec  = pow(sigmah[1][7],6.0);
  //printf("%12.8f,%12.8f,%12.8f\n",sigmah[1][9],sigmah[1][8],sigmah[1][10]);
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imol = atom->molecule[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      jmol = atom->molecule[j];

      int iglobal, jglobal;

      iglobal = atom->tag[i];
      jglobal = atom->tag[j];
      // we do not want to apply pair potential for beads in the same segment 

      //Kartik comment
      //Dictionary of arrays and the file names from where they are built 
      //domain --- tad_index_genome_1000kb_diploid.txt
      //active --- different_ab_type_interChrom.txt
      //ideal_potential --- different_ab_type_interChrom_potential_iter00.txt
      
      if ((imol<47) && (jmol<47) && (rsq < cutsq[itype][jtype])) {
      //if (rsq < cutsq[itype][jtype]) {
        
        int iscale = 3;
        // **** add intra-chrom interaction to inter-homologs ****
        // if ( (imol==jmol) || (imol%2==1 && (jmol-imol)==1) || (imol%2==0 && (imol-jmol)==1 ) ) {
        // **** ****

        if (imol==jmol) {
          if (domain[iglobal-1] == domain[jglobal-1]) {
            iscale = 1;
            //Kartik comment
            //printf("Is this executed?\n");
            //Printed this to check that this is never executed
          } else {
            //iscale = 2;
            iscale = 3;
            //printf("This will be executed?\n");
          }
          // 	} else if (abs(jmol-imol)==1){
    	    	// iscale = 3;
    	  } 

    	// else {
    	//     iscale = 1;
  	  //   }

    	else {
        //**** Major change made at this place ****//
    		iscale =3;
        //Override to use only intra interaction throughout
        //iscale = 2;
        /*if ((imol+1)/2 == (jmol+1)/2) {
    			iscale = 3;
          //Kartik comment
          //integer point arithmetic important here
          //The if loop is true if imol and jmol are homologues 1-2,3-4
          //which would mean that the else is executed for all the remaining conditions

          //printf("Is this different from imol==jmol, %d %d %d %d %12.6f %12.6f\n",iglobal,jglobal,imol,jmol,(imol+1)/2,(jmol+1)/2);
    		}
    		else{
          
    		  //iscale = 1;
          //printf("When is this executed? %d %d %d %d %12.6f %12.6f\n",iglobal,jglobal,imol,jmol,(imol+1)/2,(jmol+1)/2);
    		}*/

  	    }
        //Kartik comment 
        //The net effect of the if/else is imol=jmol results in iscale = 2 and for
        //imol!=jmol iscale=3.
        //

        alpha = ptanh[iscale][itype][jtype];
        //Override to use only intra interaction throughout
        

        //***********************************************
        //Following lines are part of the original code
        //Overide of alpha from the file is blocked to include only type specific interaction
        //June 16 2021 Note that for different cell type the active array has to be generated from a different types info file
        /*    
        
        if ( abs(jmol-imol)>1 || ((jmol-imol) == 1 && (imol%2 == 0)) || ((imol-jmol) == 1 && (jmol%2 == 0)) ) {
        //The if loop does not get executed for intra and homologs
        
          if ( active[iglobal-1] > 0 && active[jglobal-1] > 0 ) {
           // This if loop basically checks if the atoms are A and B because types 3 and 4 are denoted as = 0 in the file that was read to intialize active array
           alpha = ideal_potential[active[iglobal-1]-1][active[jglobal-1]-1];

         } else {
          	alpha = ptanh[3][itype][jtype];
            //picks out the third column from generic pair potential
          }
        }
        */ 
        //Kartik comment 
        //brief commentary on how the above loop works 
        //clearly imol=jmol means the loop is not executed and iscale=2 picked from above holds leading to the usage of second column in generic coeff. 
        //say imol = 1 and jmol = 2 which is the homolog interaction between Chr 1 and 1'. None of the three conditions are met for this. so iscale = 3 picked from above leads to the 
        //usage of the third column in generic coeff.  
        //say imol = 1 and jmol = 3 which would have lead to iscale =1 above for value of alpha. The if loop is executed for this condition and file lookups the interaction (between A and B's only) 
        //say itype = 1 and jtype = 3 then the generic third column is used
        //This execution then implies that homolog interactions such as 00 01 11, 22 23 33,  in the ab interaction are never read. 

		// printf("%d %d %d %d %d %d %12.6f\n",imol, jmol, itype, jtype, active[iglobal-1],active[jglobal-1],alpha*2);
        // debug
        //printf("ever here i? %d %d %d %d %d \n", iglobal, jglobal, domain[iglobal], domain[jglobal], iscale);
        //printf("type: %d %d %12.6f %12.6f\n", itype,jtype, ptanh[1][itype][jtype], ptanh[2][itype][jtype]);

        r = sqrt(rsq);
        if (r <= rmh[itype][jtype]) { 
            rexp = (rmh[itype][jtype]-r)*sigmah[itype][jtype];
            // for us, p is 0.5 * h
            tanr = tanh(rexp);
            utanh = alpha*(1.0+ tanr);
            // the extra negative sign is taken care in f
            fpair = factor_lj/r*alpha*sigmah[itype][jtype]*(1-tanr*tanr);
        } else {
            utanh = alpha*rmh4[itype][jtype] / rsq / rsq;
            // there is an extra r here which is to normalize position vector
            fpair = factor_lj/r*alpha*rmh4[itype][jtype]* (4.0) /rsq/rsq/r;
        }


        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          //evdwl = utanh - offset[iscale][itype][jtype];
          //printf("%.6f %d %d %d\n",offset[iscale][itype][jtype],iscale,itype,jtype);
          
          //Correct offset 
          evdwl = utanh - (alpha*offset[iscale][itype][jtype]);
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
      //Kartik edit 
      //I do not want lamina-lamina interaction
      //but want to build in chromatin-lamina interaction
      //also adding in chromatin-speckles interaction
      

      //The lamina molecule id 547 is hard coded for now but can be lifted in future.
      //So is the type id of nucleolus (8)
      //speckle molecule ids > 547
      
      //first two brackets accounts for chr-lamina (at equality) and chr-speckle interaction 
      //second two brackets adds in chr-nucleolus interaction
      //else if ( ((imol>=547) && (jmol<=46)) || ((imol<=46) && (jmol>=547)) || ((itype==8) && (jmol<=46)) || ((imol<=46) && (jtype==8)) ) {
      /// *** The above else if was generalised to suit the changing number of nucleolus because 547 is hard code molecule id of first lamina particle with 500 nucleolus particles*** ///  

      else if ((itype>=5 && jtype<5) || (itype<5 && jtype>=5)){  
        //Can include this in the above condition but did this seperately for clarity
        //the following condition changes forces only if within cutoff
        if ((rsq < cutsq[itype][jtype])){
        //Picks out LAD interaction
        alpha = ptanh[1][itype][jtype];
        //Picks out NAD interaction 
        beta =  ptanh[2][itype][jtype];
        //Picks out Speckle interaction
        gamma = ptanh[3][itype][jtype];
        
        //The way the types are defined... for a given itype and jtype 
        //two of the three quantities above are zero
        
        //printf("Param Check %d %d %12.6f %12.6f %12.6f\n",itype,jtype,alpha,beta,gamma);
        //printf("%.6f,%.6f\n",offset[2][itype][jtype],offset[3][itype][jtype]);
        
        //Identify which one of i,j is chr
        //for use with subcomp
        //subcomp[1][0] is for atom id 1 and first column=LADs
        int chrbead;
        if(imol>=47){chrbead = jglobal;}
        else{chrbead = iglobal;}
        
        
        
        // **** Tanh Potential **** //
        //r = sqrt(rsq);
        //rexp = (rmh[itype][jtype]-r)*sigmah[itype][jtype];
        //tanr = tanh(rexp);
        //utanh = subcomp[chrbead][0]*alpha*(1.0+ tanr);
        //fpair = factor_lj/r*subcomp[chrbead][0]*alpha*sigmah[itype][jtype]*(1-tanr*tanr);


        // **** Lennard-Jones Potential **** //
        alpha = 2.0*alpha;    //Because the factor of 0.5 is premutiplied in ptanh but I want to use the alpha as epsilon here
        beta  = 2.0*beta;
        gamma = 2.0*gamma;
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;

        utanh = 0.0;
        fpair = 0.0;
        //LAD portion via Mie
        if(alpha!=0.0){
        rgamA = pow(r2inv,(gamA/2.0));
        rgamR = pow(r2inv,(gamR/2.0));

        utanh += subcomp[chrbead][0]*Cmie*alpha*((sigRlad*rgamR)-(sigAlad*rgamA));
        fpair += subcomp[chrbead][0]*factor_lj*Cmie*alpha*((gamR*sigRlad*rgamR)-(gamA*sigAlad*rgamA))*r2inv;}  
        //NAD portion 
        if(beta!=0.0){
        utanh += subcomp[chrbead][1]*4.0*beta*r6inv*(sig12nad*r6inv-sig6nad);
        fpair += subcomp[chrbead][1]*factor_lj*4.0*beta*r6inv*((12.0*sig12nad*r6inv)-(6.0*sig6nad))*r2inv;}
        //Speckle portion
        if(gamma!=0.0){
        utanh += subcomp[chrbead][2]*4.0*gamma*r6inv*(sig12spec*r6inv-sig6spec);
        fpair += subcomp[chrbead][2]*factor_lj*4.0*gamma*r6inv*((12.0*sig12spec*r6inv)-(6.0*sig6spec))*r2inv;}
        //printf("%d  %.5f  %d  %d %.14f\n",iglobal,subcomp[chrbead][0],iglobal,jglobal,offset[1][itype][jtype]);
        //printf("%d  %d %.5f\n",itype, jtype,rsq);

        //Adding the WCA repulsion for chromatin beads with exactly 0.0 lamina probability 
        if ((subcomp[chrbead][0]==0.0) && ((itype==6) || (jtype==6)) && (rsq<=cutsqwca) ){
            //Applies forces and energy only till the cutoff of wca i.e. 2^1/6 sigma
            //Note that the wca portion of the energy function is already offset here 
            utanh += (4.0*alpha*r6inv*(sig12wca*r6inv-sig6wca))-(4.0*alpha*offsetwca);
            fpair += factor_lj*4.0*alpha*r6inv*((12.0*sig12wca*r6inv)-(6.0*sig6wca))*r2inv;}

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          //Kartik edit
          //offset is unique to a bead and also determined by the SPIN state prob

          // LJ interaction. Applied offset for each SPIN portion
          // The way this works is for a given itype and jtype two out of three offset terms are zero
          evdwl = utanh -((subcomp[chrbead][0]*offset[1][itype][jtype])+(subcomp[chrbead][1]*offset[2][itype][jtype])+(subcomp[chrbead][2]*offset[3][itype][jtype]));
          evdwl *= factor_lj;
          //printf("%.16f %d %d\n",offset[3][itype][jtype],itype,jtype);
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(htanh,nscale+1,n+1,n+1,"pair:htanh");
  memory->create(sigmah,n+1,n+1,"pair:sigmah");
  memory->create(rmh,n+1,n+1,"pair:rmh");
  memory->create(rmh4,n+1,n+1,"pair:rmh4");
  memory->create(ptanh,nscale+1,n+1,n+1,"pair:ptanh");
  memory->create(offset,nscale+1,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::settings(int narg, char **arg)
{
  // format: tanh/cut/multiscale cutoff nscale lowerB upperB
  //

  //Kartik edit 
  //Passed an extra subcompartment file so changed this
  //if (narg != 6) error->all(FLERR,"Illegal pair_style command");
  if (narg != 7) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  
  nscale = atoi(arg[1]);
  int nChrom = atoi(arg[2]);
    
  // pass through the domain file
  FILE *fp = fopen(arg[3],"r");
  if (fp == NULL)
    error->all(FLERR,"Cannot open the domain file");

  int ntotal;
  ntotal = atom->natoms;

  memory->create(domain,ntotal,"pair:domain");
  

  char line[1024];
  char *ptr;
  int idx;
  // Loop through the rest of the lines
  while(fgets(line,1024,fp)) {
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;

    // skip comment lines starting with #
    if (*ptr == '#') continue;

    idx = atoi(ptr)-1;
    ptr = strtok(NULL," \t\n\r\f");

    // The second site
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted domain file");
    //printf("reading: %d %d\n", idx,atoi(ptr));
    domain[idx] = atoi(ptr); 

  }
  fclose(fp); 

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  // pass through the ab state file
  fp = fopen(arg[4],"r");
  if (fp == NULL)
    error->all(FLERR,"Cannot open the ideal active id file");

  memory->create(active,ntotal,"pair:active");

  // Loop through the rest of the lines
  while(fgets(line,102400,fp)) {
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;

    // skip comment lines starting with #
    if (*ptr == '#') continue;

    idx = atoi(ptr)-1;
    ptr = strtok(NULL," \t\n\r\f");

    // The second site
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted ideal active id file");
    active[idx] = atoi(ptr); 
  }
  fclose(fp); 

  // 
  fp = fopen(arg[5],"r");
  if (fp == NULL)
    error->all(FLERR,"Cannot open the ideal chromosome potential file");
  
  memory->create(ideal_potential,nChrom,nChrom,"pair:ideal");

  int ichr, jchr;
  while(fgets(line,102400,fp)) {
    ptr = strtok(line," \t\n\r\f");
    if (!ptr) continue;
    if (*ptr == '#') continue;
    ichr = atoi(ptr);
    ptr = strtok(NULL," \t\n\r\f");
    jchr = atoi(ptr);
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr) error->all(FLERR,"Incorrectly formatted ideal chromosome potential file");
    ideal_potential[ichr][jchr] = force->numeric(FLERR, ptr) * 0.5;
    ideal_potential[jchr][ichr] = ideal_potential[ichr][jchr];
  }
  fclose(fp);

  //Kartik edit 
  //Subcompartment File
  //Reading all the 4 columns but will be using only the first one (LADs)
  //General form maybe useful in the future
  memory->create(subcomp,ntotal,4,"pair:subcomp");
  fp = fopen(arg[6],"r");
  if (fp == NULL)error->all(FLERR,"Cannot open the subcompartment data file");
  int ibead;
  while(fgets(line,102400,fp)) {
    ptr = strtok(line," \t\n\r\f");
    if (!ptr) continue;
    if (*ptr == '#') continue;
    ibead = atoi(ptr);

    for (int s=0;s<4;s++){
      ptr   = strtok(NULL," \t\n\r\f");
      subcomp[ibead][s] = force->numeric(FLERR, ptr);}
    
    //Note that the first element of subcomp is not used.
    //basically atom id was directly used instead of doing a -1

    //printf("%d\t%.5f\t%.5f\t%.5f\t%.5f\t\n",ibead,subcomp[ibead][0],subcomp[ibead][1],subcomp[ibead][2],subcomp[ibead][3]);
    if (!ptr) error->all(FLERR,"Incorrectly formatted subcompartment file");
  }
  fclose(fp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::coeff(int narg, char **arg)
{
  if (narg < 4+nscale || narg > 5+nscale) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double htanh_one[nscale+1];
  for (int s=1; s <= nscale; s++)
    htanh_one[s] = force->numeric(FLERR,arg[s+1]);
  double rmh_one = force->numeric(FLERR,arg[nscale+2]);
  double sigmah_one = force->numeric(FLERR,arg[nscale+3]);

  double cut_one = cut_global;
  if (narg == nscale+5) cut_one = force->numeric(FLERR,arg[nscale+4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      for (int s = 1; s <= nscale; s++) {
        htanh[s][i][j] = htanh_one[s];
      }
      sigmah[i][j] = sigmah_one;
      rmh[i][j] = rmh_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTanhlrCutDomainab::init_one(int i, int j)
{
  ///Defining parameters for lamina-chromatin Mie
  ///Hard Code. WARNING:Needs to be same gamR and gamA here and in compute()

  double gamA,gamR,Cmie;

  gamR     = 12.0;
  gamA     = 6.0;

  Cmie     = ((gamR/(gamR-gamA))*pow((gamR/gamA),(gamA/(gamR-gamA))));
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  for (int s=1; s <= nscale; s++)
    ptanh[s][i][j] = htanh[s][i][j] * 0.5;

  if (offset_flag) {
    double rexp = (rmh[i][j]-cut[i][j])*sigmah[i][j];
    //printf("%d,%d\n",i,j);
    for (int s=1; s <= nscale; s++) {
        //Kartik edit
        //chromatin-chromatin offset changed to correct one
        //chromatin-chromatin offset used in the if statement of compute()
        //chromatin-lamina, chromatin-nucleolus,chromatin-speckles offset used in the "else if" statement of compute()
        if(i<5 && j<5){
        //Correct offset 
        offset[s][i][j] = pow(rmh[i][j],4.0)/pow(cut[i][j],4.0);
        }
        //chr-lamina/nuc/speck interact with LJ although read from the tanhlr definition from input script
        //However each bead has its own unique offset because each bead has a 
        //unique P_LADs/NADs/Speckles. So computing the offset term without the prob term
        //will evaluate the prob term directly in the force loop
        else{

        //i=1-7   j=8   is  chr-nucleolus
        //i=1-7   j=9   is  chr-lamina
        //i=1-7   j=10  is  chr-speckles
        //every other combination of i,j is never used (for e.g. 8,8) in the compute function
        
        //multiplying 2.0 so that the value entered in input file can be interpretted as epsilon of LJ
        if( s == 1 ){
          //chromatin-LJ is through Mie and hence has a different offset compared to LJ
           offset[s][i][j] = Cmie*(2.0*ptanh[s][i][j])*(pow(sigmah[i][j]/cut[i][j],gamR)-pow(sigmah[i][j]/cut[i][j],gamA));  

        }
        else{
        offset[s][i][j] = 4.0*(2.0*ptanh[s][i][j])*(pow(sigmah[i][j]/cut[i][j],12.0)-pow(sigmah[i][j]/cut[i][j],6.0));  
        }
        //logic for easy interpretation
        //say i=1-7 j=8
        //sigmah is type dependent (and automatically changes) and cut is basically the global cutoff
        //However note that only the second column in pair_coeff will be non zero for this(i=1-7 j=8) interaction and hence 
        //offset[1][i][j]=0 and offset[3][i][j]=0
        //Thus the offset is not only type dependent but the first dimension accounts for LADs/NADs/Speckles
        }
    }
  } else {
    //offset[i][j] = 0.0;
    for (int s=1; s <= nscale; s++) {
        offset[s][i][j] = 0.0;
    }
  }

  for (int s=1; s <= nscale; s++) {
    htanh[s][j][i] = htanh[s][i][j];
    ptanh[s][j][i] = ptanh[s][i][j];
    offset[s][j][i] = offset[s][i][j];
  }
  sigmah[j][i] = sigmah[i][j];
  rmh[j][i] = rmh[i][j];
  rmh4[i][j] = pow(rmh[i][j],4.0);
  rmh4[j][i] = rmh4[i][j];
  cut[j][i] = cut[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        for (int s=1; s <= nscale; s++) {
            fwrite(&htanh[s][i][j],sizeof(double),1,fp);
        }
        fwrite(&rmh[i][j],sizeof(double),1,fp);
        fwrite(&sigmah[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          for (int s=1; s <= nscale; s++) {
            fread(&htanh[s][i][j],sizeof(double),1,fp);
          }
          fread(&rmh[i][j],sizeof(double),1,fp);
          fread(&sigmah[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        for (int s=1; s <= nscale; s++) {
            MPI_Bcast(&htanh[s][i][j],1,MPI_DOUBLE,0,world);
        }
        MPI_Bcast(&rmh[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigmah[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairTanhlrCutDomainab::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{ // This function seems to be never executed in the run
  // I checked this by printing lines over a long run 
  //
  //But I will keep all the lines same as the main compute just in case.
  //It will be quite an expensive function (because of the revaluation 
  //of the LJ bit) if executed repeatedly and the user 
  //will see some slow down. Will also leave a error message.
  double r, rexp,utanh,phitanh, tanr;

  double sigRlad,sigAlad,sig12nad,sig6nad,sig12spec,sig6spec;
  double alpha,beta,gamma;
  double r2inv,r6inv;
  double gamA,gamR,Cmie,rgamA,rgamR;

  int iglobal, jglobal;
  int imol, jmol;
  imol = atom->molecule[i];
  jmol = atom->molecule[j];

  iglobal = atom->tag[i];
  jglobal = atom->tag[j];

  gamR     = 12.0;
  gamA     = 6.0;
  
  Cmie     = ((gamR/(gamR-gamA))*pow((gamR/gamA),(gamA/(gamR-gamA))));
  sigRlad  = pow(sigmah[1][6],gamR);
  sigAlad  = pow(sigmah[1][6],gamA);
  //Chr-nucelolus via NAD 
  sig12nad = pow(sigmah[1][5],12.0);
  sig6nad  = pow(sigmah[1][5],6.0); 
  //Quantities for chr-speckle interactions. Because sigma of chr-lamina, chr-nucleolus and chr-speckle interaction maybe different.
  //Chr-Speckle
  //Type 1 is Chr Type 10 is Speckles.
  sig12spec = pow(sigmah[1][7],12.0);
  sig6spec  = pow(sigmah[1][7],6.0);

  // ptanh index starts with 1
  int iscale = 3;
  //Kartik edit 
  //The following iscale lines should only be executed for chromatin
  //blocking it out for lamina particles
  if (imol==jmol && imol<47 && jmol<47) {
    if (domain[iglobal-1] == domain[jglobal-1]) {
        iscale = 1;
    } else {
        iscale = 3;
    }
  } 

  r=sqrt(rsq);
  //Kartik edit
  if (imol<47 && jmol <47){
    if (r <= rmh[itype][jtype]) {
        rexp = (rmh[itype][jtype]-r)*sigmah[itype][jtype];
        tanr = tanh(rexp);
        utanh = ptanh[iscale][itype][jtype]*(1.0 + tanr);
        fforce = factor_lj/r*ptanh[iscale][itype][jtype]*(1.0-tanr*tanr)*sigmah[itype][jtype];
    } 
    else {
        utanh = ptanh[iscale][itype][jtype]*rmh4[itype][jtype] / rsq / rsq;
        // there is an extra r here which is to normalize position vector
        fforce = factor_lj/r*ptanh[iscale][itype][jtype]*rmh4[itype][jtype]* (4.0) /rsq/rsq/r;
    }
    phitanh = utanh - (ptanh[iscale][itype][jtype]*offset[iscale][itype][jtype]);
  }
  //chr with nuclolus,lamina and speckle
  else if ((itype>=5 && jtype<5) || (itype<5 && jtype>=5)){ 
    alpha = ptanh[1][itype][jtype];
    beta =  ptanh[2][itype][jtype];
    gamma = ptanh[3][itype][jtype];
    int chrbead;
    if(imol>=47){chrbead = jglobal;}
    else{chrbead = iglobal;}
    alpha = 2.0*alpha;    //Because the factor of 0.5 is premutiplied in ptanh but I want to use the alpha as epsilon here
    beta  = 2.0*beta;
    gamma = 2.0*gamma;
    r2inv = 1.0/rsq;
    r6inv = r2inv*r2inv*r2inv;

    utanh   = 0.0;
    fforce  = 0.0;
    //LAD portion
    if(alpha!=0.0){
    rgamA = pow(r2inv,(gamA/2.0));
    rgamR = pow(r2inv,(gamR/2.0));
    utanh  += subcomp[chrbead][0]*Cmie*alpha*((sigRlad*rgamR)-(sigAlad*rgamA));
    fforce  += subcomp[chrbead][0]*factor_lj*Cmie*alpha*((gamR*sigRlad*rgamR)-(gamA*sigAlad*rgamA))*r2inv;}  
    //NAD portion 
    if(beta!=0.0){
    utanh   += subcomp[chrbead][1]*4.0*beta*r6inv*(sig12nad*r6inv-sig6nad);
    fforce  += subcomp[chrbead][1]*factor_lj*4.0*beta*r6inv*((12.0*sig12nad*r6inv)-(6.0*sig6nad))*r2inv;}
    //Speckle portion
    if(gamma!=0.0){
    utanh   += subcomp[chrbead][2]*4.0*gamma*r6inv*(sig12spec*r6inv-sig6spec);
    fforce  += subcomp[chrbead][2]*factor_lj*4.0*gamma*r6inv*((12.0*sig12spec*r6inv)-(6.0*sig6spec))*r2inv;}

    phitanh = utanh -((subcomp[chrbead][0]*offset[1][itype][jtype])+(subcomp[chrbead][1]*offset[2][itype][jtype])+(subcomp[chrbead][2]*offset[3][itype][jtype]));
  } 

  //printing just to check if this function is executed
  printf("CUSTOM CODE WARNING:Code block within single.");
  
  return factor_lj*phitanh;

  //This function is never executed in a run
  //Look for the above error message in the log file to check if it ever got executed.
}

/* ---------------------------------------------------------------------- */
double PairTanhlrCutDomainab::memory_usage()
{
  const int n=atom->ntypes;

  double bytes = Pair::memory_usage();

  bytes += 7*((n+1)*(n+1) * sizeof(double) + (n+1)*sizeof(double *));
  bytes += 1*((n+1)*(n+1) * sizeof(int) + (n+1)*sizeof(int *));

  return bytes;
}

