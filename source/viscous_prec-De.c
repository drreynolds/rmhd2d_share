/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Daniel R. Reynolds
 * UC San Diego, Mathematics
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "viscous_prec.h"     /* precond. routine prototypes      */


static void *VPDeData;


/********************************************************************/
/* Fortran callable interface routines                              */
/********************************************************************/


/* De preconditioner dataspace allocation wrapper routine */
void VISCPREC_DE_INIT(long int *Nx, long int *Ny, double *dx, double *dy,
		      long int *Ns, long int *NGx, long int *NGy, 
		      int *NPx, int *NPy, int *iPx, int *iPy, 
		      int *NBlt, int *NBrt, int *NBtp, int *NBbt, 
		      int *XBcond, int *YBcond, int *ier)
{
  /* allocate preconditioner data */
  VPDeData = VPrecDeAlloc(*Nx, *Ny, *dx, *dy, *Ns, *NGx, *NGy, *NPx, 
			  *NPy, *iPx, *iPy,  *NBlt, *NBrt, *NBtp, *NBbt, 
			  *XBcond, *YBcond);
  if (VPDeData == NULL) *ier = -1; 
  else                  *ier = 0;

  return;
}



/* De preconditioner dataspace deallocation wrapper routine */
void VISCPREC_DE_FREE()
{
  VPrecDeFree(VPDeData);
  return;
}



/* De preconditioner setup wrapper routine */
void VISCPREC_DE_SETUP(double *uu, double *gamdt, 
		       double *Gamma, double *Kappa, 
		       double *Re, double *Pr, double *RGas, 
		       double *v1, double *v2, int *ier)
{
  /* call the C preconditioner setup routine */
  *ier = VPrecDeSetup(uu, *gamdt, *Gamma, *Kappa, *Re, 
		      *Pr, *RGas, v1, v2, VPDeData);
  return;
}



/* energy preconditioner solve wrapper routine */
void VISCPREC_DE_SOLVE(double *xx, double *bb, double *tmp, double *delta, int *ier)
{
  /* call the C preconditioner solve routine */
  *ier = VPrecDeSolve(xx, bb, tmp, *delta, VPDeData);
  return;
}



/* energy preconditioner multiply wrapper routine */
void VISCPREC_DE_MULTIPLY(double *xx, double *bb, double *tmp, int *ier)
{
  /* call the C preconditioner multiply routine */
  *ier = VPrecDeMultiply(xx, bb, tmp, VPDeData);
  return;
}



/* De preconditioner options routine */
void SET_SOL_DE_OPTS(int *iopt, int *ier)
{
  if (VPDeData == NULL) *ier = 1;
  else {
    /* cast VPDeData as the correct structure */
    ViscPrecDeData pdata;
    pdata = (ViscPrecDeData) VPDeData;
    /* set maximum number of iterations into pdata */
    if (iopt[0] != 0) pdata->sol_maxit = iopt[0];
    /* set relative change stopping criteria into pdata */
    if (iopt[1] != 0) pdata->sol_relch = iopt[1];
    /* set relaxation type into pdata */
    if (iopt[2] != 0) pdata->sol_rlxtype = iopt[2];
    /* set number of pre-relaxation sweeps into pdata */
    if (iopt[3] != 0) pdata->sol_npre = iopt[3];
    /* set number of post-relaxation sweeps into pdata */
    if (iopt[4] != 0) pdata->sol_npost = iopt[4];
    /* set print level into pdata */
    if (iopt[5] != 0) pdata->sol_printl = iopt[5];
    /* set amount of logging into pdata */
    if (iopt[6] != 0) pdata->sol_log = iopt[6];
    /* set choice of zero initial guess */
    if (iopt[7] != 0) pdata->sol_zeroguess = 1;
    *ier = 0;
  }
  return;
}



/* De preconditioner diagnostic output routine */
void VISCPREC_DE_NUMITERS(int *Niters)
{
  *Niters = VPrecDeNumIters(VPDeData);
  return;
}




/********************************************************************/
/* Internal Preconditioner Routines                                 */
/********************************************************************/


/* -----------------------------------------------------------------
 * Function : VPrecDeAlloc
 * -----------------------------------------------------------------
 * VPrecDeAlloc is called at initialization to set aside space for 
 * any internal storage that will be required by VPrecDeSetup and
 * VPrecDeSolve.
 * -------------------------------------------------------------- */
void *VPrecDeAlloc(long int Nx, long int Ny, double dx, double dy, 
		   long int Ns, long int NGx, long int NGy, int NPx, 
		   int NPy, int iPx, int iPy, int NBlt, int NBrt, 
		   int NBtp, int NBbt, int XBcond, int YBcond) 
{
  /* define necessary local variables, output variable */
  ViscPrecDeData pdata;
  int xtag, ytag;
  int ilower[2], iupper[2];
  long int Nxl, Nyl;

  /* allocate preconditioner data, cast as ViscPrecData */
  pdata = (ViscPrecDeData) malloc(sizeof *pdata);
  if (pdata == NULL) return(NULL);

  /* local domain information */
  pdata->ndim = 2;          /* 2D grid */
  pdata->Nx   = Nx;         /* num points in x-dir. */
  pdata->Ny   = Ny;         /* num points in y-dir. */
  pdata->dx   = dx;         /* mesh size in x-dir.  */
  pdata->dy   = dy;         /* mesh size in y-dir.  */
  pdata->Ns   = Ns;         /* num species          */
  pdata->NGx  = NGx;        /* num x-ghost points   */
  pdata->NGy  = NGy;        /* num y-ghost points   */
  pdata->xbc  = XBcond;     /* x-boundary condition */
  pdata->ybc  = YBcond;     /* y-boundary condition */

  /* processor layout information */
  pdata->Xprocs = NPx;      /* num procs in x-dir. */
  pdata->Yprocs = NPy;      /* num procs in y-dir. */
  pdata->iprocx = iPx;      /* x-loc in proc. grid */
  pdata->iprocy = iPy;      /* y-loc in proc. grid */
  if ((iPx==1) && (iPy==1)) {pdata->outproc=1;}
  else {pdata->outproc=0;}

  /* global domain information */
  Nxl  = 0;    /* size of left neighbor's x-grid */
  Nyl  = 0;    /* size of bottom neighbor's y-grid */
  xtag = 101;  /* x-communication tag */
  ytag = 202;  /* y-communication tag */
#ifdef PARALLEL
  pdata->comm = MPI_COMM_WORLD;
  MPI_Status status;
  if (iPx < NPx) 
    {MPI_Send(&Nx, 1, MPI_LONG, NBrt, xtag, pdata->comm);}
  if (iPx > 1) 
    {MPI_Recv(&Nxl, 1, MPI_LONG, NBlt, xtag, pdata->comm, &status);}
  if (iPy < NPy) 
    {MPI_Send(&Ny, 1, MPI_LONG, NBtp, ytag, pdata->comm);}
  if (iPy > 1) 
    {MPI_Recv(&Nyl, 1, MPI_LONG, NBbt, ytag, pdata->comm, &status);}
#else
  pdata->comm = 0;
#endif
  pdata->iXL = Nxl*(iPx-1) + 1;
  pdata->iXR = Nxl*(iPx-1) + Nx;
  pdata->iYL = Nyl*(iPy-1) + 1;
  pdata->iYR = Nyl*(iPy-1) + Ny;


  /* HYPRE-specific information */

  /*    set the general domain information */
  pdata->totIters = 0;
  pdata->Tsetup = 0.0;
  pdata->Tsolve = 0.0;

  /*    set up the grid */
  /*       create grid object */
  HYPRE_StructGridCreate(pdata->comm, pdata->ndim, &(pdata->grid));

  /*       set my grid extents as if we have one part with multiple boxes.
	   Have each processor describe it's own global extents */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_StructGridSetExtents(pdata->grid, ilower, iupper);

  /*       set grid periodicity */
  int periodicity[2] = {0, 0};
  if (XBcond == 1) {
    long int xlo_idx, xhi_idx;
    MPI_Allreduce(&(pdata->iXL), &xlo_idx, 1, MPI_LONG, MPI_MIN, pdata->comm);
    MPI_Allreduce(&(pdata->iXR), &xhi_idx, 1, MPI_LONG, MPI_MAX, pdata->comm);
    periodicity[0] = xhi_idx-xlo_idx+1;
  }
  if (YBcond == 1) {
    long int ylo_idx, yhi_idx;
    MPI_Allreduce(&(pdata->iYL), &ylo_idx, 1, MPI_LONG, MPI_MIN, pdata->comm);
    MPI_Allreduce(&(pdata->iYR), &yhi_idx, 1, MPI_LONG, MPI_MAX, pdata->comm);
    periodicity[1] = yhi_idx-ylo_idx+1;
  }
  HYPRE_StructGridSetPeriodic(pdata->grid, periodicity);


  /*       assemble the grid */
  HYPRE_StructGridAssemble(pdata->grid);
  

  /*    set up the stencil */
  /*       create e stencil */
  HYPRE_StructStencilCreate(pdata->ndim, 5, &(pdata->eStencil));
  

  /*       set stencil entries */
  int offset[2];
  /*         dependency to bottom */
  offset[0] = 0;  offset[1] = -1;
  HYPRE_StructStencilSetElement(pdata->eStencil, 0, offset);
  /*         dependency to left */
  offset[0] = -1;  offset[1] = 0;
  HYPRE_StructStencilSetElement(pdata->eStencil, 1, offset);
  /*         dependency to self */
  offset[0] = 0;  offset[1] = 0;
  HYPRE_StructStencilSetElement(pdata->eStencil, 2, offset);
  /*         dependency to right */
  offset[0] = 1;  offset[1] = 0;
  HYPRE_StructStencilSetElement(pdata->eStencil, 3, offset);
  /*         dependency to top */
  offset[0] = 0;  offset[1] = 1;
  HYPRE_StructStencilSetElement(pdata->eStencil, 4, offset);


  /*    set up the matrix and vectors */
  /*       create/initialize the matrix */
  HYPRE_StructMatrixCreate(pdata->comm, pdata->grid, pdata->eStencil, &(pdata->De));
  HYPRE_StructMatrixInitialize(pdata->De);

  /*       create/initialize the vectors */
  HYPRE_StructVectorCreate(pdata->comm, pdata->grid, &(pdata->bvec));
  HYPRE_StructVectorCreate(pdata->comm, pdata->grid, &(pdata->xvec));
  HYPRE_StructVectorInitialize(pdata->bvec);
  HYPRE_StructVectorInitialize(pdata->xvec);
  

  /********************************/
  /*  continue with general setup */

  /* set init flag to 0 at first */
  pdata->DeInit = 0;

  /* set Solver default options into pdata */
  pdata->sol_maxit     = 50;
  pdata->sol_relch     = 0;
  pdata->sol_rlxtype   = 1;
  pdata->sol_npre      = 1;
  pdata->sol_npost     = 1;
  pdata->sol_printl    = 1;
  pdata->sol_log       = 1;
  pdata->sol_zeroguess = 0;

  /* return preconditioner data, recast as void *. */
  return((void *)pdata);
}



/* -----------------------------------------------------------------
 * Function : VPrecDeFree
 * -----------------------------------------------------------------
 * VPrecDeFree frees storage allocated by VPrecDeAlloc
 * -------------------------------------------------------------- */
void VPrecDeFree(void *P_data)
{
  /* ensure that P_data is non-null, and free space as required */
  if ( P_data != NULL ) {

    /* cast P_data as the correct structure */
    ViscPrecDeData pdata;
    pdata = (ViscPrecDeData) P_data;

    /* free the HYPRE objects */ 
    HYPRE_StructMatrixDestroy(pdata->De);
    HYPRE_StructVectorDestroy(pdata->bvec);
    HYPRE_StructVectorDestroy(pdata->xvec);

    /* finally, free the pdata structure */
    free(pdata);
  }
}



/* -----------------------------------------------------------------
 * Function : VPrecDeSetup
 * -----------------------------------------------------------------
 * VPrecDeSetup sets up the viscous preconditioning matrix for the 
 * momentum equations.
 *
 * The parameters of VPrecSetup used here are as follows:
 *
 * uu      is the current state of the system
 * gamdt   is the time scaling in the Newton matrix, M = I+gamdt*J
 * Gamma   is the plasma specific heat ratio
 * kappa   is the plasma heat conductivity
 * Re      is the plasma Reynolds number
 * Pr      is the plasma Prandtl number
 * RGas    is the plasma gas constant
 * v1, v2  give already allocated arrays which may be used as 
 *         temporary storage or work space (of length Nx*Ny*Nz*8)
 * pdata   is the preconditioner data returned by VPrecAlloc.
 *
 * Return value:
 * The value returned by this VPrecDeSetup function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDeSetup(double *uu, double gamdt, double Gamma, 
		 double Kappa, double Re, double Pr, 
		 double RGas, double *v1, double *v2, void *P_data)
{
  /* local variables */
  long int Nx, Ny, NGx, NGy;
  long int ix, iy, idx, Ybl, Ybl_p, Ybl_m;
  int ilower[2], iupper[2], entries[5], IdLoc;
  int xLface, xRface, yLface, yRface;
  double dxi2, dyi2, dxsqsum, kappafact, IdVal, Tstart, Tstop;

  /* start timer */
  Tstart = MPI_Wtime();

  /* recast P_data as the correct structure */
  ViscPrecDeData pdata;
  pdata = (ViscPrecDeData) P_data;

  /* get grid information shortcuts */
  Nx  = pdata->Nx;
  Ny  = pdata->Ny;
  NGx = pdata->NGx;
  NGy = pdata->NGy;
  dxi2 = 1.0/(pdata->dx)/(pdata->dx);
  dyi2 = 1.0/(pdata->dy)/(pdata->dy);
  dxsqsum = 2.0*(dxi2+dyi2);

  /* kappafact contains parameter & time step scaling */
  kappafact = -gamdt*Gamma*Kappa/(0.7)/RGas;

  /* entries holds stencil locations */
  entries[0] = 0;
  entries[1] = 1;
  entries[2] = 2;
  entries[3] = 3;
  entries[4] = 4;

  /* IdLoc, IdVal hold identity location, value */
  IdLoc = 2;
  IdVal = 1.0;

  /* set flags determining whether proc owns external faces */
  xLface = (pdata->iprocx == 1);
  xRface = (pdata->iprocx == pdata->Xprocs);
  yLface = (pdata->iprocy == 1);
  yRface = (pdata->iprocy == pdata->Yprocs);


  /* set matrix values over grid */
  /*       internal cells */
  idx = 0;
  for (iy=0; iy<Ny; iy++) {
    Ybl = (iy+NGy) * (Nx+2*NGx);
    Ybl_p = (iy+NGy+1) * (Nx+2*NGx);
    Ybl_m = (iy+NGy-1) * (Nx+2*NGx);
    for (ix=0; ix<Nx; ix++) {
      /* set entry to bottom */
      v2[idx++] = kappafact*dyi2/uu[Ybl_m + ix + NGx];
      
      /* set entry to left */
      v2[idx++] = kappafact*dxi2/uu[Ybl + ix + NGx - 1];
      
      /* set entry to self */
      v2[idx++] = 1.0 - kappafact*dxsqsum/uu[Ybl + ix + NGx];
      
      /* set entry to right */
      v2[idx++] = kappafact*dxi2/uu[Ybl + ix + NGx + 1];
      
      /* set entry to top */
      v2[idx++] = kappafact*dyi2/uu[Ybl_p + ix + NGx];
    }
  }
  if (pdata->xbc != 1) {     /* not periodic */
    /*       ix=0 face adjustment */
    if (xLface) {
      ix=0;
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	/* 0-gradient and reflecting BCs are the same for energy
	   vals: move stencil coupling from ghost zone to self */
	v2[idx+2] += v2[idx+1];  v2[idx+1] = 0.0;
      }
    }
    /*       ix=Nx-1 face adjustment */
    if (xRface) {
      ix=Nx-1;
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	/* 0-gradient and reflecting BCs are the same for energy
	   vals: move stencil coupling from ghost zone to self */
	v2[idx+2] += v2[idx+3];  v2[idx+3] = 0.0;
      }
    }
  }
  if (pdata->ybc != 1) {     /* not periodic */
    /*       iy=0 face adjustment */
    if (yLface) {
      iy=0;
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	/* 0-gradient and reflecting BCs are the same for energy
	   vals: move stencil coupling from ghost zone to self */
	v2[idx+2] += v2[idx];  v2[idx] = 0.0;
      }
    }
    /*       iy=Ny-1 face adjustment */
    if (yRface) {
      iy=Ny-1;
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	/* 0-gradient and reflecting BCs are the same for energy
	   vals: move stencil coupling from ghost zone to self */
	v2[idx+2] += v2[idx+4];  v2[idx+4] = 0.0;
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_StructMatrixSetBoxValues(pdata->De, ilower, iupper, 5, entries, v2);

  /* assemble matrix */
  HYPRE_StructMatrixAssemble(pdata->De);

  /*       set init flag */
  pdata->DeInit = 1;

/*   /\* PRINT OUT THE MATRIX TO FILES FOR BUG-CHECKING *\/ */
/*   if ((pdata->outproc)==1) */
/*     {printf("      printing HYPRE De matrix to file \n");} */
/*   HYPRE_StructMatrixPrint("De.mat", pdata->De, 0); */
    
  /* stop timer, add to totals, and output to screen */ 
  Tstop = MPI_Wtime();
  pdata->Tsetup += Tstop - Tstart;
  if ((pdata->outproc) == 1) {
    printf("De cumulative setup time: %g\n",pdata->Tsetup);
    printf("De cumulative solve time: %g\n",pdata->Tsolve);
  }

  /* return success */ 
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDeSolve
 * -----------------------------------------------------------------
 * VPrecDeSolve solves a linear system P x = b, with the
 * preconditioner matrix P generated by VPrecSetup and solved
 * using the HYPRE library.
 *
 * The parameters of VPrecDeSolve used here are as follows:
 *
 * xx      is the rhs vector on input
 * bb      is the sol vector on output
 * tmp     is a temporary vector the size of xx
 * delta   is the desired linear solve tolerance (if used in 
 *         some iterative method)
 * pdata   is the pre-computed De preconditioner data
 *
 * The value returned by this VPrecDeSolve function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDeSolve(double *xx, double *bb, double *tmp, double delta, void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDeData pdata;
  pdata = (ViscPrecDeData) P_data;

  /* local variables */
  int Sits, Pits, ilower[2], iupper[2];
  long int ix, iy, idx, Nx, NGx, Ny, NGy;
  long int Vbl, Ybl;
  double finalresid, resid, val, bnorm, Tstart, Tstop;
  HYPRE_StructSolver solver;
  HYPRE_StructSolver preconditioner;
  int printl = ((pdata->outproc)==1) ? pdata->sol_printl : 0;

  /* start timer */
  Tstart = MPI_Wtime();

  /* check that De matrix initialized */
  if (pdata->DeInit == 0) {
    printf("VPrecDeSolve error: De matrix uninitialized!\n");
    return(1);
  }

  /* set local variables */
  Nx  = pdata->Nx;
  NGx = pdata->NGx;
  Ny  = pdata->Ny;
  NGy = pdata->NGy;

  /* convert rhs, solution vectors to HYPRE format      */
  /*    insert rhs vector entries into HYPRE vector b   */
  /*    and sol vector entries into HYPRE vector x      */
  Vbl = 7 * (Nx+2*NGx) * (Ny+2*NGy);
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  idx = 0;
  for (iy=NGy; iy<Ny+NGy; iy++) {
    Ybl = iy * (Nx+2*NGx);
    for (ix=NGx; ix<Nx+NGx; ix++)
      tmp[idx++] = bb[Vbl+Ybl+ix];
  }
  HYPRE_StructVectorSetBoxValues(pdata->bvec, ilower, iupper, tmp);
  for (idx=0; idx<Nx*Ny; idx++)  tmp[idx] = 0.0;
  HYPRE_StructVectorSetBoxValues(pdata->xvec, ilower, iupper, tmp);

  /*    assemble vectors */
  HYPRE_StructVectorAssemble(pdata->xvec);
  HYPRE_StructVectorAssemble(pdata->bvec);

/*   if (printl) printf("Writing out rhs to file e_rhs.vec\n"); */
/*   HYPRE_StructVectorPrint("e_rhs.vec", pdata->bvec, 0); */

  /* set up the solver [PCG/PFMG] */
  /*    create the solver & preconditioner */
  HYPRE_StructPCGCreate(pdata->comm, &solver);
  HYPRE_StructPFMGCreate(pdata->comm, &preconditioner);
 
  /*    set preconditioner & solver options */
  HYPRE_StructPFMGSetMaxIter(preconditioner, pdata->sol_maxit/4);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner, pdata->sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner, pdata->sol_npost);
  HYPRE_StructPCGSetPrintLevel(solver, pdata->sol_printl);
  HYPRE_StructPCGSetLogging(solver, pdata->sol_log);
  HYPRE_StructPCGSetRelChange(solver, 1);
  HYPRE_StructPCGSetMaxIter(solver, pdata->sol_maxit);
  if (delta != 0.0)  HYPRE_StructPCGSetTol(solver, delta);
  HYPRE_StructPCGSetPrecond(solver, 
			    (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
			    (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
			    preconditioner);
  HYPRE_StructPCGSetup(solver, pdata->De, pdata->bvec, pdata->xvec);
  
  /* solve the linear system */
  HYPRE_StructPCGSolve(solver, pdata->De, pdata->bvec, pdata->xvec);

  /* extract solver statistics */
  finalresid = 1.0;  Sits = 0;  Pits = 0;
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_StructPCGGetNumIterations(solver, &Sits);
  HYPRE_StructPFMGGetNumIterations(preconditioner, &Pits);
  pdata->totIters += Sits;
  if (printl) printf("      De lin resid = %.1e (tol = %.1e) its = (%i,%i) \n",
		     finalresid, delta, Sits, Pits);

/*   if (printl) printf("Writing out solution to file e_sol.vec\n"); */
/*   HYPRE_StructVectorPrint("e_sol.vec", pdata->xvec, 0); */

  /* extract values from solution vector  */
  Vbl = 7 * (Nx+2*NGx) * (Ny+2*NGy);
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_StructVectorGetBoxValues(pdata->xvec, ilower, iupper, tmp);
  idx = 0;
  for (iy=NGy; iy<Ny+NGy; iy++) {
    Ybl = iy * (Nx+2*NGx);
    for (ix=NGx; ix<Nx+NGx; ix++) 
      xx[Vbl+Ybl+ix] = tmp[idx++];
  }

  /* destroy solver and preconditioner structures */
  HYPRE_StructPCGDestroy(solver);
  HYPRE_StructPFMGDestroy(preconditioner);
      
  /* stop timer, add to totals, and output to screen */ 
  Tstop = MPI_Wtime();
  pdata->Tsolve += Tstop - Tstart;

  /* return success.  */
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDeMultiply
 * -----------------------------------------------------------------
 * VPrecDeMultiply performs the matrix-vector product P x = b, with 
 * the preconditioner matrix P generated by VPrecSetup and solved
 * using the HYPRE library.
 *
 * The parameters of VPrecDeMultiply used here are as follows:
 *
 * xx      is the x vector on input
 * bb      is the b vector on output
 * tmp     is a temporary vector the size of xx
 * pdata   is the pre-computed De preconditioner data
 *
 * The value returned by this VPrecDeMultiply function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDeMultiply(double *xx, double *bb, double *tmp, void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDeData pdata;
  pdata = (ViscPrecDeData) P_data;

  /* local variables */
  int ilower[2], iupper[2];
  long int ix, iy, idx, Nx, NGx, Ny, NGy;
  long int Vbl, Ybl;
  double val;

  /* check that De matrix initialized */
  if (pdata->DeInit == 0) {
    printf("VPrecDeMultiply error: De matrix uninitialized!\n");
    return(1);
  }

  /* set local variables */
  Nx  = pdata->Nx;
  NGx = pdata->NGx;
  Ny  = pdata->Ny;
  NGy = pdata->NGy;

  /* convert product, result vectors to HYPRE format        */
  /*    insert product vector entries into HYPRE vector x   */
  Vbl = 7 * (Nx+2*NGx) * (Ny+2*NGy);
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  idx = 0;
  for (iy=NGy; iy<Ny+NGy; iy++) {
    Ybl = iy * (Nx+2*NGx);
    for (ix=NGx; ix<Nx+NGx; ix++)
      tmp[idx++] = xx[Vbl+Ybl+ix];
  }
  HYPRE_StructVectorSetBoxValues(pdata->xvec, ilower, iupper, tmp);

  /*    assemble vectors */
  HYPRE_StructVectorAssemble(pdata->xvec);
  HYPRE_StructVectorAssemble(pdata->bvec);

  /* computing the matvec */
  HYPRE_StructMatrixMatvec(1.0, pdata->De, pdata->xvec, 1.0, pdata->bvec);

  /* extract product vector into bb */
  Vbl = 7 * (Nx+2*NGx) * (Ny+2*NGy);
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_StructVectorGetBoxValues(pdata->bvec, ilower, iupper, tmp);
  idx = 0;
  for (iy=NGy; iy<Ny+NGy; iy++) {
    Ybl = iy * (Nx+2*NGx);
    for (ix=NGx; ix<Nx+NGx; ix++) 
      bb[Vbl+Ybl+ix] = tmp[idx++];
  }
  
  /* return success.  */
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDeNumIters
 * -----------------------------------------------------------------
 * VPrecDeNumIters returns the current cumulative number of 
 * Multigrid Iterations used by the HYPRE preconditioner for 
 * the De solves.
 * -------------------------------------------------------------- */
int VPrecDeNumIters(void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDeData pdata;
  pdata = (ViscPrecDeData) P_data;
  
  /* return with the result from pdata */
  return(pdata->totIters);
}


/********************************************************************/
