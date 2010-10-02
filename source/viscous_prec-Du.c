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


static void *VPDuData;


/********************************************************************/
/* Fortran callable interface routines                              */
/********************************************************************/


/* Du preconditioner dataspace allocation wrapper routine */
void VISCPREC_DU_INIT(long int *Nx, long int *Ny, double *dx, double *dy,
		      long int *Ns, long int *NGx, long int *NGy, 
		      int *NPx, int *NPy, int *iPx, int *iPy, 
		      int *NBlt, int *NBrt, int *NBtp, int *NBbt, 
		      int *XBcond, int *YBcond, int *ier)
{
  /* allocate preconditioner data */
  VPDuData = VPrecDuAlloc(*Nx, *Ny, *dx, *dy, *Ns, *NGx, *NGy, *NPx, 
			  *NPy, *iPx, *iPy, *NBlt, *NBrt, *NBtp, *NBbt, 
			  *XBcond, *YBcond);
  if (VPDuData == NULL) *ier = -1; 
  else                  *ier = 0;

  return;
}



/* Du preconditioner dataspace deallocation wrapper routine */
void VISCPREC_DU_FREE()
{
  VPrecDuFree(VPDuData);
  return;
}



/* Du preconditioner setup wrapper routine */
void VISCPREC_DU_SETUP(double *uu, double *gamdt, double *Mu, 
		       double *Re, double *v1, double *v2, int *ier)
{
  /* call the C preconditioner setup routine */
  *ier = VPrecDuSetup(uu, *gamdt, *Mu, *Re, v1, v2, VPDuData);

  return;
}



/* momentum preconditioner solve wrapper routine */
void VISCPREC_DU_SOLVE(double *xx, double *bb, double *tmp, double *delta, int *ier)
{
  /* call the C preconditioner solve routine */
  *ier = VPrecDuSolve(xx, bb, tmp, *delta, VPDuData);

  return;
}



/* momentum preconditioner multiply wrapper routine */
void VISCPREC_DU_MULTIPLY(double *xx, double *bb, double *tmp, int *ier)
{
  /* call the C preconditioner multiply routine */
  *ier = VPrecDuMultiply(xx, bb, tmp, VPDuData);

  return;
}



/* Du preconditioner options routine */
void SET_SOL_DU_OPTS(int *iopt, int *ier)
{
  if (VPDuData == NULL) *ier = 1;
  else {
    /* cast VPDuData as the correct structure */
    ViscPrecDuData pdata;
    pdata = (ViscPrecDuData) VPDuData;
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



/* Du preconditioner diagnostic output routine */
void VISCPREC_DU_NUMITERS(int *Niters)
{
  *Niters = VPrecDuNumIters(VPDuData);
  return;
}




/********************************************************************/
/* Internal Preconditioner Routines                                 */
/********************************************************************/


/* -----------------------------------------------------------------
 * Function : VPrecDuAlloc
 * -----------------------------------------------------------------
 * VPrecDuAlloc is called at initialization to set aside space for 
 * any internal storage that will be required by VPrecDuSetup and
 * VPrecDuSolve.
 * -------------------------------------------------------------- */
void *VPrecDuAlloc(long int Nx, long int Ny, double dx, double dy, 
		   long int Ns, long int NGx, long int NGy, int NPx, int NPy, 
		   int iPx, int iPy, int NBlt, int NBrt, int NBtp, int NBbt, 
		   int XBcond, int YBcond) 
{
  /* define necessary local variables, output variable */
  ViscPrecDuData pdata;
  int xtag, ytag;
  int ilower[2], iupper[2];
  long int Nxl, Nyl;

  /* allocate preconditioner data, cast as ViscPrecData */
  pdata = (ViscPrecDuData) malloc(sizeof *pdata);
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
  pdata->ybc  = YBcond;     /* x-boundary condition */

  /* processor layout information */
  pdata->Xprocs = NPx;      /* num procs in x-dir. */
  pdata->Yprocs = NPy;      /* num procs in y-dir. */
  pdata->iprocx = iPx;      /* x-loc in proc. grid */
  pdata->iprocy = iPy;      /* x-loc in proc. grid */
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

  /*     set the matrix type: HYPRE_SSTRUCT or HYPRE_PARCSR */
  pdata->mattype = HYPRE_SSTRUCT;

  /*    initializing solver diagnostic information */
  pdata->totIters = 0;
  pdata->Tsetup = 0.0;
  pdata->Tsolve = 0.0;

  /*    set up the grid */
  /*       create grid object */
  HYPRE_SStructGridCreate(pdata->comm, pdata->ndim, 1, &(pdata->grid));

  /*       set my grid extents as if we have one part with multiple boxes.
	   Have each processor describe it's own global extents */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_SStructGridSetExtents(pdata->grid, 0, ilower, iupper);

  /*       set grid variables for this part */
#ifdef TWO_HALF_D
  int vartypes[3] = {HYPRE_SSTRUCT_VARIABLE_CELL,
		     HYPRE_SSTRUCT_VARIABLE_CELL,
		     HYPRE_SSTRUCT_VARIABLE_CELL};
  HYPRE_SStructGridSetVariables(pdata->grid, 0, 3, vartypes);
#else
  int vartypes[2] = {HYPRE_SSTRUCT_VARIABLE_CELL,
		     HYPRE_SSTRUCT_VARIABLE_CELL};
  HYPRE_SStructGridSetVariables(pdata->grid, 0, 2, vartypes);
#endif

  /*       set grid periodicity */
  /*         [this currently does not work simply in the interface,  */
  /*          so we must manually set neighbor boxes into the grid]  */
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
  HYPRE_SStructGridSetPeriodic(pdata->grid, 0, periodicity);

  /*       assemble the grid */
  HYPRE_SStructGridAssemble(pdata->grid);
  
  /*       create wx stencil */
  HYPRE_SStructStencilCreate(pdata->ndim, 14, &(pdata->wxStencil));
  /*       create wy stencil */
  HYPRE_SStructStencilCreate(pdata->ndim, 14, &(pdata->wyStencil));
#ifdef TWO_HALF_D
  /*       create wz stencil */
  HYPRE_SStructStencilCreate(pdata->ndim, 5, &(pdata->wzStencil));
#endif

  /*       set stencil entries */
  int offset[2];
  /*         wx dependency on wx to bottom */
  offset[0] = 0;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 0, offset, 0);
  /*         wx dependency on wx to left */
  offset[0] = -1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 1, offset, 0);
  /*         wx dependency on wx to self */
  offset[0] = 0;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 2, offset, 0);
  /*         wx dependency on wx to right */
  offset[0] = 1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 3, offset, 0);
  /*         wx dependency on wx to top */
  offset[0] = 0;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 4, offset, 0);

  /*         wx dependency on wy to left bottom */
  offset[0] = -1;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 5, offset, 1);
  /*         wx dependency on wy to bottom */
  offset[0] = 0;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 6, offset, 1);
  /*         wx dependency on wy to right bottom*/
  offset[0] = 1;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 7, offset, 1);
  /*         wx dependency on wy to left */
  offset[0] = -1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 8, offset, 1);
  /*         wx dependency on wy to self */
  offset[0] = 0;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 9, offset, 1);
  /*         wx dependency on wy to right */
  offset[0] = 1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 10, offset, 1);
  /*         wx dependency on wy to left top */
  offset[0] = -1;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 11, offset, 1);
  /*         wx dependency on wy to top */
  offset[0] = 0;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 12, offset, 1);
  /*         wx dependency on wy to right top */
  offset[0] = 1;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 13, offset, 1);


  /*         wy dependency on wy to bottom */
  offset[0] = 0;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 0, offset, 1);
  /*         wy dependency on wy to left */
  offset[0] = -1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 1, offset, 1);
  /*         wy dependency on wy to self */
  offset[0] = 0;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 2, offset, 1);
  /*         wy dependency on wy to right */
  offset[0] = 1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 3, offset, 1);
  /*         wy dependency on wy to top */
  offset[0] = 0;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 4, offset, 1);

  /*         wy dependency on wx to left bottom */
  offset[0] = -1;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 5, offset, 0);
  /*         wy dependency on wx to bottom */
  offset[0] = 0;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 6, offset, 0);
  /*         wy dependency on wx to right bottom */
  offset[0] = 1;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 7, offset, 0);
  /*         wy dependency on wx to left */
  offset[0] = -1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 8, offset, 0);
  /*         wy dependency on wx to self */
  offset[0] = 0;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 9, offset, 0);
  /*         wy dependency on wx to right */
  offset[0] = 1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 10, offset, 0);
  /*         wy dependency on wx to left top */
  offset[0] = -1;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 11, offset, 0);
  /*         wy dependency on wx to top */
  offset[0] = 0;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 12, offset, 0);
  /*         wy dependency on wx to right top */
  offset[0] = 1;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 13, offset, 0);


#ifdef TWO_HALF_D
  /*         wz dependency on wz to bottom */
  offset[0] = 0;  offset[1] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 0, offset, 2);
  /*         wz dependency on wz to left */
  offset[0] = -1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 1, offset, 2);
  /*         wz dependency on wz to self */
  offset[0] = 0;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 2, offset, 2);
  /*         wz dependency on wz to right */
  offset[0] = 1;  offset[1] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 3, offset, 2);
  /*         wz dependency on wz to top */
  offset[0] = 0;  offset[1] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 4, offset, 2);
#endif
  
  /*    set up the graph */
  /*       create the graph object */
  HYPRE_SStructGraphCreate(pdata->comm, pdata->grid, &(pdata->graph));
  
  /*       set graph type according to solver desired */
  HYPRE_SStructGraphSetObjectType(pdata->graph, pdata->mattype);
  
  /*       set stencils into graph */
  /*          set wx stencil */
  HYPRE_SStructGraphSetStencil(pdata->graph, 0, 0, pdata->wxStencil);
  /*          set wy stencil */
  HYPRE_SStructGraphSetStencil(pdata->graph, 0, 1, pdata->wyStencil);
#ifdef TWO_HALF_D
  /*          set wz stencil */
  HYPRE_SStructGraphSetStencil(pdata->graph, 0, 2, pdata->wzStencil);
#endif

  /*       add additional non-stencil entries into graph */
  /*       (none that I can think of) */

  /*       assemble the graph */
  HYPRE_SStructGraphAssemble(pdata->graph);

  /*    set up the matrix */
  HYPRE_SStructMatrixCreate(pdata->comm, pdata->graph, &(pdata->Du));
  HYPRE_SStructMatrixSetObjectType(pdata->Du, pdata->mattype);
  HYPRE_SStructMatrixInitialize(pdata->Du);

  /*    set up the vectors */
  HYPRE_SStructVectorCreate(pdata->comm, pdata->grid, &(pdata->bvec));
  HYPRE_SStructVectorCreate(pdata->comm, pdata->grid, &(pdata->xvec));
  HYPRE_SStructVectorSetObjectType(pdata->bvec, pdata->mattype);
  HYPRE_SStructVectorSetObjectType(pdata->xvec, pdata->mattype);
  HYPRE_SStructVectorInitialize(pdata->bvec);
  HYPRE_SStructVectorInitialize(pdata->xvec);
  
  /* allocate temporary array containing matrix values */
  /* (since the temporary vectors are only 8*Nx*Ny long) */
  pdata->matvals = (double *) malloc(9*Nx*Ny*sizeof(double));


  /********************************/
  /*  continue with general setup */

  /*    init flags to 0 at first */
  pdata->DuInit = 0;

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
 * Function : VPrecDuFree
 * -----------------------------------------------------------------
 * VPrecDuFree frees storage allocated by VPrecDuAlloc
 * -------------------------------------------------------------- */
void VPrecDuFree(void *P_data)
{
  /* ensure that P_data is non-null, and free space as required */
  if ( P_data != NULL ) {

    /* cast P_data as the correct structure */
    ViscPrecDuData pdata;
    pdata = (ViscPrecDuData) P_data;

    /* free the various components */
    HYPRE_SStructMatrixDestroy(pdata->Du);
    HYPRE_SStructVectorDestroy(pdata->xvec);
    HYPRE_SStructVectorDestroy(pdata->bvec);
    HYPRE_SStructGraphDestroy(pdata->graph);
#ifdef TWO_HALF_D
    HYPRE_SStructStencilDestroy(pdata->wzStencil);
#endif
    HYPRE_SStructStencilDestroy(pdata->wyStencil);
    HYPRE_SStructStencilDestroy(pdata->wxStencil);
    HYPRE_SStructGridDestroy(pdata->grid);
    free(pdata->matvals);

    /* finally, free the pdata structure */
    free(pdata);
  }
}



/* -----------------------------------------------------------------
 * Function : VPrecDuSetup
 * -----------------------------------------------------------------
 * VPrecDuSetup sets up the viscous preconditioning matrix for the 
 * momentum equations.
 *
 * The parameters of VPrecSetup used here are as follows:
 *
 * uu      is the current state of the system
 * gamdt   is the time scaling in the Newton matrix, M = I+gamdt*J
 * Mu      is the plasma viscosity coefficient
 * Re      is the plasma Reynolds number
 * v1, v2  give already allocated arrays which may be used as 
 *         temporary storage or work space
 * pdata   is the preconditioner data returned by VPrecAlloc.
 *
 * Return value:
 * The value returned by this VPrecDuSetup function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDuSetup(double *uu, double gamdt, double Mu, 
		 double Re, double *v1, double *v2, void *P_data)
{
  /* local variables */
  long int Nx, Ny, NGx, NGy;
  long int ix, iy, idx, Ybl;
  int ilower[2], iupper[2], selfentries[5], nextentries[9];
  int xLface, xRface, yLface, yRface, IdLoc;
  double dx, dy, dxi2, dyi2, dxdyfac, rhofact, IdVal, Tstart, Tstop;
  double uxvals[14], uyvals[14], uzvals[5];

  /* start timer */
  Tstart = MPI_Wtime();

  /* recast P_data as the correct structure */
  ViscPrecDuData pdata;
  pdata = (ViscPrecDuData) P_data;

  /* get grid information shortcuts */
  Nx  = pdata->Nx;
  Ny  = pdata->Ny;
  NGx = pdata->NGx;
  NGy = pdata->NGy;
  dx  = pdata->dx;
  dy  = pdata->dy;
  dxi2 = 1.0/dx/dx;
  dyi2 = 1.0/dy/dy;
  dxdyfac = 1.0/dx/dy/12.0;

  /* uxvals holds [unscaled] template for ux stencil */
  uxvals[0]  = dyi2;                         /* ux, bottom */
  uxvals[1]  = dxi2*4.0/3.0;                 /* ux, left */
  uxvals[2]  = -2.0*(4.0/3.0*dxi2 + dyi2);   /* ux, self */
  uxvals[3]  = dxi2*4.0/3.0;                 /* ux, right */
  uxvals[4]  = dyi2;                         /* ux, top */
  uxvals[5]  = dxdyfac;                      /* uy, left bottom */
  uxvals[6]  = 0.0;                          /* uy, bottom */
  uxvals[7]  = -dxdyfac;                     /* uy, right bottom */
  uxvals[8] = 0.0;                           /* uy, left */
  uxvals[9] = 0.0;                           /* uy, self */
  uxvals[10] = 0.0;                          /* uy, right */
  uxvals[11] = -dxdyfac;                     /* uy, left top */
  uxvals[12] = 0.0;                          /* uy, top */
  uxvals[13] = dxdyfac;                      /* uy, right top */

  /* uyvals holds [unscaled] template for uy stencil */
  uyvals[0]  = dyi2*4.0/3.0;                 /* uy, bottom */
  uyvals[1]  = dxi2;                         /* uy, left */
  uyvals[2]  = -2.0*(dxi2 + 4.0/3.0*dyi2);   /* uy, self */
  uyvals[3]  = dxi2;                         /* uy, right */
  uyvals[4]  = dyi2*4.0/3.0;                 /* uy, top */
  uyvals[5] = dxdyfac;                       /* ux, left bottom */
  uyvals[6] = 0.0;                           /* ux, bottom */
  uyvals[7] = -dxdyfac;                      /* ux, right bottom */
  uyvals[8] = 0.0;                           /* ux, left */
  uyvals[9] = 0.0;                           /* ux, self */
  uyvals[10] = 0.0;                          /* ux, right */
  uyvals[11] = -dxdyfac;                     /* ux, left top */
  uyvals[12] = 0.0;                          /* ux, top */
  uyvals[13] = dxdyfac;                      /* ux, right top */

#ifdef TWO_HALF_D
  /* uzvals holds [unscaled] template for uz stencil */
  uzvals[0]  = dyi2;                         /* uz, bottom */
  uzvals[1]  = dxi2;                         /* uz, left */
  uzvals[2]  = -2.0*(dxi2 + dyi2);           /* uz, self */
  uzvals[3]  = dxi2;                         /* uz, right */
  uzvals[4]  = dyi2;                         /* uz, top */
#endif

  /* selfentries holds stencil locations for variable couplings to self */
  selfentries[0] = 0;
  selfentries[1] = 1;
  selfentries[2] = 2;
  selfentries[3] = 3;
  selfentries[4] = 4;

  /* nextentries holds stencil locations for variable couplings to next var */
  nextentries[0] = 5;
  nextentries[1] = 6;
  nextentries[2] = 7;
  nextentries[3] = 8;
  nextentries[4] = 9;
  nextentries[5] = 10;
  nextentries[6] = 11;
  nextentries[7] = 12;
  nextentries[8] = 13;

  /* IdLoc, IdVal hold identity location, value */
  IdLoc = 2;
  IdVal = 1.0;

  /* set flags determining whether proc owns external faces */
  xLface = (pdata->iprocx == 1);
  xRface = (pdata->iprocx == pdata->Xprocs);
  yLface = (pdata->iprocy == 1);
  yRface = (pdata->iprocy == pdata->Yprocs);


  /* set ux stencil couplings to self */
  /*       internal cells */
  idx = 0;
  for (iy=0; iy<Ny; iy++) {
    Ybl = (iy+NGy) * (Nx+2*NGx);
    for (ix=0; ix<Nx; ix++) {
      rhofact = -gamdt*Mu*uu[Ybl + ix + NGx];
      pdata->matvals[idx++] = rhofact*uxvals[0];
      pdata->matvals[idx++] = rhofact*uxvals[1];
      pdata->matvals[idx++] = rhofact*uxvals[2];
      pdata->matvals[idx++] = rhofact*uxvals[3];
      pdata->matvals[idx++] = rhofact*uxvals[4];
    }
  }
  /*       ix=0 face adjustment */
  if (xLface) {
    ix=0;
    if (pdata->xbc == 2) {       /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	pdata->matvals[idx+2] -= pdata->matvals[idx+1];  
	pdata->matvals[idx+1] = 0.0;
      }
    }
    else if (pdata->xbc == 0) {  /* zero-gradient */
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	pdata->matvals[idx+2] += pdata->matvals[idx+1];  
	pdata->matvals[idx+1] = 0.0;
      }
    }
  }
  /*       ix=Nx-1 face adjustment */
  if (xRface) {
    ix=Nx-1;
    if (pdata->xbc == 2) {     /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	pdata->matvals[idx+2] -= pdata->matvals[idx+3];  
	pdata->matvals[idx+3] = 0.0;
      }
    }
  }
  else if (pdata->xbc == 0) {  /* zero-gradient */
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	pdata->matvals[idx+2] += pdata->matvals[idx+3];  
	pdata->matvals[idx+3] = 0.0;
      }
  }
  if (pdata->ybc != 1) {       /* not periodic */
    /*       iy=0 face adjustment */
    if (yLface) {
      iy=0;
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	/* reflecting and zero-gradient the same at this face */
	pdata->matvals[idx+2] += pdata->matvals[idx];  
	pdata->matvals[idx] = 0.0;
      }
    }
    /*       iy=Ny-1 face adjustment */
    if (yRface) {
      iy=Ny-1;
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	/* reflecting and zero-gradient the same at this face */
	pdata->matvals[idx+2] += pdata->matvals[idx+4];  
	pdata->matvals[idx+4] = 0.0;
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, iupper, 0, 5, 
				  selfentries, pdata->matvals);


  /* set ux stencil couplings to uy */
  /*       internal cells */
  idx = 0;
  for (iy=0; iy<Ny; iy++) {
    Ybl = (iy+NGy) * (Nx+2*NGx);
    for (ix=0; ix<Nx; ix++) {
      rhofact = -gamdt*Mu*uu[Ybl + ix + NGx];
      pdata->matvals[idx++] = rhofact*uxvals[5];
      pdata->matvals[idx++] = rhofact*uxvals[6];
      pdata->matvals[idx++] = rhofact*uxvals[7];
      pdata->matvals[idx++] = rhofact*uxvals[8];
      pdata->matvals[idx++] = rhofact*uxvals[9];
      pdata->matvals[idx++] = rhofact*uxvals[10];
      pdata->matvals[idx++] = rhofact*uxvals[11];
      pdata->matvals[idx++] = rhofact*uxvals[12];
      pdata->matvals[idx++] = rhofact*uxvals[13];
    }
  }
  if (pdata->xbc != 1) {       /* not periodic */
    /*       ix=0 face adjustment */
    if (xLface) {
      ix=0;
      for (iy=0; iy<Ny; iy++) {
	idx = 9*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+1] += pdata->matvals[idx];    
	pdata->matvals[idx]   = 0.0;
	pdata->matvals[idx+4] += pdata->matvals[idx+3];  
	pdata->matvals[idx+3] = 0.0;
	pdata->matvals[idx+7] += pdata->matvals[idx+6];  
	pdata->matvals[idx+6] = 0.0;
      }
    }
    /*       ix=Nx-1 face adjustment */
    if (xRface) {
      ix=Nx-1;
      for (iy=0; iy<Ny; iy++) {
	idx = 9*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+1] += pdata->matvals[idx+2];  
	pdata->matvals[idx+2] = 0.0;
	pdata->matvals[idx+4] += pdata->matvals[idx+5];  
	pdata->matvals[idx+5] = 0.0;
	pdata->matvals[idx+7] += pdata->matvals[idx+8];  
	pdata->matvals[idx+8] = 0.0;
      }
    }
  }
  /*       iy=0 face adjustment */
  if (yLface) {
    iy=0;
    if (pdata->ybc == 2) {       /* reflecting */
      for (ix=0; ix<Nx; ix++) {
	idx = 9*(iy*Nx + ix);
	pdata->matvals[idx+3] -= pdata->matvals[idx];    
	pdata->matvals[idx]   = 0.0;
	pdata->matvals[idx+4] -= pdata->matvals[idx+1];  
	pdata->matvals[idx+1] = 0.0;
	pdata->matvals[idx+5] -= pdata->matvals[idx+2];  
	pdata->matvals[idx+2] = 0.0;
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (ix=0; ix<Nx; ix++) {
	idx = 9*(iy*Nx + ix);
	pdata->matvals[idx+3] += pdata->matvals[idx];    
	pdata->matvals[idx]   = 0.0;
	pdata->matvals[idx+4] += pdata->matvals[idx+1];  
	pdata->matvals[idx+1] = 0.0;
	pdata->matvals[idx+5] += pdata->matvals[idx+2];  
	pdata->matvals[idx+2] = 0.0;
      }
    }
  }
  /*       iy=Ny-1 face adjustment */
  if (yRface) {
    iy=Ny-1;
    if (pdata->ybc == 2) {       /* reflecting */
      for (ix=0; ix<Nx; ix++) {
	idx = 9*(iy*Nx + ix);
	pdata->matvals[idx+3] -= pdata->matvals[idx+6];  
	pdata->matvals[idx+6] = 0.0;
	pdata->matvals[idx+4] -= pdata->matvals[idx+7];  
	pdata->matvals[idx+7] = 0.0;
	pdata->matvals[idx+5] -= pdata->matvals[idx+8];  
	pdata->matvals[idx+8] = 0.0;
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (ix=0; ix<Nx; ix++) {
	idx = 9*(iy*Nx + ix);
	pdata->matvals[idx+3] += pdata->matvals[idx+6];  
	pdata->matvals[idx+6] = 0.0;
	pdata->matvals[idx+4] += pdata->matvals[idx+7];  
	pdata->matvals[idx+7] = 0.0;
	pdata->matvals[idx+5] += pdata->matvals[idx+8];  
	pdata->matvals[idx+8] = 0.0;
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, iupper, 0, 9, 
				  nextentries, pdata->matvals);


  /* add one to matrix diagonal for identity contribution */
  for (ix=0; ix<Nx*Ny; ix++)  pdata->matvals[ix] = IdVal;
  HYPRE_SStructMatrixAddToBoxValues(pdata->Du, 0, ilower, iupper, 0, 1, 
				    &IdLoc, pdata->matvals);



  /* set uy stencil couplings to uy */
  /*       internal cells */
  idx = 0;
  for (iy=0; iy<Ny; iy++) {
    Ybl = (iy+NGy) * (Nx+2*NGx);
    for (ix=0; ix<Nx; ix++) {
      rhofact = -gamdt*Mu*uu[Ybl + ix + NGx];
      pdata->matvals[idx++] = rhofact*uyvals[0];
      pdata->matvals[idx++] = rhofact*uyvals[1];
      pdata->matvals[idx++] = rhofact*uyvals[2];
      pdata->matvals[idx++] = rhofact*uyvals[3];
      pdata->matvals[idx++] = rhofact*uyvals[4];
    }
  }
  if (pdata->xbc != 1) {      /* not periodic */
    /*       ix=0 face adjustment */
    if (xLface) {
      ix=0;
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+2] += pdata->matvals[idx+1];  
	pdata->matvals[idx+1] = 0.0;
      }
    }
    /*       ix=Nx-1 face adjustment */
    if (xRface) {
      ix=Nx-1;
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+2] += pdata->matvals[idx+3];  
	pdata->matvals[idx+3] = 0.0;
      }
    }
  }
  /*       iy=0 face adjustment */
  if (yLface) {
    iy=0;
    if (pdata->ybc == 2) {       /* reflecting */
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	pdata->matvals[idx+2] -= pdata->matvals[idx];  
	pdata->matvals[idx] = 0.0;
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	pdata->matvals[idx+2] += pdata->matvals[idx];  
	pdata->matvals[idx] = 0.0;
      }
    }
  }
  /*       iy=Ny-1 face adjustment */
  if (yRface) {
    iy=Ny-1;
    if (pdata->ybc == 2) {       /* reflecting */
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	pdata->matvals[idx+2] -= pdata->matvals[idx+4];  
	pdata->matvals[idx+4] = 0.0;
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	pdata->matvals[idx+2] += pdata->matvals[idx+4];  
	pdata->matvals[idx+4] = 0.0;
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, iupper, 1, 5, 
				  selfentries, pdata->matvals);


  /* set uy stencil couplings to ux */
  idx = 0;
  for (iy=0; iy<Ny; iy++) {
    Ybl = (iy+NGy) * (Nx+2*NGx);
    for (ix=0; ix<Nx; ix++) {
      rhofact = -gamdt*Mu*uu[Ybl + ix + NGx];
      pdata->matvals[idx++] = rhofact*uyvals[5];
      pdata->matvals[idx++] = rhofact*uyvals[6];
      pdata->matvals[idx++] = rhofact*uyvals[7];
      pdata->matvals[idx++] = rhofact*uyvals[8];
      pdata->matvals[idx++] = rhofact*uyvals[9];
      pdata->matvals[idx++] = rhofact*uyvals[10];
      pdata->matvals[idx++] = rhofact*uyvals[11];
      pdata->matvals[idx++] = rhofact*uyvals[12];
      pdata->matvals[idx++] = rhofact*uyvals[13];
    }
  }
  /*       ix=0 face adjustment */
  if (xLface) {
    ix=0;
    if (pdata->xbc == 2) {       /* reflecting */ 
      for (iy=0; iy<Ny; iy++) {
	idx = 9*(iy*Nx + ix);
	pdata->matvals[idx+1] -= pdata->matvals[idx];    
	pdata->matvals[idx]   = 0.0;
	pdata->matvals[idx+4] -= pdata->matvals[idx+3];  
	pdata->matvals[idx+3] = 0.0;
	pdata->matvals[idx+7] -= pdata->matvals[idx+6];  
	pdata->matvals[idx+6] = 0.0;
      }
    }
    else if (pdata->xbc == 0) {  /* zero gradient */
      for (iy=0; iy<Ny; iy++) {
	idx = 9*(iy*Nx + ix);
	pdata->matvals[idx+1] += pdata->matvals[idx];    
	pdata->matvals[idx]   = 0.0;
	pdata->matvals[idx+4] += pdata->matvals[idx+3];  
	pdata->matvals[idx+3] = 0.0;
	pdata->matvals[idx+7] += pdata->matvals[idx+6];  
	pdata->matvals[idx+6] = 0.0;
      }
    }
  }
  /*       ix=Nx-1 face adjustment */
  if (xRface) {
    ix=Nx-1;
    if (pdata->xbc == 2) {       /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	idx = 9*(iy*Nx + ix);
	pdata->matvals[idx+1] -= pdata->matvals[idx+2];  
	pdata->matvals[idx+2] = 0.0;
	pdata->matvals[idx+4] -= pdata->matvals[idx+5];  
	pdata->matvals[idx+5] = 0.0;
	pdata->matvals[idx+7] -= pdata->matvals[idx+8];  
	pdata->matvals[idx+8] = 0.0;
      }
    }
    else if (pdata->xbc == 0) {  /* zero gradient */
      for (iy=0; iy<Ny; iy++) {
	idx = 9*(iy*Nx + ix);
	pdata->matvals[idx+1] += pdata->matvals[idx+2];  
	pdata->matvals[idx+2] = 0.0;
	pdata->matvals[idx+4] += pdata->matvals[idx+5];  
	pdata->matvals[idx+5] = 0.0;
	pdata->matvals[idx+7] += pdata->matvals[idx+8];  
	pdata->matvals[idx+8] = 0.0;
      }
    }
  }
  /*       iy=0 face adjustment */
  if (pdata->ybc != 1) {       /* not periodic */
    if (yLface) {
      iy=0;
      for (ix=0; ix<Nx; ix++) {
	idx = 9*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+3] += pdata->matvals[idx];    
	pdata->matvals[idx]   = 0.0;
	pdata->matvals[idx+4] += pdata->matvals[idx+1];  
	pdata->matvals[idx+1] = 0.0;
	pdata->matvals[idx+5] += pdata->matvals[idx+2];  
	pdata->matvals[idx+2] = 0.0;
      }
    }
  }
  /*       iy=Ny-1 face adjustment */
  if (pdata->ybc != 1) {       /* not periodic */
    if (yRface) {
      iy=Ny-1;
      for (ix=0; ix<Nx; ix++) {
	idx = 9*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+3] += pdata->matvals[idx+6];  
	pdata->matvals[idx+6] = 0.0;
	pdata->matvals[idx+4] += pdata->matvals[idx+7];  
	pdata->matvals[idx+7] = 0.0;
	pdata->matvals[idx+5] += pdata->matvals[idx+8];  
	pdata->matvals[idx+8] = 0.0;
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, iupper, 1, 9, 
				  nextentries, pdata->matvals);


  /* add one to matrix diagonal for identity contribution */
  for (ix=0; ix<Nx*Ny; ix++)  pdata->matvals[ix] = IdVal;
  HYPRE_SStructMatrixAddToBoxValues(pdata->Du, 0, ilower, iupper, 1, 1, 
				    &IdLoc, pdata->matvals);


#ifdef TWO_HALF_D
  /* set uz stencil couplings to uz */
  /*       internal cells */
  idx = 0;
  for (iy=0; iy<Ny; iy++) {
    Ybl = (iy+NGy) * (Nx+2*NGx);
    for (ix=0; ix<Nx; ix++) {
      rhofact = -gamdt*Mu*uu[Ybl + ix + NGx];
      pdata->matvals[idx++] = rhofact*uzvals[0];
      pdata->matvals[idx++] = rhofact*uzvals[1];
      pdata->matvals[idx++] = rhofact*uzvals[2];
      pdata->matvals[idx++] = rhofact*uzvals[3];
      pdata->matvals[idx++] = rhofact*uzvals[4];
    }
  }
  /*       ix=0 face adjustment */
  if (pdata->xbc != 1) {      /* not periodic */
    if (xLface) {
      ix=0;
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+2] += pdata->matvals[idx+1];  
	pdata->matvals[idx+1] = 0.0;
      }
    }
  }
  /*       ix=Nx-1 face adjustment */
  if (pdata->xbc != 1) {      /* not periodic */
    if (xRface) {
      ix=Nx-1;
      for (iy=0; iy<Ny; iy++) {
	idx = 5*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+2] += pdata->matvals[idx+3];  
	pdata->matvals[idx+3] = 0.0;
      }
    }
  }
  /*       iy=0 face adjustment */
  if (pdata->ybc != 1) {      /* not periodic */
    if (yLface) {
      iy=0;
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+2] += pdata->matvals[idx];  
	pdata->matvals[idx] = 0.0;
      }
    }
  }
  /*       iy=Ny-1 face adjustment */
  if (pdata->ybc != 1) {      /* not periodic */
    if (yRface) {
      iy=Ny-1;
      for (ix=0; ix<Nx; ix++) {
	idx = 5*(iy*Nx + ix);
	/* reflecting and zero gradient are the same at this face */
	pdata->matvals[idx+2] += pdata->matvals[idx+4];  
	pdata->matvals[idx+4] = 0.0;
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, iupper, 2, 5, 
				  selfentries, pdata->matvals);


  /* add one to matrix diagonal for identity contribution */
  for (ix=0; ix<Nx*Ny; ix++)  pdata->matvals[ix] = IdVal;
  HYPRE_SStructMatrixAddToBoxValues(pdata->Du, 0, ilower, iupper, 2, 1, 
				    &IdLoc, pdata->matvals);
#endif


  /* assemble matrix */
  HYPRE_SStructMatrixAssemble(pdata->Du);
    
  /*       set init flag */
  pdata->DuInit = 1;

/*   /\* output matrix to file *\/ */
/*   if ((pdata->outproc)==1) */
/*     {printf("      printing HYPRE Du matrix to file \n");} */
/*   HYPRE_SStructMatrixPrint("Du.mat", pdata->Du, 0); */

  /* stop timer, add to totals, and output to screen */ 
  Tstop = MPI_Wtime();
  pdata->Tsetup += Tstop - Tstart;
  if ((pdata->outproc) == 1) {
    printf("Du cumulative setup time: %g\n",pdata->Tsetup);
    printf("Du cumulative solve time: %g\n",pdata->Tsolve);
  }

  /* return success */ 
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDuSolve
 * -----------------------------------------------------------------
 * VPrecDuSolve solves a linear system P x = b, with the
 * preconditioner matrix P generated by VPrecSetup and solved
 * using the HYPRE library.
 *
 * The parameters of VPrecDuSolve used here are as follows:
 *
 * xx      is the sol vector on output
 * bb      is the rhs vector on input
 * tmp     is a temporary vector the size of xx
 * delta   is the desired linear solve tolerance (if used in 
 *         some iterative method)
 * pdata   is the pre-computed Du preconditioner data
 *
 * The value returned by this VPrecDuSolve function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDuSolve(double *xx, double *bb, double *tmp, double delta, void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDuData pdata;
  pdata = (ViscPrecDuData) P_data;

  /* local variables */
  int Sits, Pits, ilower[2], iupper[2];
  long int iy, ix, iv, idx, Nx, NGx, Ny, NGy;
  long int Vbl, Ybl;
  double finalresid, resid, val, bnorm, Tstart, Tstop;
  HYPRE_SStructSolver solver;
  HYPRE_SStructSolver preconditioner;
  int printl = ((pdata->outproc)==1) ? pdata->sol_printl : 0;

  /*   if (printl) printf("     solving Du preconditioner system\n"); */

  /* start timer */
  Tstart = MPI_Wtime();

  /* check that Du matrix initialized */
  if (pdata->DuInit == 0) {
    printf("VPrecDuSolve error: Du matrix uninitialized!\n");
    return(1);
  }

  /* set local variables */
  Nx  = pdata->Nx;
  NGx = pdata->NGx;
  Ny  = pdata->Ny;
  NGy = pdata->NGy;

  /* convert rhs, solution vectors to HYPRE format:         */
  /*    insert rhs vector entries into HYPRE vectors bvec   */
  /*    and xvec (use the rhs as the sol initial guess)     */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  for (iv=1; iv<3; iv++) {
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy);
    idx = 0;
    for (iy=NGy; iy<Ny+NGy; iy++) {
      Ybl = iy * (Nx+2*NGx);
      for (ix=NGx; ix<Nx+NGx; ix++)
	tmp[idx++] = bb[Vbl+Ybl+ix];
    }
    HYPRE_SStructVectorSetBoxValues(pdata->bvec, 0, ilower, iupper, iv-1, tmp);
    for (idx=0; idx<Nx*Ny; idx++)  tmp[idx] = 0.0;
    HYPRE_SStructVectorSetBoxValues(pdata->xvec, 0, ilower, iupper, iv-1, tmp);
  }
#ifdef TWO_HALF_D
  iv = 3;
  Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy);
  idx = 0;
  for (iy=NGy; iy<Ny+NGy; iy++) {
    Ybl = iy * (Nx+2*NGx);
    for (ix=NGx; ix<Nx+NGx; ix++)
      tmp[idx++] = bb[Vbl+Ybl+ix];
  }
  HYPRE_SStructVectorSetBoxValues(pdata->bvec, 0, ilower, iupper, iv-1, tmp);
  for (idx=0; idx<Nx*Ny; idx++)  tmp[idx] = 0.0;
  HYPRE_SStructVectorSetBoxValues(pdata->xvec, 0, ilower, iupper, iv-1, tmp);
#endif

  /*    assemble vectors */
  HYPRE_SStructVectorAssemble(pdata->xvec);
  HYPRE_SStructVectorAssemble(pdata->bvec);

/*   if (printl) printf("Writing out rhs to file u_rhs.vec\n"); */
/*   HYPRE_SStructVectorPrint("u_rhs.vec", pdata->bvec, 0); */

  /* set up the solver [PCG/SysPFMG] */
  /*    create the solver & preconditioner */
  HYPRE_SStructPCGCreate(pdata->comm, &solver);
  HYPRE_SStructSysPFMGCreate(pdata->comm, &preconditioner);
 
  /*    set preconditioner & solver options */
  HYPRE_SStructSysPFMGSetMaxIter(preconditioner, pdata->sol_maxit/4);
  HYPRE_SStructSysPFMGSetRelChange(preconditioner, pdata->sol_relch);
  HYPRE_SStructSysPFMGSetRelaxType(preconditioner, pdata->sol_rlxtype);
  HYPRE_SStructSysPFMGSetNumPreRelax(preconditioner, pdata->sol_npre);
  HYPRE_SStructSysPFMGSetNumPostRelax(preconditioner, pdata->sol_npost);
  HYPRE_SStructPCGSetPrintLevel(solver, pdata->sol_printl);
  HYPRE_SStructPCGSetLogging(solver, pdata->sol_log);
  HYPRE_SStructPCGSetRelChange(solver, pdata->sol_relch);
  HYPRE_SStructPCGSetMaxIter(solver, pdata->sol_maxit);
  if (delta != 0.0)  HYPRE_SStructPCGSetTol(solver, delta);
  HYPRE_SStructPCGSetPrecond(solver, 
			     (HYPRE_PtrToSStructSolverFcn) HYPRE_SStructSysPFMGSolve,  
			     (HYPRE_PtrToSStructSolverFcn) HYPRE_SStructSysPFMGSetup, 
			     preconditioner);
  HYPRE_SStructPCGSetup(solver, pdata->Du, pdata->bvec, pdata->xvec);

  /* solve the linear system */
  HYPRE_SStructPCGSolve(solver, pdata->Du, pdata->bvec, pdata->xvec);

  /* extract solver statistics */
  finalresid = 1.0;  Sits = 0;  Pits = 0;
  HYPRE_SStructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_SStructPCGGetNumIterations(solver, &Sits);
  HYPRE_SStructSysPFMGGetNumIterations(preconditioner, &Pits);
  pdata->totIters += Sits;
  if (printl) printf("      Du lin resid = %.1e (tol = %.1e) its = (%i,%i) \n",
		     finalresid, delta, Sits, Pits);

/*   if (printl) printf("Writing out solution to file u_sol.vec\n"); */
/*   HYPRE_SStructVectorPrint("u_sol.vec", pdata->xvec, 0); */

  /* gather the solution vector and extract values */
  HYPRE_SStructVectorGather(pdata->xvec);
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  for (iv=1; iv<3; iv++) {
    HYPRE_SStructVectorGetBoxValues(pdata->xvec, 0, ilower, iupper, iv-1, tmp);
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy);
    idx = 0;
    for (iy=NGy; iy<Ny+NGy; iy++) {
      Ybl = iy * (Nx+2*NGx);
      for (ix=NGx; ix<Nx+NGx; ix++) 
	xx[Vbl+Ybl+ix] = tmp[idx++];
    }
  }
#ifdef TWO_HALF_D
  iv = 3;
  HYPRE_SStructVectorGetBoxValues(pdata->xvec, 0, ilower, iupper, iv-1, tmp);
  Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy);
  idx = 0;
  for (iy=NGy; iy<Ny+NGy; iy++) {
    Ybl = iy * (Nx+2*NGx);
    for (ix=NGx; ix<Nx+NGx; ix++) 
      xx[Vbl+Ybl+ix] = tmp[idx++];
  }
#endif

  /* destroy solver & preconditioner structures */
  HYPRE_SStructPCGDestroy(solver);
  HYPRE_SStructSysPFMGDestroy(preconditioner);
      
  /* stop timer, add to totals, and output to screen */ 
  Tstop = MPI_Wtime();
  pdata->Tsolve += Tstop - Tstart;

  /* return success */
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDuMultiply
 * -----------------------------------------------------------------
 * VPrecDuMultiply performs the matrix-vector product P x = b, with 
 * the preconditioner matrix P generated by VPrecSetup and 
 * multiplied using the HYPRE library.
 *
 * The parameters of VPrecDuMultiply used here are as follows:
 *
 * xx      is the x vector on input
 * bb      is the b vector on output
 * tmp     is a temporary vector the size of xx
 * pdata   is the pre-computed Du preconditioner data
 *
 * The value returned by this VPrecDuMultiply function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDuMultiply(double *xx, double *bb, double *tmp, void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDuData pdata;
  pdata = (ViscPrecDuData) P_data;

  /* local variables */
  int ilower[2], iupper[2];
  long int ix, iy, iv, idx, Nx, NGx, Ny, NGy;
  long int Vbl, Ybl;
  double val;

  /* check that Du matrix initialized */
  if (pdata->DuInit == 0) {
    printf("VPrecDuMultiply error: Du matrix uninitialized!\n");
    return(1);
  }

  /* set local variables */
  Nx  = pdata->Nx;
  NGx = pdata->NGx;
  Ny  = pdata->Ny;
  NGy = pdata->NGy;

  /* convert product, result vectors to HYPRE format        */
  /*    insert product vector entries into HYPRE vector x   */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  for (iv=1; iv<3; iv++) {
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy);
    idx = 0;
    for (iy=NGy; iy<Ny+NGy; iy++) {
      Ybl = iy * (Nx+2*NGx);
      for (ix=NGx; ix<Nx+NGx; ix++)
	tmp[idx++] = xx[Vbl+Ybl+ix];
    }
    HYPRE_SStructVectorSetBoxValues(pdata->xvec, 0, ilower, iupper, iv-1, tmp);
  }
#ifdef TWO_HALF_D
  iv = 3;
  Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy);
  idx = 0;
  for (iy=NGy; iy<Ny+NGy; iy++) {
    Ybl = iy * (Nx+2*NGx);
    for (ix=NGx; ix<Nx+NGx; ix++)
      tmp[idx++] = xx[Vbl+Ybl+ix];
  }
  HYPRE_SStructVectorSetBoxValues(pdata->xvec, 0, ilower, iupper, iv-1, tmp);
#endif

  /*    assemble vectors */
  HYPRE_SStructVectorAssemble(pdata->xvec);
  HYPRE_SStructVectorAssemble(pdata->bvec);

  /* computing the matvec */
  HYPRE_SStructMatrixMatvec(1.0, pdata->Du, pdata->xvec, 1.0, pdata->bvec);

  /* gather the solution vector before extracting values */
  HYPRE_SStructVectorGather(pdata->bvec);

  /* extract solution vector into bb */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;
  for (iv=1; iv<3; iv++) {
    HYPRE_SStructVectorGetBoxValues(pdata->bvec, 0, ilower, iupper, iv-1, tmp);
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy);
    idx = 0;
    for (iy=NGy; iy<Ny+NGy; iy++) {
      Ybl = iy * (Nx+2*NGx);
      for (ix=NGx; ix<Nx+NGx; ix++) 
	bb[Vbl+Ybl+ix] = tmp[idx++];
    }
  }
#ifdef TWO_HALF_D
  iv = 3;
  HYPRE_SStructVectorGetBoxValues(pdata->bvec, 0, ilower, iupper, iv-1, tmp);
  Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy);
  idx = 0;
  for (iy=NGy; iy<Ny+NGy; iy++) {
    Ybl = iy * (Nx+2*NGx);
    for (ix=NGx; ix<Nx+NGx; ix++) 
      bb[Vbl+Ybl+ix] = tmp[idx++];
  }
#endif

  /* return success.  */
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDuNumIters
 * -----------------------------------------------------------------
 * VPrecDuNumIters returns the current cumulative number of 
 * Multigrid Iterations used by the HYPRE preconditioner for 
 * the Du solves.
 * -------------------------------------------------------------- */
int VPrecDuNumIters(void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDuData pdata;
  pdata = (ViscPrecDuData) P_data;
  
  /* return with the result from pdata */
  return(pdata->totIters);
}


/********************************************************************/
