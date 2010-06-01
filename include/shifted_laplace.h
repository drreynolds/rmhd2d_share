/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Daniel R. Reynolds
 * UC San Diego, Mathematics
 * -----------------------------------------------------------------
 * This is the header file for shifted_laplace.c.
 * -----------------------------------------------------------------
 */

/*
 * ===========================================================================
 *
 *         Description and Usage of the Shifted Laplace Operstor Interface 
 *
 * ===========================================================================
 */

#ifndef _SHIFTED_LAPLACE_OPERATOR_H
#define _SHIFTED_LAPLACE_OPERATOR_H

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* HYPRE header files  */
#include <stdio.h>
#include <math.h>
#include "HYPRE_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"

/* MPI header file */
#ifdef PARALLEL
#include "mpi.h"
#else
typedef int MPI_Comm;
#endif


/* Definitions of interface function names */

#if defined(F77_FUNC)

#define SHIFTED_LAPLACE_INIT       F77_FUNC(shlaplaceinit,     SHLAPLACEINIT)
#define SHIFTED_LAPLACE_SETUP      F77_FUNC(shlaplacesetup,    SHLAPLACESETUP)
#define SHIFTED_LAPLACE_SOLVE      F77_FUNC(shlaplacesolve,    SHLAPLACESOLVE)
#define SHIFTED_LAPLACE_FREE       F77_FUNC(shlaplacefree,     SHLAPLACEFREE)
#define SHIFTED_LAPLACE_NUMITERS   F77_FUNC(shlaplacenumiters, SHLAPLACENUMITERS)
#define SHIFTED_LAPLACE_OPTS       F77_FUNC(shlaplaceopts,     SHLAPLACEOPTS)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define SHIFTED_LAPLACE_INIT       shlaplaceinit
#define SHIFTED_LAPLACE_SETUP      shlaplacesetup
#define SHIFTED_LAPLACE_SOLVE      shlaplacesolve
#define SHIFTED_LAPLACE_FREE       shlaplacefree
#define SHIFTED_LAPLACE_NUMITERS   shlaplacenumiters
#define SHIFTED_LAPLACE_OPTS       shlaplaceopts

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define SHIFTED_LAPLACE_INIT       SHLAPLACEINIT
#define SHIFTED_LAPLACE_SETUP      SHLAPLACESETUP
#define SHIFTED_LAPLACE_SOLVE      SHLAPLACESOLVE
#define SHIFTED_LAPLACE_FREE       SHLAPLACEFREE
#define SHIFTED_LAPLACE_NUMITERS   SHLAPLACENUMITERS
#define SHIFTED_LAPLACE_OPTS       SHLAPLACEOPTS

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define SHIFTED_LAPLACE_INIT       shlaplaceinit_
#define SHIFTED_LAPLACE_SETUP      shlaplacesetup_
#define SHIFTED_LAPLACE_SOLVE      shlaplacesolve_
#define SHIFTED_LAPLACE_FREE       shlaplacefree_
#define SHIFTED_LAPLACE_NUMITERS   shlaplacenumiters_
#define SHIFTED_LAPLACE_OPTS       shlaplaceopts_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define SHIFTED_LAPLACE_INIT       SHLAPLACEINIT_
#define SHIFTED_LAPLACE_SETUP      SHLAPLACESETUP_
#define SHIFTED_LAPLACE_SOLVE      SHLAPLACESOLVE_
#define SHIFTED_LAPLACE_FREE       SHLAPLACEFREE_
#define SHIFTED_LAPLACE_NUMITERS   SHLAPLACENUMITERS_
#define SHIFTED_LAPLACE_OPTS       SHLAPLACEOPTS_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define SHIFTED_LAPLACE_INIT       shlaplaceinit__
#define SHIFTED_LAPLACE_SETUP      shlaplacesetup__
#define SHIFTED_LAPLACE_SOLVE      shlaplacesolve__
#define SHIFTED_LAPLACE_FREE       shlaplacefree__
#define SHIFTED_LAPLACE_NUMITERS   shlaplacenumiters__
#define SHIFTED_LAPLACE_OPTS       shlaplaceopts__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define SHIFTED_LAPLACE_INIT       SHLAPLACEINIT__
#define SHIFTED_LAPLACE_SETUP      SHLAPLACESETUP__
#define SHIFTED_LAPLACE_SOLVE      SHLAPLACESOLVE__
#define SHIFTED_LAPLACE_FREE       SHLAPLACEFREE__
#define SHIFTED_LAPLACE_NUMITERS   SHLAPLACENUMITERS__
#define SHIFTED_LAPLACE_OPTS       SHLAPLACEOPTS__

#endif



/* Data structure for internal precondtioner data */
typedef struct {

  /* Domain-specific variables */
  int     ndim;    /* dimensional size of grid */
  long int  Nx;    /* mesh cells in x-direction */
  long int  Ny;    /* mesh cells in y-direction */
  long int  Nz;    /* mesh cells in z-direction */
  double    dx;    /* mesh size in x-direction */
  double    dy;    /* mesh size in y-direction */
  double    dz;    /* mesh size in z-direction */
  long int  NGx;   /* ghost cells in x-direction */
  long int  NGy;   /* ghost cells in y-direction */
  long int  NGz;   /* ghost cells in z-direction */
  long int  iXL;   /* lower x-bound for this processor in global array */
  long int  iXR;   /* upper x-bound for this processor in global array */
  long int  iYL;   /* lower y-bound for this processor in global array */
  long int  iYR;   /* upper y-bound for this processor in global array */
  long int  iZL;   /* lower z-bound for this processor in global array */
  long int  iZR;   /* upper z-bound for this processor in global array */
  int       xbc;   /* x-boundary condition (0->0-grad; 1->perdc; 3->dirichlet) */
  int       ybc;   /* y-boundary condition (0->0-grad; 1->perdc; 3->dirichlet) */
  int       zbc;   /* y-boundary condition (0->0-grad; 1->perdc; 3->dirichlet) */

  /* Problem-specific variables */
  double alpha;    /* shift: A = (alpha*I - Laplace) */

  /* MPI-specific variables */
  int Xprocs;     /* number of processors in x-direction */
  int Yprocs;     /* number of processors in y-direction */
  int Zprocs;     /* number of processors in z-direction */
  int iprocx;     /* current process location in x-direction */
  int iprocy;     /* current process location in y-direction */
  int iprocz;     /* current process location in z-direction */
  int outproc;    /* integer representing output processor */
  MPI_Comm comm;  /* MPI communicator */

  /* HYPRE Struct-specific data */
  int  stsize;                     /* stencil size */
  HYPRE_StructGrid    grid;        /* HYPRE grid object for solve */
  HYPRE_StructStencil stencil;     /* stencil object (solver) */

  /* HYPRE Solver-specific data */
  HYPRE_StructMatrix A;      /* shifted Laplace system matrix */
  int sol_zeroguess;         /* use zero initial guess */
  int sol_maxit;             /* max iterations */
  int sol_relch;             /* rel. change stopping criteria */
  int sol_rlxtype;           /* relaxation type */
  int sol_npre;              /* num. pre-relaxation sweeps */
  int sol_npost;             /* num. post-relaxation sweeps */
  int sol_printl;            /* print level */
  int sol_log;               /* amount of logging */

  /* Extra variables for solver diagnostics */
  int AInit;     /* flag denoting initialization of A matrix */
  int totIters;  /* total MG iterations solve */

} ShLaplaceOpData;






/* Prototypes of C Functions */
ShLaplaceOpData *ShLaplaceOpAlloc2D(long int Nx, long int Ny, double dx, 
				    double dy, long int NGx, long int NGy, 
				    int NPx, int NPy, int iPx, int iPy, 
				    int NbXl, int NbXr, int NbYl, int NbYr, 
				    int xbcond, int ybcond);
ShLaplaceOpData *ShLaplaceOpAlloc3D(long int Nx, long int Ny, long int Nz, 
				    double dx, double dy, double dz, 
				    long int NGx, long int NGy, long int NGz, 
				    int NPx, int NPy, int NPz, int iPx, 
				    int iPy, int iPz, int NbXl, int NbXr, 
				    int NbYl, int NbYr, int NbZl, int NbZr, 
				    int xbcond, int ybcond, int zbcond);
void ShLaplaceOpFree(ShLaplaceOpData *hdata);
int ShLaplaceOpSetup2D(double alpha, ShLaplaceOpData *hdata);
int ShLaplaceOpSetup3D(double alpha, ShLaplaceOpData *hdata);
int ShLaplaceOpSolve2D(double *xx, double *bb, double delta, 
		       ShLaplaceOpData *hdata);
int ShLaplaceOpSolve3D(double *xx, double *bb, double delta, 
		       ShLaplaceOpData *hdata);
int ShLaplaceOpIters(ShLaplaceOpData *hdata);



/* Prototypes of Fortran-callable C interface functions */
void SHIFTED_LAPLACE_INIT(int *ndim, long int *Nx, long int *Ny, long int *Nz, 
			  double *dx, double *dy, double *dz, long int *NGx, 
			  long int *NGy, long int *NGz, int *NPx, int *NPy, 
			  int *NPz, int *iPx, int *iPy, int *iPz, int *NbXl, 
			  int *NbXr, int *NbYl, int *NbYr, int *NbZl, 
			  int *NbZr, int *XBcond, int *YBcond, int *ZBcond, 
			  int *hdata, int *ier);
void SHIFTED_LAPLACE_FREE(int *Ahandle);
void SHIFTED_LAPLACE_SETUP(double *alpha, int *Ahandle, int *ier);
void SHIFTED_LAPLACE_SOLVE(double *xx, double *bb, double *delta, 
			   int *Ahandle, int *ier);
void SHIFTED_LAPLACE_NUMITERS(int *Niters, int *Ahandle);
void SHIFTED_LAPLACE_OPTS(int *iopt, int *Ahandle, int *ier);


#endif
/********************************************************************/
