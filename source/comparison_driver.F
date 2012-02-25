c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: comparison_driver.F,v $
c=================================================================




      subroutine ComparisonDriver
c-----------------------------------------------------------------
c     Description: main MHD driver, handles initialization of 
c        variables, timestepping and data I/O.  This MHD driver
c        routine solves an implicit formulation of 2D, 2.5D, 
c        3D ideal, nonlinearly-resistive, MHD equations.  All 
c        initial setup and physics calculations in this 
c        formulation rely on a previous explicit formulation of 
c        these equations by Ravi Samtaney, Princeton Plasma 
c        Physics Laboratory.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use iounits
      use tvdcoeffs
      use sundials_fcmix
      
c======= Declarations =========
      implicit none

c     local driver variables
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: ux0(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: ux_exp(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI)
      double precision :: t0, ttot, told, dt, tdesired
      integer          :: istart, maxiter, lastiter
      integer          :: output_flag, binary_flag
      integer          :: ndump, ndiag, new
      integer          :: tvdFlag
      integer          :: timestep

c     cvode timestepper/solver information
      integer          :: useMyJv, useMyPrec
      integer          :: intmethod, itermethod, taskchoice
      integer          :: prectype, gstype, maxl, mxbdf, mxsteps
      integer          :: mxwarn, stlimit, mxerrfails, mxnewtit
      integer          :: mxconvfails, cvierr, cvierrtot
      double precision :: lintol, dtinit, maxdt, mindt, nonlincoeff
      double precision :: rtol, atol1, atol2, atol3, atol4, atol5
      double precision :: atol6, atol7, atol8, tstop, dtfac
      double precision :: abstol(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: ropt(40)
      integer(kind=sund_long_int) :: nglobal, iopt(40)
      integer(kind=sund_long_int) :: xlocal, ylocal, zlocal, nspecies
      integer(kind=sund_long_int) :: xghost, yghost, zghost

c     comparison variables
      double precision :: d2n, dIn
      double precision :: rho2n, rhoIn
      double precision :: ux2n, uxIn
      double precision :: uy2n, uyIn
      double precision :: uz2n, uzIn
      double precision :: bx2n, bxIn
      double precision :: by2n, byIn
      double precision :: bz2n, bzIn
      double precision :: e2n, eIn
      integer          :: maxproc, maxnormloc(4)
      integer          :: maxlocrecv(5)

c     declare internal namelists for problem specification
      namelist /inparams/ maxiter, new, output_flag, binary_flag
      namelist /dumps/ ndump, ndiag
      namelist /gridbounds/ xl, xr, yl, yr, zl, zr
      namelist /cvode_controls/ useMyJv, useMyPrec,
     &     intmethod, itermethod, taskchoice, prectype, gstype, maxl, 
     &     mxbdf, mxsteps, mxwarn, stlimit, mxerrfails, mxnewtit,  
     &     mxconvfails, lintol, dtinit, maxdt, mindt, nonlincoeff, rtol, 
     &     atol1, atol2, atol3, atol4, atol5, atol6, atol7, atol8, 
     &     tstop, dtfac


c======= Internals ============


c------- Implicit Scheme / Nonlinear Solver Parameters -------------
c        (all parameters set in file cvcontrol.inp)
c
c     specifying default vs custom routines (0=>default, 1=>custom)
c         useMyJv   = 1  =>  use my fcvjtimes routine
c         useMyPrec = 1  =>  use my fcvpset, fcvpsol routines
c
c     required CVODE solver information
c         intmethod: integration method 
c                       1 => Adams (nonstiff), 2 => BDF (stiff)
c         itermethod: nonlinear iteration method
c                       1 => functional iteration, 2 => Newton
c         taskchoice: CVODE task indicator
c                       1 => normal mode (overshoot & interpolate)
c                       2 => return after each internal step taken
c                       3 => like 1, but must stop at tstop
c                       4 => like 2, but must stop at tstop
c                      NOTE: modes 3 and 4 imply that tstop/=0
c         gstype: Gram-Schmidt orthogonalization method
c                       1 => modified GS, 2 => classical GS
c         rtol: relative tolerance for nonlinear time step, should 
c               be on the order of spatial accuracy
c         atol1->atol8: absolute tolerances for each variable, 
c               on the order of the noise in each variable
c         tstop: physical stopping time
c         dtfac: desired output time factor (dtout = dtcfl*dtfac)
c
c     optional CVODE solver information [0 => default]
c         prectype: preconditioning type (default=none)
c                       1 => left, 2 => right, 3 => both
c         maxl: maximum Krylov subspace dimension (default=5)
c         mxbdf: maximum order of BDF method (1-5, default=5)
c         mxsteps: max time steps to reach solution (default=500)
c         mxwarn: max warnings for tn + h = tn (default=10)
c         stlimit: use stability limit detection (0=>off, 1=>on)
c         mxerrfails: max error test fails per step (default=7)
c         mxnewtit: max nonlin solver iters per step (default=3)
c         mxconvfails: max convergence fails per step (default=10)
c         lintol: ratio between lin and nonlin tols (default=0.05)
c         dtinit: initial timestep size (default is estimated)
c         maxdt: maximum absolute step size (default=infinity)
c         mindt: minimum absolute step size (default=0.0)
c         nonlincoeff: coeff in nonlin conv test (default=0.1)
c-------------------------------------------------------------------
      open(411,file='cvcontrol.inp',form='formatted')
      read(411,cvode_controls)
      close(411)

c     set up cvode options according to namelist
      iopt = 0; ropt = 0.0d0
      iopt(1)  = mxbdf
      iopt(2)  = mxsteps
      iopt(3)  = mxwarn
      iopt(14) = stlimit
      iopt(22) = mxerrfails
      iopt(23) = mxnewtit
      iopt(24) = mxconvfails
      ropt(1)  = dtinit
      ropt(2)  = maxdt
      ropt(3)  = mindt
      ropt(8)  = tstop
      ropt(9)  = nonlincoeff


c     input general MHD simulation information using namelists
      open(16,file='mhd.inp',form='formatted')
      read(16,inparams)
      read(16,dumps)
      read(16,gridbounds) 
      close(16)

c     root node outputs problem parameters
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(A,A,i6,3(A,i2))') '\nInput parameters:',
     &        '\n        maxiter =',maxiter,
     &        '\n            new =',new,
     &        '\n    output_flag =',output_flag,
     &        '\n    binary_flag =',binary_flag
         write(6,'(A,2(A,i4))') '\nDump intervals:',
     &        '\n    ndump =',ndump,
     &        '\n    ndiag =',ndiag
         write(6,'(A,6(A,es12.5))') '\nGrid bounds:',
     &        '\n    xl =',xl,',  xr =',xr,
     &        '\n    yl =',yl,',  yr =',yr,
     &        '\n    zl =',zl,',  zr =',zr
         write(6,'(A,2(A,i4),3(A,i2))') '\nGrid indices:',
     &        '\n    nx =',nx,',  ny =',ny,',  nz =',nz,
     &        ',  nvar =',nvar,',  nghost =',nghost
         write(6,'(A,A,es9.2)') '\nImplicit Scheme:',
     &        '\n    Adaptive 5th-order BDF'
         write(6,'(A,3(A,i6),11(A,es9.2))') '\nSolver parameters:',
     &        '\n           maxl =',maxl,
     &        '\n          mxbdf =',mxbdf,
     &        '\n        mxsteps =',mxsteps,
     &        '\n         lintol =',lintol,
     &        '\n    nonlincoeff =',nonlincoeff,
     &        '\n          atol1 =',atol1,
     &        '\n          atol2 =',atol2,
     &        '\n          atol3 =',atol3,
     &        '\n          atol4 =',atol4,
     &        '\n          atol5 =',atol5,
     &        '\n          atol6 =',atol6,
     &        '\n          atol7 =',atol7,
     &        '\n          atol8 =',atol8
      endif

c     initialize timestepping and data output information
      lastiter = 0
      istart   = 1

c     set up domain, grid
      call SetupDomain
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\nAfter domain setup...'
         write(6,'(A,3(A,es12.5))') '\nMesh spacing:',
     &        '\n   dx =',dx,',  dy =',dy,',  dz =',dz
      end if

c     determine local domain location in general domain
      call SetupLocalDomain
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\nAfter local domain setup...'
         write(6,'(A,3(i3,1x),A,3(2(es12.5,1x),A))') 
     &        '\nLocal subdomain for iproc (',
     &        iprocx,iprocy,iprocz,') is the box [',XLloc,XRloc,
     &        '] x [',YLloc,YRloc,'] x [',ZLloc,ZRloc,']'
      endif

c     set up TVD coefficients
c      tvdFlag=1
c      if(tvdFlag.eq.1) call SetTVDCoeffs

c     initialize state variables
      call Initialize(new,ux0,phi,dt,t0,istart,lastiter)
      ttot = t0
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\nAfter local dataspace initialization...'
      end if

c     communicate initial boundary data among processors
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) 'Communicating initial condition to neighbors...'
      end if
      call CommState(ux0)
      ux = ux0


c     set up CVODE NVector interface information using mesh variables
      nglobal = nprocs*nxlsize*nylsize*nzlsize*nvar
      xghost = nghost;  yghost = 0;  zghost = 0
#ifndef ONE_D
      yghost = nghost
#ifndef TWO_D 
      zghost = nghost
#endif
#endif

c     ensure that we send the correct integer kind through the interface
      xlocal = nxlsize
      ylocal = nylsize
      zlocal = nzlsize
      nspecies = nvar


c     initialize CVODE data & operators

c        Nvector specs
      call fnvinitmhd(1, xlocal, ylocal, zlocal, xghost, xghost, 
     &     yghost, yghost, zghost, zghost, nspecies, nglobal, cvierr)
      if (cvierr /= 0) then
         print '(A,A,i2)',' comparison_driver: fnvinitmhd',
     &        ' failed, ier =',cvierr 
         return
      end if

c        cvode dataspace
c        [we must access abstol values before sending
c         to fcvmalloc.]
      abstol(:,:,:,1) = atol1;  abstol(:,:,:,2) = atol2
      abstol(:,:,:,3) = atol3;  abstol(:,:,:,4) = atol4
      abstol(:,:,:,5) = atol5;  abstol(:,:,:,6) = atol6
      abstol(:,:,:,7) = atol7;  abstol(:,:,:,8) = atol8
      if (maxval(abstol) > 0.d0) 
     &     write(6,*) 'Setting abstol vector'
      call fcvmalloc(t0, ux, intmethod, itermethod, 2,
     &     rtol, abstol, 1, iopt, ropt, cvierr)
      if (cvierr /= 0) then
         print '(A,A,i2)',' comparison_driver: fcvmalloc',
     &        ' failed, ier =',cvierr 
         return
      end if

c        spgmr dataspace
      call fcvspgmr(prectype, gstype, maxl, lintol, cvierr)
      if (cvierr /= 0) then
         print '(A,A,i2)',' comparison_driver: fcvspgmr',
     &        ' failed, ier =',cvierr 
         return
      end if

c        J*v routine specification
      call fcvspgmrsetjac(useMyJv, cvierr)
      if (cvierr /= 0) then
         print '(A,A,i2)',' comparison_driver: fcvspgmrsetjac',
     &        ' failed, ier =',cvierr 
         return
      end if

c        preconditioner setup/solver specification
      call fcvspgmrsetprec(useMyPrec, cvierr)
      if (cvierr /= 0) then
         print '(A,A,i2)',' comparison_driver: fcvspgmrsetprec',
     &        ' failed, ier =',cvierr 
         return
      end if



c     Start explicit method time stepping
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\nAfter initialization, starting explicit:'
      end if
      do timeStep = 1,2000

c        solve for new time step
         call MHDSolve(ux,phi,dt)

c        update time
         ttot = ttot + dt

c        root outputs some timing information
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            write(6,*) 'timeStep',timeStep, dt, ttot
         endif

      enddo
      tstop = ttot
      
c     store final explicit solution in ux_exp
      ux_exp = ux



c     re-set CVODE to use rtol=1e-5, and re-initialize variables
      rtol = 1.d-5
      ux   = ux0
      ttot = t0
      call fcvreinit(ttot, ux, 2, rtol, abstol, 1, iopt, ropt, cvierr)

c     evolve to tstop using CVODE
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\n Starting implicit, rtol=1e-5:'
      end if
      call fcvode(tstop, ttot, ux, taskchoice, cvierr)

c     output difference between implicit and explicit solutions
      ux = ux - ux_exp
c        ||diff||
      d2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,:)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,:))
      dIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,:)))
      maxnormloc = maxloc(abs(ux))-(/nghost,nghost,0,0/)

c        ||rho||
      rho2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,1)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,1))
      rhoIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,1)))

c        ||ux||
      ux2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,2)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,2))
      uxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,2)))

c        ||uy||
      uy2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,3)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,3))
      uyIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,3)))

c        ||uz||
      uz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,4)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,4))
      uzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,4)))

c        ||bx||
      bx2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,5)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,5))
      bxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,5)))

c        ||by||
      by2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,6)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,6))
      byIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,6)))

c        ||bz||
      bz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,7)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,7))
      bzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,7)))

c        ||e||
      e2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,8)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,8))
      eIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,8)))

c        scale wrms norms appropriately
      d2n   = dsqrt(d2n/nxlocal/nylocal/nzlocal/nvar)
      rho2n = dsqrt(rho2n/nxlocal/nylocal/nzlocal)
      ux2n  = dsqrt(ux2n/nxlocal/nylocal/nzlocal)
      uy2n  = dsqrt(uy2n/nxlocal/nylocal/nzlocal)
      uz2n  = dsqrt(uz2n/nxlocal/nylocal/nzlocal)
      bx2n  = dsqrt(bx2n/nxlocal/nylocal/nzlocal)
      by2n  = dsqrt(by2n/nxlocal/nylocal/nzlocal)
      bz2n  = dsqrt(bz2n/nxlocal/nylocal/nzlocal)
      e2n   = dsqrt(e2n/nxlocal/nylocal/nzlocal)
      maxlocrecv = (/ 0, maxnormloc(1), maxnormloc(2), 
     &     maxnormloc(3), maxnormloc(4) /)
      
c        root node outputs nonlinear residual information
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(2(A,es11.4))')
     &        '  ||diff(u)||  =',d2n,  ',   max|diff(u)|  =',dIn
         write(6,'(A,i4,A,4(i4),A)')'   max residual on proc',
     &        maxlocrecv(1),',\n     at (x,y,z,v) = (',
     &        maxlocrecv(2:5),')'
         write(6,'(A)') 
     &        '  -------------------------------------------------'
         write(6,'(2(A,es11.4))')
     &        '  ||diff_rho|| =',rho2n,',   max|diff_rho| =',rhoIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_ux||  =',ux2n, ',   max|diff_ux|  =',uxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uy||  =',uy2n, ',   max|diff_uy|  =',uyIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uz||  =',uz2n, ',   max|diff_uz|  =',uzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bx||  =',bx2n, ',   max|diff_Bx|  =',bxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_By||  =',by2n, ',   max|diff_By|  =',byIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bz||  =',bz2n, ',   max|diff_Bz|  =',bzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_e||   =',e2n,  ',   max|diff_e|   =',eIn
      end if
      


c     re-set CVODE to use rtol=1e-6, and re-initialize variables
      rtol = 1.d-6
      ux   = ux0
      ttot = t0
      call fcvreinit(ttot, ux, 2, rtol, abstol, 1, iopt, ropt, cvierr)

c     evolve to tstop using CVODE
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\n Starting implicit, rtol=1e-6:'
      end if
      call fcvode(tstop, ttot, ux, taskchoice, cvierr)

c     output difference between implicit and explicit solutions
      ux = ux - ux_exp
c        ||diff||
      d2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,:)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,:))
      dIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,:)))
      maxnormloc = maxloc(abs(ux))-(/nghost,nghost,0,0/)

c        ||rho||
      rho2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,1)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,1))
      rhoIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,1)))

c        ||ux||
      ux2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,2)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,2))
      uxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,2)))

c        ||uy||
      uy2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,3)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,3))
      uyIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,3)))

c        ||uz||
      uz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,4)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,4))
      uzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,4)))

c        ||bx||
      bx2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,5)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,5))
      bxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,5)))

c        ||by||
      by2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,6)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,6))
      byIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,6)))

c        ||bz||
      bz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,7)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,7))
      bzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,7)))

c        ||e||
      e2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,8)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,8))
      eIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,8)))

c        scale wrms norms appropriately
      d2n   = dsqrt(d2n/nxlocal/nylocal/nzlocal/nvar)
      rho2n = dsqrt(rho2n/nxlocal/nylocal/nzlocal)
      ux2n  = dsqrt(ux2n/nxlocal/nylocal/nzlocal)
      uy2n  = dsqrt(uy2n/nxlocal/nylocal/nzlocal)
      uz2n  = dsqrt(uz2n/nxlocal/nylocal/nzlocal)
      bx2n  = dsqrt(bx2n/nxlocal/nylocal/nzlocal)
      by2n  = dsqrt(by2n/nxlocal/nylocal/nzlocal)
      bz2n  = dsqrt(bz2n/nxlocal/nylocal/nzlocal)
      e2n   = dsqrt(e2n/nxlocal/nylocal/nzlocal)
      maxlocrecv = (/ 0, maxnormloc(1), maxnormloc(2), 
     &     maxnormloc(3), maxnormloc(4) /)
      
c        root node outputs nonlinear residual information
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(2(A,es11.4))')
     &        '  ||diff(u)||  =',d2n,  ',   max|diff(u)|  =',dIn
         write(6,'(A,i4,A,4(i4),A)')'   max residual on proc',
     &        maxlocrecv(1),',\n     at (x,y,z,v) = (',
     &        maxlocrecv(2:5),')'
         write(6,'(A)') 
     &        '  -------------------------------------------------'
         write(6,'(2(A,es11.4))')
     &        '  ||diff_rho|| =',rho2n,',   max|diff_rho| =',rhoIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_ux||  =',ux2n, ',   max|diff_ux|  =',uxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uy||  =',uy2n, ',   max|diff_uy|  =',uyIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uz||  =',uz2n, ',   max|diff_uz|  =',uzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bx||  =',bx2n, ',   max|diff_Bx|  =',bxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_By||  =',by2n, ',   max|diff_By|  =',byIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bz||  =',bz2n, ',   max|diff_Bz|  =',bzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_e||   =',e2n,  ',   max|diff_e|   =',eIn
      end if
      


c     re-set CVODE to use rtol=1e-7, and re-initialize variables
      rtol = 1.d-7
      ux   = ux0
      ttot = t0
      call fcvreinit(ttot, ux, 2, rtol, abstol, 1, iopt, ropt, cvierr)

c     evolve to tstop using CVODE
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\n Starting implicit, rtol=1e-7:'
      end if
      call fcvode(tstop, ttot, ux, taskchoice, cvierr)

c     output difference between implicit and explicit solutions
      ux = ux - ux_exp
c        ||diff||
      d2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,:)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,:))
      dIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,:)))
      maxnormloc = maxloc(abs(ux))-(/nghost,nghost,0,0/)

c        ||rho||
      rho2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,1)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,1))
      rhoIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,1)))

c        ||ux||
      ux2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,2)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,2))
      uxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,2)))

c        ||uy||
      uy2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,3)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,3))
      uyIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,3)))

c        ||uz||
      uz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,4)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,4))
      uzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,4)))

c        ||bx||
      bx2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,5)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,5))
      bxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,5)))

c        ||by||
      by2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,6)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,6))
      byIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,6)))

c        ||bz||
      bz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,7)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,7))
      bzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,7)))

c        ||e||
      e2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,8)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,8))
      eIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,8)))

c        scale wrms norms appropriately
      d2n   = dsqrt(d2n/nxlocal/nylocal/nzlocal/nvar)
      rho2n = dsqrt(rho2n/nxlocal/nylocal/nzlocal)
      ux2n  = dsqrt(ux2n/nxlocal/nylocal/nzlocal)
      uy2n  = dsqrt(uy2n/nxlocal/nylocal/nzlocal)
      uz2n  = dsqrt(uz2n/nxlocal/nylocal/nzlocal)
      bx2n  = dsqrt(bx2n/nxlocal/nylocal/nzlocal)
      by2n  = dsqrt(by2n/nxlocal/nylocal/nzlocal)
      bz2n  = dsqrt(bz2n/nxlocal/nylocal/nzlocal)
      e2n   = dsqrt(e2n/nxlocal/nylocal/nzlocal)
      maxlocrecv = (/ 0, maxnormloc(1), maxnormloc(2), 
     &     maxnormloc(3), maxnormloc(4) /)
      
c        root node outputs nonlinear residual information
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(2(A,es11.4))')
     &        '  ||diff(u)||  =',d2n,  ',   max|diff(u)|  =',dIn
         write(6,'(A,i4,A,4(i4),A)')'   max residual on proc',
     &        maxlocrecv(1),',\n     at (x,y,z,v) = (',
     &        maxlocrecv(2:5),')'
         write(6,'(A)') 
     &        '  -------------------------------------------------'
         write(6,'(2(A,es11.4))')
     &        '  ||diff_rho|| =',rho2n,',   max|diff_rho| =',rhoIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_ux||  =',ux2n, ',   max|diff_ux|  =',uxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uy||  =',uy2n, ',   max|diff_uy|  =',uyIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uz||  =',uz2n, ',   max|diff_uz|  =',uzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bx||  =',bx2n, ',   max|diff_Bx|  =',bxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_By||  =',by2n, ',   max|diff_By|  =',byIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bz||  =',bz2n, ',   max|diff_Bz|  =',bzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_e||   =',e2n,  ',   max|diff_e|   =',eIn
      end if
      



c     re-set CVODE to use rtol=1e-8, and re-initialize variables
      rtol = 1.d-8
      ux   = ux0
      ttot = t0
      call fcvreinit(ttot, ux, 2, rtol, abstol, 1, iopt, ropt, cvierr)

c     evolve to tstop using CVODE
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\n Starting implicit, rtol=1e-8:'
      end if
      call fcvode(tstop, ttot, ux, taskchoice, cvierr)

c     output difference between implicit and explicit solutions
      ux = ux - ux_exp
c        ||diff||
      d2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,:)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,:))
      dIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,:)))
      maxnormloc = maxloc(abs(ux))-(/nghost,nghost,0,0/)

c        ||rho||
      rho2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,1)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,1))
      rhoIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,1)))

c        ||ux||
      ux2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,2)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,2))
      uxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,2)))

c        ||uy||
      uy2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,3)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,3))
      uyIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,3)))

c        ||uz||
      uz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,4)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,4))
      uzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,4)))

c        ||bx||
      bx2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,5)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,5))
      bxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,5)))

c        ||by||
      by2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,6)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,6))
      byIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,6)))

c        ||bz||
      bz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,7)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,7))
      bzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,7)))

c        ||e||
      e2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,8)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,8))
      eIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,8)))

c        scale wrms norms appropriately
      d2n   = dsqrt(d2n/nxlocal/nylocal/nzlocal/nvar)
      rho2n = dsqrt(rho2n/nxlocal/nylocal/nzlocal)
      ux2n  = dsqrt(ux2n/nxlocal/nylocal/nzlocal)
      uy2n  = dsqrt(uy2n/nxlocal/nylocal/nzlocal)
      uz2n  = dsqrt(uz2n/nxlocal/nylocal/nzlocal)
      bx2n  = dsqrt(bx2n/nxlocal/nylocal/nzlocal)
      by2n  = dsqrt(by2n/nxlocal/nylocal/nzlocal)
      bz2n  = dsqrt(bz2n/nxlocal/nylocal/nzlocal)
      e2n   = dsqrt(e2n/nxlocal/nylocal/nzlocal)
      maxlocrecv = (/ 0, maxnormloc(1), maxnormloc(2), 
     &     maxnormloc(3), maxnormloc(4) /)
      
c        root node outputs nonlinear residual information
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(2(A,es11.4))')
     &        '  ||diff(u)||  =',d2n,  ',   max|diff(u)|  =',dIn
         write(6,'(A,i4,A,4(i4),A)')'   max residual on proc',
     &        maxlocrecv(1),',\n     at (x,y,z,v) = (',
     &        maxlocrecv(2:5),')'
         write(6,'(A)') 
     &        '  -------------------------------------------------'
         write(6,'(2(A,es11.4))')
     &        '  ||diff_rho|| =',rho2n,',   max|diff_rho| =',rhoIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_ux||  =',ux2n, ',   max|diff_ux|  =',uxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uy||  =',uy2n, ',   max|diff_uy|  =',uyIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uz||  =',uz2n, ',   max|diff_uz|  =',uzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bx||  =',bx2n, ',   max|diff_Bx|  =',bxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_By||  =',by2n, ',   max|diff_By|  =',byIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bz||  =',bz2n, ',   max|diff_Bz|  =',bzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_e||   =',e2n,  ',   max|diff_e|   =',eIn
      end if
      



c     re-set CVODE to use rtol=1e-9, and re-initialize variables
      rtol = 1.d-9
      ux   = ux0
      ttot = t0
      call fcvreinit(ttot, ux, 2, rtol, abstol, 1, iopt, ropt, cvierr)

c     evolve to tstop using CVODE
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\n Starting implicit, rtol=1e-9:'
      end if
      call fcvode(tstop, ttot, ux, taskchoice, cvierr)

c     output difference between implicit and explicit solutions
      ux = ux - ux_exp
c        ||diff||
      d2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,:)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,:))
      dIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,:)))
      maxnormloc = maxloc(abs(ux))-(/nghost,nghost,0,0/)

c        ||rho||
      rho2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,1)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,1))
      rhoIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,1)))

c        ||ux||
      ux2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,2)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,2))
      uxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,2)))

c        ||uy||
      uy2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,3)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,3))
      uyIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,3)))

c        ||uz||
      uz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,4)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,4))
      uzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,4)))

c        ||bx||
      bx2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,5)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,5))
      bxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,5)))

c        ||by||
      by2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,6)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,6))
      byIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,6)))

c        ||bz||
      bz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,7)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,7))
      bzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,7)))

c        ||e||
      e2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,8)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,8))
      eIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,8)))

c        scale wrms norms appropriately
      d2n   = dsqrt(d2n/nxlocal/nylocal/nzlocal/nvar)
      rho2n = dsqrt(rho2n/nxlocal/nylocal/nzlocal)
      ux2n  = dsqrt(ux2n/nxlocal/nylocal/nzlocal)
      uy2n  = dsqrt(uy2n/nxlocal/nylocal/nzlocal)
      uz2n  = dsqrt(uz2n/nxlocal/nylocal/nzlocal)
      bx2n  = dsqrt(bx2n/nxlocal/nylocal/nzlocal)
      by2n  = dsqrt(by2n/nxlocal/nylocal/nzlocal)
      bz2n  = dsqrt(bz2n/nxlocal/nylocal/nzlocal)
      e2n   = dsqrt(e2n/nxlocal/nylocal/nzlocal)
      maxlocrecv = (/ 0, maxnormloc(1), maxnormloc(2), 
     &     maxnormloc(3), maxnormloc(4) /)
      
c        root node outputs nonlinear residual information
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(2(A,es11.4))')
     &        '  ||diff(u)||  =',d2n,  ',   max|diff(u)|  =',dIn
         write(6,'(A,i4,A,4(i4),A)')'   max residual on proc',
     &        maxlocrecv(1),',\n     at (x,y,z,v) = (',
     &        maxlocrecv(2:5),')'
         write(6,'(A)') 
     &        '  -------------------------------------------------'
         write(6,'(2(A,es11.4))')
     &        '  ||diff_rho|| =',rho2n,',   max|diff_rho| =',rhoIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_ux||  =',ux2n, ',   max|diff_ux|  =',uxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uy||  =',uy2n, ',   max|diff_uy|  =',uyIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uz||  =',uz2n, ',   max|diff_uz|  =',uzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bx||  =',bx2n, ',   max|diff_Bx|  =',bxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_By||  =',by2n, ',   max|diff_By|  =',byIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bz||  =',bz2n, ',   max|diff_Bz|  =',bzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_e||   =',e2n,  ',   max|diff_e|   =',eIn
      end if
      



c     re-set CVODE to use rtol=1e-10, and re-initialize variables
      rtol = 1.d-10
      ux   = ux0
      ttot = t0
      call fcvreinit(ttot, ux, 2, rtol, abstol, 1, iopt, ropt, cvierr)

c     evolve to tstop using CVODE
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\n Starting implicit, rtol=1e-10:'
      end if
      call fcvode(tstop, ttot, ux, taskchoice, cvierr)

c     output difference between implicit and explicit solutions
      ux = ux - ux_exp
c        ||diff||
      d2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,:)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,:))
      dIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,:)))
      maxnormloc = maxloc(abs(ux))-(/nghost,nghost,0,0/)

c        ||rho||
      rho2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,1)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,1))
      rhoIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,1)))

c        ||ux||
      ux2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,2)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,2))
      uxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,2)))

c        ||uy||
      uy2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,3)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,3))
      uyIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,3)))

c        ||uz||
      uz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,4)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,4))
      uzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,4)))

c        ||bx||
      bx2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,5)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,5))
      bxIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,5)))

c        ||by||
      by2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,6)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,6))
      byIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,6)))

c        ||bz||
      bz2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,7)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,7))
      bzIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,7)))

c        ||e||
      e2n = sum(ux(1:nxlocal,1:nylocal,1:nzlocal,8)
     &     *ux(1:nxlocal,1:nylocal,1:nzlocal,8))
      eIn = maxval(abs(ux(1:nxlocal,1:nylocal,1:nzlocal,8)))

c        scale wrms norms appropriately
      d2n   = dsqrt(d2n/nxlocal/nylocal/nzlocal/nvar)
      rho2n = dsqrt(rho2n/nxlocal/nylocal/nzlocal)
      ux2n  = dsqrt(ux2n/nxlocal/nylocal/nzlocal)
      uy2n  = dsqrt(uy2n/nxlocal/nylocal/nzlocal)
      uz2n  = dsqrt(uz2n/nxlocal/nylocal/nzlocal)
      bx2n  = dsqrt(bx2n/nxlocal/nylocal/nzlocal)
      by2n  = dsqrt(by2n/nxlocal/nylocal/nzlocal)
      bz2n  = dsqrt(bz2n/nxlocal/nylocal/nzlocal)
      e2n   = dsqrt(e2n/nxlocal/nylocal/nzlocal)
      maxlocrecv = (/ 0, maxnormloc(1), maxnormloc(2), 
     &     maxnormloc(3), maxnormloc(4) /)
      
c        root node outputs nonlinear residual information
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(2(A,es11.4))')
     &        '  ||diff(u)||  =',d2n,  ',   max|diff(u)|  =',dIn
         write(6,'(A,i4,A,4(i4),A)')'   max residual on proc',
     &        maxlocrecv(1),',\n     at (x,y,z,v) = (',
     &        maxlocrecv(2:5),')'
         write(6,'(A)') 
     &        '  -------------------------------------------------'
         write(6,'(2(A,es11.4))')
     &        '  ||diff_rho|| =',rho2n,',   max|diff_rho| =',rhoIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_ux||  =',ux2n, ',   max|diff_ux|  =',uxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uy||  =',uy2n, ',   max|diff_uy|  =',uyIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_uz||  =',uz2n, ',   max|diff_uz|  =',uzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bx||  =',bx2n, ',   max|diff_Bx|  =',bxIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_By||  =',by2n, ',   max|diff_By|  =',byIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_Bz||  =',bz2n, ',   max|diff_Bz|  =',bzIn
         write(6,'(2(A,es11.4))')
     &        '  ||diff_e||   =',e2n,  ',   max|diff_e|   =',eIn
      end if



c     free up solver dataspace memory
      call fcvfree

c     dump all remaining data, output diagnostics, close I/O files
      call WrapUp(ux,phi,dt,ttot,timestep)


      return
      end subroutine ComparisonDriver
c=================================================================




      subroutine MHDMain
c-----------------------------------------------------------------
c     Description: main MHD routine, initializes more local grid
c        variables and allocates directional meshpoints.
c
c     Note: this is virtually unchanged from Ravi's code, it has 
c        only been cleaned up for commenting and clarity.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common

c======= Declarations =========
      implicit none

c======= Internals ============

c     mesh variables for a dyanmic grid
#ifdef DYNAMIC
      IXLO=1-nghost
      IYLO=1-nghost
      IZLO=1-nghost
      IXHI=nxlocal+nghost
      IYHI=nylocal+nghost
      IZHI=nzlocal+nghost
      
c     mesh variable adjustments for 2-D, dynamic grid
#ifdef TWO_D
      IZLO=1; IZHI=nzlocal
#endif

c     mesh variables for general dynamic grid
      inlo=min(IXLO,IYLO,IZLO)
      inhi=max(IXHI,IYHI,IZHI)
      nxlsize=nxlocal
      nylsize=nylocal
      nzlsize=nzlocal
#endif

c     Allocate directional meshpoints
      allocate(xc(IXLO:IXHI))
      allocate(yc(IYLO:IYHI))
      allocate(zc(IZLO:IZHI))

c     call the MHD driver routine
      call ComparisonDriver

c     free space used by directional meshpoints
      deallocate(xc,yc,zc)

      return
      end subroutine MHDMain
c=================================================================
