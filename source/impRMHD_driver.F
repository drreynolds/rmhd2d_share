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
c     $Log: impRMHD_driver.F,v $
c=================================================================




      subroutine ImpRMHDDriver
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
      use boundary_conds
      use iounits
      use tvdcoeffs
      use timelevel_data
      use sundials_fcmix
      use profiling
#ifdef PARALLEL
      use mpistuff
#endif
      
c======= Declarations =========
      implicit none

#ifdef DYNAMIC
      double precision, allocatable :: ux(:,:,:,:)
      double precision, allocatable :: phi(:,:)
#else
      double precision :: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: phi(ixlo:ixhi,iylo:iyhi)
#endif
      double precision :: ttot, dtexp, dtfac, tstop
      integer          :: istart, maxiter, lastiter
      integer          :: output_flag, binary_flag
      integer          :: ndump, ndiag, new
      integer          :: tvdFlag
      integer          :: timestep
      logical          :: outnode

c     kinsol nonlinear solver information
      integer(kind=sund_long_int) :: iopt(40)
      integer          :: totnewton, totkrylov, totpset, totpsol
      integer          :: globalstrat, kier
      double precision :: ropt(40)
#ifdef DYNAMIC
      double precision, dimension(:,:,:,:), allocatable :: uscale,fscale
#else
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uscale, fscale
#endif

c     declare internal namelists for problem specification
      namelist /inparams/ maxiter, new, output_flag, binary_flag
      namelist /dumps/ ndump, ndiag
      namelist /gridbounds/ xl, xr, yl, yr, zl, zr


c======= Internals ============

c     see if this is the output node
      outnode = .false.
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1)  outnode = .true.

c     if necessary, allocate local data arrays
#ifdef DYNAMIC
      allocate(ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(phi(ixlo:ixhi,iylo:iyhi))
      allocate(uscale(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(fscale(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#endif      


c     input general MHD simulation information using namelists
      open(16,file='mhd.inp',form='formatted')
      read(16,inparams)
      read(16,dumps)
      read(16,gridbounds) 
      close(16)

c     root node outputs problem parameters
      if(outnode) then
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
         write(6,'(A)') '\nBoundary conditions:'
         if (xbc == BCperiodic)  then
            write(6,'(A)') '    x: periodic'
         else if (xbc == BCreflecting) then
            write(6,'(A)') '    x: reflecting'
         else 
            write(6,'(A)') '    x: zero-gradient'
         endif
         if (ybc == BCperiodic) then
            write(6,'(A)') '    y: periodic'
         else if (ybc == BCreflecting) then
            write(6,'(A)') '    y: reflecting'
         else 
            write(6,'(A)') '    y: zero-gradient'
         endif
         if (zbc == BCperiodic) then
            write(6,'(A)') '    z: periodic'
         else if (zbc == BCreflecting) then
            write(6,'(A)') '    z: reflecting'
         else 
            write(6,'(A)') '    z: zero-gradient'
         endif
         write(6,'(A,3(A,i4),2(A,i2))') '\nGrid indices:',
     &        '\n    nx =',nx,',  ny =',ny,',  nz =',nz,
     &        ',  nvar =',nvar,',  nghost =',nghost
         write(6,'(A,3(A,i4))') '\nParallelism information:',
     &        '\n    xprocs =',xprocs,',  yprocs =',yprocs,
     &        ',  zprocs =',zprocs
      endif

c     initialize timestepping and data output information
      lastiter = 0
      istart   = 1

c     set up domain, grid
      call SetupDomain
      if(outnode) then
         write(6,*) '\nAfter domain setup...'
         write(6,'(A,3(A,es12.5))') '\nMesh spacing:',
     &        '\n   dx =',dx,',  dy =',dy,',  dz =',dz
      end if

c     determine local domain location in general domain
      call SetupLocalDomain
      if(outnode) then
         write(6,*) '\nAfter local domain setup...'
         write(6,'(A,3(i3,1x),A,3(2(es12.5,1x),A))') 
     &        '\nLocal subdomain for iproc (',
     &        iprocx,iprocy,iprocz,') is the box [',XLloc,XRloc,
     &        '] x [',YLloc,YRloc,'] x [',ZLloc,ZRloc,']'
      endif

c$$$c     set up TVD coefficients
c$$$      tvdFlag=1
c$$$      if(tvdFlag.eq.1) call SetTVDCoeffs
  
c     initialize state variables, Root node outputs progress
      call Initialize(new,ux,phi,dtexp,ttot,istart,lastiter)
      if(outnode)  write(6,*) '\nFinished local initialization'
      
c     communicate initial boundary data among processors
      if(outnode) 
     &     write(6,*) 'Communicating initial condition to neighbors...'
      call CommState(ux)

c     initialize KINSOL solver structures
      if(outnode) 
     &     write(6,*) 'Initializing KINSOL solver...'
      call KinsolInit(uscale, fscale, iopt, ropt, 
     &     dtfac, tstop, globalstrat)

c     save rhs at initial state
      call computerhs(ux,dFnew)

c     if parallel, wait until all processors are ready
#ifdef PARALLEL
      call mpi_barrier(comm3d, ierr)
#endif


c     output initial state and diagnostics (energies, mass, etc), 
c        start the diags & I/O profiler
      call prof_start(6)
      if(output_flag==1) then
         if(binary_flag==0) then
            call WriteGnuplotFile(ux,0)
c            call WriteAVSFile(ux,0)
         endif
      endif
c     Root node outputs progress
      call Diagnostics(ux,phi,ttot)
c        stop the diags & I/O profiler
      call prof_stop(6)

      if(outnode) write(6,*) '\n Finished initialization,',
     &     ' starting time steps'
      totnewton = 0
      totkrylov = 0
      totpset   = 0
      totpsol   = 0

c     start MPI timer
      call prof_start(1)

c     Start time stepping
      do timeStep = 1,maxiter

c        get timestep size as a scaled explicit time step
         call newdt(dtexp,ux)
         dt = min(dtfac, tstop-ttot)
c$$$         dt = min(dtfac*dtexp, tstop-ttot)
         if(outnode) then
            write(6,*) 'Explicit time step =',dtexp
            write(6,*) 'Desired time step  =',dt
         endif

c        set current timelevel information into saved variables
         dFold = dFnew
         uold  = ux

c        output initial residual information
         call prof_start(6)
         if(outnode) write(6,'(A)')  
     &        '\n  ------------- before kinsol solution ------------'
         call fkfunDiags(ux,fscale)
         if(outnode) write(6,'(A)') 
     &        '  -------------------------------------------------'
         call prof_stop(6)

c        call Kinsol to solve for new time step
         call fkinsol(ux, globalstrat, uscale, fscale, kier)

c        wait for all processes to catch up
#ifdef PARALLEL
         call mpi_barrier(comm3d, ierr)
#endif

c        output updated residual information after solve
         if(outnode) write(6,'(A)')  
     &        '\n  ------------- after kinsol solution -------------'
         call fkfunDiags(ux,fscale)
         if(outnode) write(6,'(A)') 
     &        '  -------------------------------------------------'

c        error if failure, solver diagnostics if success
         call kinerr(kier,iopt,ropt)
         totnewton = totnewton + iopt(3)
         totkrylov = totkrylov + iopt(14)
         totpset   = totpset   + iopt(12)
         totpsol   = totpsol   + iopt(13)
         if ( kier == 1 ) then
            if(outnode) write(6,'(A,i7,2(A,es12.5))') 'timeStep =',
     &           timeStep,',   dt =',dt,',   ttot =',ttot            
            exit
         endif

c        various div B cleaning routines
#ifdef DIVBCLEAN
         call CommState(ux)
         call DivSource(ux,uold,dt)
cccccccc TURNED OFF DIVBDIFFUSIONFILTER PER RAVI'S REQUEST, 12/3/03
c         call DivBDiffusionFilter(ux)
#endif

c        do a final update of boundary information [?is this necessary?]
         call CommState(ux)
         
c        output updated residual information after div B cleaning
#ifdef DIVBCLEAN
         call prof_start(6)
         if(outnode) write(6,'(A)') 
     &        '  ------------- after (divB) cleaning -------------'
         call fkfunDiags(ux,fscale)
         if(outnode) write(6,'(A)') 
     &        '  -------------------------------------------------'
         call prof_stop(6)
#endif

c        update time, periodically dump data
         call prof_start(6)
         ttot = ttot + dt
         if(mod(timeStep,ndump)==0) then
            if(output_flag==1) then
               if(binary_flag==0) then
c                  call WriteTextFileParallel(ux,timeStep)
                  call WriteGnuplotFile(ux,timeStep)
c                  call WriteAVSFile(ux,timeStep)
               else
                  call WriteBinaryFileParallel(ux,timeStep)
               endif
            endif
            if (output_flag==2) then
c              dump restart information only
               if(outnode) then
                  open(itime,file='time.txt',form='formatted')
                  write(itime,*) dt
                  write(itime,*) ttot
                  write(itime,*) timestep-1
                  close(itime)
               endif
               call DumpOutFluidParallel(ux,phi)
            endif
         endif

c        root outputs timing information
         if(outnode) write(6,   '(A,i7,2(A,es12.5))') 'timeStep =',
     &        timeStep,',   dt =',dt,',   ttot =',ttot
         

c        periodically output diagnostic information
         if(mod(timeStep,ndiag)==0) then
            call Diagnostics(ux,phi,ttot)
c$$$            call WaveSolutionError(ux,ttot)
         endif
         call prof_stop(6)

c        exit time stepping loop if we have exceeded tstop
         if ( ttot >= tstop ) then
            if (outnode) write(6,'(A,es9.2,A)') 
     &           '\n Physical stop time, t =',tstop,
     &           ' reached.  Stopping simulation.'
            call Diagnostics(ux,phi,ttot)
c$$$            call WaveSolutionError(ux,ttot)
            exit
         endif

      enddo

c     stop MPI timer
      call prof_stop(1)

c     output final solution
      call prof_start(6)
      if(output_flag == 1) then
         if(binary_flag == 0) then
c            call WriteTextFileParallel(ux,timeStep)
            call WriteGnuplotFile(ux,timeStep)
c            call WriteAVSFile(ux,timeStep)
         else
            call WriteBinaryFileParallel(ux,timeStep)
         endif
      endif
      call prof_stop(6)


c     output final solver diagnostics
      if (outnode) write(6,'(A,A,2(A,i9),13(A,es12.5),A)') 
     &     '\n  ---------------------------------------',
     &     '\n  Runtime statistics:\n',
     &     '\n    Total Newton Iterations =',totnewton, 
     &     '\n    Total Krylov Iterations =',totkrylov,
     &     '\n    Total Solution Time     =',prof_time(1),
     &     '\n    Total RHS Time          =',prof_time(2),
     &     '\n    Total Comm Time         =',prof_time(3),
     &     '\n    Total Diags & I/O Time  =',prof_time(6),
     &     '\n    Total Cons->Prim Time   =',prof_time(7),
     &     '\n    Total Flux Init Time    =',prof_time(8),
     &     '\n    Total Divergence Time   =',prof_time(15),
     &     '\n    Avg RHS Time            =',prof_time(2)/prof_count(2),
     &     '\n    Avg Comm Time           =',prof_time(3)/prof_count(3),
     &     '\n    Avg Diags & I/O Time    =',prof_time(6)/prof_count(6),
     &     '\n    Avg Cons->Prim Time     =',prof_time(7)/prof_count(7),
     &     '\n    Avg Flux Init Time      =',prof_time(8)/prof_count(8),
     &   '\n    Avg Divergence Time     =',prof_time(15)/prof_count(15),
     &     '\n  ---------------------------------------\n'

c     dump all remaining data, output diagnostics, close I/O files
      call WrapUp(ux,phi,dt,ttot,timestep)

c     free up allocated memory
      call KinsolFree
#ifdef DYNAMIC
      deallocate(ux,phi,uscale,fscale)
#endif      


      return
      end subroutine ImpRMHDDriver
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
      ixlo=1-nghost
      iylo=1-nghost
      izlo=1-nghost
      ixhi=nxlocal+nghost
      iyhi=nylocal+nghost
      izhi=nzlocal+nghost
      
c     mesh variable adjustments for 2-D, dynamic grid
#ifdef TWO_D
      izlo=1; izhi=nzlocal
#endif

c     mesh variables for general dynamic grid
      inlo=min(ixlo,iylo,izlo)
      inhi=max(ixhi,iyhi,izhi)
      nxlsize=nxlocal
      nylsize=nylocal
      nzlsize=nzlocal
#endif

c     Allocate directional meshpoints
      allocate(xc(ixlo:ixhi))
      allocate(yc(iylo:iyhi))
      allocate(zc(izlo:izhi))

c     call the MHD driver routine
      call ImpRMHDDriver

c     free space used by directional meshpoints
      deallocate(xc,yc,zc)

      return
      end subroutine MHDMain
c=================================================================
