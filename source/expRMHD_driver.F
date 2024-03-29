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
c     $Log: expRMHD_driver.F,v $
c=================================================================



      subroutine MHDMain
c-----------------------------------------------------------------
c     Description: main MHD routine, initializes grid variables
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
      call MHDDriver
      end
c=================================================================



      subroutine MHDDriver
c-----------------------------------------------------------------
c     Description: main MHD driver, handles initialization of 
c     variables, timestepping and data I/O.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use iounits
c      use timemod
c      use tvdcoeffs
      use profiling
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI)
      double precision :: dt, ttot, tstop
      integer          :: istart, iter, maxiter, lastiter
      integer          :: output_flag, binary_flag
      integer          :: ndump, ndiag
      integer          :: dump,new
      integer          :: ipar
      integer          :: tvdFlag
      integer          :: timestep


c======= Internals ============

c     declare internal namelists for problem specification
      namelist /inparams/ maxiter, new, output_flag, binary_flag
      namelist /dumps/ ndump, ndiag
      namelist /gridbounds/ xl, xr, yl, yr, zl, zr

c     set up physical stopping time
      tstop = 4.0d2

c     input general MHD simulation information using namelists
      open(16,file='mhd.inp',form='formatted')
      read(16,inparams)
      read(16,dumps)
      read(16,gridbounds) 
      close(16)

c     root node opens output streams and outputs problem parameters
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(6,inparams)
         write(6,dumps)
         write(6,gridbounds)
      endif

c     initialize timestepping and data output information
      lastiter = 0
      istart   = 1
      dump     = 1

c     set up domain, grid
      call SetupDomain
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) 'After domain setup'
      end if

c     determine local domain location in general domain
      call SetupLocalDomain
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\nAfter local domain setup...'
      endif

c     set up TVD coefficients
      tvdFlag=1
      if(tvdFlag.eq.1) call SetTVDCoeffs
      
c     initialize state variables
      call Initialize(new,ux,phi,dt,ttot,istart,lastiter)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) 'After local initialization'
      end if

#ifdef PARALLEL
      call MPI_Barrier(comm3d, ierr)
#endif

c     communicate initial boundary data among processors
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)

c     if parallel, wait until all processors are ready
#ifdef PARALLEL
      call MPI_Barrier(comm3d, ierr)
#endif

c     output initial diagnostics (energies, mass, etc)
      call Diagnostics(ux,phi,ttot)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) 'After final initialization'
      end if

c     start timing information for speed check
c      call Timing

c     Output initial state 
      if(output_flag==1) then
         if(binary_flag==0) then
            call WriteGnuplotFile(ux,0)
c            call WriteAVSFile(ux,0)
         endif
      endif

c     start MPI timer
      call prof_start(1)

c     Start time stepping
      do timeStep = istart,istart+maxiter-1,1

c        solve for new time step
         call MHDSolve(ux,phi,dt)

c        root outputs timing information
         if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
            write(6,'(A,i8,2(A,es12.5))') 'timeStep =',timeStep, 
     &           ',   dt =',dt,',   ttot =',ttot
         endif

c        update time, periodically dump data
         ttot = ttot + dt
         if(mod(timeStep,ndump).eq.0) then
            if(output_flag == 1) then
               if(binary_flag == 0) then
c                  call WriteTextFileParallel(ux,timeStep)
                  call WriteGnuplotFile(ux,timeStep)
c                  call WriteAVSFile(ux,timeStep)
               else
                  call WriteBinaryFileParallel(ux,timeStep)
               endif
            endif
            if (output_flag == 2) then
c              dump restart information only
               if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
                  open(itime,file='time.txt',form='formatted')
                  write(itime,*) dt
                  write(itime,*) ttot
                  write(itime,*) timestep-1
                  close(itime)
               endif
               call DumpOutFluidParallel(ux,phi)
            endif
         endif

c        periodically output diagnostic information
         if(mod(timeStep,ndiag).eq.0) then
            call Diagnostics(ux,phi,ttot)
c            call WaveSolutionError(ux,ttot)
         endif

c        update time step to satisfy CFL condition
         call newdt(dt,ux)

c        exit time stepping loop if tstop reached
         if ( ttot > tstop ) then
            if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
               write(6,'(A,es9.2,A)') '\n Physical stop time, t =',
     &              tstop,' reached.  Stopping simulation.'
            endif
            call Diagnostics(ux,phi,ttot)
c            call WaveSolutionError(ux,ttot)
            exit
         endif

      enddo
     
c     stop MPI timer
      call prof_stop(1)

c     output final solution
      if(output_flag == 1) then
         if(binary_flag == 0) then
c            call WriteTextFileParallel(ux,timeStep)
            call WriteGnuplotFile(ux,timeStep)
c            call WriteAVSFile(ux,timeStep)
         else
            call WriteBinaryFileParallel(ux,timeStep)
         endif
      endif

c     output final solver diagnostics
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) 
     &     write(6,'(A,A,A,es12.5,A)') 
     &     '\n  ---------------------------------------',
     &     '\n  Runtime statistics:\n',
     &     '\n    Total Solution Time     =',prof_time(1),
     &     '\n  ---------------------------------------\n'

c     dump all remaining data, output diagnostics, close I/O files
      call WrapUp(ux,phi,dt,ttot,timestep)

c     stop timing information for speed check
c      call Timing

      return
      end
c=================================================================
