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
c     $Log: bdryexchange.F,v $
c-----------------------------------------------------------------
      
      integer function iprocnum(ipx, ipy, ipz)
c     
c     precondition:
c               ipx, ipy are coordinates of process in process
c               grid.
c     postcondition:
c               return value is a single integer denoting process
c               -- unique for each process.
c       
      use mesh_parms
      integer ipx, ipy, ipz
      
      iprocnum = (ipz-1)*XPROCS*YPROCS+(ipy-1)*XPROCS + ipx
      
      return
      end function iprocnum
c     ----------------------------------------------------------------
      
      subroutine FluidProcCoord(fluidprocid)
c     
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif
      integer::fluidprocid
c     integer:: ijkp(3)
c     ijkp(1)=iprocx; ijkp(2)=iprocy; ijkp(3)=iprocz;
      fluidprocid=iproc_idx-1
      return
      end subroutine FluidProcCoord
c     ----------------------------------------------------------------
      
      subroutine pxpypz(iproc, ipx, ipy, ipz)
c     
c     precondition:
c		iproc is a single integer denoting process
c			-- unique for each process.
c     postcondition:
c		ipx, ipy are coordinates of process in process.
c			grid.
c		(i.e., this subroutine is the inverse of function
c			iprocnum above.)
c
      use mesh_parms
#ifdef PARALLEL
      use mpistuff
      integer NDIM
      PARAMETER(NDIM=3)
      integer coords(NDIM)
#endif
      integer ipx, ipy, ipz,  iproc
      
      
c     compute ipx, ipy, ipz
      
c02/21/98 for reorder
#ifdef PARALLEL
      call MPI_Cart_Coords(comm3D, iproc-1, NDIM, coords, ierr)
      call ErrorHandler(ierr,ERROR_CARTCOORDS)
      ipx = coords(1) + 1
      ipy = coords(2) + 1
      ipz = coords(3) + 1
      call MPI_Cart_Shift( comm3D, 0, 1, left, right, ierr)
      call ErrorHandler(ierr,ERROR_CARTSHIFT)
      call MPI_Cart_Shift( comm3D, 1, 1, bottom, top, ierr)
      call ErrorHandler(ierr,ERROR_CARTSHIFT)
      call MPI_cart_Shift( comm3D, 2, 1, behind, forward, ierr)
      call ErrorHandler(ierr,ERROR_CARTSHIFT)
#else
      ipx = mod(iproc-1, XPROCS)+1
      ipy = mod( (iproc - ipx) / XPROCS, YPROCS)+1
      ipz = mod((iproc-ipx-(ipy-1)*XPROCS)/XPROCS/YPROCS,ZPROCS) + 1
#endif
      
      return
      end subroutine pxpypz
      
c     ================================================================
c     
c     subroutines to map from global to local indices
c     
c     
      subroutine ilocal(ii, ipx, iiloc)
c     
c     precondition:
c		ii is a global X index.
c     postcondition:
c		ipx, iiloc are the corresponding X index into the
c			process grid and local X index.
c     
      use mesh_parms
      use mesh_common
      integer ii, ipx, iiloc
      
      
      ipx = (ii-1) / NXlsize + 1
      iiloc = mod((ii-1), NXlsize) + 1
      return
      end subroutine ilocal
      
c     ----------------------------------------------------------------
c     
      subroutine jlocal(jj, ipy, jjloc)
c     
c     precondition:
c		jj is a global Y index.
c     postcondition:
c		ipy, jjloc are the corresponding Y index into the
c			process grid and local Y index.
c     
      use mesh_parms
      use mesh_common
      integer jj, ipy, jjloc
      
      ipy = (jj-1) / NYlsize + 1
      jjloc = mod((jj-1), NYlsize) + 1
      
      return
      end subroutine jlocal
      
c     ----------------------------------------------------------------
      
      subroutine klocal(kk, ipz, kkloc)
c     
c     precondition:
c		kk is a global Z index.
c     postcondition:
c		ipy, kkloc are the corresponding Z index into the
c			process grid and local Z index.
c     
      use mesh_parms
      use mesh_common
      integer kk, ipz, kkloc
      
      ipz = (kk-1) / NZlsize + 1
      kkloc = mod((kk-1), NZlsize) + 1
      
      return
      end subroutine klocal
      
c     ================================================================
c     
c     subroutines to map from local to global indices
c     
c     ================================================================
      
      subroutine iglobal(ipx, iiloc, ii)
c     
c     precondition:
c		ipx, iiloc are an X index into the process grid and
c			a local X index.
c     postcondition:
c		ii is the corresponding global X index.
c     
      
      use mesh_parms
      use mesh_common
      integer ipx, iiloc, ii
      ii = (ipx-1)*NXlsize + iiloc
      return
      end subroutine iglobal
      
c     ----------------------------------------------------------------
      
      subroutine jglobal(ipy, jjloc, jj)
c     
c     precondition:
c		ipy, jjloc are an Y index into the process grid and
c			a local Y index.
c     postcondition:
c		jj is the corresponding global Y index.
c
      use mesh_parms
      use mesh_common
      integer ipy, jjloc, jj
      jj = (ipy-1)*NYlsize + jjloc
      return
      end subroutine jglobal
      
c     ----------------------------------------------------------------
      
      subroutine kglobal(ipz, kkloc, kk)
c     
c     precondition:
c		ipz, kkloc are an Z index into the process grid and
c			a local Z index.
c     postcondition:
c		kk is the corresponding global Z index.
c     
      use mesh_parms
      use mesh_common
      integer ipz, kkloc, kk
      kk = (ipz-1)*NZlsize + kkloc
      return
      end subroutine kglobal
      
c     ================================================================
      
      subroutine mesh_update_bdry_async(grid_array,nfields)
c     
c     precondition:
c               grid_array holds the local data to be communicated
c               nfields holds the number of distinct fluid variables
c     postcondition:
c               grid_array holds the local data, post communication
c
c     Note: there are a number of #ifdef statements, meant to test 
c     for directional periodicity, depending on whether processes 
c     are on the domain boundary or not.
c     
      use mesh
      use mesh_common
      use boundary_conds
#ifdef PARALLEL
      use mpistuff
#endif
      
      integer :: nfields
      double precision, dimension(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nfields) 
     &     :: grid_array
      
#ifdef PARALLEL
      double precision, dimension(NGHOST,IYLO:IYHI,IZLO:IZHI,nfields) 
     &     :: xbuffer_send, xbuffer_recv
#ifndef ONE_D
      double precision, dimension(IXLO:IXHI,NGHOST,IZLO:IZHI,nfields) 
     &     :: ybuffer_send, ybuffer_recv
#ifndef TWO_D
      double precision, dimension(IXLO:IXHI,IYLO:IYHI,NGHOST,nfields) 
     &     :: zbuffer_send, zbuffer_recv
#endif
#endif

      integer:: XBUFFSIZE, YBUFFSIZE, ZBUFFSIZE
      integer:: ii,jj,mm,kk,l,idest
      integer:: iprocnum
      integer:: msg_id_send_x_low
      integer:: msg_id_send_x_hi
      integer:: msg_id_recv_x_low
      integer:: msg_id_recv_x_hi
      
      integer:: msg_id_send_y_low
      integer:: msg_id_send_y_hi
      integer:: msg_id_recv_y_low
      integer:: msg_id_recv_y_hi
      
      integer:: msg_id_send_z_low
      integer:: msg_id_send_z_hi
      integer:: msg_id_recv_z_low
      integer:: msg_id_recv_z_hi
c     
c
c     set buffer sizes
      XBUFFSIZE = nfields*(NGHOST)*(IYHI-IYLO+1)*(IZHI-IZLO+1)
      YBUFFSIZE = nfields*(IXHI-IXLO+1)*(NGHOST)*(IZHI-IZLO+1)
      ZBUFFSIZE = nfields*(IXHI-IXLO+1)*(IYHI-IYLO+1)*(NGHOST) 

      
c     -------X DIRECTION COMMUNICATION, done for all simulations
c
c     
c     Update x-low boundaries
c
      if ((iprocx > 1) .or. (xbc == BCperiodic)) then
         call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        left, MSG_XCH_XHI_TAG, comm3D, msg_id_recv_x_hi, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
      endif


      if ((iprocx < XPROCS) .or. (xbc == BCperiodic)) then
         do l = 1,nfields
            do kk = IZLO,IZHI
               do jj = 1-NGHOST,NYlocal+NGHOST
                  do mm = 1,NGHOST
                     xbuffer_send(mm,jj,kk,l) = 
     &                    grid_array(NXlsize+1-mm,jj,kk,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        right, MSG_XCH_XHI_TAG, comm3D, msg_id_send_x_hi, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
      endif


      if ((iprocx > 1) .or. (xbc == BCperiodic)) then
         call MPI_Wait(msg_id_recv_x_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do kk = IZLO,IZHI
               do jj = 1-NGHOST,NYlocal+NGHOST
                  do mm = 1,NGHOST
                     grid_array(1-mm,jj,kk,l) = 
     &                    xbuffer_recv(mm,jj,kk,l)
                  enddo
               enddo
            enddo
         enddo
      endif
      
      if ((iprocx < XPROCS) .or. (xbc == BCperiodic)) then
         call MPI_Wait(msg_id_send_x_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
      endif
c
c      
c     update x-high boundaries
c
      if ((iprocx < XPROCS) .or. (xbc == BCperiodic)) then
         call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        right, MSG_XCH_XLOW_TAG, comm3D, msg_id_recv_x_low, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
      endif
      

      if ((iprocx > 1) .or. (xbc == BCperiodic)) then
         do l = 1,nfields
            do kk = IZLO,IZHI
               do jj = 1-NGHOST,NYlocal+NGHOST
                  do mm = 1,NGHOST
                     xbuffer_send(mm,jj,kk,l) = 
     &                    grid_array(mm,jj,kk,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        left, MSG_XCH_XLOW_TAG, comm3D, msg_id_send_x_low, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
      endif

      
      if ((iprocx < XPROCS) .or. (xbc == BCperiodic)) then
         call MPI_Wait(msg_id_recv_x_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do kk = IZLO,IZHI
               do jj = 1-NGHOST,NYlocal+NGHOST
                  do mm = 1,NGHOST
                     grid_array(NXlsize+mm,jj,kk,l) = 
     &                    xbuffer_recv(mm,jj,kk,l)
                  enddo
               enddo
            enddo
         enddo
      endif

      
      if ((iprocx > 1) .or. (xbc == BCperiodic)) then
         call MPI_Wait(msg_id_send_x_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
      endif
c
c     
c     
c     -------Y DIRECTION COMMUNICATION, done for 2D, 2.5D, 3D runs
#ifndef ONE_D
c     
c     
c     update y-low boundaries
      if ((iprocy > 1) .or. (ybc == BCperiodic)) then
         call  MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        bottom, MSG_XCH_YHI_TAG, comm3D, msg_id_recv_y_hi, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
      endif
      

      if ((iprocy < YPROCS) .or. (ybc == BCperiodic)) then
         do l = 1,nfields
            do kk = IZLO,IZHI
               do mm = 1,NGHOST
                  do ii = 1-NGHOST,NXlocal+NGHOST
                     ybuffer_send(ii,mm,kk,l) = 
     &                    grid_array(ii,NYlsize+1-mm,kk,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        top, MSG_XCH_YHI_TAG, comm3D, msg_id_send_y_hi, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
      endif

      
      if ((iprocy > 1) .or. (ybc == BCperiodic)) then
         call MPI_Wait(msg_id_recv_y_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do kk = IZLO,IZHI
               do mm = 1,NGHOST
                  do ii = 1-NGHOST,NXlocal+NGHOST
                     grid_array(ii,1-mm,kk,l) = 
     &                    ybuffer_recv(ii,mm,kk,l)
                  enddo
               enddo
            enddo
         enddo
      endif
      

      if ((iprocy < YPROCS) .or. (ybc == BCperiodic)) then
         call MPI_Wait(msg_id_send_y_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
      endif
c
c      
c     update y-high boundaries
c
      if ((iprocy < YPROCS) .or. (ybc == BCperiodic)) then
         call MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        top, MSG_XCH_YLOW_TAG, comm3D, msg_id_recv_y_low, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
      endif

      
      if ((iprocy > 1) .or. (ybc == BCperiodic)) then
         do l = 1,nfields
            do kk = IZLO,IZHI
               do mm = 1,NGHOST
                  do ii = 1-NGHOST,NXlocal+NGHOST
                     ybuffer_send(ii,mm,kk,l) = 
     &                    grid_array(ii,mm,kk,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        bottom, MSG_XCH_YLOW_TAG, comm3D, msg_id_send_y_low, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
      endif

      
      if ((iprocy < YPROCS) .or. (ybc == BCperiodic)) then
         call MPI_Wait(msg_id_recv_y_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do kk = IZLO,IZHI
               do mm = 1,NGHOST
                  do ii = 1-NGHOST,NXlocal+NGHOST
                     grid_array(ii,NYlsize+mm,kk,l) = 
     &                    ybuffer_recv(ii,mm,kk,l)
                  enddo
               enddo
            enddo
         enddo
      endif

      
      if ((iprocy > 1) .or. (ybc == BCperiodic)) then
         call MPI_Wait(msg_id_send_y_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
      endif
c
c     
c     
c     -------Z DIRECTION COMMUNICATION, done for 3D simulations only
#ifndef TWO_D
c
c     
c     update z-low boundaries
      if ((iprocz > 1) .or. (zbc == BCperiodic)) then
         call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        behind, MSG_XCH_ZHI_TAG, comm3D, msg_id_recv_z_hi, ierr) 
         call ErrorHandler(ierr,ERROR_RECV)
      endif

      
      if ((iprocz < ZPROCS) .or. (zbc == BCperiodic)) then
         do l = 1,nfields
            do mm = 1,NGHOST
               do jj = 1-NGHOST,NYlocal+NGHOST
                  do ii = 1-NGHOST,NXlocal+NGHOST
                     zbuffer_send(ii,jj,mm,l) = 
     &                    grid_array(ii,jj,NZlsize+1-mm,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        forward, MSG_XCH_ZHI_TAG, comm3D, msg_id_send_z_hi, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
      endif
      

      if ((iprocz > 1) .or. (zbc == BCperiodic)) then
         call MPI_Wait(msg_id_recv_z_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do mm = 1,NGHOST
               do jj = 1-NGHOST,NYlocal+NGHOST
                  do ii = 1-NGHOST,NXlocal+NGHOST
                     grid_array(ii,jj,1-mm,l) = 
     &                    zbuffer_recv(ii,jj,mm,l)
                  enddo
               enddo
            enddo
         enddo
      endif
      

      if ((iprocz < ZPROCS) .or. (zbc == BCperiodic)) then
         call MPI_Wait(msg_id_send_z_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
      endif
c      
c     update z-high boundaries
c
      if ((iprocz < ZPROCS) .or. (zbc == BCperiodic)) then
         call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        forward, MSG_XCH_ZLOW_TAG, comm3D,msg_id_recv_z_low, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
      endif

      
      if ((iprocz > 1) .or. (zbc == BCperiodic)) then
         do l = 1,nfields
            do mm = 1,NGHOST
               do jj = 1-NGHOST,NYlocal+NGHOST
                  do ii = 1-NGHOST,NXlocal+NGHOST
                     zbuffer_send(ii,jj,mm,l) = 
     &                    grid_array(ii,jj,mm,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &	      behind, MSG_XCH_ZLOW_TAG, comm3D, msg_id_send_z_low, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
      endif

      
      if ((iprocz < ZPROCS) .or. (zbc == BCperiodic)) then
         call MPI_Wait(msg_id_recv_z_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do mm = 1,NGHOST
               do jj = 1-NGHOST,NYlocal+NGHOST
                  do ii = 1-NGHOST,NXlocal+NGHOST
                     grid_array(ii,jj,NZlsize+mm,l) = 
     &                    zbuffer_recv(ii,jj,mm,l)
                  enddo
               enddo
            enddo
         enddo
      endif

      
      if ((iprocz > 1) .or. (zbc == BCperiodic)) then
         call MPI_Wait(msg_id_send_z_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
      endif

#endif
c end TWO_D
#endif
c end ONE_D
#endif
c end PARALLEL
	
      return
      end subroutine mesh_update_bdry_async
c-----------------------------------------------------------------------
      
      subroutine ErrorHandler(mpierr,errortype)
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif
      integer::mpierr,errortype
#ifdef PARALLEL
      if(mpierr.ne.MPI_SUCCESS) then
         write(0,*) 'FLUID: MPI RETURN VALUE',iproc_idx,mpierr,errortype
      endif
#endif
      return
      end subroutine ErrorHandler
c-----------------------------------------------------------------------

