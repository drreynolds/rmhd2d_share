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
c     $Log: newdt.F,v $
c=================================================================



      subroutine newdt (dt,ux)
c-----------------------------------------------------------------
c     Description: compute stable time step according to CFL cond.
c
c     Arguments:
c          dt - (out) new timestep
c          ux - (in) conserved variables
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
#ifdef PARALLEL
      use mpistuff
#endif  

c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI, nvar)
      double precision :: dt
      double precision :: rho
      double precision :: xspmax, yspmax, zspmax
      double precision :: u, v, w, press
      double precision :: bi,bj,bk
      double precision :: dt_local, dtdiff
      double precision :: rhoInv
      double precision :: Asq
      double precision :: axsq, aysq, azsq
      double precision :: csndsq
      double precision :: cfastx, cfasty, cfastz
      double precision :: cslowx, cslowy, cslowz
      double precision :: csnd, wavespeedx, wavespeedy, wavespeedz
      double precision :: tmp, tmpx, tmpy, tmpz, dm, etaPrime
      double precision :: cfl
      integer :: i,j,k


c======= Internals ============

c     compute maximum velocity and adjust time step to satisfy
c     cfl criterion
      cfl = 0.2
      dt_local = -1.0e0
      dtdiff = 1.0d16
      xspmax = 1.D99
      yspmax = 1.D99
      zspmax = 1.D99
      dm = dmin1(dx,dy,dz)
      do k = 1,nzlocal
         do j = 1,nylocal
            do i = 1,nxlocal
               rho = ux(i,j,k,1)
               u = ux(i,j,k,2)/rho
               v = ux(i,j,k,3)/rho
               w = ux(i,j,k,4)/rho
               bi = ux(i,j,k,5)
               bj = ux(i,j,k,6)
               bk = ux(i,j,k,7)
               rhoInv = 1.D0/ux(i,j,k,1)
#ifdef ENTROPY_FORMULATION
               press = rho**gamma*dexp(ux(i,j,k,8)/ux(i,j,k,1))
#else
               press = (gamma-1.D0)*(ux(i,j,k,8)-0.5D0*rho*
     &              (u*u+v*v+w*w)-0.5*(bi*bi+bj*bj+bk*bk))
#endif

               axsq = bi**2*rhoInv
               aysq = bj**2*rhoInv
               azsq = bk**2*rhoInv
c               Asq = (bi**2+bj**2+bk**2)*rhoInv
               Asq = axsq+aysq+azsq
               csndsq = gamma*press*rhoInv
               tmpx = dsqrt(dmax1(((csndsq+Asq)**2-4.0*csndsq*axsq),
     &              0.D0))
               tmpy = dsqrt(dmax1(((csndsq+Asq)**2-4.0*csndsq*aysq),
     &              0.D0))
c               tmpz = dsqrt(dmax1(((csndsq+Asq)**2-4.0*csndsq*azsq),
c    &	            0.D0))
               tmp = csndsq+Asq
               cfastx = dsqrt(0.5*(tmp+tmpx))
c               cslowx = dsqrt(0.5*(tmp-tmpx))
               cfasty = dsqrt(0.5*(tmp+tmpy))
c               cslowy = dsqrt(0.5*(tmp-tmpy))
c               cfastz = dsqrt(0.5*(tmp+tmpz))
c               cslowz = dsqrt(0.5*(tmp-tmpz))
               csnd = dsqrt(csndsq)
c               waveSpeed = dmax1(cfast,cslow,dsqrt(axsq))
               waveSpeedx = dmax1(cfastx,dsqrt(axsq))
               waveSpeedy = dmax1(cfasty,dsqrt(aysq))
c               waveSpeedz = dmax1(cfastz,dsqrt(azsq))
               xspmax = dmin1(xspmax,dx/(dabs(u)+waveSpeedx))
               yspmax = dmin1(yspmax,dy/(dabs(v)+waveSpeedy))
c               zspmax = dmin1(zspmax,dz/(dabs(w)+waveSpeedz))
#ifndef TWO_D
               dt_local = dmin1(xspmax,yspmax,zspmax)*cfl
#else
               dt_local = dmin1(xspmax,yspmax)*cfl
#endif

#ifndef NOVISCFLUX
c$$$               if(bk.le.-1.D-04) then		 
c$$$                  etaPrime = eta + etaFactor*eta*
c$$$     &			 (1.D0-dexp(-(bk/0.075D0/7.8D0)**2))
c$$$               else
                  etaPrime = eta
c               endif
               dtdiff = dmin1(dm*dm*cfl/etaPrime,dtdiff)
               dtdiff = dmin1(dm*dm*cfl/mu,dtdiff)
#else
               dtdiff = dt_local
#endif
            enddo
         enddo
      enddo
c      write(6,*) 'NEWDT',dt_local,dtdiff
      dt_local=dmin1(dt_local,dtdiff)
      
c     if parallel, get best overall dt
#ifdef PARALLEL
      call MPI_Allreduce(dt_local, dt, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_MIN, comm3D, ierr)
      Call ErrorHandler(ierr,ERROR_ALLREDUCE)

c     if not parallel, set output variable
#else   
      dt = dt_local
#endif  
      return
      end
c=================================================================
