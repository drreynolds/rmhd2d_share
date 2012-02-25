C -*- Mode: Fortran; -*-
c----------------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c----------------------------------------------------------------------
c     $Log$
c----------------------------------------------------------------------
c Subroutine to set up the initial conditions for 
c Rayleigh-Taylor instability test problem
c----------------------------------------------------------------------
      subroutine InitialConditions(ux,phi)
c----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
      use iounits
      implicit none

c=======Declarations=========
c
c     input variables
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
c
c     local variables
      integer          :: i,j,k
      double precision:: densRatio
      double precision:: mol1,mol2
      double precision:: psint
      double precision:: const
      double precision:: x,y,z
      integer jg
      double precision:: tmp, erf
      double precision:: atwood
      double precision:: y0
      double precision:: drhodx,drhody
      double precision:: Froude, Mach
      double precision:: xloc,yloc,zloc,pi
      double precision::amp
      integer, parameter:: N=120
      double precision::rk1(N), rk2(N)
c     
      namelist /molwts/ wmol
      namelist /gammas/ gamma
      namelist /diffusiveproperties/ mu, eta, etaFactor, kappa
      namelist /rtinit/densRatio
c

c=======Internals=========
c     
c     input general problem parameters
      open(16,file='prop.inp',form='formatted')
      read(16,molwts)
      read(16,gammas) 
      read(16,diffusiveproperties) 
      close(16)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(A,6(A,es12.5))') '\nPlasma properties:',
     &        '\n       wmol =',wmol,
     &        '\n      gamma =',gamma,
     &        '\n         mu =',mu,
     &        '\n        eta =',eta,
     &        '\n  etaFactor =',etaFactor,
     &        '\n      kappa =',kappa
      end if
      rgas = 1.0d0/wmol
c     
c     input Rayleigh-Taylor parameters
      open(17,file='rt_init.inp',form='formatted')
      read(17,rtinit)
      close(17)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(A,1(A,es12.5))') 
     &        '\nRayleigh-Taylor Parameters',
     &        '\n         densratio =',densRatio
      end if
c     
c     
      mol1=wmol
      mol2=densRatio
c      Mach=0.1D0
      Froude=1.D04
c
      const=1.0
      amp=0.5D0
      atwood=(mol2-mol1)/(mol2+mol1)
      y0=0.5D0*(yl+yr)
c
      open(33,file='rm5.coef',form='formatted',status='old')
      do k=1,N,1
         read(33,*) rk1(k), rk2(k)
         write(6,*)'Rk12=', rk1(k), rk2(k)
      enddo
      close(33)
c     initialize state to constants
      ux(:,:,:,1) = 1.d0
      ux(:,:,:,2) = 0.d0
      ux(:,:,:,3) = 0.d0
      ux(:,:,:,4) = 0.d0
      ux(:,:,:,5) = 0.d0
      ux(:,:,:,6) = 0.d0
      ux(:,:,:,7) = 0.d0
      ux(:,:,:,8) = 1.D0
c     
c     set pi
      pi = 4.D0*datan(1.D0)
c     
c
c     iterate over the domain to set up x-velocity
      do k=izlo,izhi,1
         zloc = ZLloc + (k-0.5d0)*dz

         do j=iylo,iyhi,1
            yloc = YLloc + (j-0.5D0)*dy

            do i=ixlo,ixhi,1
               xloc = XLloc + (i-0.5D0)*dx

               tmp=tanh(const*psint(xloc,yloc,zloc,amp,rk1,rk2))
c               tmp=erf(const*psint(x,y,z)/delta)
               ux(i,j,k,1)=0.5*(mol1+mol2)+0.5*(mol2-mol1)*tmp
c               ux(i,j,k,1)=2.D0*(
c     &              1.D0+atwood*dtanh(psint(x,y,z,amp,rk1,rk2)))
            enddo
         enddo
      enddo
c
      do k=izlo,izhi,1
         zloc = ZLloc + (k-0.5d0)*dz

         do j=iylo,iyhi,1
            yloc = YLloc + (j-0.5D0)*dy

            do i=ixlo,ixhi,1
               xloc = XLloc + (i-0.5D0)*dx
               ux(i,j,k,8)=1.0-0.5/Froude*((mol1+mol2)*(yloc-y0)+
     &              (mol2-mol1)/const*(dlog(dcosh(const*yloc))-
     &              dlog(dcosh(const*y0))))
            enddo
         enddo
      enddo
c
      return
c     
      end
c
c-------------------------------------------------------------------------
      double precision function psint(x,y,z,amp,rk1,rk2)
      double precision:: x,y,z
      double precision:: y0
      double precision:: amp
      double precision:: sum
      double precision:: sigma,scal,cosx,sinx
c      double precision,external:: dcos, dsin, dexp
      integer, parameter:: N=120
      double precision:: rk1(N), rk2(N)
      integer::k
c      psint=y-sum
      sum = 0.D0
      sigma=0.01D0

      do k=1,N,1
      cosx = dcos(k * x)
      sinx = dsin(k * x)
      scal = dexp(-sigma*(k-2)*(k-2))

      sum= sum+scal*((0.5-rk1(k))*cosx + (0.5-rk2(k))*sinx)
      enddo
c
c      write(6,*) x,y,tanh(y-sum)
c      write(12,*) x,sum
      psint=y-0.5*sum
c      psint=y
      return
      end
