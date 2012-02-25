C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c-----------------------------------------------------------------------
c     $Log: entsource.F,v $
c-----------------------------------------------------------------------
c
      subroutine EntropySource(esrc,vx,dudx,dbdx)
c-----------------------------------------------------------------------
c     Description:  Calculation of entropy source
c     Input: vx - primitive variables
c            direction - 1,2,3 (=x,y,z)
c            dt - time step (not used)
c            dudx - grad(velocity)
c            dbdx - grad(magnetic field)
c            dudx and dbdx are precomputed.
c            indx,indy,indz - specify which direction (sort of 
c            redundant to specify both direction and these)
c     Output: esrc - entropy source
c     Issues: Works for constant eta 
c-----------------------------------------------------------------------
      use mesh
      use properties
c     
c 
      double precision:: esrc(nxlsize,nylsize,nzlsize)
      double precision:: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: dudx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,3,3)
      double precision:: dbdx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,3,3)
c
#ifdef ENTROPY_FORMULATION
c
      double precision:: Jcur(3)
      double precision:: temp(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision:: rho,entropy,kappaLapT
      integer:: i,j,k
c     
      double precision:: dxi,dyi,dzi
      double precision:: d0,d1,d2
      integer :: l

c     This subroutine determines the vector of diffusive fluxes.
c     v: Vector of variables - rho, u,v,w,Bx,By,Bz,p
c     rho=density, u=x-velocity, v=y-velocity
c     w=z-velocity
c     Bx,By,Bz=x,y,z components of the magnetic field
c     p=thermal pressure
c
      dxi=1.D0/dx
      dyi=1.D0/dy
      dzi=1.D0/dz
c     Fourth order central second derivative coefficients
      d0=-2.5D0
      d1=4.D0/3.D0
      d2=-1.D0/12.D0
c
c     
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
c     T=p/rho/rgas = rho^(gamma-1)*exp(s)/rgas
              temp(i,j,k)=vx(i,j,k,1)**(gamma-1.D0)
     &              *dexp(vx(i,j,k,8))/rgas
           enddo
         enddo
      enddo
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               rho=vx(i,j,k,1)
               entropy=vx(i,j,k,8)
               Jcur(1)=dbdx(i,j,k,3,2)-dbdx(i,j,k,2,3)
               Jcur(2)=dbdx(i,j,k,1,3)-dbdx(i,j,k,3,1)
               Jcur(3)=dbdx(i,j,k,1,2)-dbdx(i,j,k,2,1)
               kappaLapT=d2*(temp(i-2,j,k)+temp(i+2,j,k)
     &              +temp(i,j-2,k)+temp(i,j+2,k))
     &              +d1*(temp(i-1,j,k)+temp(i+1,j,k)
     &              +temp(i,j-1,k)+temp(i,j+1,k))
     &              +2.0D0*d0*temp(i,j,k)
#ifndef TWO_D
               kappaLapT=kappaLapT+
     &                  d2*(temp(i,j,k+2)+temp(i,j,k-2))
     &                 +d1*(temp(i,j,k-1)+temp(i+1,j,k-1))
     &                 +d0*temp(i,j,k)
#endif
               kappaLapT=kappa*kappaLapT
               esrc(i,j,k)=(gamma-1.D0)/rho**(gamma-1.D0)*
     &                     dexp(-entropy)*(
     &                     mu*sum(dudx(i,j,k,:,:)*dudx(i,j,k,:,:))+
     &                     eta*sum(Jcur(:)*Jcur(:))+
     &                     kappaLapT)
            enddo
         enddo
      enddo
c
#endif
      return
      end subroutine EntropySource
c
