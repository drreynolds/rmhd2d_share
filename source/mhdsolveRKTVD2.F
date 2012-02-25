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
c    $Log: mhdsolverk2.F,v $
c-----------------------------------------------------------------
      subroutine MHDSolve(ux,phi,dt)
c-----------------------------------------------------------------
c     Description: Second order Runge-Kutta time integration 
c                  of the semi-discrete form of the MHD equations
c                  using RK TVD method.
c     Input:  dt -  time step
c             ux - conserved variables at time ttot
c     Output: ux - conserved variables at time ttot+dt
c             phi - auxillary variable used for projection
c-----------------------------------------------------------------
      use mesh
      use properties
#ifdef PARALLEL
      use mpistuff
#endif
c     
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: phi(IXLO:IXHI,IYLO:IYHI)
      double precision:: dt
c     
c     Local variables
c     ux0 = storage for ux at t=ttot
      double precision:: ux0(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: ux1(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: ux2(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     vx - cell centered primitive variables
      double precision:: vx(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     vrp{x,y,z}: Face centered primitive variables for faces {x,y,z}
      double precision:: vrpx(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: vrpy(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: vrpz(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     {f,g,h}inv: Hyperbolic fluxes in the {x,y,z} directions
      double precision:: finv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: ginv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: hinv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     
c     {f,g,h}vis: Parabolic or diffusive fluxes in the {x,y,z} directions
      double precision:: fvis(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: gvis(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: hvis(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c
c     dudx, dbdx: Grad of U and B
      double precision:: dudx(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,3,3)
      double precision:: dbdx(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,3,3)
c
c     divb: Cell centered divergence computed from face centered B
      double precision:: divb(NXLSIZE,NYLSIZE,NZLSIZE)
c     
      integer:: indx,indy,indz
      integer:: i,j,k,l
cc
c     Flag for projection: 0=no 1=yes
        integer:: doProjection
c     Flag for non-conservative source prop to divB: 0=no 1=yes
        integer:: doDivSource
c     Flag for divB Filter: 0=no 1=yes
        integer:: doDivFilter
c     Flag for flux Interpolated CD: 0=no 1=yes
        integer:: doFluxInterpolateCD


c     Exchange ghost values and set BCs
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c
      doProjection=0
      doDivSource=0
      doDivFilter=0
      doFluxInterpolateCD=0
c     
      ux0=ux
      ginv=0.0d0
      hinv=0.0d0
      finv=0.0d0
      gvis=0.0d0
      hvis=0.0d0
      fvis=0.0d0
      dudx=0.0d0
      dbdx=0.0d0
c     
c     
c     First stage in RK2
      call ConservedToPrimitive(ux,vx)
#ifndef NOVISCFLUX
      dudx(:,:,:,1,1:3)=vx(:,:,:,2:4)
      call GradientVector(dudx)
      dbdx(:,:,:,1,1:3)=vx(:,:,:,5:7)
      call GradientVector(dbdx)
#endif
c     
      inlo=ixlo; inhi=ixhi
#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
c      call InviscidFluxRP(finv,vrpx,ux,vx,1,dt)
      call InviscidFluxRP(finv,vx,1,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(finv,vx,1,dt)
#else
c     using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(finv,vx,1,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
      call InviscidFluxTCD(finv,vx,1,dt)
#elif CD4
      call InviscidFluxCD4(finv,vx,1,dt)
#elif ZIP
      call InviscidFluxZip(finv,vx,1,dt)
#else
      call InviscidFluxCD2(finv,vx,1,dt)
#endif
#endif
      inlo=ixlo; inhi=ixhi
      indx=1; indy=0; indz=0;
#ifndef NOVISCFLUX
#ifdef TCD
      call ViscousFluxTCD(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
#elif CD4
      call ViscousFluxCD4(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
#else
      call ViscousFluxCD2(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
#endif
#endif
c     
c     
#ifndef ONE_D
      inlo=iylo; inhi=iyhi
#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
c      call InviscidFluxRP(ginv,vrpy,ux,vx,2,dt)
      call InviscidFluxRP(ginv,vx,2,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(ginv,vx,2,dt)
#else
c     using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(ginv,vx,2,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
      call InviscidFluxTCD(ginv,vx,2,dt)
#elif CD4
      call InviscidFluxCD4(ginv,vx,2,dt)
#elif ZIP
      call InviscidFluxZip(ginv,vx,2,dt)
#else
      call InviscidFluxCD2(ginv,vx,2,dt)
#endif
#endif
      inlo=iylo; inhi=iyhi
      indx=0; indy=1; indz=0;
#ifndef NOVISCFLUX
#ifdef TCD
      call ViscousFluxTCD(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
#elif CD4
      call ViscousFluxCD4(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
#else
      call ViscousFluxCD2(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
#endif
#endif
c     
#ifndef TWO_D
      inlo=izlo; inhi=izhi
#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
c      call InviscidFluxRP(hinv,vrpz,ux,vx,3,dt)
      call InviscidFluxRP(hinv,vx,3,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(hinv,vx,3,dt)
#else
c     using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(hinv,vx,3,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
      call InviscidFluxTCD(hinv,vx,3,dt)
#elif CD4
      call InviscidFluxCD4(hinv,vx,3,dt)
#elif ZIP
      call InviscidFluxZip(hinv,vx,3,dt)
#else
      call InviscidFluxCD2(hinv,vx,3,dt)
#endif
#endif
      inlo=izlo; inhi=izhi
      indx=0; indy=0; indz=1;
#ifndef NOVISCFLUX
#ifdef TCD
      call ViscousFluxTCD(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#elif CD4
      call ViscousFluxCD4(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#else
      call ViscousFluxCD2(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif
#endif
#endif
c     
c     
c     Flux Interpolated Central Difference
      if(doFluxInterpolateCD.eq.1) then
         call fluxInterpolateCD2(finv,ginv,hinv)
         call fluxInterpolateCD2(fvis,gvis,hvis)
      endif
c
c     Update ux using fluxes
      call RK(ux,ux0,finv,ginv,hinv,fvis,gvis,hvis,dt)
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c     
#ifdef DIVBCLEAN
c     Update ux using non-conservative divergence source
      call DivSource(ux,ux0,dt/2.0d0)
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
#endif
      ux1=ux
c     End of first stage of RK2

c     Second stage of the RK scheme
      finv=0.0d0; ginv=0.0d0; hinv=0.0d0; dudx=0.0d0
      fvis=0.0d0; gvis=0.0d0; hvis=0.0d0; dbdx=0.0d0
c
      call ConservedToPrimitive(ux,vx)
c
#ifndef NOVISCFLUX
      dudx(:,:,:,1,1:3)=vx(:,:,:,2:4)
      call GradientVector(dudx)
      dbdx(:,:,:,1,1:3)=vx(:,:,:,5:7)
      call GradientVector(dbdx)
#endif
c     
      divb=0.0D0
      inlo=ixlo; inhi=ixhi
#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
c      call InviscidFluxRP(finv,vrpx,ux,vx,1,dt)
      call InviscidFluxRP(finv,vx,1,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(finv,vx,1,dt)
#else
c     using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(finv,vx,1,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
      call InviscidFluxTCD(finv,vx,1,dt)
#elif CD4
      call InviscidFluxCD4(finv,vx,1,dt)
#elif ZIP
      call InviscidFluxZip(finv,vx,1,dt)
#else
      call InviscidFluxCD2(finv,vx,1,dt)
#endif
#endif
c     Compute divb using x-face normal B
      inlo=ixlo; inhi=ixhi
      indx=1; indy=0; indz=0;
#ifndef NOVISCFLUX
#ifdef TCD
      call ViscousFluxTCD(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
#elif CD4
      call ViscousFluxCD4(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
#else
      call ViscousFluxCD2(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
#endif
#endif
#ifdef DIVBCLEAN
      call BFaceDivergenceF(divb,vrpx,1)
#endif
c     
#ifndef ONE_D
      inlo=iylo; inhi=iyhi
#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
c      call InviscidFluxRP(ginv,vrpy,ux,vx,2,dt)
      call InviscidFluxRP(ginv,vx,2,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(ginv,vx,2,dt)
#else
c     using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(ginv,vx,2,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
      call InviscidFluxTCD(ginv,vx,2,dt)
#elif CD4
      call InviscidFluxCD4(ginv,vx,2,dt)
#elif ZIP
      call InviscidFluxZip(ginv,vx,2,dt)
#else
      call InviscidFluxCD2(ginv,vx,2,dt)
#endif
#endif
c     Augment divb using y-face normal B
#ifdef DIVBCLEAN
      call BfaceDivergenceF(divb,vrpy,2)
#endif
      inlo=iylo; inhi=iyhi
      indx=0; indy=1; indz=0;
#ifndef NOVISCFLUX
#ifdef TCD
      call ViscousFluxTCD(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
#elif CD4
      call ViscousFluxCD4(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
#else
      call ViscousFluxCD2(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
#endif
#endif
c     
#ifndef TWO_D
      inlo=izlo; inhi=izhi
#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
c      call InviscidFluxRP(hinv,vrpz,ux,vx,3,dt)
      call InviscidFluxRP(hinv,vx,3,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(hinv,vx,3,dt)
#else
c     using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(hinv,vx,3,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
      call InviscidFluxTCD(hinv,vx,3,dt)
#elif CD4
      call InviscidFluxCD4(hinv,vx,3,dt)
#elif ZIP
      call InviscidFluxZip(hinv,vx,3,dt)
#else
      call InviscidFluxCD2(hinv,vx,3,dt)
#endif
#endif
c     Augment divb using y-face normal B
      inlo=ixlo; inhi=ixhi
      indx=0; indy=0; indz=1;
#ifndef NOVISCFLUX
#ifdef TCD
      call ViscousFluxTCD(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#elif CD4
      call ViscousFluxCD4(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#else
      call ViscousFluxCD2(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif
#ifdef DIVBCLEAN
      call BfaceDivergenceF(divb,vrpz,3)
#endif
#endif
#endif
c     write(6,*) 'DIVB MAX/MIN= ',maxval(divb),minval(divb)
      
c     Projection of B
#ifdef DIVBCLEAN
      call Projection(phi,divb)
      call mesh_update_bdry_async(phi,1)
c     Correct fluxes based on projected face normal B
      call CorrectFlux(finv,phi,vrpx,1)
      call CorrectFlux(ginv,phi,vrpy,2)
c      call CorrectFlux(hinv,phi,vrpz,3)
#endif
      
c     Flux Interpolated Central Difference
      if(doFluxInterpolateCD.eq.1) then
         call fluxInterpolateCD2(finv,ginv,hinv)
         call fluxInterpolateCD2(fvis,gvis,hvis)
      endif

c     Update ux based on corrected fluxes
      Call RK(ux,ux1,finv,ginv,hinv,fvis,gvis,hvis,dt)
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
      
      ux2=ux
      ux=0.5D0*(ux0+ux2)
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c
c     Adjust ux based on divB update
#ifdef DIVBCLEAN
      call DivSource(ux,ux0,dt)
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
      
cccccccTURNED OFF DIVBDIFFUSIONFILTER PER RAVI'S REQUEST, 12/3/03 
cc     Diffuse divergence of B based on cell-centered B
c     call DivBDiffusionFilter(ux)
c     call mesh_update_bdry_async(ux,nvar)
c     call SetBoundaryValues(ux)
#endif
c
      
      return
      end
c     
#ifdef TWO_D
c----------------------------------------------------------------------
      subroutine RK(ux,ux0,finv,ginv,hinv,fvis,gvis,hvis,dt)
      use mesh
      use mesh_common
      double precision:: dt
      double precision:: ux(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     ux0(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     finv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     ginv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     hinv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     fvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     gvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     hvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      integer:: i,j,k,l
      double precision:: dtodx,dtody
c     
C     $DOACROSS LOCAL(i,j,k)
c
      dtodx=dt/dx
      dtody=dt/dy
c     
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1,NXlocal,1
c
c     Density
               ux(i,j,k,1)=ux0(i,j,k,1)
     &              -dtodx*(
     &              finv(i,j,k,1)
     &              -finv(i-1,j,k,1))
     &              -dtody*(ginv(i,j,k,1)  
     &              -ginv(i,j-1,k,1))
c     momentum
               ux(i,j,k,2)=ux0(i,j,k,2)
     &              -dtodx*(
     &              finv(i,j,k,2)  
     &              -finv(i-1,j,k,2))
     &              -dtody*(
     &              ginv(i,j,k,4)  
     &              -ginv(i,j-1,k,4))
     &              +dtodx*(
     &              fvis(i,j,k,2)  
     &              -fvis(i-1,j,k,2))
     &              +dtody*(
     &              gvis(i,j,k,4)  
     &              -gvis(i,j-1,k,4))
c     
               ux(i,j,k,3)=ux0(i,j,k,3)
     &              -dtodx*(
     &              finv(i,j,k,3)  
     &              -finv(i-1,j,k,3))
     &              -dtody*(
     &              ginv(i,j,k,2)  
     &              -ginv(i,j-1,k,2))
     &              +dtodx*(
     &              fvis(i,j,k,3)  
     &              -fvis(i-1,j,k,3))
     &              +dtody*(
     &              gvis(i,j,k,2)  
     &              -gvis(i,j-1,k,2))
c
#ifdef TWO_HALF_D
               ux(i,j,k,4)=ux0(i,j,k,4)
     &              -dtodx*(
     &              finv(i,j,k,4)  
     &              -finv(i-1,j,k,4))
     &              -dtody*(
     &              ginv(i,j,k,3)  
     &              -ginv(i,j-1,k,3))
     &              +dtodx*(
     &              fvis(i,j,k,4)  
     &              -fvis(i-1,j,k,4))
     &              +dtody*(
     &              gvis(i,j,k,3)  
     &              -gvis(i,j-1,k,3))
#endif
c     X Magnetic Field
               ux(i,j,k,5)=ux0(i,j,k,5)
     &              -dtodx*(
     &              finv(i,j,k,5)  
     &              -finv(i-1,j,k,5))
     &              -dtody*(
     &              ginv(i,j,k,7)  
     &              -ginv(i,j-1,k,7))
     &              +dtodx*(
     &              fvis(i,j,k,5)  
     &              -fvis(i-1,j,k,5))
     &              +dtody*(
     &              gvis(i,j,k,7)  
     &              -gvis(i,j-1,k,7))
c     
c     Y Magnetic Field
               ux(i,j,k,6)=ux0(i,j,k,6)
     &              -dtodx*(
     &              finv(i,j,k,6)  
     &              -finv(i-1,j,k,6))
     &              -dtody*(
     &              ginv(i,j,k,5)  
     &              -ginv(i,j-1,k,5))
     &              +dtodx*(
     &              fvis(i,j,k,6)  
     &              -fvis(i-1,j,k,6))
     &              +dtody*(
     &              gvis(i,j,k,5)  
     &              -gvis(i,j-1,k,5))
c     
#ifdef TWO_HALF_D
c     Z Magnetic Field
               ux(i,j,k,7)=ux0(i,j,k,7)
     &              -dtodx*(
     &              finv(i,j,k,7)  
     &              -finv(i-1,j,k,7))
     &              -dtody*(
     &              ginv(i,j,k,6)  
     &              -ginv(i,j-1,k,6))
     &              +dtodx*(
     &              fvis(i,j,k,7)  
     &              -fvis(i-1,j,k,7))
     &              +dtody*(
     &              gvis(i,j,k,6)  
     &              -gvis(i,j-1,k,6))
#endif
c     
c     Energy
               ux(i,j,k,8)=ux0(i,j,k,8)
     &              -dtodx*(
     &              finv(i,j,k,8)  
     &              -finv(i-1,j,k,8))
     &              -dtody*(
     &              ginv(i,j,k,8)  
     &              -ginv(i,j-1,k,8))
     &              +dtodx*(
     &              fvis(i,j,k,8)  
     &              -fvis(i-1,j,k,8))
     &              +dtody*(
     &              gvis(i,j,k,8)  
     &              -gvis(i,j-1,k,8))
c
            enddo
         enddo
      enddo
c     
c
#ifdef RTSOURCE
c     Rayleigh-Taylor sources due to gravity (Froude=1e4)
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1,NXlocal,1
               ux(i,j,k,3) = ux(i,j,k,3) - dt*ux(i,j,k,1)*1.D-4
               ux(i,j,k,8) = ux(i,j,k,8) - dt*ux(i,j,k,3)*1.D-4
            enddo
         enddo
      enddo
#endif
c
c
      return
      end
c     
c-----------------------------------------------------------------------
#else 
c----------------------------------------------------------------------
      subroutine RK(ux,ux0,finv,ginv,hinv,fvis,gvis,hvis,dt)
      use mesh
      use mesh_common
      double precision:: dt
      double precision:: ux(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     ux0(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     finv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     ginv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     hinv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
       double precision,intent(in):: 
      &     fvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
       double precision,intent(in):: 
      &     gvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
       double precision,intent(in):: 
      &     hvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      integer:: i,j,k
      double precision:: dtodx,dtody,dtodz
c     
C     $DOACROSS LOCAL(i,j,k)
      dtodx=dt/dx
      dtody=dt/dy
      dtodz=dt/dz
c     
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1,NXlocal,1
               ux(i,j,k,1)=ux0(i,j,k,1)
     &              -dtodx*(
     &              finv(i,j,k,1)
     &              -finv(i-1,j,k,1))
     &              -dtody*(ginv(i,j,k,1)  
     &              -ginv(i,j-1,k,1))
     &              -dtodz*(hinv(i,j,k,1)  
     &              -hinv(i,j,k-1,1))
c     momentum
               ux(i,j,k,2)=ux0(i,j,k,2)
     &              -dtodx*(
     &              finv(i,j,k,2)  
     &              -finv(i-1,j,k,2))
     &              -dtody*(
     &              ginv(i,j,k,4)  
     &              -ginv(i,j-1,k,4))
     &              -dtodz*(hinv(i,j,k,3)  
     &              -hinv(i,j,k-1,3))
c     
               ux(i,j,k,3)=ux0(i,j,k,3)
     &              -dtodx*(
     &              finv(i,j,k,3)  
     &              -finv(i-1,j,k,3))
     &              -dtody*(
     &              ginv(i,j,k,2)  
     &              -ginv(i,j-1,k,2))
     &              -dtodz*(hinv(i,j,k,4)  
     &              -hinv(i,j,k-1,4))
c     
               ux(i,j,k,4)=ux0(i,j,k,4)
     &              -dtodx*(
     &              finv(i,j,k,4)  
     &              -finv(i-1,j,k,4))
     &              -dtody*(
     &              ginv(i,j,k,3)  
     &              -ginv(i,j-1,k,3))
     &              -dtodz*(hinv(i,j,k,2)  
     &              -hinv(i,j,k-1,2))
c     
c     X Magnetic Field
               ux(i,j,k,5)=ux0(i,j,k,5)
     &              -dtodx*(
     &              finv(i,j,k,5)  
     &              -finv(i-1,j,k,5))
     &              -dtody*(
     &              ginv(i,j,k,7)  
     &              -ginv(i,j-1,k,7))
     &              -dtodz*(hinv(i,j,k,6)  
     &              -hinv(i,j,k-1,6))
c     
c     Y Magnetic Field
               ux(i,j,k,6)=ux0(i,j,k,6)
     &              -dtodx*(
     &              finv(i,j,k,6)  
     &              -finv(i-1,j,k,6))
     &              -dtody*(
     &              ginv(i,j,k,5)  
     &              -ginv(i,j-1,k,5))
     &              -dtodz*(hinv(i,j,k,7)  
     &              -hinv(i,j,k-1,7))
c     
c     Z Magnetic Field
               ux(i,j,k,7)=ux0(i,j,k,7)
     &              -dtodx*(
     &              finv(i,j,k,7)  
     &              -finv(i-1,j,k,7))
     &              -dtody*(
     &              ginv(i,j,k,6)  
     &              -ginv(i,j-1,k,6))
     &              -dtodz*(hinv(i,j,k,5)  
     &              -hinv(i,j,k-1,5))
c     
c     Energy
               ux(i,j,k,8)=ux0(i,j,k,8)
     &              -dtodx*(
     &              finv(i,j,k,8)  
     &              -finv(i-1,j,k,8))
     &              -dtody*(
     &              ginv(i,j,k,8)  
     &              -ginv(i,j-1,k,8))
     &              -dtodz*(hinv(i,j,k,8)  
     &              -hinv(i,j,k-1,8))
c     
            enddo
         enddo
      enddo
c     
c
#ifdef RTSOURCE
c     Rayleigh-Taylor sources due to gravity (Froude=1e4)
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1,NXlocal,1
               ux(i,j,k,3) = ux(i,j,k,3) - dt*ux(i,j,k,1)*1.D-4
               ux(i,j,k,8) = ux(i,j,k,8) - dt*ux(i,j,k,3)*1.D-4
            enddo
         enddo
      enddo
#endif
c
c
      return
      end
c     
c-----------------------------------------------------------------------
#endif
c-----------------------------------------------------------------------
      
c----------------------------------------------------------------------
      subroutine fluxInterpolateCD2(finv,ginv,hinv)
      use mesh
      use mesh_common
      double precision:: 
     &     finv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: 
     &     ginv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: 
     &     hinv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      integer:: i,j,k,l 
      double precision:: 
     &     fmod(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI)
      double precision:: 
     &     gmod(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI)
c     
      fmod=finv(:,:,:,6)
      gmod=ginv(:,:,:,7)
c     
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal-1,1
               fmod(i,j,k)=0.125D0*(
     &              2.D0*finv(i,j,k,6)
     &              +finv(i+1,j,k,6)
     &              +finv(i-1,j,k,6)
     &              -ginv(i+1,j-1,k,7)
     &              -ginv(i+1,j,k,7)
     &              -ginv(i,j-1,k,7)
     &              -ginv(i,j,k,7))
               
            enddo
         enddo
      enddo
c     
      do k=1,NZlocal,1
         do j=1,NYlocal-1,1
            do i=1, NXlocal,1
               gmod(i,j,k)=0.125D0*(
     &              2.D0*ginv(i,j,k,7)
     &              +ginv(i,j+1,k,7)
     &              +ginv(i,j-1,k,7)
     &              -finv(i-1,j,k,6)
     &              -finv(i-1,j+1,k,6)
     &              -finv(i,j,k,6)
     &              -finv(i,j+1,k,6))
            enddo
         enddo
      enddo
c     
      finv(:,:,:,6)=fmod(:,:,:)
      ginv(:,:,:,7)=gmod(:,:,:)
      return
      end
c----------------------------------------------------------------------
