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
c     $Log: properties.h,v $
c=================================================================


      module properties
c-----------------------------------------------------------------
c     Description: contains physical constants for the plasma
c
c     Contains:
c             rgas - (dbl) universal gas constant
c            gamma - (dbl) ratio of specific heats
c             wmol - (dbl) molecular weight
c               mu - (dbl) viscosity coefficient
c              eta - (dbl) permittivity constant
c        etaFactor - (dbl) ??
c            kappa - (dbl) thermal conductivity
c-----------------------------------------------------------------
      save

      double precision :: rgas
      double precision :: gamma
      double precision :: wmol
      double precision :: mu
      double precision :: eta
      double precision :: etaFactor
      double precision :: kappa

      end module properties
c=================================================================
