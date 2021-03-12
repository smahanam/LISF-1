!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: set_albedo_attribs
!  \label{set_albedo_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_albedo_attribs(n,source)

  ! !USES:
  use ESMF
  use LDT_albedoMod
  use LDT_coreMod, only  :  LDT_rc, LDT_config
  
  implicit none

  integer,         intent(in) :: n
  character(len=*),intent(in) :: source
  character*20                :: albInterval
  integer                     :: rc
  
! !ARGUMENTS: 

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[source]
!     Albedo dataset source
!   \end{description}
!EOP      
!

   select case( source )

    case( "NCEP_LIS", "NCEP_Native", "CLSMF2.5" )
      LDT_albedo_struc(n)%albedo%num_bins = 1
      LDT_albedo_struc(n)%albedo%num_times = 12

    case( "NCEP_NativeQtr" )
      LDT_albedo_struc(n)%albedo%num_bins = 1
      LDT_albedo_struc(n)%albedo%num_times = 4

    case( "CONSTANT" )
      LDT_albedo_struc(n)%albedo%num_bins = 1
      LDT_albedo_struc(n)%albedo%num_times = 12

   case( "MCD43GFv6-CLSM", "MCD43GF")
      call ESMF_ConfigGetAttribute(LDT_config,albInterval,label = "Albedo climatology interval:", rc=rc)
      
            if(albInterval == "monthly") LDT_albedo_struc(:)%albedo%num_times = 12
            if(albInterval == "8day"   ) LDT_albedo_struc(:)%albedo%num_times = 46
            if(albInterval == "5day"   ) LDT_albedo_struc(:)%albedo%num_times = 73
            if(albInterval == "daily"  ) LDT_albedo_struc(:)%albedo%num_times = 365
            LDT_albedo_struc(n)%albedo%num_bins = 1
      
    case default
      write(*,*) "[ERR] Albedo source not recognized: ",trim(source)
      write(*,*) " Please select: NCEP_LIS, NCEP_Native, NCEP_NativeQtr, " 
      write(*,*) "                CLSMF2.5, or CONSTANT"
      write(*,*) " Program stopping ..."
      stop
!      call LDT_endrun
   end select

end subroutine set_albedo_attribs
