!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

!BOP
! 
! !ROUTINE: noah33_getirrigationstates
! \label{noah33_getirrigationstates}
! 
! !INTERFACE:

subroutine noah33_getirrigationstates(n,irrigState)  

  ! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_TKFRZ, LIS_CONST_LATVAP
  use noah33_lsmMod

  implicit none 
    
  ! !DESCRIPTION:        
  !
  ! Calculate water requirement and apply the amount to precipitation.
  !
  ! Irrigate when root zone soil moisture falls below 50 % of 
  ! the field capacity (reference soil moiture) at 6 am LST.  
  ! The root zone is actual maximum root depth rather than NOAH root zone.
  ! Method of irrigation is by precipitation between 6-10 am LST.
  !
  ! Irrigation amount is scaled to grid total crop fraction when intensity
  ! is less than the fraction.  Irrigation is expanded to non-crop, non-forest,
  ! non-baresoil/urban tiles if intensity exceeds grid total crop fraction.
  ! In latter case, scaled irrigation is applied to grassland first, 
  ! then further applied over the rest of tiles equally if the intensity 
  ! exceeds grassland fraction as well.   
  !
  ! Optionally efficiency correction is applied to account for field loss. 
  !
  ! Optionally outputs amount of water put into the system to a text file. 
  !
  ! This version includes modifications to irr4 as follows:
  ! 1) Use location specific growing season threshold (40% of GFRAC range)
  ! 2) Allow irrigation in non-crop/non-forest tiles when irrigation 
  !    intensity exceeds total crop fraction
  !
  ! REVISION HISTORY:
  !
  ! Aug 2008: Hiroko Kato; Initial code
  ! Nov 2012: Sujay Kumar, Incorporated into LIS
  ! Jun 2014: Ben Zaitchik; Added flood scheme
  !EOP
  
  ! Sprinkler parameters
  real, parameter      :: otimess = 6.0 ! local trigger check start time [hour]
  real, parameter      :: irrhrs = 4.   ! duration of irrigation hours 
  ! Drip parameters 
  real, parameter      :: otimeds = 6.0 ! local trigger check start time [hour]
  real, parameter      :: irrhrd = 12.0   ! duration of irrigation hours 
  ! Flood parameters
  real, parameter      :: otimefs = 6.0 ! local trigger check start time [hour]
  real, parameter      :: irrhrf = 1.0   ! duration of irrigation hours 
!!!real, parameter      :: ffreq = 0.0 ! frequency of flood irrig [days] set to 0.0 to use thresh instead
  
  real, parameter      :: efcor = 0.0      ! Efficiency Correction (%)
  integer, parameter   :: nsoil = 4

  integer              :: n
  type(ESMF_State)     :: irrigState
  
  if(LIS_rc%irrigation_type.eq."Concurrent") then
     CALL UPDATE_IRRIG_STATE_CONCURRENT (n,irrigState)
  else
     CALL UPDATE_IRRIG_STATE_SINGLE (n,irrigState)
  endif
contains 
      
  ! _______________________________________________________
  
  SUBROUTINE UPDATE_IRRIG_STATE_SINGLE (n,irrigState)

    implicit none
    
    integer              :: n
    integer              :: rc
    integer              :: t,k,gid,vegt,l
    type(ESMF_State)     :: irrigState
    type(ESMF_Field)     :: irrigRateField,irrigFracField
    type(ESMF_Field)     :: irrigRootDepthField,irrigScaleField    
    real,  pointer       :: irrigRate(:), irrigFrac(:)
    real,  pointer       :: irrigRootDepth(:), irrigScale(:)
    integer              :: chhr, lhr
    real                 :: asmc, tsmcwlt, tsmcref, ma, otimes, otimee, irrhr
    real                 :: sldpth(nsoil)
    real                 :: rdpth(nsoil)
    real                 :: zdpth(nsoil)
    real                 :: water(nsoil)
    real                 :: twater, twater1, twater2
    real                 :: ippix, crootd
    real                 :: smcmax, smcref, smcwlt, psisat, dksat
    real                 :: smcref1, smcwlt1,shdfac
    real                 :: bexp, smhigh, smlow
    integer              :: lroot,veg_index1,veg_index2
    real                 :: gsthresh, ltime
    logical              :: irrig_check_frozen_soil
    real                 :: timestep, shift_otimes, shift_otimee
    
    call ESMF_StateGet(irrigState, "Irrigation rate",irrigRateField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation rate')    
    call ESMF_FieldGet(irrigRateField, localDE=0,farrayPtr=irrigRate,rc=rc)
    
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation rate')
    
    call ESMF_StateGet(irrigState, "Irrigation frac",&
         irrigFracField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation frac')    
    call ESMF_FieldGet(irrigFracField, localDE=0,&
         farrayPtr=irrigFrac,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation frac')
    
    call ESMF_StateGet(irrigState, "Irrigation max root depth",&
         irrigRootDepthField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation max root depth')    
    call ESMF_FieldGet(irrigRootDepthField, localDE=0,&
         farrayPtr=irrigRootDepth,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation root depth')
    
    call ESMF_StateGet(irrigState, "Irrigation scale",&
         irrigScaleField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation scale')    
    call ESMF_FieldGet(irrigScaleField, localDE=0,&
         farrayPtr=irrigScale,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation scale')  
    
    irrigRate = 0.0  
    
    !----------------------------------------------------------------------
    ! Set start and end times for selected irrigation type
    !----------------------------------------------------------------------
    if(LIS_rc%irrigation_type.eq."Sprinkler") then
       otimes = otimess
       irrhr = irrhrs
       otimee = otimess + irrhrs
    elseif(LIS_rc%irrigation_type.eq."Drip") then
       otimes = otimeds
       irrhr = irrhrd
       otimee = otimeds + irrhrd
    elseif(LIS_rc%irrigation_type.eq."Flood") then
       otimes = otimefs
       irrhr = irrhrf
       otimee = otimefs + irrhrf
    endif
        
    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

       timestep = LIS_rc%ts
       
       ! Adjust bounds by timestep to account for the fact that LIS_rc%hr, etc.
       ! will represents the END of the integration timestep window
       
       shift_otimes = otimes + (timestep/3600.)
       shift_otimee = otimee + (timestep/3600.)
       
       irrig_check_frozen_soil = .false.
       
       if((noah33_struc(n)%noah(t)%smc(1) - &
            noah33_struc(n)%noah(t)%sh2o(1)).gt.0001) then 
          irrig_check_frozen_soil = .true. 
       elseif((noah33_struc(n)%noah(t)%smc(2) - &
            noah33_struc(n)%noah(t)%sh2o(2)).gt.0001) then 
          irrig_check_frozen_soil = .true. 
       elseif((noah33_struc(n)%noah(t)%smc(3) - &
            noah33_struc(n)%noah(t)%sh2o(3)).gt.0001) then 
          irrig_check_frozen_soil = .true. 
       elseif((noah33_struc(n)%noah(t)%smc(4) - &
            noah33_struc(n)%noah(t)%sh2o(4)).gt.0001) then 
          irrig_check_frozen_soil = .true. 
       elseif(noah33_struc(n)%noah(t)%stc(2).le.LIS_CONST_TKFRZ) then
          irrig_check_frozen_soil = .true. 
       elseif(noah33_struc(n)%noah(t)%stc(3).le.LIS_CONST_TKFRZ) then
          irrig_check_frozen_soil = .true. 
       elseif(noah33_struc(n)%noah(t)%stc(4).le.LIS_CONST_TKFRZ) then
          irrig_check_frozen_soil = .true. 
       endif
       
       if(.not.irrig_check_frozen_soil) then 
          twater  = 0.0
          water   = 0.0
          asmc    = 0.0
          tsmcwlt = 0.0
          tsmcref = 0.0
          ma      = 0.0
          crootd  = 0.0
          lroot   = 0
          
          sldpth(1) = 0.1         ! Soil layer thicknesses (m)
          sldpth(2) = 0.3
          sldpth(3) = 0.6
          sldpth(4) = 1.0
          zdpth(1) = sldpth(1)         ! Soil layer depth from top (m)
          zdpth(2) = sldpth(1) + sldpth(2)
          zdpth(3) = sldpth(1) + sldpth(2) + sldpth(3)
          zdpth(4) = sldpth(1) + sldpth(2) + sldpth(3) + sldpth(4)
          
          smcmax =  noah33_struc(n)%noah(t)%smcmax
          psisat =  noah33_struc(n)%noah(t)%psisat
          dksat  =  noah33_struc(n)%noah(t)%dksat
          bexp   =  noah33_struc(n)%noah(t)%bexp
          smlow = 0.5
          smhigh = 6.0
          smcref1 = smcmax*(5.79e-9/dksat)**(1.0/(2.0*bexp+3.0))
          smcref = smcref1 + (smcmax-smcref1) / smhigh
          smcwlt1 = smcmax * (200.0/psisat)**(-1.0/bexp)
          smcwlt = smcwlt1 - smlow * smcwlt1
          
          gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
          chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
          if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
               (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
               chhr = chhr -1
          lhr = LIS_rc%hr +chhr
          if(lhr.ge.24) lhr = lhr-24
          if(lhr.lt.0) lhr = lhr+24
          
          ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0
          
          shdfac =  noah33_struc(n)%noah(t)%shdfac
          
          ! Calculate vegetation and root depth parameters

          ! If we are outside of the irrigation window, set rate to 0
          if ((ltime.ge.shift_otimee).or.(ltime.lt.shift_otimes)) then
             irrigRate(t) = 0.0
          endif
          
          if((ltime.ge.shift_otimes).and.(ltime.lt.shift_otimee)) then 
             vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
             !----------------------------------------------------------------------       
             !    Proceed if it is non-forest, non-baresoil, non-urban
             !----------------------------------------------------------------------       
             if(LIS_rc%lcscheme.eq."UMD") then !UMD
                veg_index1 = 6
                veg_index2 = 11
             elseif(LIS_rc%lcscheme.eq."UMD+MIRCAIrrig") then !UMD+MIRCAIrrig (Temporary, KRA)
                veg_index1 = 6
                veg_index2 = 16
!             elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP") then 
                ! TEMPORARILY ADDED HERE (HKB) for single crop tile option
                ! once multiple crop tile is unabled, +MIRCA should be separated
                ! into a case where veg_index1 and veg_index2 are set for croptypes
             elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP".or.&
                  LIS_rc%lcscheme.eq."IGBPNCEP+MIRCA") then  
                veg_index1 = 6
                veg_index2 = 14
             elseif(LIS_rc%lcscheme.eq."USGS") then !UMD
                veg_index1 = 2
                veg_index2 = 10
             else
                write(LIS_logunit,*) '[ERR] The landcover scheme ',trim(LIS_rc%lcscheme)
                write(LIS_logunit,*) '[ERR] is not supported for the Noah.3.3 irrigation module.'
                call LIS_endrun()
             endif
             
             if(vegt.ge.veg_index1.and.vegt.le.veg_index2&
                  .and.vegt.ne.LIS_rc%bareclass.and.&
                  vegt.ne.LIS_rc%urbanclass) then 
                if(irrigFrac(t).gt.0) then 
                   ippix = irrigFrac(t)*0.01
                   
                   ! Determine the amount of irrigation to apply if irrigated tile
                   if( IrrigScale(t).gt.0.0 ) then ! irrigated tile
                      !                if(ippix.gt.0.0) then  ! irrigated tile
                      
                      gsthresh = noah33_struc(n)%noah(t)%shdmin + & 
                           0.40 * (noah33_struc(n)%noah(t)%shdmax - &
                           noah33_struc(n)%noah(t)%shdmin)
                      
                      if(shdfac .ge. gsthresh) then 
                         crootd = irrigRootdepth(t)*shdfac
                         if(crootd.gt.0.and.crootd.lt.zdpth(1)) then 
                            lroot = 1
                            rdpth(1) = crootd
                         elseif(crootd .ge. zdpth(1).and.crootd .lt. zdpth(2) ) then
                            lroot = 2
                            rdpth(1) = sldpth(1)
                            rdpth(2) = crootd - zdpth(1)
                         elseif ( crootd.ge.zdpth(2).and.crootd .lt. zdpth(3) ) then
                            lroot = 3
                            rdpth(1) = sldpth(1)
                            rdpth(2) = sldpth(2)
                            rdpth(3) = crootd - zdpth(2)
                         elseif ( crootd.ge.zdpth(3).and.crootd .lt. zdpth(4) ) then
                            lroot = 4
                            rdpth(1) = sldpth(1)
                            rdpth(2) = sldpth(2)
                            rdpth(3) = sldpth(3)
                            rdpth(4) = crootd - zdpth(3)
                            !                      else
                            !                         print*,'error getting root depth'
                            !                         stop
                         endif
                         
                         !!!!! SPRINKLER IRRIGATION
                         if(LIS_rc%irrigation_type.eq."Sprinkler") then
                            !----------------------------------------------------------------------       
                            !    Set the irrigation rate at start time; keep the value till next day
                            !    If local time at the tile fall in the irrigation check
                            !    hour then check the root zone average soil moisture
                            !----------------------------------------------------------------------       
                            if(ltime.eq.shift_otimes) then 
                               irrigRate(t) = 0.0
                               !-------------------------------------------------------------
                               !     Compute the root zone accumlative soil moisture [mm], 
                               !     field capacity [mm], and wilting point [mm] 
                               !-------------------------------------------------------------
                               if(lroot.gt.0) then

                                  call rz_soil_moisture_conditions  &
                                       (lroot, noah33_struc(n)%noah(t)%smc(1:lroot),     &
                                       rdpth(1:lroot), smcwlt, smcref, asmc, tsmcwlt, tsmcref)

                                  !---------------------------------------------------------------
                                  !     Get the root zone moisture availability to the plant
                                  !--------------------------------------------------------------- 
                                  ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)
                                  if(ma.le.LIS_rc%irrigation_thresh) then 
                                     !-----------------------------------------------------------------------------
                                     !     Compute irrigation rate
                                     !-----------------------------------------------------------------------------
                                     irrigRate(t) = sprinkler_irrig_water (lroot, noah33_struc(n)%noah(t)%smc(1:lroot), &
                                          rdpth(1:lroot), smcref, irrigScale(t)) /(irrhr*3600.0)
                                  endif
                               endif
                            endif
                !!!!! DRIP IRRIGATION
                         elseif(LIS_rc%irrigation_type.eq."Drip") then

                            ! Need to get crop coefficient so that we can caculate unstressed Transp
                            !       RC=RSMIN/(XLAI*RCS*RCT*RCQ)
                            !       PCIRR=(RR+DELTA)/(RR*(1.+RC*CH)+DELTA)
                            ! CALL TRANSP (with PCIRR)
                            ! 
                            ! SM : 
                            ! (1) added below 2 variables to noah33dec structure in noah33_module.F90
                            !     real :: PC    ! PLANT COEFFICIENT (PC: UNITLESS FRACTION, 0-1) PC  WHERE PC*ETP = ACTUAL TRANSP
                            !     real :: PCIRR ! PLANT COEFFICIENT (PC: UNITLESS FRACTION, 0-1) with no soil moisutre stress-i.e., perfect irrigation.
                            ! (2) total transpriation with no soil moisture stress (perfect irrigation) from SUBROUTINE TRANSP would be
                            !     noah33_struc(n)%noah(t)%tveg * noah33_struc(n)%noah(t)%PCIRR / noah33_struc(n)%noah(t)%PC [with appropriate unit conversion]  
                            !     and converted TVEG units [W / m2] to irrigrate units [kg / m /s]
       
                            twater = noah33_struc(n)%noah(t)%tveg * noah33_struc(n)%noah(t)%PCIRR / noah33_struc(n)%noah(t)%PC /LIS_CONST_LATVAP

                            !-----------------------------------------------------------------------------
                            !     Apply efficiency correction
                            !-----------------------------------------------------------------------------
                            twater2 = twater
                            twater = twater*(100.0/(100.0-efcor))
                            !-----------------------------------------------------------------------------
                            !     Compute irrigation rate
                            !-----------------------------------------------------------------------------
                            irrigRate(t) = twater  ! for drip calculation, twater is a rate [kg m-2/s]
                            noah33_struc(n)%noah(t)%smc(1) = &
                                 noah33_struc(n)%noah(t)%smc(1) + (twater-twater2)/(sldpth(1)*1000.0) !! check this with Sujay
                            
                !!!!! FLOOD IRRIGATION
                         elseif(LIS_rc%irrigation_type.eq."Flood") then
                            !-------------------------------------------------------------
                            !     Compute the root zone accumlative soil moisture [mm], 
                            !     field capacity [mm], and wilting point [mm] 
                            !-------------------------------------------------------------
                            if(lroot.gt.0) then 

                               call rz_soil_moisture_conditions  &
                                    (lroot, noah33_struc(n)%noah(t)%smc(1:lroot),     &
                                    rdpth(1:lroot), smcwlt, smcref, asmc, tsmcwlt, tsmcref)

                               !---------------------------------------------------------------
                               !     Get the root zone moisture availability to the plant
                               !--------------------------------------------------------------- 
                               !                            ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)   ! Original
                               ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)/IrrigScale(t) ! BZ added IrrigScale
                               
                               if( ma .le. LIS_rc%irrigation_thresh ) then

                                  !-----------------------------------------------------------------------------
                                  !     Compute irrigation rate
                                  !-----------------------------------------------------------------------------

                                  irrigRate(t) = flood_irrig_water (LIS_rc%irrigation_mxsoildpth, &
                                       noah33_struc(n)%noah(t)%smc(1:LIS_rc%irrigation_mxsoildpth), &
                                       sldpth(1:LIS_rc%irrigation_mxsoildpth), smcmax, irrigScale(t))/LIS_rc%ts
                                  
                                  ! BZ modification 4/2/2015 to account for ippix and all soil layers:
                                  do l = 1, LIS_rc%irrigation_mxsoildpth
                                     noah33_struc(n)%noah(t)%smc(l) = IrrigScale(t)*smcmax + &
                                          (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(l)
                                  end do
                                  !                              noah33_struc(n)%noah(t)%smc(1) = IrrigScale(t)*smcmax + &
                                  !                                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(1)
                                  !                              noah33_struc(n)%noah(t)%smc(2) = IrrigScale(t)*smcmax + &
                                  !                                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(2)
                                  !                              noah33_struc(n)%noah(t)%smc(3) = IrrigScale(t)*smcmax + &
                                  !                                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(3)
                                  !                              noah33_struc(n)%noah(t)%smc(4) = IrrigScale(t)*smcmax + &
                                  !                                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(4)
                               endif
                            endif                            
                         endif                         
                      endif
                   end if
                end if
             end if
          end if
       endif
    enddo
    
  END SUBROUTINE UPDATE_IRRIG_STATE_SINGLE

  ! ---------------------------------------------------------------------------------------------------
  
  SUBROUTINE UPDATE_IRRIG_STATE_CONCURRENT (n,irrigState)

    implicit none
    
    integer              :: n
    integer              :: rc
    integer              :: t,k,gid,vegt,l, iType
    type(ESMF_State)     :: irrigState
    type(ESMF_Field)     :: irrigRateField,irrigFracField
    type(ESMF_Field)     :: irrigRootDepthField,irrigScaleField
    
    real,  pointer       :: irrigRate(:), irrigFrac(:,:)
    real,  pointer       :: irrigRootDepth(:), irrigScale(:,:)
    integer              :: chhr, lhr
    real                 :: asmc, tsmcwlt, tsmcref, ma, otimes, otimee, irrhr
    real                 :: sldpth(nsoil)
    real                 :: rdpth(nsoil)
    real                 :: zdpth(nsoil)
    real                 :: water(nsoil)
    real                 :: twater, twater1, twater2
    real                 :: ippix, crootd
    real                 :: smcmax, smcref, smcwlt, psisat, dksat
    real                 :: smcref1, smcwlt1,shdfac
    real                 :: bexp, smhigh, smlow
    integer              :: lroot,veg_index1,veg_index2
    real                 :: gsthresh, ltime
    logical              :: irrig_check_frozen_soil
    real, dimension (3)  :: irrigRate_frac
    real                 :: timestep, shift_otimes, shift_otimee

    call ESMF_StateGet(irrigState, "Irrigation rate",irrigRateField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation rate')    
    call ESMF_FieldGet(irrigRateField, localDE=0,farrayPtr=irrigRate,rc=rc)
    
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation rate')
    
    call ESMF_StateGet(irrigState, "Irrigation frac",&
         irrigFracField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation frac')    
    call ESMF_FieldGet(irrigFracField, localDE=0,&
         farrayPtr=irrigFrac,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation frac')
    
    call ESMF_StateGet(irrigState, "Irrigation max root depth",&
         irrigRootDepthField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation max root depth')    
    call ESMF_FieldGet(irrigRootDepthField, localDE=0,&
         farrayPtr=irrigRootDepth,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation root depth')
    
    call ESMF_StateGet(irrigState, "Irrigation scale",&
         irrigScaleField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation scale')    
    call ESMF_FieldGet(irrigScaleField, localDE=0,&
         farrayPtr=irrigScale,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation scale')  
    
    irrigRate = 0.0  
    
    !----------------------------------------------------------------------
    ! Set start and end times for selected irrigation type
    !----------------------------------------------------------------------
        
    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

       timestep = LIS_rc%ts
       
       ! Adjust bounds by timestep to account for the fact that LIS_rc%hr, etc.
       ! will represents the END of the integration timestep window
       
       shift_otimes = otimes + (timestep/3600.)
       shift_otimee = otimee + (timestep/3600.)
       
       irrig_check_frozen_soil = .false.
       
       if((noah33_struc(n)%noah(t)%smc(1) - &
            noah33_struc(n)%noah(t)%sh2o(1)).gt.0001) then 
          irrig_check_frozen_soil = .true. 
       elseif((noah33_struc(n)%noah(t)%smc(2) - &
            noah33_struc(n)%noah(t)%sh2o(2)).gt.0001) then 
          irrig_check_frozen_soil = .true. 
       elseif((noah33_struc(n)%noah(t)%smc(3) - &
            noah33_struc(n)%noah(t)%sh2o(3)).gt.0001) then 
          irrig_check_frozen_soil = .true. 
       elseif((noah33_struc(n)%noah(t)%smc(4) - &
            noah33_struc(n)%noah(t)%sh2o(4)).gt.0001) then 
          irrig_check_frozen_soil = .true. 
       elseif(noah33_struc(n)%noah(t)%stc(2).le.LIS_CONST_TKFRZ) then
          irrig_check_frozen_soil = .true. 
       elseif(noah33_struc(n)%noah(t)%stc(3).le.LIS_CONST_TKFRZ) then
          irrig_check_frozen_soil = .true. 
       elseif(noah33_struc(n)%noah(t)%stc(4).le.LIS_CONST_TKFRZ) then
          irrig_check_frozen_soil = .true. 
       endif
       
       if(.not.irrig_check_frozen_soil) then 

          twater  = 0.0
          water   = 0.0
          asmc    = 0.0
          tsmcwlt = 0.0
          tsmcref = 0.0
          ma      = 0.0
          crootd  = 0.0
          lroot   = 0
          
          sldpth(1) = 0.1         ! Soil layer thicknesses (m)
          sldpth(2) = 0.3
          sldpth(3) = 0.6
          sldpth(4) = 1.0
          zdpth(1) = sldpth(1)         ! Soil layer depth from top (m)
          zdpth(2) = sldpth(1) + sldpth(2)
          zdpth(3) = sldpth(1) + sldpth(2) + sldpth(3)
          zdpth(4) = sldpth(1) + sldpth(2) + sldpth(3) + sldpth(4)
          
          smcmax =  noah33_struc(n)%noah(t)%smcmax
          psisat =  noah33_struc(n)%noah(t)%psisat
          dksat  =  noah33_struc(n)%noah(t)%dksat
          bexp   =  noah33_struc(n)%noah(t)%bexp
          smlow = 0.5
          smhigh = 6.0
          smcref1 = smcmax*(5.79e-9/dksat)**(1.0/(2.0*bexp+3.0))
          smcref = smcref1 + (smcmax-smcref1) / smhigh
          smcwlt1 = smcmax * (200.0/psisat)**(-1.0/bexp)
          smcwlt = smcwlt1 - smlow * smcwlt1
          
          gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
          chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
          if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
               (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
               chhr = chhr -1
          lhr = LIS_rc%hr +chhr
          if(lhr.ge.24) lhr = lhr-24
          if(lhr.lt.0) lhr = lhr+24
          
          ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0
          
          shdfac =  noah33_struc(n)%noah(t)%shdfac

          ! Calculate vegetation and root depth parameters
          ! If we are outside of the irrigation window, set rate to 0
          if ((ltime.ge.shift_otimee).or.(ltime.lt.shift_otimes)) then
             irrigRate(t) = 0.0
          endif
          
          if((ltime.ge.shift_otimes).and.(ltime.lt.shift_otimee)) then 
             vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
             !----------------------------------------------------------------------       
             !    Proceed if it is non-forest, non-baresoil, non-urban
             !----------------------------------------------------------------------       
             if(LIS_rc%lcscheme.eq."UMD") then !UMD
                veg_index1 = 6
                veg_index2 = 11
             elseif(LIS_rc%lcscheme.eq."UMD+MIRCAIrrig") then !UMD+MIRCAIrrig (Temporary, KRA)
                veg_index1 = 6
                veg_index2 = 16
!             elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP") then 
                ! TEMPORARILY ADDED HERE (HKB) for single crop tile option
                ! once multiple crop tile is unabled, +MIRCA should be separated
                ! into a case where veg_index1 and veg_index2 are set for croptypes
             elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP".or.&
                  LIS_rc%lcscheme.eq."IGBPNCEP+MIRCA") then  
                veg_index1 = 6
                veg_index2 = 14
             elseif(LIS_rc%lcscheme.eq."USGS") then !UMD
                veg_index1 = 2
                veg_index2 = 10
             else
                write(LIS_logunit,*) '[ERR] The landcover scheme ',trim(LIS_rc%lcscheme)
                write(LIS_logunit,*) '[ERR] is not supported for the Noah.3.3 irrigation module.'
                call LIS_endrun()
             endif
             
             if(vegt.ge.veg_index1.and.vegt.le.veg_index2&
                  .and.vegt.ne.LIS_rc%bareclass.and.&
                  vegt.ne.LIS_rc%urbanclass) then
 
                gsthresh = noah33_struc(n)%noah(t)%shdmin + & 
                     0.40 * (noah33_struc(n)%noah(t)%shdmax - &
                     noah33_struc(n)%noah(t)%shdmin)
                
                if(shdfac .ge. gsthresh) then 
                   crootd = irrigRootdepth(t)*shdfac
                   if(crootd.gt.0.and.crootd.lt.zdpth(1)) then 
                      lroot = 1
                      rdpth(1) = crootd
                   elseif(crootd .ge. zdpth(1).and.crootd .lt. zdpth(2) ) then
                      lroot = 2
                      rdpth(1) = sldpth(1)
                      rdpth(2) = crootd - zdpth(1)
                   elseif ( crootd.ge.zdpth(2).and.crootd .lt. zdpth(3) ) then
                      lroot = 3
                      rdpth(1) = sldpth(1)
                      rdpth(2) = sldpth(2)
                      rdpth(3) = crootd - zdpth(2)
                   elseif ( crootd.ge.zdpth(3).and.crootd .lt. zdpth(4) ) then
                      lroot = 4
                      rdpth(1) = sldpth(1)
                      rdpth(2) = sldpth(2)
                      rdpth(3) = sldpth(3)
                      rdpth(4) = crootd - zdpth(3)
                      !                      else
                      !                         print*,'error getting root depth'
                      !                         stop
                   endif

                   if(lroot.gt.0) &                                  
                        call rz_soil_moisture_conditions  &
                        (lroot, noah33_struc(n)%noah(t)%smc(1:lroot),     &
                        rdpth(1:lroot), smcwlt, smcref, asmc, tsmcwlt, tsmcref)

                   ! -----------------------------
                   ! Loop through irrigation types
                   ! -----------------------------
                   
                   irrigRate_frac = 0.
                   
                   do iType = 1,3  
                   
                      if(iType == 1) then              ! Sprinkler
                         otimes = otimess
                         irrhr = irrhrs
                         otimee = otimess + irrhrs
                      elseif(iType == 2) then          ! Drip
                         otimes = otimeds
                         irrhr = irrhrd
                         otimee = otimeds + irrhrd
                      elseif(iType == 3) then          ! Flood
                         otimes = otimefs
                         irrhr = irrhrf
                         otimee = otimefs + irrhrf
                      endif

                      if(irrigFrac(t, IType).gt.0) then 
                         ippix = irrigFrac(t, iType)*0.01
                         
                         ! Determine the amount of irrigation to apply if irrigated tile
                         if( IrrigScale(t,iType).gt.0.0 ) then ! irrigated tile
                            !                if(ippix.gt.0.0) then  ! irrigated tile                             
                            !!!!! SPRINKLER IRRIGATION
                            if(iType == 1) then
                               !----------------------------------------------------------------------       
                               !    Set the irrigation rate at start time; keep the value till next day
                               !    If local time at the tile fall in the irrigation check
                               !    hour then check the root zone average soil moisture
                               !----------------------------------------------------------------------       
                               if(ltime.eq.shift_otimes) then 
                                  irrigRate_frac(iType) = 0.0
                                  !-------------------------------------------------------------
                                  !     Compute the root zone accumlative soil moisture [mm], 
                                  !     field capacity [mm], and wilting point [mm] 
                                  !-------------------------------------------------------------
                                  if(lroot.gt.0) then 
 
                                     !---------------------------------------------------------------
                                     !     Get the root zone moisture availability to the plant
                                     !--------------------------------------------------------------- 
                                     ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)
                                     if(ma.le.LIS_rc%irrigation_thresh) then 
                                        !-----------------------------------------------------------------------------
                                        !     Compute irrigation rate
                                        !-----------------------------------------------------------------------------
                                        irrigRate_frac(Itype) = sprinkler_irrig_water (lroot, noah33_struc(n)%noah(t)%smc(1:lroot), &
                                          rdpth(1:lroot), smcref, irrigScale(t,iType))/(irrhr*3600.0)
                                     endif
                                  endif
                               endif
!!!!! DRIP IRRIGATION 
                            elseif(iType == 2) then
                               ! Need to get crop coefficient so that we can caculate unstressed Transp
                               !       RC=RSMIN/(XLAI*RCS*RCT*RCQ)
                               !       PCIRR=(RR+DELTA)/(RR*(1.+RC*CH)+DELTA)
                               ! CALL TRANSP (with PCIRR)
                               ! SM : 
                               ! (1) added below 2 variables to noah33dec structure in noah33_module.F90
                               !     real :: PC    ! PLANT COEFFICIENT (PC: UNITLESS FRACTION, 0-1) PC  WHERE PC*ETP = ACTUAL TRANSP
                               !     real :: PCIRR ! PLANT COEFFICIENT (PC: UNITLESS FRACTION, 0-1) with no soil moisutre stress-i.e., perfect irrigation.
                               ! (2) total transpriation with no soil moisture stress (perfect irrigation) from SUBROUTINE TRANSP would be
                               !     noah33_struc(n)%noah(t)%tveg * noah33_struc(n)%noah(t)%PCIRR / noah33_struc(n)%noah(t)%PC [w
                               !     and converted TVEG units [W / m2] to irrigrate units [kg / m /s]

                               twater = noah33_struc(n)%noah(t)%tveg * noah33_struc(n)%noah(t)%PCIRR / noah33_struc(n)%noah(t)%PC / LIS_CONST_LATVAP

                               !-----------------------------------------------------------------------------
                               !     Apply efficiency correction
                               !-----------------------------------------------------------------------------
                               twater2 = twater
                               twater = twater*(100.0/(100.0-efcor))
                               !-----------------------------------------------------------------------------
                               !     Compute irrigation rate
                               !-----------------------------------------------------------------------------
                               irrigRate_frac(iType) = twater  ! for drip calculation, twater is a rate [kg m-2/s]
                               noah33_struc(n)%noah(t)%smc(1) = &
                                    noah33_struc(n)%noah(t)%smc(1) + (twater-twater2)/(sldpth(1)*1000.0) !! check this with Sujay
                               
!!!!! FLOOD IRRIGATION
                            elseif(iType == 3) then
                               !-------------------------------------------------------------
                               !     Compute the root zone accumlative soil moisture [mm], 
                               !     field capacity [mm], and wilting point [mm] 
                               !-------------------------------------------------------------
                               if(lroot.gt.0) then 
 
                                  !---------------------------------------------------------------
                                  !     Get the root zone moisture availability to the plant
                                  !--------------------------------------------------------------- 
                                  !                            ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)   ! Original
                                  ma = (asmc-tsmcwlt) /(tsmcref - tsmcwlt)/IrrigScale(t,iType) ! BZ added IrrigScale
                                  
                                  if( ma .le. LIS_rc%irrigation_thresh ) then

                                     !-----------------------------------------------------------------------------
                                     !     Compute irrigation rate
                                     !-----------------------------------------------------------------------------
                                     irrigRate_frac(iType) = flood_irrig_water (LIS_rc%irrigation_mxsoildpth, &
                                          noah33_struc(n)%noah(t)%smc(1:LIS_rc%irrigation_mxsoildpth), &
                                          sldpth(1:LIS_rc%irrigation_mxsoildpth), smcmax, irrigScale(t,iType))/LIS_rc%ts
                                     
                                     ! BZ modification 4/2/2015 to account for ippix and all soil layers:
                                     do l = 1, LIS_rc%irrigation_mxsoildpth
                                        noah33_struc(n)%noah(t)%smc(l) = IrrigScale(t,iType)*smcmax + &
                                             (1-IrrigScale(t,iType))*noah33_struc(n)%noah(t)%smc(l)
                                     end do
                                     !                              noah33_struc(n)%noah(t)%smc(1) = IrrigScale(t)*smcmax + &
                                     !                                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(1)
                                     !                              noah33_struc(n)%noah(t)%smc(2) = IrrigScale(t)*smcmax + &
                                     !                                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(2)
                                     !                              noah33_struc(n)%noah(t)%smc(3) = IrrigScale(t)*smcmax + &
                                     !                                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(3)
                                     !                              noah33_struc(n)%noah(t)%smc(4) = IrrigScale(t)*smcmax + &
                                     !                                             (1-IrrigScale(t))*noah33_struc(n)%noah(t)%smc(4)
                                  endif
                               endif
                            endif
                         endif
                      end if
                   end do
                   do iType = 1, 3
                      irrigRate (t) = irrigRate(t) + irrigRate_frac(iType) * irrigFrac(t, iType)
                   end do
                end if
             end if
          endif
       endif
    enddo
    
  END SUBROUTINE UPDATE_IRRIG_STATE_CONCURRENT

  ! ********************************************************************

  SUBROUTINE rz_soil_moisture_conditions (lroot, smc, rdpth, smcwlt, smcref, &
       asmc, tsmcwlt, tsmcref)

    !-------------------------------------------------------------
    !     Compute the root zone accumlative soil moisture [mm], 
    !     field capacity [mm], and wilting point [mm] 
    !-------------------------------------------------------------

    implicit none
    integer, intent (in)                 :: lroot
    real, dimension (lroot), intent (in) :: smc, rdpth
    real, intent (in)                    :: smcwlt, smcref
    real, intent (out)                   :: asmc, tsmcwlt, tsmcref
    integer                              :: k

    asmc    = 0.0
    tsmcwlt = 0.0
    tsmcref = 0.0

    do k=1,lroot
       asmc = asmc + smc(k)*rdpth(k)*1000.0
       tsmcwlt = tsmcwlt + smcwlt * rdpth(k)*1000.0
       tsmcref = tsmcref + smcref * rdpth(k)*1000.0
    enddo
    
  END SUBROUTINE rz_soil_moisture_conditions

  ! ********************************************************************

  REAL FUNCTION sprinkler_irrig_water (lroot,  smc, rdpth, smcref, irrigScale)

    implicit none

    integer, intent (in)                 :: lroot
    real, dimension (lroot), intent (in) :: smc, rdpth
    real, intent (in)                    :: smcref, irrigScale
    real                                 :: water
    integer                              :: k

    water = 0.
    do k=1,lroot
       water = water + (smcref- smc(k))*rdpth(k)*1000.0
    enddo
                                         
    !-----------------------------------------------------------------------------
    !     Scale the irrigation intensity to the crop % when intensity < crop%.
    !     Expand irrigation for non-crop, non-forest when intensity > crop %
    !     in preference order of grassland first then rest.
    !     *scale is pre-computed for each tile in getirrpmapetc module in a way
    !     that is transparent to every tile irrigated or non-irrigated
    !-----------------------------------------------------------------------------

    water = water * irrigScale
    
    !-----------------------------------------------------------------------------
    !     Apply efficiency correction
    !-----------------------------------------------------------------------------

    sprinkler_irrig_water = water * (100.0/(100.0-efcor))

  END FUNCTION sprinkler_irrig_water

  ! ********************************************************************

  REAL FUNCTION drip_irrig_water

    implicit none

    real :: twater, twater2
 
    ! Need to get crop coefficient so that we can caculate unstressed Transp
    !       RC=RSMIN/(XLAI*RCS*RCT*RCQ)
    !       PCIRR=(RR+DELTA)/(RR*(1.+RC*CH)+DELTA)
    ! CALL TRANSP (with PCIRR)
    
    ! Then add enough water to get from actual Transp to unstressed Transp

    twater = 0.0
    !-----------------------------------------------------------------------------
    !     Apply efficiency correction
    !-----------------------------------------------------------------------------
    twater2 = twater
    drip_irrig_water = twater*(100.0/(100.0-efcor))

  END FUNCTION drip_irrig_water

  ! ********************************************************************
 
  REAL FUNCTION flood_irrig_water (mxsoildpth, smc, soildep, smcmax, irrigScale)

    implicit none

    integer, intent (in)                      :: mxsoildpth
    real, dimension (mxsoildpth), intent (in) :: smc, soildep
    real, intent (in)                         :: smcmax, irrigScale
    real                                      :: water
    integer                                   :: l

    water = 0.

    do l = 1, mxsoildpth
       if( l == 1 ) then
          water = (smcmax - smc(l))*soildep(l)*1000.0
       else
          ! BZ modification 4/2/2015 to saturate entire column and apply ippix 
          water = water + (smcmax - smc(l))*soildep(l)*1000.0
          !   twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(2))*sldpth(2)*1000.0
          !   twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(3))*sldpth(3)*1000.0
          !   twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(4))*sldpth(4)*1000.0
       endif
    end do

    !-----------------------------------------------------------------------------
    !     Scale the irrigation intensity to the crop % when intensity < crop%.
    !     Expand irrigation for non-crop, non-forest when intensity > crop %
    !     in preference order of grassland first then rest.
    !     *scale is pre-computed for each tile in getirrpmapetc module in a way
    !     that is transparent to every tile irrigated or non-irrigated
    !-----------------------------------------------------------------------------

    water = water * irrigScale
    
    !-----------------------------------------------------------------------------
    !     Apply efficiency correction
    !-----------------------------------------------------------------------------

    flood_irrig_water = water * (100.0/(100.0-efcor))
    
  END FUNCTION flood_irrig_water

end subroutine noah33_getirrigationstates

