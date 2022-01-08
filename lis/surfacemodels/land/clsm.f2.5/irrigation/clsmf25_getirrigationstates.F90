!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

!BOP
! 
! !ROUTINE: clsmf25_getirrigationstates
! \label{clsmf25_getirrigationstates}
! 
! !INTERFACE:
subroutine clsmf25_getirrigationstates(nest,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use clsmf25_lsmMod
  use LIS_vegDataMod
  use LIS_irrigationMod
  use LIS_constantsMod, ONLY : pi => LIS_CONST_PI
  use IRRIGATION_MODULE
  
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
! Aug 2008: Hiroko Kato; Initial code for Noah LSM
! Feb 2014: Sujay Kumar; Implemented in LIS based on the work of
!           John Bolten and student. 
! Jul 2014: Ben Zaitchik; added flood routine

!EOP
  implicit none
  ! moved to lis.config
  ! Sprinkler parameters
  !real, parameter      :: otimess = 6.0 ! local trigger check start time [hour]
  !real, parameter      :: irrhrs = 4.   ! duration of irrigation hours 
  ! Drip parameters (not currently implemented)
  !real, parameter      :: otimeds = 6.0 ! local trigger check start time [hour]
  !real, parameter      :: irrhrd = 12.0   ! duration of irrigation hours 
 ! Flood parameters
  !real, parameter      :: otimefs = 6.0 ! local trigger check start time [hour]
  !real, parameter      :: irrhrf = 1.0   ! duration of irrigation hours 
  !!!real, parameter      :: ffreq = 0.0 ! frequency of flood irrig [days] set to 0.0 to use thresh instead
  
  ! real, parameter      :: efcor = 76.0      ! Efficiency Correction (%)

  integer,intent(in)   :: nest
  integer              :: rc
  integer              :: TileNo,tid,gid,vegt
  type(ESMF_State)     :: irrigState
  type(ESMF_Field)     :: irrigRateField,irrigFracField,irrigTypeField
  type(ESMF_Field)     :: irrigRootDepthField,irrigScaleField
  
  real,  pointer       :: irrigRate(:), irrigFrac(:), irrigType(:)
  real,  pointer       :: irrigRootDepth(:), irrigScale(:)
  real                 :: ltime, T1, T2, rdpth, laifac, laithresh, smcwlt, smcref, smcmax
  integer              :: veg_index1,veg_index2, season
  real,  allocatable   :: laimax(:,:),laimin(:,:)
  logical              :: season_active
  type(irrigation_model) :: IM
  
  if(clsmf25_struc(nest)%modelStart) then 
     clsmf25_struc(nest)%modelStart = .false. 

     allocate(laimax(LIS_rc%lnc(nest),LIS_rc%lnr(nest)))
     allocate(laimin(LIS_rc%lnc(nest),LIS_rc%lnr(nest)))

     call LIS_read_laimax(nest,laimax)
     call LIS_read_laimin(nest,laimin)

     do TileNo=1,LIS_rc%npatch(nest,LIS_rc%lsm_index)
        if (laimax(LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,       &
             LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row).ne.-9999.00) then
           clsmf25_struc(nest)%cat_param(TileNo)%laimax =            &
                laimax(LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,     &
                LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row)
        else
           clsmf25_struc(nest)%cat_param(TileNo)%laimax = 1.0
        endif
        if (laimin(LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,        &
             LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row).ne.-9999.00) then
           clsmf25_struc(nest)%cat_param(TileNo)%laimin =            &
                laimin(LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%col,     &
                LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%row)
        else
           clsmf25_struc(nest)%cat_param(TileNo)%laimin = 0.0
        endif
     enddo
     
     deallocate(laimax)
     deallocate(laimin)

  endif

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

  call ESMF_StateGet(irrigState, "Irrigation type",&
       irrigTypeField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Irrigation type')    
  call ESMF_FieldGet(irrigTypeField, localDE=0,&
       farrayPtr=irrigType,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation type')
  
  ! Set vegetation type index to be irrigated   
  select case ( LIS_rc%lcscheme )
  case( "UMD" )
     veg_index1 = 6
     veg_index2 = 11
  case( "IGBP", "IGBPNCEP", "MODIS" )
     veg_index1 = 6
     veg_index2 = 14
  case( "IGBPNCEP+MIRCA")
     veg_index1 = 6
     veg_index2 = 46
  case( "USGS" )
     veg_index1 = 2
     veg_index2 = 10
  case default
     write(LIS_logunit,*) "The landcover scheme, ",trim(LIS_rc%lcscheme),","
     write(LIS_logunit,*) "is not supported for irrigation. Stopping program ... "
     call LIS_endrun()
  end select

  !---------------------------------------------------------------
  ! Main tile loop
  !---------------------------------------------------------------
  
  do TileNo = 1,LIS_rc%npatch(nest,LIS_rc%lsm_index)
     
     gid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%index
     tid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%tile_id

     ! Set time
     ! --------
     ltime = real(LIS_rc%hr) + real(LIS_rc%mn)/60. + real(LIS_rc%ss)/3600. + &
          12.* (LIS_domain(nest)%grid(gid)%lon/pi)
     IF (ltime >= 24.) ltime = ltime - 24.
     IF (ltime <   0.) ltime = ltime + 24.
     T1 = CEILING (ltime)     - real(LIS_rc%ts)/3600.
     T2 = FLOOR   (ltime + 1) + real(LIS_rc%ts)/3600.

     if((ltime >= T1).and.(ltime < T2))then
        ltime = real(NINT(ltime))
     end if
     
     ! Process only irrigated tiles
     IRRF: if(irrigFrac(TileNo).gt.0) then
        ! Proceed if it is non-forest, non-baresoil, non-urban
        vegt = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%vegt
        VEGIF:  if(vegt.ge.veg_index1.and.vegt.le.veg_index2&
             .and.vegt.ne.LIS_rc%bareclass.and.&
             vegt.ne.LIS_rc%urbanclass) then
           
           ! Calculate vegetation and root depth parameters
           if(clsmf25_struc(nest)%cat_param(TileNo)%laimax.ne.&
                clsmf25_struc(nest)%cat_param(TileNo)%laimin) then 
              laifac = (LIS_lai(nest)%tlai(tid) -                     &
                   clsmf25_struc(nest)%cat_param(TileNo)%laimin)/  &
                   (clsmf25_struc(nest)%cat_param(TileNo)%laimax - &
                   clsmf25_struc(nest)%cat_param(TileNo)%laimin)
           else
              laifac = 0.0
           endif
           
           rdpth = irrigRootdepth(TileNo)*laifac
           season_active = .false.
           
           CALENDAR: if ( LIS_irrig_struc(nest)%cropcalendar .eq. "none" ) then
              
              ! -----------------------------
              ! Vegetation Trigger
              ! -----------------------------

              laithresh =   clsmf25_struc(nest)%cat_param(TileNo)%laimin + & 
                   LIS_irrig_struc(nest)%veg_thresh *         &
                   (clsmf25_struc(nest)%cat_param(TileNo)%laimax - &
                   clsmf25_struc(nest)%cat_param(TileNo)%laimin)

              if(LIS_lai(nest)%tlai(tid) .ge. laithresh) season_active = .true.

           else
              
              ! -----------------------------
              ! Crop calendar
              ! -----------------------------

              NOF_SEASONS: do season = 1, SIZE(LIS_irrig_struc(nest)%plantDay(TileNo,:))
                 IF(IM%IS_WITHIN_SEASON(LIS_rc%doy,NINT(LIS_irrig_struc(nest)%PLANTDAY(TileNo, season)), &
                      NINT(LIS_irrig_struc(nest)%harvestDay(TileNo, season)))) season_active = .true.
              END do NOF_SEASONS
              
           endif CALENDAR
           
           ! Run irrigation model if the crop growing season is active
           ! ----------------------------------------------------------
           
           CROP_GROWING_SEASON: if (season_active) then
              
              smcmax = clsmf25_struc(nest)%cat_param(TileNo)%poros
              smcwlt = clsmf25_struc(nest)%cat_param(TileNo)%wpwet*              &  
                      clsmf25_struc(nest)%cat_param(TileNo)%poros
              smcref = (clsmf25_struc(nest)%cat_param(TileNo)%wpwet +            & 
                      0.333* (1.-clsmf25_struc(nest)%cat_param(TileNo)%wpwet))*  &
                      clsmf25_struc(nest)%cat_param(TileNo)%poros
              
              ! Sprinkler
              if (irrigType(TileNo) == 1) call IM%irrig_by_type (nest,ltime,     &
                   smcwlt,                                                       &
                   smcmax,                                                       &
                   smcref,                                                       &
                   (/clsmf25_struc(nest)%cat_diagn(TileNo)%rzmc/),               &
                   (/rdpth/),IrrigScale(TileNo),SRATE = irrigRate(TileNo))
              
              ! Drip
              if (irrigType(TileNo) == 2) call IM%irrig_by_type (nest,ltime,     &
                   smcwlt,                                                       &
                   smcmax,                                                       &
                   smcref,                                                       &
                   (/clsmf25_struc(nest)%cat_diagn(TileNo)%rzmc/),               &
                   (/rdpth/),IrrigScale(TileNo),DRATE = irrigRate(TileNo))
              
              ! Flood
              if (irrigType(TileNo) == 3) call IM%irrig_by_type (nest,ltime,     &
                   smcwlt,                                                       &
                   smcmax,                                                       &
                   smcref,                                                       &
                   (/clsmf25_struc(nest)%cat_diagn(TileNo)%sfmc/),               &
                   (/clsmf25_struc(nest)%cat_param(TileNo)%dzsf/1000./),         &
                   IrrigScale(TileNo),FRATE = irrigRate(TileNo))
              
           else
              !  Outside the season
              irrigRate(TileNo) = 0.
              
           endif CROP_GROWING_SEASON
                      
        endif VEGIF
     endif IRRF
  end do
  
end subroutine clsmf25_getirrigationstates
