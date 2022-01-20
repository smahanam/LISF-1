!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! BOP
! 
! !ROUTINE: irrigation_model
! \label{irrigation_model}

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
! May 2019: Jessica Erlingis; Incorporate W. Nie's updates into LIS
!                             and add optional flag for groundwater abstraction
! Feb 2020: Jessica Erlingis; Fix sprinkler irrigation winodw 
! Dec 2020: Hiroko Beaudoing; Updated things based on old LIS/Noah and 
!                             incorporated Sarith's concurrent irrigation types
!                             and Wanshu/Ben's modifications.


!EOP

MODULE IRRIGATION_MODULE

  use ESMF
  use LIS_coreMod
  use LIS_irrigationMod
  use LIS_logMod
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod, only : LIS_FORC_State
  
  IMPLICIT NONE

  PRIVATE

  real, parameter :: NOF_paccum_DAYS = 3, paccum_cutoff = 25.4, schedule_irate = 25.4/2. ! 1 inch in 48 hrs in [mm/day]
  
  type, public :: irrig_state
     
     real,  pointer :: irrigRate(:)
     real,  pointer :: irrigFrac(:)
     real,  pointer :: irrigType(:)
     real,  pointer :: irrigScale(:)
     real,  pointer :: irrigRootDepth(:)
     real,  pointer :: irriggwratio(:)
     real,  pointer :: prcp(:)
     real,  pointer :: aprcp(:)
     
  end type irrig_state
     
  type, public, extends (irrig_state) :: irrigation_model
     
   contains
     
     ! public
     procedure, public :: get_irrigstate
     procedure, public :: update_irrigrate
          
  end type irrigation_model

contains

  SUBROUTINE get_irrigstate (IM,nest,irrigState)

    implicit none
    
    class (irrigation_model), intent(inout) :: IM
    integer, intent(in)                     :: nest
    type(ESMF_State)                        :: irrigState
    type(ESMF_Field)                        :: irrigRateField,irrigFracField,irrigTypeField,prcpField
    type(ESMF_Field)                        :: irrigRootDepthField,irrigScaleField,irriggwratioField
    integer                                 :: rc
    real                                    :: Nsteps_paccum
    
    call ESMF_StateGet(irrigState, "Irrigation rate",irrigRateField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation rate')    
    call ESMF_FieldGet(irrigRateField, localDE=0,farrayPtr=IM%irrigRate,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation rate')
    
    call ESMF_StateGet(irrigState, "Irrigation frac",&
         irrigFracField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation frac')    
    call ESMF_FieldGet(irrigFracField, localDE=0,&
         farrayPtr=IM%irrigFrac,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation frac')
    
    call ESMF_StateGet(irrigState, "Irrigation max root depth",&
         irrigRootDepthField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation max root depth')    
    call ESMF_FieldGet(irrigRootDepthField, localDE=0,&
         farrayPtr=IM%irrigRootDepth,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation root depth')
    
    call ESMF_StateGet(irrigState, "Irrigation scale",&
         irrigScaleField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation scale')    
    call ESMF_FieldGet(irrigScaleField, localDE=0,&
         farrayPtr=IM%irrigScale,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation scale')  

    call ESMF_StateGet(irrigState, "Irrigation type",&
         irrigTypeField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Irrigation type')    
    call ESMF_FieldGet(irrigTypeField, localDE=0,&
         farrayPtr=IM%irrigType,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation type')

    if (LIS_irrig_struc(nest)%irrigation_GWabstraction == 1) then
       call ESMF_StateGet(irrigState, "Groundwater irrigation ratio",&
            irriggwratioField,rc=rc)
       call LIS_verify(rc,'ESMF_StateGet failed for Groundwater irrigation ratio')
       call ESMF_FieldGet(irriggwratioField, localDE=0,&
            farrayPtr=IM%irriggwratio,rc=rc)
       call LIS_verify(rc,'ESMF_FieldGet failed for Groundwater irrigation ratio')
    endif
    
    call ESMF_StateGet(LIS_FORC_State(nest),       &
         trim(LIS_FORC_Rainf%varname(1)),prcpField,&
         rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for prcp')
    call ESMF_FieldGet(prcpField, localDE=0,&
         farrayPtr=IM%prcp,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for rainfall')

    call ESMF_StateGet(irrigState, "Accumulated Precipitation",&
         irrigTypeField,rc=rc)
    call LIS_verify(rc,'ESMF_StateGet failed for Accumulated Precipitation')    
    call ESMF_FieldGet(irrigTypeField, localDE=0,&
         farrayPtr=IM%aprcp,rc=rc)
    call LIS_verify(rc,'ESMF_FieldGet failed for Accumulated Precipitation')

    !  set number of time steps within the preceding accumulative precip period "running average"
    ! -------------------------------------------------------------------------------------------

    Nsteps_paccum = NOF_paccum_DAYS * 86400. / LIS_rc%ts
    IM%aprcp = ((Nsteps_paccum - 1.)*IM%aprcp + IM%prcp) /  Nsteps_paccum ! mean precip rate during NOF_paccum_DAYS [kg/m^2/s]
    
  END SUBROUTINE get_irrigstate
  
  ! ----------------------------------------------------------------------------

  SUBROUTINE update_irrigrate (IM, nest, TileNo, croptype, longitude, veg_trigger, veg_thresh, SMWP, SMSAT, &
       SMREF, SMCNT, RDPTH)
    
    ! INPUTS:
    ! -------
    ! NEST         : LIS grid nest identifier
    ! TileNo       : Tile ID
    ! croptype     : 1-26
    ! LONGITUDE    : Tile longitude
    ! VEG_TRIGGER  : Current vegetation trigger value
    ! VEG_THRESH   : vegetation threshold to turn the trigger on
    ! SMWP         : soil moisture content at wilting point [m^3/m^3]
    ! SMSAT        : soil moisture content at saturation [m^3/m^3]
    ! SMREF        : ~soil field capacity - the upper limit of water content that soil can hold for plants [m^3/m^3]
    ! SMCNT(layers): soil moisture content in soil layers where plant roots are active [m^3/m^3]
    ! RDPTH(layers): thicknesses of active soil layers [m]

    ! OUTPUT
    ! ------
    ! irrigRate    : irrigation rate  [kg m-2 s-1] - internal state
    
    implicit none
    class (irrigation_model), intent(inout) :: IM
    integer, intent (in)                    :: nest, TileNo, croptype
    real, intent (in)                       :: longitude, veg_trigger, veg_thresh, SMWP, SMSAT, SMREF, SMCNT(:), RDPTH(:)
 
    ! locals
    real     :: HC, T1, T2, ma, asmc, tsmcwlt,tsmcref
    logical  :: season_active
    INTEGER  :: layer, season

    asmc    = 0.0
    tsmcwlt = 0.0
    tsmcref = 0.0
    ma      = 0.0

    !-------------------------------------------------------------
    !     Compute the root zone accumlative soil moisture [mm], 
    !     field capacity [mm], and wilting point [mm] 
    !-------------------------------------------------------------
    
    SOIL_LAYERS :do layer = 1, SIZE (RDPTH)
       asmc = asmc + SMCNT(layer)*RDPTH(layer)*1000.
       tsmcwlt = tsmcwlt + smwp * rdpth(layer)*1000.0
       tsmcref = tsmcref + smref* rdpth(layer)*1000.0
    end do SOIL_LAYERS
    
    !---------------------------------------------------------------
    !     Get the root zone moisture availability to the plant
    !---------------------------------------------------------------
    
    if(tsmcref > tsmcwlt)then
       ma = (asmc - tsmcwlt) /(tsmcref - tsmcwlt)
    else
       ma = -1.
    endif

    ! --------
    ! Set time
    ! --------
    
    HC = real(LIS_rc%hr) + real(LIS_rc%mn)/60. + real(LIS_rc%ss)/3600. + &
         12.* longitude/180.
    IF (HC >= 24.) HC = HC - 24.
    IF (HC <   0.) HC = HC + 24.
    T1 = CEILING (HC)     - real(LIS_rc%ts)/3600.
    T2 = FLOOR   (HC + 1) + real(LIS_rc%ts)/3600.
    
    if((HC >= T1).and.(HC < T2))then
       HC = real(NINT(HC))
    end if
    
    season_active = .false.
           
    CALENDAR: if ( LIS_irrig_struc(nest)%cropcalendar .eq. "none" ) then
              
       ! -----------------------------
       ! Vegetation Trigger
       ! -----------------------------
              
       if(veg_trigger >= veg_thresh) season_active = .true.

    else
              
       ! -----------------------------
       ! Crop calendar
       ! -----------------------------
       
       NOF_SEASONS: do season = 1, SIZE(LIS_irrig_struc(nest)%plantDay(TileNo,:))
          IF(IS_WITHIN_SEASON(LIS_rc%doy,NINT(LIS_irrig_struc(nest)%PLANTDAY(TileNo, season)), &
               NINT(LIS_irrig_struc(nest)%harvestDay(TileNo, season)))) season_active = .true.
       END do NOF_SEASONS
       
    endif CALENDAR
    
    ! Run irrigation model if the crop growing season is active
    ! ----------------------------------------------------------
    
    CROP_GROWING_SEASON: if (season_active) then

       PADDY: if(croptype == 3) then
          ! PADDY
          T1 = LIS_irrig_struc(nest)%flood_start
          T2 = LIS_irrig_struc(nest)%flood_start + LIS_irrig_struc(nest)%flood_duration
          
          if ((HC >= T1).AND.(HC < T2)) then
             ! within flood irrigation period
             if (T1 == HC) then
                ! set rate at start time
                do layer = 1,SIZE (SMCNT)
                   IM%irrigRate(TileNo) = IM%irrigRate(TileNo) + (smsat -smcnt(layer))*rdpth(layer)*1000.0*100.0/ &
                        (100.0-LIS_irrig_struc(nest)%flood_efcor)
                enddo
                IM%irrigRate(TileNo) = IM%irrigRate(TileNo)*IM%IrrigScale(TileNo)/(T2 - T1)/3600.
             endif
          else
             IM%irrigRate(TileNo) = 0
          endif
       else
          
          SM_DEFICIT_OR_SCHEDULE: if (LIS_irrig_struc(nest)%irrigation_schedule ==0) then

             SOILM: if(ma >= 0) then                      
                if (IM%irrigType(TileNo) == 1.) call irrig_by_type (nest,HC, ma, smref,SMCNT, RDPTH, IM%IrrigScale(TileNo), SRATE = IM%irrigRate(TileNo))
                if (IM%irrigType(TileNo) == 2.) call irrig_by_type (nest,HC, ma, smref,SMCNT, RDPTH, IM%IrrigScale(TileNo), DRATE = IM%irrigRate(TileNo))
                if (IM%irrigType(TileNo) == 3.) call irrig_by_type (nest,HC, ma, smref,SMCNT, RDPTH, IM%IrrigScale(TileNo), FRATE = IM%irrigRate(TileNo))                
             endif SOILM

          else
             
             CHECK_ACUUM_PRECIP: if (IM%aprcp(TileNo)*86400*NOF_paccum_DAYS < paccum_cutoff) then                
                if (IM%irrigType(TileNo) == 1.) call irrig_on_schedule (nest,HC, IM%IrrigScale(TileNo), SRATE = IM%irrigRate(TileNo))
                if (IM%irrigType(TileNo) == 2.) call irrig_on_schedule (nest,HC, IM%IrrigScale(TileNo), DRATE = IM%irrigRate(TileNo))
                if (IM%irrigType(TileNo) == 3.) call irrig_on_schedule (nest,HC, IM%IrrigScale(TileNo), FRATE = IM%irrigRate(TileNo))                
             endif CHECK_ACUUM_PRECIP             

          endif SM_DEFICIT_OR_SCHEDULE
       endif PADDY
    else
       !  Outside the season
       IM%irrigRate(TileNo) = 0.
       
    endif CROP_GROWING_SEASON
    
  END SUBROUTINE update_irrigrate
  
  ! ----------------------------------------------------------------------------
  
  SUBROUTINE irrig_on_schedule (nest, HC, IrrigScale, SRATE, DRATE, FRATE)
    
    implicit none

    INTEGER, intent (in)                    :: nest    
    REAL, intent (in)                       :: HC, IrrigScale
    REAL, optional, intent (inout)          :: SRATE, DRATE, FRATE
    REAL                                    :: H1, H2, IT
    
    if(present (SRATE)) then
       ! SPRINKLER IRRIGATION
       H1 = LIS_irrig_struc(nest)%sprinkler_start
       H2 = LIS_irrig_struc(nest)%sprinkler_start + LIS_irrig_struc(nest)%sprinkler_duration
       IT = LIS_irrig_struc(nest)%sprinkler_thresh
       
       if ((HC >= H1).AND.(HC < H2)) then
          ! Irrirate is set to a predifined value
          ! rates for the day and maintains the same rate through out the irrigation
          ! duration (H1 <= HC < H2).
          SRATE = schedule_irate*IrrigScale*100.0/(100.0-LIS_irrig_struc(nest)%sprinkler_efcor)/(H2 - H1)/3600.
       else
          SRATE = 0.
       endif
    endif
       
    if(present (DRATE)) then
       ! DRIP IRRIGATION
       H1 = LIS_irrig_struc(nest)%drip_start
       H2 = LIS_irrig_struc(nest)%drip_start + LIS_irrig_struc(nest)%drip_duration
       IT = LIS_irrig_struc(nest)%drip_thresh
       
       if ((HC >= H1).AND.(HC < H2)) then
          ! Notice drip uses the same soil moisture threshold of sprinkler but with 0.% efficiency correction.
          DRATE = schedule_irate*IrrigScale*100.0/(100.0-LIS_irrig_struc(nest)%drip_efcor)/(H2 - H1)/3600.
       else
          DRATE = 0.
       endif
    endif
    
    if(present (FRATE)) then
       ! FLOOD IRRIGATION
       H1 = LIS_irrig_struc(nest)%flood_start
       H2 = LIS_irrig_struc(nest)%flood_start + LIS_irrig_struc(nest)%flood_duration
       IT = LIS_irrig_struc(nest)%flood_thresh
       
       if ((HC >= H1).AND.(HC < H2)) then
          FRATE = schedule_irate*IrrigScale*100.0/(100.0-LIS_irrig_struc(nest)%flood_efcor)/(H2 - H1)/3600.
       else
          FRATE = 0.
       endif
    endif
     
  END SUBROUTINE irrig_on_schedule
  
  ! ----------------------------------------------------------------------------
  
  SUBROUTINE irrig_by_type (nest, HC, ma, SMREF, SMCNT, RDPTH, IrrigScale, SRATE, DRATE, FRATE)
    
    implicit none

    INTEGER, intent (in)                    :: nest    
    REAL, intent (in)                       :: HC, ma, SMREF, SMCNT(:), RDPTH(:), IrrigScale
    REAL, optional, intent (inout)          :: SRATE, DRATE, FRATE
    REAL                                    :: H1, H2, IT
    
    if(present (SRATE)) then
       ! SPRINKLER IRRIGATION
       H1 = LIS_irrig_struc(nest)%sprinkler_start
       H2 = LIS_irrig_struc(nest)%sprinkler_start + LIS_irrig_struc(nest)%sprinkler_duration
       IT = LIS_irrig_struc(nest)%sprinkler_thresh
       
       if ((HC >= H1).AND.(HC < H2)) then
          ! The model uses rootzone soil moisture state at H1 to compute irrigation
          ! rates for the day and maintains the same rate through out the irrigation
          ! duration (H1 <= HC < H2).
          if((ma <= IT).AND.(H1 == HC)) &
               SRATE = crop_water_deficit (SMCNT, RDPTH, SMREF,LIS_irrig_struc(nest)%sprinkler_efcor)*IrrigScale/(H2 - H1)/3600.
       else
          SRATE = 0.
       endif
    endif
       
    if(present (DRATE)) then
       ! DRIP IRRIGATION
       H1 = LIS_irrig_struc(nest)%drip_start
       H2 = LIS_irrig_struc(nest)%drip_start + LIS_irrig_struc(nest)%drip_duration
       IT = LIS_irrig_struc(nest)%drip_thresh
       
       if ((HC >= H1).AND.(HC < H2)) then
          ! use SMCNT at H1 during H1 <= HC < H2 to compute irrigrate.
          ! Notice drip uses the same soil moisture threshold of sprinkler but with 0.% efficiency correction.
          if((ma <= IT).AND.(H1 == HC)) &
               DRATE = crop_water_deficit (SMCNT, RDPTH, SMREF,LIS_irrig_struc(nest)%drip_efcor)*IrrigScale/(H2 - H1)/3600.
       else
          DRATE = 0.
       endif
    endif
    
    if(present (FRATE)) then
       ! FLOOD IRRIGATION
       H1 = LIS_irrig_struc(nest)%flood_start
       H2 = LIS_irrig_struc(nest)%flood_start + LIS_irrig_struc(nest)%flood_duration
       IT = LIS_irrig_struc(nest)%flood_thresh
       
       if ((HC >= H1).AND.(HC < H2)) then
          ! use SMCNT at H1 during H1 <= HC < H2 to compute irrigrate.
          if((ma <= IT).AND.(H1 == HC)) &
               FRATE = crop_water_deficit (SMCNT, RDPTH, SMREF,LIS_irrig_struc(nest)%flood_efcor)*IrrigScale/(H2 - H1)/3600.
       else
          FRATE = 0.
       endif
    endif
     
  END SUBROUTINE irrig_by_type
  
  ! ----------------------------------------------------------------------------

  REAL FUNCTION crop_water_deficit (SMCNT, RDPTH, SMREF, efcor)

    implicit none

    real, intent (in)                    :: SMCNT(:), RDPTH(:), SMREF, efcor
    real                                 :: twater
    integer                              :: layer

    crop_water_deficit = 0.
    
    !---------------------------------------------------------------
    !     Get the moisture availability to the plant
    !---------------------------------------------------------------    

    crop_water_deficit = 0.
    
    do layer = 1,SIZE (SMCNT)
       crop_water_deficit = crop_water_deficit + (smref -smcnt(layer))*rdpth(layer)*1000.0
    enddo
    
    crop_water_deficit = crop_water_deficit*100.0/(100.0-efcor)
    
  END FUNCTION crop_water_deficit
  
  ! ----------------------------------------------------------------------------

  logical FUNCTION IS_WITHIN_SEASON (DOY,DP, DH)

    implicit none

    integer, intent (in)                 :: DOY,DP, DH

    IS_WITHIN_SEASON = .false.
    if(DH > DP) then
       if((DOY >= DP).AND.(DOY <=  DH)) IS_WITHIN_SEASON = .true.
    elseif (DH < DP) then
       if((DOY >= DP).AND.(DOY <= 366)) IS_WITHIN_SEASON = .true.
       if((DOY >=  1).AND.(DOY <=  DH)) IS_WITHIN_SEASON = .true.
    endif
      
  end FUNCTION IS_WITHIN_SEASON
  
END MODULE IRRIGATION_MODULE
