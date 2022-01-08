MODULE IRRIGATION_MODULE

  use ESMF
  use LIS_coreMod
  use LIS_irrigationMod
  
  IMPLICIT NONE

  ! This module computes irrigation rates by 3 different methods: sprinkler, flood and drip.

  PRIVATE
  
  type, public :: irrigation_model
     
   contains
     
     ! public
     procedure, public :: update_irrigrate
          
  end type irrigation_model

contains

  ! ----------------------------------------------------------------------------
  SUBROUTINE update_irrigrate (this, nest, TileNo, longitude, veg_trigger, veg_thresh, SMWP, SMSAT, &
       SMREF, SMCNT, RDPTH, IrrigScale,irrigType, irrigRate)
    
    ! INPUTS:
    ! -------
    ! NEST         : LIS grid nest identifier
    ! TileNo       : Tile ID
    ! LONGITUDE    : Tile longitude
    ! VEG_TRIGGER  : Current vegetation trigger value
    ! VEG_THRESH   : vegetation threshold to turn the trigger on
    ! SMWP         : soil moisture content at wilting point [m^3/m^3]
    ! SMSAT        : soil moisture content at saturation [m^3/m^3]
    ! SMREF        : ~soil field capacity - the upper limit of water content that soil can hold for plants [m^3/m^3]
    ! SMCNT(layers): soil moisture content in soil layers where plant roots are active [m^3/m^3]
    ! RDPTH(layers): thicknesses of active soil layers [m]
    ! IrrigScale   : IrrigScale parameter
    ! irrigType    : Irrigation Type : (1) Sprinkler; (2) Drip; and (3) Flood

    ! OUTPUT
    ! ------
    ! irrigRate    : irrigation rate  [kg m-2 s-1]
    
    implicit none
    class (irrigation_model), intent(inout) :: this
    integer, intent (in)                    :: nest, TileNo
    real, intent (in)                       :: longitude, veg_trigger, veg_thresh, SMWP, SMSAT, SMREF, SMCNT(:), RDPTH(:), IrrigScale, irrigType 
    REAL, intent (inout)                    :: irrigRate
 
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
       
       SOILM: if(ma >= 0) then       

          if (irrigType == 1.) call irrig_by_type (nest,HC, ma, smref,SMCNT, RDPTH, IrrigScale, SRATE = irrigRate)
          if (irrigType == 2.) call irrig_by_type (nest,HC, ma, smref,SMCNT, RDPTH, IrrigScale, DRATE = irrigRate)
          if (irrigType == 3.) call irrig_by_type (nest,HC, ma, smref,SMCNT, RDPTH, IrrigScale, FRATE = irrigRate)
 
       endif SOILM
    else
       !  Outside the season
       irrigRate = 0.
       
    endif CROP_GROWING_SEASON
    
  END SUBROUTINE update_irrigrate
  
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
