MODULE IRRIGATION_MODULE

  use LIS_irrigationMod
  
  IMPLICIT NONE

  ! This module computes irrigation rates by 3 different methods: sprinkler, flood and drip.

  PRIVATE
  
  type, public :: irrigation_model
     
   contains
     
     ! public
     procedure, public :: irrig_by_type
     procedure, public :: IS_WITHIN_SEASON
     
     ! private
     procedure, private :: cwd => crop_water_deficit
     
  end type irrigation_model

contains

  ! ----------------------------------------------------------------------------

  SUBROUTINE irrig_by_type (this,nest, HC, SMWP, SMSAT, SMREF, SMCNT, RDPTH, IrrigScale, SRATE, DRATE, FRATE)

    ! INPUTS:
    ! -------
    ! NEST         : LIS grid nest identifier
    ! HC           : local time [hours]
    ! SMWP         : soil moisture content at wilting point [m^3/m^3]
    ! SMSAT        : soil moisture content at saturation [m^3/m^3]
    ! SMREF        : ~soil field capacity - the upper limit of water content that soil can hold for plants [m^3/m^3]
    ! SMCNT(layers): soil moisture content in soil layers where plant roots are active [m^3/m^3]
    ! RDPTH(layers): thicknesses of active soil layers [m]
    ! IrrigScale   : IrrigScale parameter

    ! OUTPUTS
    ! -------
    ! if(present (SRATE)) : Sprinkler irrigation rate  [kg m-2 s-1]
    ! if(present (DRATE)) : Drip irrigation rate  [kg m-2 s-1]
    ! if(present (FRATE)) : Flood irrigation rate  [kg m-2 s-1]
    
    implicit none
    class (irrigation_model), intent(inout) :: this
    INTEGER, intent (in)                    :: nest    
    REAL, intent (in)                       :: HC, SMWP, SMSAT, SMREF, SMCNT(:), RDPTH(:), IrrigScale
    REAL, optional, intent (inout)          :: SRATE, DRATE, FRATE 
    REAL                                    :: H1, H2, IT, ma, asmc, tsmcwlt,tsmcref
    INTEGER                                 :: layer

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

    SOILM: if(ma >= 0) then
    
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
                  SRATE = this%cwd (SMCNT, RDPTH, SMREF,LIS_irrig_struc(nest)%sprinkler_efcor)*IrrigScale/(H2 - H1)/3600.
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
                  DRATE = this%cwd(SMCNT, RDPTH, SMREF,LIS_irrig_struc(nest)%drip_efcor)*IrrigScale/(H2 - H1)/3600.
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
                  FRATE = this%cwd (SMCNT, RDPTH, SMREF,LIS_irrig_struc(nest)%flood_efcor)*IrrigScale/(H2 - H1)/3600.
          else
             FRATE = 0.
          endif
       endif
    endif SOILM
    
  END SUBROUTINE irrig_by_type
  
  ! ----------------------------------------------------------------------------

  REAL FUNCTION crop_water_deficit (this, SMCNT, RDPTH, SMREF, efcor)

    implicit none
    class(irrigation_model),intent(inout):: this
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

  logical FUNCTION IS_WITHIN_SEASON (this, DOY,DP, DH)

    implicit none
    class(irrigation_model),intent(inout):: this
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
