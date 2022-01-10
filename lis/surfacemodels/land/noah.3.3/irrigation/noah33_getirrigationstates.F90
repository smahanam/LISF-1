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
! !ROUTINE: noah33_getirrigationstates
! \label{noah33_getirrigationstates}
! 
! !INTERFACE:
subroutine noah33_getirrigationstates(nest,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use noah33_lsmMod
  use LIS_irrigationMod
  use IRRIGATION_MODULE
  
! !DESCRIPTION:        
!
! Calculate water requirement and apply the amount for Sprinkler, Drip, or
! Flood irrigation types accordingly. The water requirement is checked once
! a day at specified local time. Irrigation check time, duration, and threshold
! are specifield for each irrigation type in the lis.config. 
!
! Sprinker: water requirement is based on the actual root zone soil moisture
!           availability. The irrigation threshold determines when to irrigate
!           (i.e. when root zone soil moisture falls below a percentage of 
!           the fields capacity). Apply the amount to precipitation  
!           during the specified hours, eg 6-10 am LST.
! Drip: water requirement is based on the transpiration stress.  Apply the 
!       amount to the top layer soil moisture during the specified hours.
! Flood: water requirement is based on the amount to bring the entire column
!        to saturation and applied to all soil layers at once.
!
! The root zone is actual maximum root depth rather than NOAH's vegetation 
! parameter "NROOT". The root depth grows following the GVF seasonality to 
! reflect the water demand increase/decrease for crops. (note: need age of
! trees for orchard).
! 
! Irrigation amount is scaled to grid total crop fraction when intensity
! is less than the fraction.  Irrigation is expanded to non-crop, non-forest,
! non-baresoil/urban tiles if intensity exceeds grid total crop fraction.
! In the latter case, scaled irrigation is applied to grassland first, 
! then further applied over the rest of tiles equally if the intensity 
! exceeds grassland fraction as well. The irrigation scale and alighning 
! intensity with the crops or land cover types, as well as irrigation type
! determination for the tile are done before this routine is called.   
!
! Optionally efficiency correction is applied to account for field loss. 
!
! This version includes modifications and updates as follows:
! 1) Change growing season threshold (i.e. 40% of GFRAC range) in lis.config
! 2) Use Crop Plant/Harvest dates in addition to GVF for trigger
! 3) Irrigation check time, duration, threshold, and efficiency are set in 
!    lis.config and can be different for irrigation types 
!
! Note: need to add rice paddy in flood 
!       add option to use crop calender for growing season check
!
! REVISION HISTORY:
!
! Aug 2008: Hiroko Kato; Initial code
! Nov 2012: Sujay Kumar, Incorporated into LIS
! Jun 2014: Ben Zaitchik; Added flood scheme
! Feb 2020: Jessica Erlingis; Fix sprinkler irrigation winodw 
! Dec 2020: Hiroko Beaudoing; Updated things based on old LIS/Noah and 
!                             incorporated Sarith's concurrent irrigation types
!                             and Wanshu/Ben's modifications.
!                             
!EOP
  implicit none

  integer, intent(in)  :: nest
  integer              :: rc
  integer              :: TileNo,tid,gid,vegt
  type(ESMF_State)     :: irrigState
  real                 :: sldpth(noah33_struc(nest)%nslay)
  real                 :: rdpth(noah33_struc(nest)%nslay)
  real                 :: zdpth(noah33_struc(nest)%nslay)
  real                 :: crootd
  integer              :: lroot,veg_index1,veg_index2
  real                 :: gsthresh
  logical              :: irrig_check_frozen_soil
  type(irrigation_model) :: IM 

  ! _______________________________________________________

  call IM%get_irrigstate (irrigState)

  ! Set vegetation type index to be irrigated
  ! -----------------------------------------
  
  if(LIS_rc%lcscheme.eq."UMD") then !UMD
     veg_index1 = 6
     veg_index2 = 11
  elseif(LIS_rc%lcscheme.eq."UMD+MIRCA") then !UMD+MIRCA (Temporary, KRA)
     veg_index1 = 6
     veg_index2 = 40
  elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP") then
     veg_index1 = 6
     veg_index2 = 14
  elseif(LIS_rc%lcscheme.eq."IGBPNCEP+MIRCA") then  
     veg_index1 = 6
     veg_index2 = 46
  elseif(LIS_rc%lcscheme.eq."USGS") then 
     veg_index1 = 2
     veg_index2 = 10
  elseif(LIS_rc%lcscheme.eq."UMDCROPMAP") then 
     veg_index1 = 5
     veg_index2 = 32
  else
     write(LIS_logunit,*) '[ERR] The landcover scheme ',trim(LIS_rc%lcscheme)
     write(LIS_logunit,*) '[ERR] is not supported for the Noah.3.3 irrigation module.'
     call LIS_endrun()
  endif
  
  ! Set global soil  parameters
  sldpth(1) = noah33_struc(nest)%lyrthk(1)         ! Soil layer thicknesses (m)
  sldpth(2) = noah33_struc(nest)%lyrthk(2)
  sldpth(3) = noah33_struc(nest)%lyrthk(3)
  sldpth(4) = noah33_struc(nest)%lyrthk(4)
  zdpth(1) = sldpth(1)         ! Soil layer depth from top (m)
  zdpth(2) = sldpth(1) + sldpth(2)
  zdpth(3) = sldpth(1) + sldpth(2) + sldpth(3)
  zdpth(4) = sldpth(1) + sldpth(2) + sldpth(3) + sldpth(4)
  
  !---------------------------------------------------------------
  ! Main tile loop
  !---------------------------------------------------------------
  
  TILE_LOOP: do TileNo = 1,LIS_rc%npatch(nest,LIS_rc%lsm_index)
     
     gid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%index
     tid = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%tile_id
     
     ! frozen tile check
     irrig_check_frozen_soil = .false.
     
     if((noah33_struc(nest)%noah(TileNo)%smc(1) - &
          noah33_struc(nest)%noah(TileNo)%sh2o(1)).gt.0001) then 
        irrig_check_frozen_soil = .true. 
     elseif((noah33_struc(nest)%noah(TileNo)%smc(2) - &
          noah33_struc(nest)%noah(TileNo)%sh2o(2)).gt.0001) then 
        irrig_check_frozen_soil = .true. 
     elseif((noah33_struc(nest)%noah(TileNo)%smc(3) - &
          noah33_struc(nest)%noah(TileNo)%sh2o(3)).gt.0001) then 
        irrig_check_frozen_soil = .true. 
     elseif((noah33_struc(nest)%noah(TileNo)%smc(4) - &
          noah33_struc(nest)%noah(TileNo)%sh2o(4)).gt.0001) then 
        irrig_check_frozen_soil = .true. 
     elseif(noah33_struc(nest)%noah(TileNo)%stc(2).le.LIS_CONST_TKFRZ) then
        irrig_check_frozen_soil = .true. 
     elseif(noah33_struc(nest)%noah(TileNo)%stc(3).le.LIS_CONST_TKFRZ) then
        irrig_check_frozen_soil = .true. 
     elseif(noah33_struc(nest)%noah(TileNo)%stc(4).le.LIS_CONST_TKFRZ) then
        irrig_check_frozen_soil = .true. 
     endif
     
     ! Process only non-frozen tiles
     FROZEN: if(.not.irrig_check_frozen_soil) then
        
        ! Process only irrigated tiles
        IRRF: if(IM%irrigFrac(TileNo).gt.0) then
           
           ! Determine the amount of irrigation to apply if irrigated tile
           IRRS: if( IM%IrrigScale(TileNo).gt.0.0 ) then ! irrigated tile
           
              ! Proceed if it is non-forest, non-baresoil, non-urban
              vegt = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%vegt  
              VEGIF: if(vegt.ge.veg_index1.and.vegt.le.veg_index2     &
                   .and.vegt.ne.LIS_rc%bareclass                      &
                   .and.vegt.ne.LIS_rc%urbanclass) then

                 ! Calculate active root depth
                 ! ---------------------------
                 
                 crootd = IM%irrigRootdepth(TileNo)*noah33_struc(nest)%noah(TileNo)%shdfac
                 lroot  = 0
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
                 endif

                 if (lroot == 0) then
                    write(LIS_logunit,*) '[ERR] lroot should be > 0!'
                    call LIS_endrun()
                 endif

                 ! compute vegetation threshold for the trigger
                 ! --------------------------------------------
                 
                 gsthresh = noah33_struc(nest)%noah(TileNo)%shdmin +   & 
                      LIS_irrig_struc(nest)%veg_thresh * (noah33_struc(nest)%noah(TileNo)%shdmax - &
                      noah33_struc(nest)%noah(TileNo)%shdmin)

                 ! get irrigation rates from the irrigation model
                 ! ----------------------------------------------
           
                 call IM%update_irrigrate (nest,TileNo, LIS_domain(nest)%grid(gid)%lon, &
                      noah33_struc(nest)%noah(TileNo)%shdfac,gsthresh,                  &
                      noah33_struc(nest)%noah(TileNo)%smcwlt,                           &
                      noah33_struc(nest)%noah(TileNo)%smcmax,                           &
                      noah33_struc(nest)%noah(TileNo)%smcref,                           &
                      noah33_struc(nest)%noah(TileNo)%smc(:lroot),                      &
                      rdpth(:lroot))
                                   
              endif VEGIF
           endif IRRS
        endif IRRF
     endif FROZEN
  end do TILE_LOOP

  ! Update land surface model's moisture state
  ! ------------------------------------------
    
end subroutine noah33_getirrigationstates
