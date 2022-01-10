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
  use IRRIGATION_MODULE
  

  implicit none

  integer,intent(in)   :: nest
  integer              :: rc
  integer              :: TileNo,tid,gid,vegt
  type(ESMF_State)     :: irrigState
  real                 :: rdpth, laifac, laithresh, smcwlt, smcref, smcmax
  integer              :: veg_index1,veg_index2
  real,  allocatable   :: laimax(:,:),laimin(:,:)
  type(irrigation_model) :: IM

  ! _______________________________________________________

  call IM%get_irrigstate (irrigState)
  
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
     
     ! Process only irrigated tiles
     IRRF: if(IM%irrigFrac(TileNo).gt.0) then
        ! Proceed if it is non-forest, non-baresoil, non-urban
        vegt = LIS_surface(nest,LIS_rc%lsm_index)%tile(TileNo)%vegt
        VEGIF:  if(vegt.ge.veg_index1.and.vegt.le.veg_index2&
             .and.vegt.ne.LIS_rc%bareclass.and.&
             vegt.ne.LIS_rc%urbanclass) then

           ! Calculate active root depth
           ! ---------------------------           

           if(clsmf25_struc(nest)%cat_param(TileNo)%laimax.ne.&
                clsmf25_struc(nest)%cat_param(TileNo)%laimin) then 
              laifac = (LIS_lai(nest)%tlai(tid) -                     &
                   clsmf25_struc(nest)%cat_param(TileNo)%laimin)/  &
                   (clsmf25_struc(nest)%cat_param(TileNo)%laimax - &
                   clsmf25_struc(nest)%cat_param(TileNo)%laimin)
           else
              laifac = 0.0
           endif
           
           rdpth = IM%irrigRootdepth(TileNo)*laifac
           
           ! compute vegetation threshold for the trigger
           ! --------------------------------------------
           
           laithresh =   clsmf25_struc(nest)%cat_param(TileNo)%laimin + & 
                LIS_irrig_struc(nest)%veg_thresh *         &
                (clsmf25_struc(nest)%cat_param(TileNo)%laimax - &
                clsmf25_struc(nest)%cat_param(TileNo)%laimin)
       
           smcmax = clsmf25_struc(nest)%cat_param(TileNo)%poros
           smcwlt = clsmf25_struc(nest)%cat_param(TileNo)%wpwet*           &  
                clsmf25_struc(nest)%cat_param(TileNo)%poros
           smcref = (clsmf25_struc(nest)%cat_param(TileNo)%wpwet +         & 
                0.333* (1.-clsmf25_struc(nest)%cat_param(TileNo)%wpwet))*  &
                clsmf25_struc(nest)%cat_param(TileNo)%poros

           ! get irrigation rates from the irrigation model
           ! ----------------------------------------------
           
           call IM%update_irrigrate (nest,TileNo, LIS_domain(nest)%grid(gid)%lon,     &
                LIS_lai(nest)%tlai(tid),laithresh,                                    &
                smcwlt,smcmax,smcref,                                                 &
                (/clsmf25_struc(nest)%cat_diagn(TileNo)%rzmc/),                       &
                (/rdpth/))
                      
        endif VEGIF
     endif IRRF
  end do
  
  ! Update land surface model's moisture state
  ! ------------------------------------------  

end subroutine clsmf25_getirrigationstates
