#define VERIFY_(A) if(A /=0)then;print *,'ERROR code',A,'at',__LINE__;call exit(3);endif
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_MODISNative_lc
!  \label{read_MODISNative_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  Aug 2013: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!  23  Apr 2014: KR Arsenault;  Added new optimized interpolation code
!
! !INTERFACE:

module mod_ESA_lc

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
            LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use CLSM_util, only : LDT_RegridRaster,   &
       NC_VarID, LDT_g5map,                 &
       c_data => G5_BCSDIR,                 &
       NX => nc_g5_rst,                     &
       NY => nr_g5_rst, histogram 
  use LDT_numericalMethodsMod, only : LDT_quicksort

  implicit none
  include 'netcdf.inc'	

  contains

    subroutine read_ESA_lc(n, num_types, fgrd, maskarray )

    implicit none
    
    ! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(inout) :: num_types
    real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
    real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
    !
    ! !DESCRIPTION:
    !  This subroutine reads the MODIS landcover data and returns the 
    !  distribution of vegetation in each grid cell, in a lat/lon
    !  projection.  Also, the landmask is either generated and/or 
    !  read in this routine.
    !
    !  The arguments are:
    !  \begin{description}
    !   \item[n]
    !     index of nest
    !   \item[fgrd]
    !     fraction of grid covered by each vegetation type
    !   \item[maskarray]
    !     landmask for the region of interest
    !   \end{description}
    !EOP      
    !
    ! IGBP-NCEP landcover version:
    integer, parameter :: IN_cols_igbpncep = 43200*3 
    integer, parameter :: IN_rows_igbpncep = 21600*3
    real,    parameter :: IN_xres = 1.0/360.0
    real,    parameter :: IN_yres = 1.0/360.0
    !   character*1 :: read_igbpncep_veg(IN_cols_igbpncep,IN_rows_igbpncep)
    character*1, allocatable :: read_igbpncep_veg(:,:)
    
    integer :: ftn, ierr, ios1
    logical :: file_exists
    integer :: i, t, c, r, line
    integer :: input_cols, input_rows
    integer :: glpnc, glpnr               ! Parameter (global) total columns and rows
    integer :: subpnc, subpnr             ! Parameter subsetted columns and rows
    integer :: mi                         ! Total number of input param grid array points
    integer :: mo                         ! Total number of output LIS grid array points
    real    :: param_gridDesc(20)         ! Input parameter grid desc array
    real    :: subparam_gridDesc(20)      ! Subsetted parameter grid desc array
    integer, allocatable  :: lat_line(:,:), lon_line(:,:)
    real,    allocatable  :: gi(:)        ! Input parameter 1d grid
    logical*1,allocatable :: li(:)        ! Input logical mask (to match gi)
    real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
    logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)
    real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%nt) ! Output lis 1d grid
    logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%nt) ! Output logical mask (to match go)
    
    real, allocatable :: subset_veg(:,:)  ! Read input parameter
    real, allocatable :: input_vegtype(:,:)
    real      :: vegcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
    real      :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
    
    !__________________________________________________________________
    
    num_types = 20 
    !- Check if land cover file exists:
    inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists )
    if(.not. file_exists) then
       write(LDT_logunit,*)"[ERR] The landcover map: ",trim(LDT_rc%vfile(n))," does not exist."
       write(LDT_logunit,*)"Program stopping ..."
       call LDT_endrun
    endif
    
    ftn = LDT_getNextUnitNumber()
    write(LDT_logunit,*) "[INFO] Reading landcover file: ",trim(LDT_rc%vfile(n))
    
    !- Open LDT land cover file:
    select case ( LDT_rc%lc_type(n) )
       
    case ( "IGBPNCEP" )
       
       !- Assign additional land cover types, including generic water points: 
       LDT_rc%wetlandclass = 11
       LDT_rc%urbanclass   = 13
       LDT_rc%snowclass    = 15
       LDT_rc%glacierclass = 15
       LDT_rc%bareclass    = 16
       LDT_rc%waterclass   = 17
       
       input_cols = IN_cols_igbpncep
       input_rows = IN_rows_igbpncep
       allocate( input_vegtype(input_cols, input_rows) )
       input_vegtype = LDT_rc%waterclass
       
       !- Set parameter grid array inputs: 
       param_gridDesc(1)  = 0.          ! Latlon
       param_gridDesc(2)  = input_cols
       param_gridDesc(3)  = input_rows
       param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat
       param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon
       param_gridDesc(6)  = 128
       param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat
       param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon
       param_gridDesc(9)  = IN_xres     ! dx: 0.0083333
       param_gridDesc(10) = IN_yres     ! dy: 0.0083333
       param_gridDesc(20) = 64
       
       ! Open file:
       open(ftn, file=LDT_rc%vfile(n), status='old', form='unformatted',&
            access ='direct', recl=(input_cols*input_rows), iostat=ios1)
       allocate( read_igbpncep_veg(IN_cols_igbpncep,IN_rows_igbpncep) )
       
       ! Veg types are stored as 8-bit unsigned integers:
       read(ftn,rec=1) read_igbpncep_veg
       
       ! Reverse-Y and Convert 8-bit unsigned integers:
       i = 0
       do r = input_rows, 1, -1
          i = i + 1  
          do c = 1, input_cols
             input_vegtype(c,i) = Real(ichar(read_igbpncep_veg(c,r)))  
          end do
       end do
       deallocate( read_igbpncep_veg )
       
       ! --------------------------------
       !    case ( "IGBP" )
       !    case ( "UMD" )
       !    case ( "PFT" )
       
       !- Other landcover classifications associated with MODIS landcover:
    case default  ! Non-supported options
       write(LDT_logunit,*) " The native MODIS map with land classification: ",&
            trim(LDT_rc%lc_type(n)),", is not yet supported."
       write(LDT_logunit,*) " -- Please select: IGBPNCEP "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
       
    end select
    
    write(LDT_logunit,*) "[INFO] Done reading ", trim(LDT_rc%vfile(n))
    
    ! -------------------------------------------------------------------
    vegcnt   = 0.
    vegtype  = float(LDT_rc%waterclass)
    maskarray= 0.0
    fgrd     = 0.0
    
    ! -------------------------------------------------------------------
    !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
    ! -------------------------------------------------------------------
    
    !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
    subparam_gridDesc = 0.
    call LDT_RunDomainPts( n, LDT_rc%lc_proj, param_gridDesc(:), &
         glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )
    
    allocate( subset_veg(subpnc, subpnr) )
    subset_veg = LDT_rc%waterclass
    
    !- Subset parameter read-in array:
    line = 0
    do r = 1, subpnr
       do c = 1, subpnc
          subset_veg(c,r) = input_vegtype(lon_line(c,r),lat_line(c,r))
       enddo
    enddo
   
    ! -------------------------------------------------------------------
    !     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
    ! -------------------------------------------------------------------
    mi = subpnc*subpnr
    allocate( gi(mi), li(mi) )
    gi = float(LDT_rc%waterclass)
    li = .false.
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    lo1 = .false.;  lo2 = .false.
    
    !- Assign 2-D array to 1-D for aggregation routines:
    i = 0
    do r = 1, subpnr
       do c = 1, subpnc;  i = i + 1
          gi(i) = subset_veg(c,r)
          if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
       enddo
    enddo
    
    !- Aggregation/Spatial Transform Section:
    select case ( LDT_rc%lc_gridtransform(n) )
       
       !- (a) Estimate NON-TILED dominant land cover types (vegtype):
    case( "neighbor", "mode" )
       
       !- Transform parameter from original grid to LIS output grid:
       call LDT_transform_paramgrid(n, LDT_rc%lc_gridtransform(n), &
            subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )
       
       !- Convert 1D vegcnt to 2D grid arrays:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             vegtype(c,r) = go1(i)
          enddo
       enddo
       
       !- (b) Estimate TILED land cover files (vegcnt):
    case( "tile" )
       
       !- Transform parameter from original grid to LIS output grid:
       call LDT_transform_paramgrid(n, LDT_rc%lc_gridtransform(n), &
            subparam_gridDesc, mi, LDT_rc%nt, gi, li, mo, go2, lo2 )
       
       !- Convert 1D vegcnt to 2D grid arrays:
       i = 0
       do r = 1, LDT_rc%lnr(n) 
          do c = 1, LDT_rc%lnc(n)  
             i = i + 1
             do t = 1, LDT_rc%nt
                vegcnt(c,r,t) = go2(i,t)
             end do
          enddo
       enddo
       
    end select  ! End vegtype/cnt aggregation method
    deallocate( gi, li )
    
    ! ........................................................................
    
    !- Bring 2-D Vegtype to 3-D Vegcnt tile space:
    if ( LDT_rc%lc_gridtransform(n) == "none"     .or. &    ! -- NON-TILED SURFACES
         LDT_rc%lc_gridtransform(n) == "neighbor" .or. &
         LDT_rc%lc_gridtransform(n) == "mode" ) then  
       
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             if ( vegtype(c,r) .le. 0 ) then
                vegtype(c,r) = float(LDT_rc%waterclass)
             endif
             if ( (nint(vegtype(c,r)) .ne. LDT_rc%waterclass ) .and. &
                  (nint(vegtype(c,r)) .ne. LDT_rc%udef)) then
                vegcnt(c,r,NINT(vegtype(c,r))) = 1.0
             endif
          enddo
       end do
    endif   ! End NON-TILED vegetation option
    
    !- Estimate fraction of grid (fgrid) represented by vegetation type::
    call param_index_fgrdcalc( n, LDT_rc%lc_proj, LDT_rc%lc_gridtransform(n), &
         LDT_rc%waterclass, LDT_rc%nt, vegcnt, fgrd )
    
    ! -------------------------------------------------------------------
    !    CREATE OR READ-IN (OR IMPOSE) LAND MASK FILE AND CREATE
    !    SURFACE MAP
    ! -------------------------------------------------------------------
    
    !- "READ-IN" land mask file, if user-specified:
    if( LDT_rc%mask_type(n) == "readin" ) then
       
       call read_maskfile( n, vegtype, fgrd, maskarray )
       
       !- "CREATE" land mask and surface type fields (user-specified):
    elseif( LDT_rc%mask_type(n) == "create" ) then
       
       call create_maskfile( n, LDT_rc%nt, LDT_rc%lc_gridtransform(n), &
            vegtype, vegcnt, maskarray )
       
    end if
    deallocate( input_vegtype )
    deallocate( subset_veg )
    
    call LDT_releaseUnitNumber(ftn)
  end subroutine read_ESA_lc
  
  ! ---------------------------------------------------------------------
  !
  SUBROUTINE ESA2MOSAIC (ityp)
    
    implicit none

    integer, intent (inout), dimension (:)    :: ityp
    integer  , parameter   :: nc_esa = 129600, nr_esa = 64800
    integer*2, allocatable, target, dimension (:,:) :: esa_veg
    integer*2, pointer    , dimension (:,:) :: subset
    integer  , allocatable, dimension (:)   :: tile_id
    integer :: i,j, k, status, ncid, maxcat, dx,dy, esa_type, tid, cid
    integer :: mos1, mos2
    real    :: mfrac, sfrac, tfrac, tem (6)
    integer, allocatable, dimension (:) :: density, loc_int
    real   , allocatable, dimension (:) :: loc_val
    logical, allocatable, dimension (:) :: unq_mask
    real   , allocatable       :: veg (:,:)
    integer :: NBINS, NPLUS

    ! Reading ESA vegetation types
    !-----------------------------

    allocate (esa_veg (1:nc_esa, 1: nr_esa))

    status    = NF_OPEN (trim(c_data)//'/ESA_GlobalCover.nc', NF_NOWRITE, ncid)  ; VERIFY_(STATUS) 

    do j = 1,nr_esa
       status  = NF_GET_VARA_INT2 (ncid,3,(/1,j/),(/nc_esa,1/),esa_veg(:,j)) ; VERIFY_(STATUS) 
    end do
    status = NF_CLOSE(ncid)

!
! Reading number of tiles
! -----------------------

    maxcat = LDT_g5map%NT_GEOS 

!
! Loop through tile_id raster
! ___________________________

    allocate (tile_id (1:NX))   
    allocate(veg(1:maxcat,1:6))
    veg = 0.

    dx = nc_esa / NX
    dy = nr_esa / NY
    
    do j=1,NY

       ! read a row

       tile_id(:) = LDT_g5map%rst(:,j)

       do i = 1,NX
          if((tile_id (i) >= 1).and.(tile_id(i)  <= maxcat)) then
             if (associated (subset)) NULLIFY (subset)
             subset => esa_veg((i-1)*dx +1 :i*dx, (j-1)*dy +1:j*dy)
             
             NPLUS = count(subset >= 1 .and. subset <= 230)
             
             if(NPLUS > 0)  then
                allocate (loc_int (1:NPLUS))
                allocate (unq_mask(1:NPLUS))
                loc_int = pack(subset,mask = (subset >= 1 .and. subset <= 230))
                call LDT_quicksort (loc_int)
                unq_mask = .true.
                
                do k = 2,NPLUS 
                   unq_mask(k) = .not.(loc_int(k) == loc_int(k-1))
                end do
                NBINS = count(unq_mask)
                
                allocate(loc_val (1:NBINS))
                allocate(density (1:NBINS))
                loc_val = 1.*pack(loc_int,mask =unq_mask)
                call histogram (size(subset,1)*size(subset,2), NBINS, density, loc_val, real(subset))   
                
                do k = 1, nbins

                   if (density (k) > 0) then
                      esa_type = int (loc_val(k))
                      !                 if (esa_type ==  10)  veg (tile_id(i),10) = 1.* density (k)   ; lakes inland water
                      if (esa_type ==  10)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)   ! inconsistent mask/veg  
                      if (esa_type ==  11)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type ==  14)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type ==  20)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type ==  30)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type ==  40)  veg (tile_id(i), 1) = veg (tile_id(i), 1) + 1.* density (k)
                      if (esa_type ==  50)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 1.* density (k)
                      if (esa_type ==  60)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 1.* density (k)
                      if (esa_type ==  70)  veg (tile_id(i), 3) = veg (tile_id(i), 3) + 1.* density (k)
                      if (esa_type ==  90)  veg (tile_id(i), 3) = veg (tile_id(i), 3) + 1.* density (k)
                      if (esa_type == 100)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 0.5* density (k)
                      if (esa_type == 100)  veg (tile_id(i), 3) = veg (tile_id(i), 3) + 0.5* density (k)
                      if (esa_type == 110)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 0.3* density (k)
                      if (esa_type == 110)  veg (tile_id(i), 5) = veg (tile_id(i), 5) + 0.3* density (k)
                      if (esa_type == 110)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 0.4* density (k)
                      if (esa_type == 120)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 0.2* density (k)
                      if (esa_type == 120)  veg (tile_id(i), 5) = veg (tile_id(i), 5) + 0.2* density (k)
                      if (esa_type == 120)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 0.6* density (k)
                      if (esa_type == 130)  veg (tile_id(i), 5) = veg (tile_id(i), 5) + 1.* density (k)
                      if (esa_type == 140)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      
                      if((j > NINT(real(NY)*(40./180.))).and.(j < NINT(real(NY)*(140./180.)))) then
                         if (esa_type == 150)  veg (tile_id(i),5) = veg (tile_id(i),5) + 0.5* density (k)
                         if (esa_type == 150)  veg (tile_id(i),4) = veg (tile_id(i),4) + 0.5* density (k)
                      else
                         if (esa_type == 150)  veg (tile_id(i),6) = veg (tile_id(i),6) + 0.5* density (k)
                         if (esa_type == 150)  veg (tile_id(i),4) = veg (tile_id(i),4) + 0.5* density (k)
                      end if
                      
                      if((j > NINT(real(NY)*(70./180.))).and.(j < NINT(real(NY)*(110./180.)))) then 
                         if (esa_type == 160)  veg (tile_id(i), 1) = veg (tile_id(i), 1) + 1.* density (k) 
                      else
                         if (esa_type == 160)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 1.* density (k)
                      end if
                      
                      if (esa_type == 170)  veg (tile_id(i), 1) = veg (tile_id(i), 1) + 1.* density (k)
                      if (esa_type == 180)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type == 190)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type == 200)  veg (tile_id(i), 5) = veg (tile_id(i), 5) + 1.* density (k)
                      if (esa_type == 210)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)  ! inconsistent mask/veg  
                      if (esa_type == 220)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)  ! inconsistent mask/veg  
                      !                         if (esa_type == 210)  veg (tile_id(i),11) = 1.* density (k)      ; ocean
                      !                         if (esa_type == 220)  veg (tile_id(i), 9) = 1.* density (k)     ; ice     
                   endif
                enddo
                deallocate (loc_int,unq_mask,loc_val,density)
             endif
          endif
       end do
    end do

    deallocate (tile_id)   

! Canopy height and ASCAT roughness length
!
!    call ascat_r0 (NX,NY,gfile, z0)
!
!    if(jpl_height) then
!       call jpl_canoph (NX,NY,gfile, z2)
!    else
!       allocate (z2(1:maxcat))       
!    endif
!
! Now create mosaic_veg_fracs file
! --------------------------------

    do k = 1, maxcat

       read (11,'(i8,i8,5(2x,f9.4))') tid,cid
       tem = 0.
       tem(1:6)=veg (k,1:6)

       if(sum(tem).gt.0)then

          mfrac = -10.
          sfrac = -10.
          mos1  = 100
          mos2  = 100         

          do i = 1,6
             if(mfrac.le.tem(i))then
                sfrac = mfrac
                mos2  = mos1
                mfrac = tem(i)
                mos1  = i
             elseif(sfrac.le.tem(i)) then
                if(tem(i).lt.mfrac)then
                   sfrac = tem(i)
                   mos2  = i
                endif
             endif             
          end do

          mfrac = max (mfrac,0.)
          sfrac = max (sfrac,0.)
          tfrac = mfrac + sfrac
          mfrac = mfrac / tfrac
          sfrac = sfrac / tfrac

          if (mos1 == 100) then
             mos1 = 4
             mos2 = 4
             mfrac= 1.
             sfrac= 0. 
          endif

          if (sfrac == 0.) mos2 = mos1 ! No secondary type
!          if(.not.jpl_height) z2(k) = VGZ2(mos1)
          ityp (k) = mos1
       endif
    end do

    deallocate (veg)
   
  END SUBROUTINE ESA2MOSAIC


end module mod_ESA_lc

