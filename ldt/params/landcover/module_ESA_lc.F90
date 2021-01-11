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
  use LDT_coreMod,     only : LDT_rc, LDT_domain
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
            LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use CLSM_util, only : LDT_RegridRaster,   &
       NC_VarID, LDT_g5map,                 &
       c_data => G5_BCSDIR,                 &
       NX => nc_g5_rst,                     &
       NY => nr_g5_rst, histogram , write_clsm_files, init_geos2lis_mapping,G52LIS, read_clsm_maskfile
  use LDT_numericalMethodsMod, only : LDT_quicksort
  use map_utils
  use LDT_ClimateBCsReader,    ONLY : ClimateBCsReader
  
  implicit none
  include 'netcdf.inc'

  private
  public read_ESA_lc, ESA2MOSAIC,  ESA2CLM45, ESA2CLM4
  
  integer,  parameter   :: nc_esa = 129600, nr_esa = 64800
  character*100         :: ESA_FILE = '/ESA_GlobalCover.nc'
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
  SUBROUTINE ESA2MOSAIC(n, num_types, fgrd, maskarray) 
    
    implicit none
    type (ClimateBCsReader):: bcr
    integer, intent(in) :: n
    integer, intent(in) :: num_types
    real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
    real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real    :: rlat(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real    :: rlon(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real    :: param_grid(20)
    integer*2, allocatable, target, dimension (:,:) :: esa_veg
    integer*2, pointer    , dimension (:,:) :: subset
    integer  , allocatable, dimension (:)   :: tile_id
    integer :: i,j, k, status, ncid, maxcat, dx,dy, esa_type, tid, cid, fmos
    integer :: mos1, mos2
    real    :: mfrac, sfrac, tfrac, tem (6)
    integer, allocatable, dimension (:) :: density, loc_int
    real   , allocatable, dimension (:) :: loc_val
    logical, allocatable, dimension (:) :: unq_mask
    real   , allocatable       :: veg (:,:)
    integer :: NBINS, NPLUS, c,r,gr,gc, glpnc, glpnr
    real, allocatable,  dimension (:)       :: ityp, ityp_lis
    real, pointer, dimension (:,:)          :: vegtype
    integer, allocatable,  dimension (:,:)  :: vegtype_lc

    ! Initialize the global mask and GEOS5 to LIS mapping
    ! ---------------------------------------------------

    call read_clsm_maskfile(n, maskarray)
    if (associated(LDT_g5map) .eqv. .false.) call init_geos2lis_mapping (n) 
   ! if (.not.LDT_g5map%init) call init_geos2lis_mapping (n) 

    LDT_rc%waterclass   = 7
    LDT_rc%bareclass    = 5
    LDT_rc%urbanclass   = 5
    LDT_rc%glacierclass = 5
    LDT_rc%wetlandclass = 0
    LDT_rc%snowclass    = 0

    !- Double-check landcover spatial transform option:
    if( LDT_rc%lc_gridtransform(n)=="tile" .or. LDT_rc%lc_gridtransform(n)=="mode")then
       write(LDT_logunit,*) " (in read_CLSMF25_lc) :: The 'tile' or 'mode' spatial transform option"
       write(LDT_logunit,*) "          has been selected, but these options are not"
       write(LDT_logunit,*) "          currently supported for the Catchment F2.5 LSM. "
       write(LDT_logunit,*) "          Please select option 'none' and run LDT again. "
       write(LDT_logunit,*) " Program stopping ..."
       call LDT_endrun
    end if

    ! Reading ESA vegetation types
    !-----------------------------

    allocate (esa_veg (1:nc_esa, 1: nr_esa))

    status    = NF_OPEN (trim(c_data)//trim(ESA_FILE), NF_NOWRITE, ncid)  ; VERIFY_(STATUS) 

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
    allocate(ityp(1:maxcat))
    allocate(ityp_lis(1:LDT_g5map%NT_LIS))
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

    if (write_clsm_files) then
       fmos = LDT_getNextUnitNumber()
       open (fmos, file = 'LDT_clsm/mosaic_veg_typs_fracs', form = 'formatted', action = 'write')
    endif

    do k = 1, maxcat

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
          ityp (k) = real(mos1)
          if (write_clsm_files) write(fmos,'(i8,i8,2(2x,i3),2(2x,f6.2))') k,LDT_g5map%catid_index(k),mos1, mos2,100.*mfrac,100.*sfrac
       endif
    end do

    if (write_clsm_files) then
       close (fmos, status = 'keep')
       call LDT_releaseUnitNumber(fmos)
    endif

    !- Determine global/complete parameter domain number of points:
    param_grid(:) = LDT_rc%mask_gridDesc(n,:)
    glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
    glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1
    allocate( vegtype(glpnc,glpnr), stat=status );  VERIFY_(STATUS)
    allocate( vegtype_lc(LDT_rc%lnc(n),LDT_rc%lnr(n)), stat=status );  VERIFY_(STATUS)
    vegtype = LDT_rc%waterclass
    vegtype_lc = NINT(LDT_rc%udef)
    ityp_lis = G52LIS (ityp)

! - For now - complete domains:
    k = 0
    do r = 1,  LDT_rc%lnr(n) !glpnr
       do c = 1, LDT_rc%lnc(n) ! glpnc
         call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
                           rlat(c,r),rlon(c,r))
         gr = nint((rlat(c,r)-param_grid(4))/param_grid(10))+1
         gc = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1
          
          if( LDT_rc%global_mask(gc,gr) >= 1. ) then
             k = k + 1
             if( ityp_lis(k) > 0 ) then
                vegtype(c,r) = ityp_lis(k)
             elseif( ityp_lis(k) == 0 ) then
                vegtype(c,r) = 5.   ! Assign landcover for glacier points
             endif
          endif          
      end do
   enddo    

!!- Final fgrd output fields:
 !- Build in subsetting domain:
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)

         call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
                           rlat(c,r),rlon(c,r))
         gr = nint((rlat(c,r)-param_grid(4))/param_grid(10))+1
         gc = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1

         if( vegtype(gc,gr) > 0. .and. vegtype(gc,gr).ne.LDT_rc%waterclass ) then 
            fgrd(c,r,int(vegtype(gc,gr))) = 1.0
            vegtype_lc (c,r) = vegtype(gc,gr)
         else
            fgrd(c,r,LDT_rc%waterclass) = 1.0
         endif
      enddo
   enddo

   call BCR%update_sibinputs(n, ITYP = vegtype_lc)
   deallocate (veg, vegtype, vegtype_lc, ityp, ityp_lis)
  ! nullify (LDT_g5map)
   
  END SUBROUTINE ESA2MOSAIC

!
! ---------------------------------------------------------------------
!

  SUBROUTINE ESA2CLM45 (nest, nc, nr, gfile)

    implicit none

    integer  , intent (in) :: nc, nr, nest
    character (*)          :: gfile
    
    integer  , parameter   :: N_lon_clm = 7200, N_lat_clm = 3600, lsmpft = 25
    integer*2, allocatable, target, dimension (:,:) :: esa_veg
    integer*2, pointer    , dimension (:,:) :: subset
    integer  , allocatable, dimension (:)   :: tile_id, i_esa2clm, j_esa2clm
    integer :: i,j, k,n, status, ncid, varid, maxcat, dx,dy, esa_type, tid, cid, ii, jj   
    real    :: dx_clm, dy_clm, x_min_clm (N_lon_clm), y_min_clm (N_lat_clm), clm_fracs(lsmpft)
    real    :: minlon,maxlon,minlat,maxlat,tile_lat, scale, ftot
    integer :: cpt1, cpt2, cst1, cst2  ! CLM-carbon types
    real    :: cpf1, cpf2, csf1, csf2  ! CLM-carbon fractions
    DOUBLE PRECISION,  allocatable, dimension (:) :: lon_esa, lat_esa
    DOUBLE PRECISION       :: EDGEN, EDGEE, EDGES, EDGEW 
    
    REAL, ALLOCATABLE, DIMENSION (:,:,:) :: PCTPFT 
    integer, allocatable, dimension (:) :: density, loc_int
    real   , allocatable, dimension (:) :: loc_val
    logical, allocatable, dimension (:) :: unq_mask
    integer :: NBINS, NPLUS, FCLM
    integer, allocatable, dimension (:,:) :: clm_veg
    integer :: esa_clm_veg (2)
    real    :: esa_clm_frac(2)
    real, dimension (:,:), allocatable      :: maskarray

    ! Initialize the global mask and GEOS5 to LIS mapping
    ! ---------------------------------------------------

    allocate (maskarray(1: LDT_rc%lnc(1),1: LDT_rc%lnr(1)))
    call read_clsm_maskfile(nest, maskarray)

    if (associated(LDT_g5map) .eqv. .false.) call init_geos2lis_mapping (nest) 

    ! These 2 values are assumed as same as they are in surfdata_0.23x0.31_simyr2000_c100406.nc

    EDGEW = -180.
    EDGES = -90. 

    ! Reading CLM pft data file
    !--------------------------

    ALLOCATE (PCTPFT      (1:N_lon_clm, 1:N_lat_clm, 1:lsmpft))
     
    status  = NF_OPEN (trim(c_data)//'/CLM45/mksrf_24pftNT_landuse_rc2000_c121207.nc', NF_NOWRITE, ncid) ; VERIFY_(STATUS)  
    status  = NF_INQ_VARID (ncid,'PCT_PFT',VarID) ; VERIFY_(STATUS)

    do k = 1, 25 ! Natural vegetation
       status  = NF_GET_VARA_REAL (ncid,VarID,(/1,1,k/),(/N_lon_clm, N_lat_clm, 1/),PCTPFT(:,:,k)) ; VERIFY_(STATUS)
    end do

    status = NF_CLOSE(ncid)

    ! CLM 4_5 description (25)                                CLM45-carbon description (27)                                    
    ! ------------------------                                ----------------------------- 

    ! 'BARE'   1  	bare                                     (does not have bare soil)
    ! 'NLEt'   2 	needleleaf evergreen temperate tree    1
    ! 'NLEB'   3 	needleleaf evergreen boreal tree       2
    ! 'NLDB'   4  	needleleaf deciduous boreal tree       3
    ! 'BLET'   5 	broadleaf evergreen tropical tree      4
    ! 'BLEt'   6 	broadleaf evergreen temperate tree     5
    ! 'BLDT'   7 	broadleaf deciduous tropical tree      6
    ! 'BLDt'   8 	broadleaf deciduous temperate tree     7
    ! 'BLDB'   9 	broadleaf deciduous boreal tree        8
    ! 'BLEtS' 10 	broadleaf evergreen temperate shrub    9
    ! 'BLDtS' 11 	broadleaf deciduous temperate shrub   10  broadleaf deciduous temperate shrub [moisture +  deciduous]
    ! 'BLDtSm'  	broadleaf deciduous temperate shrub   11  broadleaf deciduous temperate shrub [moisture stress only]
    ! 'BLDBS' 12 	broadleaf deciduous boreal shrub      12
    ! 'AC3G'  13 	arctic c3 grass                       13
    ! 'CC3G'  14 	cool c3 grass                         14  cool c3 grass [moisture +  deciduous]
    ! 'CC3Gm'           cool c3 grass                         15  cool c3 grass [moisture stress only]
    ! 'WC4G'  15 	warm c4 grass                         16  warm c4 grass [moisture +  deciduous]
    ! 'WC4Gm'   	warm c4 grass                         17  warm c4 grass [moisture stress only]
    ! 'C3CROP' 16       c3_crop                               18
    ! 'C3IRR'  17       c3_irrigated                          19
    ! 'CORN'   18       corn                                  20
    ! 'ICORN'  19       irrigated corn                        21
    ! 'STCER'  20       spring temperate cereal               22
    ! 'ISTCER' 21       irrigated spring temperate cereal     23
    ! 'WTCER'  22       winter temperate cereal               24
    ! 'IWTCER' 23       irrigated winter temperate cereal     25
    ! 'SOYB'   24       soybean                               26
    ! 'ISOYB'  25       irrigated soybean                     27
    
!**    ! 'CROP'  16 	crop                                  18  crop [moisture +  deciduous]
!**    ! 'CROPm'   	crop                                  19  crop [moisture stress only]
!**    !         17        water

    dx_clm = 360./N_lon_clm
    dy_clm = 180./N_lat_clm

    do i = 1, N_lon_clm 
       x_min_clm (i) = (i-1)*dx_clm + EDGEW  
    end do

    do i = 1,  N_lat_clm
       y_min_clm (i) = (i-1)*dy_clm  + EDGES
    end do

    ! This data set is DE
    !PCTPFT (1:N_lon_clm/2             ,:,:) =  REAL (PCT_PFT_DBL(N_lon_clm/2 + 1: N_lon_clm,:,:))
    !PCTPFT (N_lon_clm/2 + 1: N_lon_clm,:,:) =  REAL (PCT_PFT_DBL(1:N_lon_clm/2             ,:,:))

    !DEALLOCATE (PCT_PFT_DBL)

    ! Find primary and secondary types in the CLM data file
    ! -----------------------------------------------------

    ! allocate (clm_veg (1:N_lon_clm,1:N_lat_clm,1:2))
    !
    ! do j = 1, N_lat_clm 
    !    do i = 1, N_lon_clm  
    !       if(maxval(PCT_PFT(i,j,:)) > 0.) then
    !          clm_fracs = PCT_PFT(i,j,:)
    !          if (maxval (clm_fracs) == 100.) then 
    !             clm_veg(i,j,:) = maxloc (clm_fracs)
    !          else 
    !             clm_veg(i,j,0)             = maxloc (clm_fracs)		
    !             clm_fracs (clm_veg(i,j,0)) = 0.
    !             clm_veg(i,j,1)             = maxloc (clm_fracs)	 
    !          endif
    !       else 
    !          clm_veg(i,j,:) = 17
    !       endif
    !    end do
    ! end do
    
    ! Reading ESA vegetation types
    !-----------------------------

    allocate (esa_veg (1:nc_esa, 1: nr_esa))
    allocate (lon_esa (1:nc_esa))
    allocate (lat_esa (1:nr_esa))

    status    = NF_OPEN (trim(c_data)//trim(ESA_FILE), NF_NOWRITE, ncid) ; VERIFY_(STATUS)     

    status  = NF_GET_VARA_DOUBLE (ncid,1,(/1/),(/nr_esa/),lat_esa)
    status  = NF_GET_VARA_DOUBLE (ncid,2,(/1/),(/nc_esa/),lon_esa)

    do j = 1,nr_esa
       status  = NF_GET_VARA_INT2 (ncid,3,(/1,j/),(/nc_esa,1/),esa_veg(:,j)) ; VERIFY_(STATUS)   
    end do

    status = NF_CLOSE(ncid)

    ! Find I,J of overlying CLM grid cells for each ESA pixel 
    !--------------------------------------------------------
    allocate (i_esa2clm (1:nc_esa))
    allocate (j_esa2clm (1:nr_esa))

    do i = 1, N_lon_clm 
       where ((real(lon_esa) >=  x_min_clm(i)).and.(real(lon_esa) <  (x_min_clm(i) + dx_clm))) i_esa2clm= i
    end do

    i_esa2clm(129545:nc_esa) = 1

    do j = 1, N_lat_clm 
       where ((real(lat_esa) >=  y_min_clm(j)).and.(real(lat_esa) <  (y_min_clm(j) + dy_clm))) j_esa2clm= j
    end do

    !
    ! Reading number of tiles
    ! -----------------------

    maxcat = LDT_g5map%NT_GEOS 
    
    !
    ! Loop through tile_id raster
    ! ___________________________


    allocate (tile_id (1:nc             ))   
    allocate (clm_veg (1:maxcat,1:lsmpft))
    clm_veg = 0.

    dx = nc_esa / nc
    dy = nr_esa / nr

    open (10,file=trim(gfile)//'.rst',status='old',action='read',  &
          form='unformatted',convert='little_endian')

    do j=1,nr
       
       ! read a row
       
       read(10)tile_id(:)
       
       do i = 1,nc

          ii = i_esa2clm ((i-1)*dx + dx/2)
          jj = j_esa2clm ((j-1)*dy + dy/2)

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
                do n = 2,NPLUS 
                   unq_mask(n) = .not.(loc_int(n) == loc_int(n-1))
                end do
                NBINS = count(unq_mask)
                
                allocate(loc_val (1:NBINS))
                allocate(density (1:NBINS))
                loc_val = 1.*pack(loc_int,mask =unq_mask)
                call histogram (size(subset,1)*size(subset,2), NBINS, density, loc_val, real(subset))   
                
                do k = 1, nbins
                   
                   if (density (k) > 0) then
                      
                      esa_type = int (loc_val(k))
                      
                      ! if (esa_type ==  10)  clm_veg (tile_id(i), 17) = 1.* density(k)   ! lakes inland water
                      
                      if ((esa_type ==  11).or. (esa_type ==  14).or.(esa_type ==  20).or. (esa_type == 190)) then

                         ! ESA type  11: Post-flooding or irrigated croplands 
                         ! ESA type  14: Rainfed croplands 
                         ! ESA type  20: Mosaic Cropland (50-70%) / Vegetation (grassland, shrubland, forest) (20-50%) 
                         ! ESA type 190:	Artificial surfaces and associated areas (urban areas >50%) 

                         if(sum(PCTPFT(ii,jj,16:25)) > 0.) then
                            do n = 16,25                                
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,16:25))
                            end do
                         else
                            clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   
                         endif
                      endif
                      
                      ! if (esa_type == 200)  clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 1.* density(k)   ! ESA type 200:	Bare areas
                      ! if (esa_type == 210)  clm_veg (tile_id(i), 17) = clm_veg (tile_id(i), 17) + 1.* density(k)   ! ocean
                      ! if (esa_type == 220)  clm_veg (tile_id(i), 17) = clm_veg (tile_id(i), 17) + 1.* density(k)   ! ice  
                      ! gkw: bare soil excluded! only considering vegetated land                   
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  30) then 
                         ! ESA type  30: Mosaic Vegetation (grassland, shrubland, forest) (50-70%) / Cropland (20-50%) 

                         if(sum(PCTPFT(ii,jj,16:25)) > 0.) then
                            do n = 16,25                                
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.5*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,16:25))
                            end do
                         elseif(sum(PCTPFT(ii,jj,2:15)) > 0.) then
                            do n = 2,  15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.5* density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:15))
                            enddo
                         else
                            clm_veg (tile_id(i),  16) = clm_veg (tile_id(i),  16) + 1.0* density(k)
                         endif

                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  40) then
                         ! ESA type  40:	Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,5:6)) > 0.) then
                            do n = 5, 6 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,5:6))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  5) = clm_veg (tile_id(i),  5) + 1.0* density(k)
                            else
                               clm_veg (tile_id(i),  6) = clm_veg (tile_id(i),  6) + 1.0* density(k)
                            endif
                         endif
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if ((esa_type ==  50) .or. (esa_type ==  60)) then    
                         ! ESA type  50:	Closed (>40%) broadleaved deciduous forest (>5m) 
                         ! ESA type  60:	Open (15-40%) broadleaved deciduous forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,7:9)) > 0.) then
                            do n = 7, 9 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:9))
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 1.0* density(k)
                            else
                               if(abs(y_min_clm(jj) + 0.5*dy_clm) <  60.) clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 1.0* density(k)
                               if(abs(y_min_clm(jj) + 0.5*dy_clm) >= 60.) clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 1.0* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  70) then    
                         ! ESA type  70:	Closed (>40%) needleleaved evergreen forest (>5m)
                         
                         if(sum(PCTPFT(ii,jj,2:3)) > 0.) then
                            do n = 2, 3 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:3))	
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  90) then 
                         !ESA type  90:	Open (15-40%) needleleaved deciduous or evergreen forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,2:4)) > 0.) then
                            do n = 2, 4 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n)   + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:4))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 100) then 
                         !  ESA type 100:	Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m) 
                         
                         if((sum(PCTPFT(ii,jj,2:4)) + sum(PCTPFT(ii,jj,7:9))) > 0.) then
                            do n = 2, 9 
                               if((n /= 5) .and. (n /= 6)) clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/(sum(PCTPFT(ii,jj,2:4)) + sum(PCTPFT(ii,jj,7:9)))	
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.5* density(k)
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 0.5* density(k)			
                            elseif (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.5* density(k)
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 0.5* density(k)
                            else 
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.5* density(k)
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 0.5* density(k)			
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 110) then   
                         ! ESA type 110:	Mosaic Forest/Shrubland (50-70%) / Grassland (20-50%) 
                         
                         if(sum(PCTPFT(ii,jj,7:12))  > 0.) then
                            do n = 7, 12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) +  0.6*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:12)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.3* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.3* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.3* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.3* density(k)
                            else
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.3* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.3* density(k)			
                            end if
                         end if
                         
                         if(sum(PCTPFT(ii,jj,13:15))  > 0.) then
                            do n =13, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.4*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 0.4* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.4* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.4* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------

                      if (esa_type == 120) then   
                         ! ESA type 120:	Mosaic Grassland (50-70%) / Forest/Shrubland (20-50%) 
                         
                         if(sum(PCTPFT(ii,jj,7:12)) > 0.) then
                            do n = 7, 12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) +  0.4*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:12)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.2* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.2* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.2* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.2* density(k)
                            else 
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.2* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.2* density(k)			
                            end if
                         end if
                         
                         if(sum(PCTPFT(ii,jj,13:15)) > 0.) then
                            do n =13, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.6*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15)) 
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 0.6* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.6* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.6* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 130) then
                         ! 	Closed to open (>15%) shrubland (<5m) 
                         
                         if(sum(PCTPFT(ii,jj,10:12)) > 0.) then
                            do n = 10,12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:12))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  11) = clm_veg (tile_id(i),  11) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  12) = clm_veg (tile_id(i),  12) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 140) then
                         ! ESA type 140:	Closed to open (>15%) grassland 
                         
                         if(sum(PCTPFT(ii,jj,13:15)) > 0.) then
                            do n = 13,15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i),  15) = clm_veg (tile_id(i),  15) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then	
                               clm_veg (tile_id(i),  14) = clm_veg (tile_id(i),  14) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i),  13) = clm_veg (tile_id(i),  13) + 1.0* density(k)	
                            end if
                         end if
                      end if
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 150) then
                         ! ESA type 150:	Sparse (<15%) vegetation (woody vegetation, shrubs, grassland) 
                         
                         if(sum(PCTPFT(ii,jj,10:15)) > 0.) then
                            do n = 10, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.0*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:15)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.5* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.5* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.5* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.5* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.5* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.5* density(k)			
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if((esa_type == 160) .or. (esa_type == 170)) then  
                         ! ESA type 160:	Closed (>40%) broadleaved forest regularly flooded - Fresh water      ! ESA type 170:	Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded
                         
                         if(sum(PCTPFT(ii,jj,5:9)) > 0.) then
                            do n = 5,9 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,5:9)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  5) = clm_veg (tile_id(i),  5) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 180) then
                         ! ESA type 180:	Closed to open (>15%) vegetation (grassland, shrubland, woody vegetation) on regularly flooded or waterlogged soil - Fresh, brackish or saline water 
                         
                         if(sum(PCTPFT(ii,jj,10:15)) > 0.) then
                            do n = 10,15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:15))	
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 1.0* density(k)	
                            end if
                         end if
                      endif
                   endif
                enddo
                deallocate (loc_int,unq_mask,loc_val,density)
             endif
          end if
       enddo
    end do
    
    
    deallocate (tile_id, PCTPFT,esa_veg,lon_esa,lat_esa,i_esa2clm,j_esa2clm)  
    close (10,status='keep')    

    ! Now create CLM-carbon_veg_fracs file
    ! ------------------------------------

    if (write_clsm_files) then
       fclm = LDT_getNextUnitNumber()
       open (fclm,file='LDT_clsm/CLM4.5_veg_typs_fracs',  &
            form='formatted',status='unknown')
    endif

    do k = 1, maxcat

       tile_lat = LDT_g5map%lat(k)
       scale = (ABS (tile_lat) - 32.)/10.
       scale = min (max(scale,0.),1.)

       esa_clm_veg = 0
       esa_clm_frac= 0.

       clm_fracs = clm_veg (k,:)
             
       if (sum (clm_fracs) == 0.) then ! gkw: no vegetated land found; set to BLDtS
          esa_clm_veg (1) = 11              ! broadleaf deciduous shrub 
          esa_clm_frac(1) = 100.
       else
          esa_clm_veg (1) = maxloc(clm_fracs,1)
          esa_clm_frac(1) = maxval(clm_fracs) 
       endif

       clm_fracs (esa_clm_veg (1)) = 0.

       if (sum (clm_fracs) == 0.) then ! gkw: no vegetated secondary type found, set to primary with zero fraction
          esa_clm_veg (2) = esa_clm_veg (1)
          esa_clm_frac(1) = 100.
          esa_clm_frac(2) = 0.
       else
          esa_clm_veg (2) = maxloc(clm_fracs,1)
          esa_clm_frac(1) = 100.*clm_veg (k,esa_clm_veg (1))/(clm_veg (k,esa_clm_veg (1)) + clm_veg (k,esa_clm_veg (2)))
          esa_clm_frac(2) = 100. - esa_clm_frac(1)
       end if

! Now splitting CLM types for CLM-carbon model
! --------------------------------------------
 
! CLM types 2- 10,12,13 are not being splitted.
! .............................................
     
       if ((esa_clm_veg (1) >= 2).and.(esa_clm_veg (1) <= 10)) then
          CPT1 = esa_clm_veg (1) - 1
          CPT2 = esa_clm_veg (1) - 1
          CPF1 = esa_clm_frac(1) 
          CPF2 = 0.
       endif

       if ((esa_clm_veg (2) >= 2).and.(esa_clm_veg (2) <= 10)) then
          CST1 = esa_clm_veg (2) - 1
          CST2 = esa_clm_veg (2) - 1
          CSF1 = esa_clm_frac(2) 
          CSF2 = 0.
       endif

! .............................................

       if ((esa_clm_veg (1) >= 12).and.(esa_clm_veg (1) <= 13)) then
          CPT1 = esa_clm_veg (1)
          CPT2 = esa_clm_veg (1)
          CPF1 = esa_clm_frac(1) 
          CPF2 = 0.
       endif

       if ((esa_clm_veg (2) >= 12).and.(esa_clm_veg (2) <= 13)) then
          CST1 = esa_clm_veg (2)
          CST2 = esa_clm_veg (2)
          CSF1 = esa_clm_frac(2) 
          CSF2 = 0.
       endif

! CLM4_5 crop types - we don't split

       if ((esa_clm_veg (1) >= 16).and.(esa_clm_veg (1) <= 25)) then
          CPT1 = esa_clm_veg (1) + 2
          CPT2 = esa_clm_veg (1) + 2
          CPF1 = esa_clm_frac(1) 
          CPF2 = 0.
       endif

       if ((esa_clm_veg (2) >= 16).and.(esa_clm_veg (2) <= 25)) then
          CST1 = esa_clm_veg (2) + 2
          CST2 = esa_clm_veg (2) + 2
          CSF1 = esa_clm_frac(2) 
          CSF2 = 0.
       endif

! Now splitting (broadleaf deciduous temperate shrub )
! .............

       if (esa_clm_veg (1) == 11) then
          CPT1 = 10
          CPT2 = 11
          CPF1 = esa_clm_frac(1) * scale
          CPF2 = esa_clm_frac(1) * (1. - scale)
       endif

       if (esa_clm_veg (2) == 11) then
          CST1 = 10
          CST2 = 11
          CSF1 = esa_clm_frac(2) * scale       
          CSF2 = esa_clm_frac(2) * (1. - scale) 
       endif

! ............. (cool c3 grass)

       if (esa_clm_veg (1) == 14) then
          CPT1 = 14
          CPT2 = 15
          CPF1 = esa_clm_frac(1) * scale        
          CPF2 = esa_clm_frac(1) * (1. - scale) 
       endif

       if (esa_clm_veg (2) == 14) then
          CST1 = 14
          CST2 = 15
          CSF1 = esa_clm_frac(2) * scale        
          CSF2 = esa_clm_frac(2) * (1. - scale) 
       endif

! ............. warm c4 grass

       if (esa_clm_veg (1) == 15) then
          CPT1 = 16
          CPT2 = 17
          CPF1 = esa_clm_frac(1) * scale        
          CPF2 = esa_clm_frac(1) * (1. - scale) 
       endif

       if (esa_clm_veg (2) == 15) then
          CST1 = 16
          CST2 = 17
          CSF1 = esa_clm_frac(2) * scale        
          CSF2 = esa_clm_frac(2) * (1. - scale)
       endif
! .............
! CLM_4.5 : we don't splot crop type anymore 16 has become 16-25 and they are now 18-27 in catchment-CN
!       if (esa_clm_veg (1) == 16) then
!          CPT1 = 18
!          CPT2 = 19
!          CPF1 = esa_clm_frac(1) * scale        
!          CPF2 = esa_clm_frac(1) * (1. - scale) 
!       endif
!
!       if (esa_clm_veg (2) == 16) then
!          CST1 = 18
!          CST2 = 19
!          CSF1 = esa_clm_frac(2) * scale        
!          CSF2 = esa_clm_frac(2) * (1. - scale) 
!       endif

       ! fractions must sum to 1
       ! -----------------------
       ftot = cpf1 + cpf2 + csf1 + csf2

       if(ftot /= 100.) then
          cpf1 = 100. * cpf1 / ftot  
          cpf2 = 100. * cpf2 / ftot 
          csf1 = 100. * csf1 / ftot 
          csf2 = 100. * csf2 / ftot 
       endif
    
      if (write_clsm_files)  write (fclm,'(2I8,4I3,4f7.2,2I3,2f7.2)')     &
            tid,cid,cpt1, cpt2, cst1, cst2, cpf1, cpf2, csf1, csf2, &
            esa_clm_veg (1), esa_clm_veg (2), esa_clm_frac(1), esa_clm_frac(2)
    end do
    if (write_clsm_files) then
       close (fclm, status = 'keep')
       call LDT_releaseUnitNumber(fclm)
    endif

  END SUBROUTINE ESA2CLM45

! ------------------------------------------------------------------------------------------------

  SUBROUTINE ESA2CLM4 (nest, nc, nr, gfile)

    implicit none

    integer  , intent (in) :: nc, nr, nest
    character (*)          :: gfile
    
    integer  , parameter   :: N_lon_clm = 1152, N_lat_clm = 768, lsmpft = 17
    integer*2, allocatable, target, dimension (:,:) :: esa_veg
    integer*2, pointer    , dimension (:,:) :: subset
    integer  , allocatable, dimension (:)   :: tile_id, i_esa2clm, j_esa2clm
    integer :: i,j, k,n, status, ncid, varid, maxcat, dx,dy, esa_type, tid, cid, ii, jj   
    real    :: dx_clm, dy_clm, x_min_clm (N_lon_clm), y_min_clm (N_lat_clm), clm_fracs(lsmpft)
    real    :: minlon,maxlon,minlat,maxlat,tile_lat, scale, ftot
    integer :: cpt1, cpt2, cst1, cst2  ! CLM-carbon types
    real    :: cpf1, cpf2, csf1, csf2  ! CLM-carbon fractions
    DOUBLE PRECISION,  allocatable, dimension (:) :: lon_esa, lat_esa
    DOUBLE PRECISION       :: EDGEN, EDGEE, EDGES, EDGEW 
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: PCT_PFT_DBL
    REAL, ALLOCATABLE, DIMENSION (:,:,:) :: PCTPFT 
    integer, allocatable, dimension (:) :: density, loc_int
    real   , allocatable, dimension (:) :: loc_val
    logical, allocatable, dimension (:) :: unq_mask
    integer :: NBINS, NPLUS, fclm
    integer, allocatable, dimension (:,:) :: clm_veg
    integer :: esa_clm_veg (2)
    real    :: esa_clm_frac(2)
    real, dimension (:,:), allocatable      :: maskarray

    ! Initialize the global mask and GEOS5 to LIS mapping
    ! ---------------------------------------------------

!    allocate (maskarray(1: LDT_rc%lnc(1),1: LDT_rc%lnr(1)))
    call read_clsm_maskfile(nest, maskarray)

    if (associated(LDT_g5map) .eqv. .false.) call init_geos2lis_mapping (nest) 

    ! Reading CLM pft data file
    !--------------------------

    ALLOCATE (PCTPFT      (1:N_lon_clm, 1:N_lat_clm, 1:lsmpft))
    ALLOCATE (PCT_PFT_DBL (1:N_lon_clm, 1:N_lat_clm, 1:lsmpft))
    status  = NF_OPEN (trim(c_data)//'/surfdata_0.23x0.31_simyr2000_c100406.nc', NF_NOWRITE, ncid) ; VERIFY_(STATUS)   
    status  = NF_GET_VARA_DOUBLE (ncid,1,(/1/),(/1/),EDGEN) ; VERIFY_(STATUS)
    status  = NF_GET_VARA_DOUBLE (ncid,2,(/1/),(/1/),EDGEE) ; VERIFY_(STATUS)
    status  = NF_GET_VARA_DOUBLE (ncid,3,(/1/),(/1/),EDGES) ; VERIFY_(STATUS)
    status  = NF_GET_VARA_DOUBLE (ncid,4,(/1/),(/1/),EDGEW) ; VERIFY_(STATUS)
    status  = NF_INQ_VARID (ncid,'PCT_PFT',VarID) ; VERIFY_(STATUS)

    do k = 1, lsmpft
       status  = NF_GET_VARA_DOUBLE (ncid,VarID,(/1,1,k/),(/N_lon_clm, N_lat_clm, 1/),PCT_PFT_DBL(:,:,k)) ; VERIFY_(STATUS)
    end do

    status = NF_CLOSE(ncid)

    ! change type 6 to 10 for Australia only gkw: to remove CLM artificial tree line, and stay true to ESA
    ! ----------------------------------------------------------------------------------------------------

    PCT_PFT_DBL(360:494,215:341,11) = PCT_PFT_DBL(360:494,215:341,11) + PCT_PFT_DBL(360:494,215:341, 7)
    PCT_PFT_DBL(360:494,215:341, 7) = 0.

    ! CLM description (17)                                       CLM-carbon description (19)                                    
    ! --------------------                                        -------------------------- 

    ! 'BARE'   1  	bare                                     (does not have bare soil)
    ! 'NLEt'   2 	needleleaf evergreen temperate tree    1
    ! 'NLEB'   3 	needleleaf evergreen boreal tree       2
    ! 'NLDB'   4  	needleleaf deciduous boreal tree       3
    ! 'BLET'   5 	broadleaf evergreen tropical tree      4
    ! 'BLEt'   6 	broadleaf evergreen temperate tree     5
    ! 'BLDT'   7 	broadleaf deciduous tropical tree      6
    ! 'BLDt'   8 	broadleaf deciduous temperate tree     7
    ! 'BLDB'   9 	broadleaf deciduous boreal tree        8
    ! 'BLEtS' 10 	broadleaf evergreen temperate shrub    9
    ! 'BLDtS' 11 	broadleaf deciduous temperate shrub   10  broadleaf deciduous temperate shrub [moisture +  deciduous]
    ! 'BLDtSm'  	broadleaf deciduous temperate shrub   11  broadleaf deciduous temperate shrub [moisture stress only]
    ! 'BLDBS' 12 	broadleaf deciduous boreal shrub      12
    ! 'AC3G'  13 	arctic c3 grass                       13
    ! 'CC3G'  14 	cool c3 grass                         14  cool c3 grass [moisture +  deciduous]
    ! 'CC3Gm'           cool c3 grass                         15  cool c3 grass [moisture stress only]
    ! 'WC4G'  15 	warm c4 grass                         16
    ! 'WC4Gm'   	warm c4 grass                         17
    ! 'CROP'  16 	crop                                  18  crop [moisture +  deciduous]
    ! 'CROPm'   	crop                                  19  crop [moisture stress only]
    !         17        water

    dx_clm = 360./N_lon_clm
    dy_clm = 180./N_lat_clm

    do i = 1, N_lon_clm 
       x_min_clm (i) = (i-1)*dx_clm + EDGEW - 180.
    end do

    do i = 1,  N_lat_clm
       y_min_clm (i) = (i-1)*dy_clm  + EDGES
    end do

    PCTPFT (1:N_lon_clm/2             ,:,:) =  REAL (PCT_PFT_DBL(N_lon_clm/2 + 1: N_lon_clm,:,:))
    PCTPFT (N_lon_clm/2 + 1: N_lon_clm,:,:) =  REAL (PCT_PFT_DBL(1:N_lon_clm/2             ,:,:))

    DEALLOCATE (PCT_PFT_DBL)

    ! Find primary and secondary types in the CLM data file
    ! -----------------------------------------------------

    ! allocate (clm_veg (1:N_lon_clm,1:N_lat_clm,1:2))
    !
    ! do j = 1, N_lat_clm 
    !    do i = 1, N_lon_clm  
    !       if(maxval(PCT_PFT(i,j,:)) > 0.) then
    !          clm_fracs = PCT_PFT(i,j,:)
    !          if (maxval (clm_fracs) == 100.) then 
    !             clm_veg(i,j,:) = maxloc (clm_fracs)
    !          else 
    !             clm_veg(i,j,0)             = maxloc (clm_fracs)		
    !             clm_fracs (clm_veg(i,j,0)) = 0.
    !             clm_veg(i,j,1)             = maxloc (clm_fracs)	 
    !          endif
    !       else 
    !          clm_veg(i,j,:) = 17
    !       endif
    !    end do
    ! end do
    
    ! Reading ESA vegetation types
    !-----------------------------

    allocate (esa_veg (1:nc_esa, 1: nr_esa))
    allocate (lon_esa (1:nc_esa))
    allocate (lat_esa (1:nr_esa))

    status    = NF_OPEN (trim(c_data)//trim(ESA_FILE), NF_NOWRITE, ncid) ; VERIFY_(STATUS)   

    status  = NF_GET_VARA_DOUBLE (ncid,1,(/1/),(/nr_esa/),lat_esa) ; VERIFY_(STATUS)   
    status  = NF_GET_VARA_DOUBLE (ncid,2,(/1/),(/nc_esa/),lon_esa) ; VERIFY_(STATUS)   

    do j = 1,nr_esa
       status  = NF_GET_VARA_INT2 (ncid,3,(/1,j/),(/nc_esa,1/),esa_veg(:,j))
    end do

    status = NF_CLOSE(ncid)

    ! Find I,J of overlying CLM grid cells for each ESA pixel 
    !--------------------------------------------------------
    allocate (i_esa2clm (1:nc_esa))
    allocate (j_esa2clm (1:nr_esa))

    do i = 1, N_lon_clm 
       where ((real(lon_esa) >=  x_min_clm(i)).and.(real(lon_esa) <  (x_min_clm(i) + dx_clm))) i_esa2clm= i
    end do

    i_esa2clm(129545:nc_esa) = 1

    do j = 1, N_lat_clm 
       where ((real(lat_esa) >=  y_min_clm(j)).and.(real(lat_esa) <  (y_min_clm(j) + dy_clm))) j_esa2clm= j
    end do

    !
    ! Reading number of tiles
    ! -----------------------

    maxcat = LDT_g5map%NT_GEOS 
    
    ! Loop through tile_id raster
    ! ___________________________


    allocate (tile_id (1:nc             ))   
    allocate (clm_veg (1:maxcat,1:lsmpft))
    clm_veg = 0.

    dx = nc_esa / nc
    dy = nr_esa / nr

    open (10,file=trim(gfile)//'.rst',status='old',action='read',  &
          form='unformatted',convert='little_endian')

    do j=1,nr
       
       ! read a row
       
       read(10)tile_id(:)
       
       do i = 1,nc

          ii = i_esa2clm ((i-1)*dx + dx/2)
          jj = j_esa2clm ((j-1)*dy + dy/2)

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
                do n = 2,NPLUS 
                   unq_mask(n) = .not.(loc_int(n) == loc_int(n-1))
                end do
                NBINS = count(unq_mask)
                
                allocate(loc_val (1:NBINS))
                allocate(density (1:NBINS))
                loc_val = 1.*pack(loc_int,mask =unq_mask)
                call histogram (size(subset,1)*size(subset,2), NBINS, density, loc_val, real(subset))   
                
                do k = 1, nbins
                   
                   if (density (k) > 0) then
                      
                      esa_type = int (loc_val(k))
                      
                      ! if (esa_type ==  10)  clm_veg (tile_id(i), 17) = 1.* density(k)   ! lakes inland water
                      
                      if (esa_type ==  11)  clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   ! ESA type  11: Post-flooding or irrigated croplands 
                      if (esa_type ==  14)  clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   ! ESA type  14: Rainfed croplands 
                      if (esa_type ==  20)  clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   ! ESA type  20: Mosaic Cropland (50-70%) / Vegetation (grassland, shrubland, forest) (20-50%) 
                      if (esa_type == 190)  clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   ! ESA type 190:	Artificial surfaces and associated areas (urban areas >50%) 
                      
                      ! if (esa_type == 200)  clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 1.* density(k)   ! ESA type 200:	Bare areas
                      ! if (esa_type == 210)  clm_veg (tile_id(i), 17) = clm_veg (tile_id(i), 17) + 1.* density(k)   ! ocean
                      ! if (esa_type == 220)  clm_veg (tile_id(i), 17) = clm_veg (tile_id(i), 17) + 1.* density(k)   ! ice  
                      ! gkw: bare soil excluded! only considering vegetated land                   
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  30) then 
                         ! ESA type  30: Mosaic Vegetation (grassland, shrubland, forest) (50-70%) / Cropland (20-50%) 
                         clm_veg (tile_id(i),  16) = clm_veg (tile_id(i),  16) + 0.5* density(k)
                         if(sum(PCTPFT(ii,jj,2:15)) > 0.) then
                            do n = 2,  15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.5* density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:15))
                            enddo
                         else 
                            clm_veg (tile_id(i),  16) = clm_veg (tile_id(i),  16) + 1.0* density(k)
                         endif
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  40) then
                         ! ESA type  40:	Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,5:6)) > 0.) then
                            do n = 5, 6 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,5:6))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  5) = clm_veg (tile_id(i),  5) + 1.0* density(k)
                            else
                               clm_veg (tile_id(i),  6) = clm_veg (tile_id(i),  6) + 1.0* density(k)
                            endif
                         endif
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if ((esa_type ==  50) .or. (esa_type ==  60)) then    
                         ! ESA type  50:	Closed (>40%) broadleaved deciduous forest (>5m) 
                         ! ESA type  60:	Open (15-40%) broadleaved deciduous forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,7:9)) > 0.) then
                            do n = 7, 9 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:9))
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 1.0* density(k)
                            else
                               if(abs(y_min_clm(jj) + 0.5*dy_clm) <  60.) clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 1.0* density(k)
                               if(abs(y_min_clm(jj) + 0.5*dy_clm) >= 60.) clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 1.0* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  70) then    
                         ! ESA type  70:	Closed (>40%) needleleaved evergreen forest (>5m)
                         
                         if(sum(PCTPFT(ii,jj,2:3)) > 0.) then
                            do n = 2, 3 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:3))	
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  90) then 
                         !ESA type  90:	Open (15-40%) needleleaved deciduous or evergreen forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,2:4)) > 0.) then
                            do n = 2, 4 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n)   + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:4))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 100) then 
                         !  ESA type 100:	Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m) 
                         
                         if((sum(PCTPFT(ii,jj,2:4)) + sum(PCTPFT(ii,jj,7:9))) > 0.) then
                            do n = 2, 9 
                               if((n /= 5) .and. (n /= 6)) clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/(sum(PCTPFT(ii,jj,2:4)) + sum(PCTPFT(ii,jj,7:9)))	
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.5* density(k)
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 0.5* density(k)			
                            elseif (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.5* density(k)
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 0.5* density(k)
                            else 
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.5* density(k)
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 0.5* density(k)			
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 110) then   
                         ! ESA type 110:	Mosaic Forest/Shrubland (50-70%) / Grassland (20-50%) 
                         
                         if(sum(PCTPFT(ii,jj,7:12))  > 0.) then
                            do n = 7, 12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) +  0.6*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:12)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.3* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.3* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.3* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.3* density(k)
                            else
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.3* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.3* density(k)			
                            end if
                         end if
                         
                         if(sum(PCTPFT(ii,jj,13:15))  > 0.) then
                            do n =13, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.4*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 0.4* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.4* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.4* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------

                      if (esa_type == 120) then   
                         ! ESA type 120:	Mosaic Grassland (50-70%) / Forest/Shrubland (20-50%) 
                         
                         if(sum(PCTPFT(ii,jj,7:12)) > 0.) then
                            do n = 7, 12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) +  0.4*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:12)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.2* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.2* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.2* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.2* density(k)
                            else 
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.2* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.2* density(k)			
                            end if
                         end if
                         
                         if(sum(PCTPFT(ii,jj,13:15)) > 0.) then
                            do n =13, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.6*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15)) 
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 0.6* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.6* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.6* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 130) then
                         ! 	Closed to open (>15%) shrubland (<5m) 
                         
                         if(sum(PCTPFT(ii,jj,10:12)) > 0.) then
                            do n = 10,12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:12))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  11) = clm_veg (tile_id(i),  11) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  12) = clm_veg (tile_id(i),  12) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 140) then
                         ! ESA type 140:	Closed to open (>15%) grassland 
                         
                         if(sum(PCTPFT(ii,jj,13:15)) > 0.) then
                            do n = 13,15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i),  15) = clm_veg (tile_id(i),  15) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then	
                               clm_veg (tile_id(i),  14) = clm_veg (tile_id(i),  14) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i),  13) = clm_veg (tile_id(i),  13) + 1.0* density(k)	
                            end if
                         end if
                      end if
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 150) then
                         ! ESA type 150:	Sparse (<15%) vegetation (woody vegetation, shrubs, grassland) 
                         
                         if(sum(PCTPFT(ii,jj,10:15)) > 0.) then
                            do n = 10, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.0*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:15)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.5* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.5* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.5* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.5* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.5* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.5* density(k)			
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if((esa_type == 160) .or. (esa_type == 170)) then  
                         ! ESA type 160:	Closed (>40%) broadleaved forest regularly flooded - Fresh water      ! ESA type 170:	Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded
                         
                         if(sum(PCTPFT(ii,jj,5:9)) > 0.) then
                            do n = 5,9 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,5:9)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  5) = clm_veg (tile_id(i),  5) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 180) then
                         ! ESA type 180:	Closed to open (>15%) vegetation (grassland, shrubland, woody vegetation) on regularly flooded or waterlogged soil - Fresh, brackish or saline water 
                         
                         if(sum(PCTPFT(ii,jj,10:15)) > 0.) then
                            do n = 10,15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:15))	
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 1.0* density(k)	
                            end if
                         end if
                      endif
                   endif
                enddo
                deallocate (loc_int,unq_mask,loc_val,density)
             endif
          end if
       enddo
    end do
    
    
    deallocate (tile_id, PCTPFT,esa_veg,lon_esa,lat_esa,i_esa2clm,j_esa2clm)  
    close (10,status='keep')    

    ! Now create CLM-carbon_veg_fracs file
    ! ------------------------------------

    if (write_clsm_files) then
       fclm = LDT_getNextUnitNumber()
       open (fclm,file='LDT_clsm/CLM_veg_typs_fracs',  &
            form='formatted',status='unknown')
    endif

    do k = 1, maxcat

       tile_lat = LDT_g5map%lat(k)
       scale = (ABS (tile_lat) - 32.)/10.
       scale = min (max(scale,0.),1.)

       esa_clm_veg = 0
       esa_clm_frac= 0.

       clm_fracs = clm_veg (k,:)
             
       if (sum (clm_fracs) == 0.) then ! gkw: no vegetated land found; set to BLDtS
          esa_clm_veg (1) = 11              ! broadleaf deciduous shrub 
          esa_clm_frac(1) = 100.
       else
          esa_clm_veg (1) = maxloc(clm_fracs,1)
          esa_clm_frac(1) = maxval(clm_fracs) 
       endif

       clm_fracs (esa_clm_veg (1)) = 0.

       if (sum (clm_fracs) == 0.) then ! gkw: no vegetated secondary type found, set to primary with zero fraction
          esa_clm_veg (2) = esa_clm_veg (1)
          esa_clm_frac(1) = 100.
          esa_clm_frac(2) = 0.
       else
          esa_clm_veg (2) = maxloc(clm_fracs,1)
          esa_clm_frac(1) = 100.*clm_veg (k,esa_clm_veg (1))/(clm_veg (k,esa_clm_veg (1)) + clm_veg (k,esa_clm_veg (2)))
          esa_clm_frac(2) = 100. - esa_clm_frac(1)
       end if

! Now splitting CLM types for CLM-carbon model
! --------------------------------------------
 
! CLM types 2- 10,12,13 are not being splitted.
! .............................................
     
       if ((esa_clm_veg (1) >= 2).and.(esa_clm_veg (1) <= 10)) then
          CPT1 = esa_clm_veg (1) - 1
          CPT2 = esa_clm_veg (1) - 1
          CPF1 = esa_clm_frac(1) 
          CPF2 = 0.
       endif

       if ((esa_clm_veg (2) >= 2).and.(esa_clm_veg (2) <= 10)) then
          CST1 = esa_clm_veg (2) - 1
          CST2 = esa_clm_veg (2) - 1
          CSF1 = esa_clm_frac(2) 
          CSF2 = 0.
       endif

! .............................................

       if ((esa_clm_veg (1) >= 12).and.(esa_clm_veg (1) <= 13)) then
          CPT1 = esa_clm_veg (1)
          CPT2 = esa_clm_veg (1)
          CPF1 = esa_clm_frac(1) 
          CPF2 = 0.
       endif

       if ((esa_clm_veg (2) >= 12).and.(esa_clm_veg (2) <= 13)) then
          CST1 = esa_clm_veg (2)
          CST2 = esa_clm_veg (2)
          CSF1 = esa_clm_frac(2) 
          CSF2 = 0.
       endif

! Now splitting
! .............

       if (esa_clm_veg (1) == 11) then
          CPT1 = 10
          CPT2 = 11
          CPF1 = esa_clm_frac(1) * scale
          CPF2 = esa_clm_frac(1) * (1. - scale)
       endif

       if (esa_clm_veg (2) == 11) then
          CST1 = 10
          CST2 = 11
          CSF1 = esa_clm_frac(2) * scale       
          CSF2 = esa_clm_frac(2) * (1. - scale) 
       endif

! .............

       if (esa_clm_veg (1) == 14) then
          CPT1 = 14
          CPT2 = 15
          CPF1 = esa_clm_frac(1) * scale        
          CPF2 = esa_clm_frac(1) * (1. - scale) 
       endif

       if (esa_clm_veg (2) == 14) then
          CST1 = 14
          CST2 = 15
          CSF1 = esa_clm_frac(2) * scale        
          CSF2 = esa_clm_frac(2) * (1. - scale) 
       endif

! .............

       if (esa_clm_veg (1) == 15) then
          CPT1 = 16
          CPT2 = 17
          CPF1 = esa_clm_frac(1) * scale        
          CPF2 = esa_clm_frac(1) * (1. - scale) 
       endif

       if (esa_clm_veg (2) == 15) then
          CST1 = 16
          CST2 = 17
          CSF1 = esa_clm_frac(2) * scale        
          CSF2 = esa_clm_frac(2) * (1. - scale)
       endif
! .............

       if (esa_clm_veg (1) == 16) then
          CPT1 = 18
          CPT2 = 19
          CPF1 = esa_clm_frac(1) * scale        
          CPF2 = esa_clm_frac(1) * (1. - scale) 
       endif

       if (esa_clm_veg (2) == 16) then
          CST1 = 18
          CST2 = 19
          CSF1 = esa_clm_frac(2) * scale        
          CSF2 = esa_clm_frac(2) * (1. - scale) 
       endif

       ! fractions must sum to 1
       ! -----------------------
       ftot = cpf1 + cpf2 + csf1 + csf2

       if(ftot /= 100.) then
          cpf1 = 100. * cpf1 / ftot  
          cpf2 = 100. * cpf2 / ftot 
          csf1 = 100. * csf1 / ftot 
          csf2 = 100. * csf2 / ftot 
       endif
       
       if (write_clsm_files) write (fclm,'(2I8,4I3,4f7.2,2I3,2f7.2)')     &
            tid,cid,cpt1, cpt2, cst1, cst2, cpf1, cpf2, csf1, csf2, &
            esa_clm_veg (1), esa_clm_veg (2), esa_clm_frac(1), esa_clm_frac(2)

    end do

    if (write_clsm_files) then
       close (fclm, status = 'keep')
    endif

  END SUBROUTINE ESA2CLM4
!! ---------------------------------------------------------------
!
!!BOP
!!
!! !ROUTINE: read_clsm_maskfile
!!  \label{read_clsm_maskfile}
!!
!! !REVISION HISTORY:
!!  03 Sept 2004: Sujay Kumar; Initial Specification
!!  01 June 2012: KR Arsenault; Restructured to simply read in a mask file
!!
!! !INTERFACE:
! subroutine read_clsm_maskfile( n, localmask )
!
!! !USES:
!  use LDT_coreMod, only : LDT_rc, LDT_localPet
!  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
!                          LDT_releaseUnitNumber, LDT_endrun
!  use LDT_gridmappingMod
!
!  implicit none
!
!! !ARGUMENTS: 
!  integer, intent(in)  :: n
!  real,    intent(out) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
!! !DESCRIPTION:
!!  This subroutine reads the landmask data and returns the 
!!   mask and surface type arrays.
!!
!!  The arguments are:
!!  \begin{description}
!!   \item[n]
!!    index of nest
!!   \item[localmask]
!!    landmask for the region of interest
!!   \end{description}
!!
!!EOP      
!  integer :: ftn, ios1
!  logical :: file_exists
!  integer :: c, r, t, line
!  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
!  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
!  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
!
!  integer, allocatable :: lat_line(:,:)
!  integer, allocatable :: lon_line(:,:)
!  real,    allocatable :: read_inputparm(:,:)  ! Read input parameter
!!_________________________________________________________________________________
!
!   LDT_rc%nmaskpts = 0.
!   localmask = 0.
!
!!- Check for and open landmask file:
!   inquire(file=trim(LDT_rc%mfile(n)), exist=file_exists)
!   if( file_exists ) then 
!      write(LDT_logunit,*)'[INFO] Reading CLSM mask file:',trim(LDT_rc%mfile(n)), & 
!                          ' (',LDT_localPet,')'
!
!! -------------------------------------------------------------------
!!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
!! -------------------------------------------------------------------
! !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
!    subparam_gridDesc = 0.
!
!    call LDT_RunDomainPts( n, LDT_rc%mask_proj, LDT_rc%mask_gridDesc(n,:), &
!                  glpnc, glpnr, subpnc, subpnr,  &
!                  subparam_gridDesc, lat_line, lon_line )
!
!    allocate( read_inputparm(subpnc, subpnr) )
!    read_inputparm = 0.
!
! ! -------------------------------------------------------------------
!
! ! -- Open land/water mask file:
!      ftn = LDT_getNextUnitNumber()
!      open(ftn, file=trim(LDT_rc%mfile(n)),form='unformatted', recl=4, &
!           access='direct', iostat=ios1)
!
! ! == (1) READ IN GLOBAL/ENTIRE MASK PARAMETER DATA: == 
! !     (Currently important for reading in other CLSM F2.5 parameters)
!
!    ! Global mask field for current CLSM F2.5 needs:
!      allocate(LDT_rc%global_mask(glpnc,glpnr))  !  Temporary ...
!      LDT_rc%global_mask = LDT_rc%udef
!      line = 0
!      do r = 1, glpnr 
!         do c = 1, glpnc
!            line = line + 1
!            read(ftn,rec=line)  LDT_rc%global_mask(c,r)
!         enddo
!      enddo
! !  ( to be removed later after full CLSM-F2.5 preprocesser
! !    is implemented into LDT)
! ! =======================================================
!
! ! == (2) READ IN LANDMASK DATA: == 
!      line = 0
!      do r = 1, subpnr
!         do c = 1, subpnc
!            line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
!            read(ftn,rec=line) read_inputparm(c,r)
!         enddo
!      enddo
!      localmask = read_inputparm
!
! ! == (3) INCLUDE WATER POINTS, IF SELECTED: == 
!
!      if( LDT_rc%inc_water_pts ) then
!         write(*,*) " FOR CLSM F2.5, PLEASE DO NOT INCLUDE WATER PTS "
!         write(*,*) "  AT THIS TIME .... stopping."
!         call LDT_endrun
!
!         do r = 1, glpnr
!            do c = 1, glpnc
!               if( LDT_rc%global_mask(c,r) == 0 .or. &
!                   LDT_rc%global_mask(c,r) == LDT_rc%waterclass ) then
!                 LDT_rc%global_mask(c,r) = 1 
!               endif
!            enddo
!         enddo
!         do r = 1, LDT_rc%lnr(n)
!            do c = 1, LDT_rc%lnc(n)
!               if( localmask(c,r) == 0 .or. &
!                   localmask(c,r) == LDT_rc%waterclass ) then
!                 localmask(c,r) = 1 
!               endif
!            end do
!         end do
!      end if
!
!   !- Generate total number of accounted mask points:
!! - Eventually will handle subsetted mask:
!!      do r=1,LDT_rc%lnr(n)
!!         do c=1,LDT_rc%lnc(n)
!!             if( localmask(c,r) >= 1 ) then
!! - Set for global mask for now:
!      do r=1,glpnr
!         do c=1,glpnc
!            if( LDT_rc%global_mask(c,r) >= 1 ) then
!                LDT_rc%nmaskpts(n) = LDT_rc%nmaskpts(n) + 1
!            endif
!         end do 
!      end do
!
!      call LDT_releaseUnitNumber(ftn)
!
!   else
!      write(LDT_logunit,*) "Landmask map: ",trim(LDT_rc%mfile(n)), " does not exist"
!      write(LDT_logunit,*) "program stopping ..."
!      call LDT_endrun
!   endif
!
!end subroutine read_clsm_maskfile
!
end module mod_ESA_lc

