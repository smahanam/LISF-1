#define VERIFY_(A) if(A /=0)then;print *,'ERROR code',A,'at',__LINE__;call exit(3);endif

!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_HWSD-STATSGO2_texture
! \label{read_HWSD-STATSGO2_texture}
!
! !REVISION HISTORY:
!  22 Sept 2019: Sarith Mahanama and Kristi Arsenault ; Initial Specification
!
! !INTERFACE:

module mod_HWSD_STATSGO2_texture

! !USES:
  use LDT_coreMod,  only : LDT_rc, LDT_domain ! EMK
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_numericalMethodsMod, only : LDT_quicksort
  use LDT_constantsMod, ONLY:      &
       RADIUS => LDT_CONST_REARTH, &
       PI => LDT_CONST_PI
  use catch_util, only : LDT_RegridRaster, NC_VarID, GEOS2LIS
  use get_DeLannoy_SoilClass, ONLY :  &
       mineral_perc, GDL_center_pix, &
       n_SoilClasses => n_DeLannoy_classes, & 
       soil_class => DeLannoy_class

  implicit none
  include 'netcdf.inc'	

  private
  public  read_HWSD_STATSGO2_texture

  character*300 :: c_data = '/discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001/'
  
contains

  subroutine read_HWSD_STATSGO2_texture ( n, num_bins, fgrd, texture_layers )

    implicit none
    
    integer,intent(in)    :: n
    integer,intent(in)    :: num_bins   ! Number of soil types
    real,   intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
    real,   intent(inout) :: texture_layers(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
    !
    ! !DESCRIPTION:
    !  This module is good to  process NGDC-HWSD-STATSGO merged soil properties 
    !     and produce clay%, sand%, silt% and OC% and produce De Lannoy et al. (2014)
    !     soil types on LIS grid cells.
    !     De Lannoy, G. J., R. Koster, R. H. Reichle, S. Mahanama, and Q. Liu. 2014. 
    !          "An updated treatment of soil texture and associated hydraulic properties in a global 
    !          land modeling system." J. Adv. Model. Earth Sys, 6: 957-979 [10.1002/2014MS000330]
    !
    !  The 30 sec soil data starts at: latitude -90. and longitude -180.
    !
    !  The arguments are:
    !  \begin{description}
    !  \item[n]
    !   index of the nest
    !  \item[fgrd]
    !   output field with the retrieved soil texture
    !  \end{description}
    !EOP
    
    integer, parameter :: input_cols = 360*120     ! 43200
    integer, parameter :: input_rows = 180*120     ! 21600
    real,    parameter :: input_xres = 1.0/120.0
    real,    parameter :: input_yres = 1.0/120.0
    
    integer   :: c, r, i, t 
    integer   :: ftn, length, nrec
    logical   :: file_exists
    integer   :: mi                     ! Total number of input param grid fgrd points
    integer   :: mo                     ! Total number of output LIS grid fgrd points
    integer   :: glpnc, glpnr           ! Parameter (global) total columns and rows
    integer   :: subpnc, subpnr         ! Parameter subsetted columns and rows
    real      :: param_gridDesc(20)     ! Input parameter grid desc fgrd
    real      :: subparam_gridDesc(20)  ! Input parameter grid desc fgrd
    integer, allocatable  :: lat_line(:,:), lon_line(:,:)
    real,     allocatable :: gi(:)      ! Input parameter 1d grid
    logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)
    real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output lis 1d grid
    logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output logical mask (to match go)
    real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins) ! Output lis 1d grid
    logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins) ! Output logical mask (to match go)
    
    real      :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real      :: gridcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
    integer*1, allocatable :: soil_int(:)
    real,      allocatable :: input_soiltext(:,:)
    ! ___________________________________________________________________
    
    temp = LDT_rc%udef
    gridcnt = 0.
    fgrd = 0.
    texture_layers = 0.
    
    !- Set parameter grid fgrd inputs:
    param_gridDesc(1)  = 0.             ! Latlon
    param_gridDesc(2)  = input_cols
    param_gridDesc(3)  = input_rows
    param_gridDesc(4)  = -90.0  + (input_yres/2) ! LL lat
    param_gridDesc(5)  = -180.0 + (input_xres/2) ! LL lon
    param_gridDesc(6)  = 128
    param_gridDesc(7)  =  90.0 - (input_yres/2)  ! UR lat
    param_gridDesc(8)  = 180.0 - (input_xres/2)  ! UR lon
    param_gridDesc(9)  = input_yres     ! dy: 0.0083333
    param_gridDesc(10) = input_xres     ! dx: 0.0083333
    param_gridDesc(20) = 64
    
    inquire(file=trim(LDT_rc%txtfile(n)), exist=file_exists)
    if(.not.file_exists) then 
       write(LDT_logunit,*) "[ERR] Native STATSGO+FAO Soil texture map ",&
            trim(LDT_rc%txtfile(n))," not found."
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
    endif
    write(unit=LDT_logunit,fmt=*) "[INFO] Reading Native STATSGO+FAO texture file: ",&
         trim(LDT_rc%txtfile(n))
    
    !- Open texture file:
    ftn = LDT_getNextUnitNumber()
    open(ftn,file=LDT_rc%txtfile(n),form='unformatted',status='old',&
         access='direct',recl=input_cols)  ! along longitude
    
    ! -------------------------------------------------------------------
    !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
    ! -------------------------------------------------------------------
    
    !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
    subparam_gridDesc = 0.
    call LDT_RunDomainPts( n, LDT_rc%soiltext_proj, param_gridDesc(:), &
         glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )
    
    ! _____________________________________
    
    ! nrec: from north to south (max 21600)
    ! irec: from west to east   (max 43200)
    allocate( soil_int(input_cols) )
    allocate( input_soiltext(subpnc, subpnr) )
!    input_soiltext = water_class
    
    !- Reverse-Y and read in values:
    nrec = 0
    do r = subpnr, 1, -1
       nrec = input_rows - lat_line(1,r) + 1
       read(ftn,rec=nrec) soil_int
       do c = 1, subpnc
          input_soiltext(c,r) = soil_int(lon_line(c,1))
       end do
    end do
    deallocate( soil_int )
    deallocate( lat_line, lon_line )
    
    ! -------------------------------------------------------------------
    !    AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
    ! -------------------------------------------------------------------
    mi = subpnc*subpnr
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    if( mi .ne. mo .and. LDT_rc%soiltext_gridtransform(n) == "none" ) then
       write(LDT_logunit,*) "[ERR] Spatial transform, 'none', is selected, but number of"
       write(LDT_logunit,*) "  input and output points do not match. Select other spatial"
       write(LDT_logunit,*) "  option (e.g., mode, neighbor, tile, etc.)."
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
    endif
    allocate( li(mi), gi(mi) )
!    gi = float(LDT_rc%waterclass)
    li = .false.
    lo1 = .false.;  lo2 = .false.
    
    !- Assign 2-D input soil text to 1-D for aggregation routines:
    i = 0
    do r = 1, subpnr
       do c = 1, subpnc;  i = i + 1
          gi(i) = input_soiltext(c,r)
          if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
       enddo
    enddo
    
    !- Apply the spatial transform option:
    select case( LDT_rc%soiltext_gridtransform(n) ) 
       
       !- (a) Single-layer selection:
    case( "none", "mode", "neighbor" )
       
       !- Transform parameter from original grid to LIS output grid:
       call LDT_transform_paramgrid(n, LDT_rc%soiltext_gridtransform(n), &
            subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )
       
       !- Convert 1D count to 2D grid fgrds:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             temp(c,r) = go1(i)
          enddo
       enddo
       
       !- (b) Estimate TILED soiltexture files (gridcnt):
    case( "tile" ) 
       
       !- Calculate total counts for each soil type in each coarse gridcell:
       call LDT_transform_paramgrid(n, LDT_rc%soiltext_gridtransform(n), &
            subparam_gridDesc, mi, num_bins, gi, li, mo, go2, lo2 )
       
       !- Convert 1D gridcnt to 2D grid fgrds:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             do t = 1, num_bins
                gridcnt(c,r,t) = go2(i,t)
             end do
          enddo
       enddo
       
    end select  ! End grid cnt aggregation method
    deallocate( gi, li )
    deallocate( input_soiltext )
        
    !- Bring 2-D Array to 3-D Soil tile space:
    if( LDT_rc%soiltext_gridtransform(n) == "none" .or. &  ! Non-tiled surfaces
         LDT_rc%soiltext_gridtransform(n) == "mode" .or. & 
         LDT_rc%soiltext_gridtransform(n) == "neighbor" ) then  
       
       !- Assign soil texture types of less than 0 to an actual texture value:
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             if( nint(temp(c,r)) .le. 0 ) then
                ! EMK:  Soil data is missing for south of 60S (Antarctica).  
                ! A simple fix is to use "Other" (for Land-Ice)
                if (LDT_domain(n)%lat(c + (r-1)*LDT_rc%lnc(n)) < -60.) then
                   temp(c,r) = 16   ! STATSGO -- Other (Land-Ice)
                else
                   temp(c,r) = 12   ! STATSGO -- Clay 
                end if
             endif
!             if( (nint(temp(c,r)) .ne. water_class  ) .and. &
!                  (nint(temp(c,r)) .ne. LDT_rc%udef) ) then
!                gridcnt(c,r,NINT(temp(c,r))) = 1.0
!             endif
          enddo
       enddo
    end if
    
    !- Estimate final grid fraction:
!    call param_index_fgrdcalc( n, LDT_rc%soiltext_proj, LDT_rc%soiltext_gridtransform(n), &
!         water_class, num_bins, gridcnt, fgrd )
    
    call LDT_releaseUnitNumber(ftn)
    write(unit=LDT_logunit,fmt=*) "[INFO] Done reading Native STATSGO+FAO texture file"
    
    
  end subroutine read_HWSD_STATSGO2_texture

  !----------------------------------------------------------------------  

  SUBROUTINE soil_para_hwsd (nx,ny,gfiler)
    
    ! Processing NGDC-HWSD-STATSGO merged soil properties with Woesten Soil
    ! Parameters and produces tau_param.dat and soil_param.dat files
    
    implicit none	    
    integer, intent (in) :: nx, ny 
    character(*)  :: gfiler
    real, dimension (:), allocatable ::           &
         a_sand,a_clay,a_silt,a_oc,a_bee,a_psis, &
         a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf, &
         atau_2cm,btau_2cm 
    integer, dimension (100,3) :: table_map
    integer, dimension (3) :: nsoil_pcarbon 
    type (mineral_perc) :: min_percs
    
    integer :: n,maxcat,i,j,k,ktop,ncid,i_highd,j_highd,nx_adj,ny_adj
    integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,d_undef,   &
         i1,i2,icount
    character*100 :: fname,fout
    character*10 :: string
    character*2 :: VV,HH
    
    logical, allocatable, dimension(:,:) :: land_pixels
    integer, allocatable, dimension (:,:) :: &
         net_data1,net_data2,net_data3,net_data4,net_data5,net_data6 ,net_data7 
    integer (kind=2) , allocatable, target, dimension (:,:) :: SOIL_HIGH,  &
         sand_top,clay_top,oc_top,sand_sub,clay_sub,oc_sub, grav_grid
    integer (kind=2), pointer, dimension (:,:) :: Raster, &  
         Raster1,Raster2,Raster3,Raster4,Raster5,Raster6
    integer,          allocatable, dimension (:) :: tileid_vec,arrayA,arrayB          
    integer (kind=2), allocatable, dimension (:) ::  &
         data_vec1, data_vec2,data_vec3, data_vec4,data_vec5, data_vec6
    REAL, ALLOCATABLE, dimension (:) :: soildepth, grav_vec,soc_vec,poc_vec,&
         ncells_top,ncells_top_pro,ncells_sub_pro 
    integer(kind=2) , allocatable, dimension (:) :: ss_clay,    &
         ss_sand,ss_clay_all,ss_sand_all,ss_oc_all
    REAL, ALLOCATABLE :: count_soil(:)
    integer, allocatable, target, dimension (:,:) :: tile_id
    integer, pointer :: iRaster(:,:)
    integer :: tindex, pfafindex,fac,o_cl,o_clp,fac_surf,vtype
    real,dimension(4) :: cFamily
    real   ,dimension(5) :: cF_lim
    logical :: first_entry = .true.
    logical :: regrid,write_file
    INTEGER, allocatable, dimension (:) :: soil_class_top,soil_class_com
    REAL :: sf,factor,wp_wetness,fac_count
    logical                            :: file_exists
    REAL, ALLOCATABLE, DIMENSION (:,:) :: parms4file
    REAL, PARAMETER  :: p_poros = 0.93, p_bee = 3.5, p_psis = -0.03, p_ks = 2.8e-5, pmap_thresh = 0.3
    REAL, DIMENSION (:), POINTER      :: PMAP
    REAL :: d_poros, d_bee, d_psis, d_ks

    ! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
    !
    ! NOTE: "!$" is for conditional compilation
    !
    logical :: running_omp = .false.
    !
    !$ integer :: omp_get_thread_num, omp_get_num_threads
    !
    integer :: n_threads=1, li, ui, t_count
    !
    integer, dimension(:), allocatable :: low_ind, upp_ind
    !
    ! ------------------------------------------------------------------
    
    ! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
    !
    ! FIND OUT WHETHER -omp FLAG HAS BEEN SET DURING COMPILATION
    !
    !$ running_omp = .true.         ! conditional compilation
    !
    ! ECHO BASIC OMP VARIABLES
    !
    !$OMP PARALLEL DEFAULT(NONE) SHARED(running_omp,n_threads) 
    !
    !$OMP SINGLE
    !
    !$ n_threads = omp_get_num_threads()
    !
    !$ write (*,*) 'running_omp = ', running_omp
    !$ write (*,*)
    !$ write (*,*) 'parallel OpenMP with ', n_threads, 'threads'
    !$ write (*,*)
    !$OMP ENDSINGLE
    !
    !$OMP CRITICAL
    !$ write (*,*) 'thread ', omp_get_thread_num(), ' alive'
    !$OMP ENDCRITICAL
    !
    !$OMP BARRIER
    !
    !$OMP ENDPARALLEL
    
    if (first_entry) then
       nullify(iraster) ; first_entry = .false.
    endif
    cF_lim(1) =  0.
    cF_lim(2) =  0.4      ! 0.365    ! 0.3
    cF_lim(3) =  0.64     ! 0.585    ! 4.0	   
    cF_lim(4) =  15./1.72 ! 9.885    ! 8.5
    cF_lim(5) =  100.0
    
    nsoil_pcarbon(1) = 84 ! 84
    nsoil_pcarbon(2) = nsoil_pcarbon(1) + 84 ! 84
    nsoil_pcarbon(3) = nsoil_pcarbon(2) + 84 ! 57
    
    fname='clsm/catchment.def'
    !
    ! Reading number of cathment-tiles from catchment.def file
    ! 
    open (10,file=fname,status='old',action='read',form='formatted')
    read(10,*) maxcat
    
    close (10,status='keep')
    
    fname =trim(c_data)//'SOIL-DATA/GSWP2_soildepth_H11V13.nc'
    status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
    status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
    status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
    status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
    status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
    status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
    status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
    status = NF_CLOSE(ncid); VERIFY_(STATUS)
    
    allocate(soildepth(1:maxcat))
    allocate(soil_high(1:i_highd,1:j_highd))  
    allocate(count_soil(1:maxcat))  
    allocate(tile_id(1:nx,1:ny))
    allocate(net_data1 (1:nc_10,1:nr_10))
    
    fname=trim(gfiler)//'.rst'
              
    ! Reading tile-id raster file
    
    open (10,file=fname,status='old',action='read',  &
           form='unformatted',convert='little_endian')
    
    do j=1,ny
       read(10)tile_id(:,j)
    end do
    
    close (10,status='keep')
    
    ! reading soil depth data
    
    soil_high = -9999
    do jx = 1,18
       do ix = 1,36
          write (vv,'(i2.2)')jx
          write (hh,'(i2.2)')ix 
          fname = trim(c_data)//'SOIL-DATA/GSWP2_soildepth_H'//hh//'V'//vv//'.nc'
          status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
          if(status == 0) then
             status = NF_GET_att_INT  (ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
             status = NF_GET_att_INT  (ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
             status = NF_GET_att_INT  (ncid,4,'UNDEF',d_undef); VERIFY_(STATUS)
             status = NF_GET_att_REAL (ncid,4,'ScaleFactor',sf); VERIFY_(STATUS)
             status = NF_GET_VARA_INT (ncid, 4,(/1,1/),(/nc_10,nr_10/),net_data1); VERIFY_(STATUS)
             
             do j = jLL,jLL + nr_10 -1 
                do i = iLL, iLL + nc_10 -1 
                   if(net_data1(i-iLL +1 ,j - jLL +1) /= d_undef) &
                        soil_high(i,j) = net_data1(i-iLL +1 ,j - jLL +1)
                enddo
             enddo
             status = NF_CLOSE(ncid)
          endif
       end do
    end do
    
    deallocate (net_data1)
    
    ! Regridding 
    
    nx_adj = nx
    ny_adj = ny
    
    regrid = nx/=i_highd .or. ny/=j_highd
    
    if(regrid) then
       if(nx > i_highd) then 
          allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
          call LDT_RegridRaster(soil_high,raster)	
          iRaster => tile_id
          if(ny < j_highd) then
             print *,'nx > i_highd and ny < j_highd'
             stop 
          endif
       else
          if( .not.associated(iraster) ) then 
             allocate(iraster(i_highd,j_highd),stat=STATUS); VERIFY_(STATUS)
          endif
          call LDT_RegridRaster (tile_id,iraster)	
          raster  => soil_high
          nx_adj = i_highd
          ny_adj = j_highd
          
          if(ny > j_highd) then
             print *,'nx < i_highd and ny > j_highd'
             stop 
          endif
       endif
    else
       raster  => soil_high
       iRaster => tile_id
    end if
    
    ! Interpolation or aggregation on to catchment-tiles
    
    soildepth =0.
    count_soil = 0.
 
    do j=1,ny_adj
       do i=1,nx_adj
          if((iRaster(i,j).gt.0).and.(iRaster(i,j).le.maxcat)) then
             if ((raster(i,j).gt.0)) then
                soildepth(iRaster(i,j)) = &
                     soildepth(iRaster(i,j)) + sf*raster(i,j)
                count_soil(iRaster(i,j)) = &
                     count_soil(iRaster(i,j)) + 1. 
             endif
          endif
       end do
    end do
    
    DO n =1,maxcat
       if(count_soil(n)/=0.) soildepth(n)=soildepth(n)/count_soil(n)	
       soildepth(n) = max(soildepth(n),1334.)
       !                  soildepth(n) = soildepth(n) + 2000.
       !                  soildepth(n) = min(soildepth(n),8000.) 
    END DO
    
    deallocate (SOIL_HIGH)
    deallocate (count_soil)
    NULLIFY(Raster)
    
    ! Reading NGDC-HWSD-STATSGO merged Soil Properties
    
    fname =trim(c_data)//'SOIL-DATA/SoilProperties_H11V13.nc'
    status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
    status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
    status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
    status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
    status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
    status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
    status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
    status = NF_CLOSE(ncid)
    
    regrid = nx/=i_highd .or. ny/=j_highd
    allocate(net_data1 (1:nc_10,1:nr_10))
    allocate(net_data2 (1:nc_10,1:nr_10))
    allocate(net_data3 (1:nc_10,1:nr_10))
    allocate(net_data4 (1:nc_10,1:nr_10))
    allocate(net_data5 (1:nc_10,1:nr_10))
    allocate(net_data6 (1:nc_10,1:nr_10))
    allocate(net_data7 (1:nc_10,1:nr_10))
    
    allocate(sand_top (1:i_highd,1:j_highd))  
    allocate(clay_top (1:i_highd,1:j_highd))  
    allocate(oc_top   (1:i_highd,1:j_highd))  
    allocate(sand_sub (1:i_highd,1:j_highd))  
    allocate(clay_sub (1:i_highd,1:j_highd))  
    allocate(oc_sub   (1:i_highd,1:j_highd))  
    allocate(grav_grid(1:i_highd,1:j_highd))   
    
    sand_top = -9999.
    clay_top = -9999.
    oc_top   = -9999.
    sand_sub = -9999.
    clay_sub = -9999.
    oc_sub   = -9999.
    grav_grid= -9999.
    
    do jx = 1,18
       do ix = 1,36
          write (vv,'(i2.2)')jx
          write (hh,'(i2.2)')ix 
          fname = trim(c_data)//'SOIL-DATA/SoilProperties_H'//hh//'V'//vv//'.nc'
          status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
          if(status == 0) then
             status = NF_GET_att_INT  (ncid, NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
             status = NF_GET_att_INT  (ncid, NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
             status = NF_GET_att_INT  (ncid, 4,'UNDEF',d_undef); VERIFY_(STATUS)
             status = NF_GET_att_REAL (ncid, 4,'ScaleFactor',sf); VERIFY_(STATUS)
             status = NF_GET_VARA_INT (ncid, 4,(/1,1/),(/nc_10,nr_10/),net_data1); VERIFY_(STATUS)
             status = NF_GET_VARA_INT (ncid, 5,(/1,1/),(/nc_10,nr_10/),net_data2); VERIFY_(STATUS)
             status = NF_GET_VARA_INT (ncid, 6,(/1,1/),(/nc_10,nr_10/),net_data3); VERIFY_(STATUS)
             status = NF_GET_VARA_INT (ncid, 7,(/1,1/),(/nc_10,nr_10/),net_data4); VERIFY_(STATUS)
             status = NF_GET_VARA_INT (ncid, 8,(/1,1/),(/nc_10,nr_10/),net_data5); VERIFY_(STATUS)
             status = NF_GET_VARA_INT (ncid, 9,(/1,1/),(/nc_10,nr_10/),net_data6); VERIFY_(STATUS)
             status = NF_GET_VARA_INT (ncid,10,(/1,1/),(/nc_10,nr_10/),net_data7); VERIFY_(STATUS)
             do j = jLL,jLL + nr_10 -1 
                do i = iLL, iLL + nc_10 -1 
                   if(net_data1(i-iLL +1 ,j - jLL +1) /= d_undef) &
                        clay_top(i,j) = net_data1(i-iLL +1 ,j - jLL +1)
                   if(net_data2(i-iLL +1 ,j - jLL +1) /= d_undef) &
                        sand_top(i,j) = net_data2(i-iLL +1 ,j - jLL +1)
                   if(net_data3(i-iLL +1 ,j - jLL +1) /= d_undef) &
                        oc_top  (i,j) = net_data3(i-iLL +1 ,j - jLL +1)
                   if(net_data4(i-iLL +1 ,j - jLL +1) /= d_undef) &
                        clay_sub(i,j) = net_data4(i-iLL +1 ,j - jLL +1)
                   if(net_data5(i-iLL +1 ,j - jLL +1) /= d_undef) &
                        sand_sub(i,j) = net_data5(i-iLL +1 ,j - jLL +1)
                   if(net_data6(i-iLL +1 ,j - jLL +1) /= d_undef) &
                        oc_sub  (i,j) = net_data6(i-iLL +1 ,j - jLL +1)
                   if(net_data7(i-iLL +1 ,j - jLL +1) /= d_undef) &
                        grav_grid(i,j) = net_data7(i-iLL +1 ,j - jLL +1)
                enddo
             enddo
             status = NF_CLOSE(ncid)
          endif
       end do
    end do
    
    deallocate (net_data1)
    deallocate (net_data2)
    deallocate (net_data3)
    deallocate (net_data4)
    deallocate (net_data5)
    deallocate (net_data6)
    deallocate (net_data7)
    
    ! now regridding
    
    nx_adj = nx
    ny_adj = ny
    
    regrid = nx/=i_highd .or. ny/=j_highd
    
    if(regrid) then
       if(nx > i_highd) then 
          allocate(raster1(nx,ny),stat=STATUS); VERIFY_(STATUS)
          call LDT_RegridRaster(clay_top,raster1)	
          
          allocate(raster2(nx,ny),stat=STATUS); VERIFY_(STATUS)
          call LDT_RegridRaster(sand_top,raster2)	
          
          allocate(raster3(nx,ny),stat=STATUS); VERIFY_(STATUS)
          call LDT_RegridRaster(oc_top,  raster3)	
          
          allocate(raster4(nx,ny),stat=STATUS); VERIFY_(STATUS)
          call LDT_RegridRaster(clay_sub,raster4)	
          
          allocate(raster5(nx,ny),stat=STATUS); VERIFY_(STATUS)
          call LDT_RegridRaster(sand_sub,raster5)	
          
          allocate(raster6(nx,ny),stat=STATUS); VERIFY_(STATUS)
          call LDT_RegridRaster(oc_sub,  raster6)	
          
          allocate(raster (nx,ny),stat=STATUS); VERIFY_(STATUS)
          call LDT_RegridRaster(grav_grid,raster)
          
          iRaster => tile_id
          
          if(ny < j_highd) then
             print *,'nx > i_highd and ny < j_highd'
             stop 
          endif
       else
          nx_adj = i_highd
          ny_adj = j_highd
          if( .not.associated(iraster) ) then
             allocate(iRaster(i_highd,j_highd),stat=STATUS); VERIFY_(STATUS)
          endif
          call LDT_RegridRaster(tile_id,iRaster)	
          
          raster1 => clay_top
          raster2 => sand_top
          raster3 => oc_top
          raster4 => clay_sub
          raster5 => sand_sub
          raster6 => oc_sub
          raster  => grav_grid
          
          if(ny > j_highd) then
             print *,'nx < i_highd and ny > j_highd'
             stop 
          endif
       endif
    else
       iRaster => tile_id
       raster1 => clay_top
       raster2 => sand_top
       raster3 => oc_top
       raster4 => clay_sub
       raster5 => sand_sub
       raster6 => oc_sub
       raster  => grav_grid
    end if
    
    ! Deallocate large arrays
    
    allocate(land_pixels(1:size(iRaster,1),1:size(iRaster,2)))
    land_pixels = (iRaster >=1).and.(iRaster<=maxcat)
    i1 = count(land_pixels)   
    deallocate (land_pixels) 
    
    allocate (tileid_vec(1:i1))
    allocate (data_vec1 (1:i1))
    allocate (data_vec2 (1:i1))
    allocate (data_vec3 (1:i1))
    allocate (data_vec4 (1:i1))
    allocate (data_vec5 (1:i1))
    allocate (data_vec6 (1:i1))
    allocate (grav_vec  (1:maxcat))
    allocate (soc_vec   (1:maxcat))
    allocate (poc_vec   (1:maxcat))
    allocate (ncells_top  (1:maxcat))
    allocate (ncells_top_pro  (1:maxcat))
    allocate (ncells_sub_pro  (1:maxcat))
    allocate(count_soil(1:maxcat))  
    count_soil = 0.
    grav_vec   = 0.
    soc_vec    = 0.
    poc_vec    = 0.
    ncells_top = 0.
    ncells_top_pro = 0.
    ncells_sub_pro = 0.
    
    n =1
    do j=1,ny_adj
       do i=1,nx_adj
          if((iRaster(i,j).ge.1).and.(iRaster(i,j).le.maxcat)) then
             
             tileid_vec (n) =  iRaster(i,j)
             data_vec1  (n) =  Raster1(i,j)
             data_vec2  (n) =  Raster2(i,j)
             data_vec3  (n) =  Raster3(i,j)
             data_vec4  (n) =  Raster4(i,j)
             data_vec5  (n) =  Raster5(i,j)
             data_vec6  (n) =  Raster6(i,j)
             
             if ((raster(i,j).gt.0)) then 
                grav_vec(iRaster(i,j)) = &
                     grav_vec(iRaster(i,j)) + sf*raster(i,j)
                count_soil(iRaster(i,j)) = &
                     count_soil(iRaster(i,j)) + 1. 
             endif
             n = n + 1
          endif
       end do
    end do
    
    DO n =1,maxcat
       if(count_soil(n)/=0.) grav_vec(n)=grav_vec(n)/count_soil(n)	
    END DO
    
    deallocate (grav_grid)
    deallocate (count_soil)
    NULLIFY(Raster)
    
    NULLIFY(Raster1,Raster2,Raster3,Raster4,Raster5,Raster6)
    deallocate (clay_top,sand_top,oc_top,clay_sub,sand_sub,oc_sub) 
    deallocate (tile_id)
    
    allocate (arrayA    (1:i1))
    allocate (arrayB    (1:i1))
    
    arrayA = tileid_vec
    arrayB = data_vec1
    call LDT_quicksort (arrayA, arrayB)
    data_vec1 = arrayB
    
    arrayA = tileid_vec
    arrayB = data_vec2
    call LDT_quicksort (arrayA, arrayB)
    data_vec2 = arrayB
    
    arrayA = tileid_vec
    arrayB = data_vec3
    call LDT_quicksort (arrayA, arrayB)
    data_vec3 = arrayB
    
    arrayA = tileid_vec
    arrayB = data_vec4
    call LDT_quicksort (arrayA, arrayB)
    data_vec4 = arrayB
    
    arrayA = tileid_vec
    arrayB = data_vec5
    call LDT_quicksort (arrayA, arrayB)
    data_vec5 = arrayB
    
    arrayA = tileid_vec
    arrayB = data_vec6
    call LDT_quicksort (arrayA, arrayB)
    data_vec6 = arrayB
    tileid_vec= arrayA
    deallocate (arrayA, arrayB)
    
    ! Reading Woesten Soil Parameters and CLSM tau parameters

    allocate(a_sand  (1:n_SoilClasses))
    allocate(a_clay  (1:n_SoilClasses))
    allocate(a_silt  (1:n_SoilClasses))
    allocate(a_oc    (1:n_SoilClasses))
    allocate(a_bee   (1:n_SoilClasses))
    allocate(a_psis  (1:n_SoilClasses))
    allocate(a_poros (1:n_SoilClasses))
    allocate(a_wp    (1:n_SoilClasses))
    allocate(a_aksat (1:n_SoilClasses))
    allocate(atau    (1:n_SoilClasses))
    allocate(btau    (1:n_SoilClasses))
    allocate(atau_2cm(1:n_SoilClasses))
    allocate(btau_2cm(1:n_SoilClasses))
    allocate(a_wpsurf(1:n_SoilClasses))
    allocate(a_porosurf(1:n_SoilClasses))
    
    fname = trim(c_data)//'SoilClasses-SoilHyd-TauParam.dat'
    table_map = 0
    open (11, file=trim(fname), form='formatted',status='old', &
         action = 'read')
    read (11,'(a)')fout
    do n =1,n_SoilClasses 
       read (11,'(4f7.3,4f8.4,e13.5,2f12.7,2f8.4,4f12.7)')a_sand(n),a_clay(n),a_silt(n),a_oc(n),a_bee(n),a_psis(n), &
            a_poros(n),a_wp(n),a_aksat(n),atau(n),btau(n),a_wpsurf(n),a_porosurf(n),atau_2cm(n),btau_2cm(n)
       
       min_percs%clay_perc = a_clay(n)
       min_percs%silt_perc = a_silt(n)
       min_percs%sand_perc = a_sand(n)
       if(n <= nsoil_pcarbon(1))                              table_map(soil_class (min_percs),1) = n  
       if((n > nsoil_pcarbon(1)).and.(n <= nsoil_pcarbon(2))) table_map(soil_class (min_percs),2) = n  
       if((n > nsoil_pcarbon(2)).and.(n <= nsoil_pcarbon(3))) table_map(soil_class (min_percs),3) = n 
       
    end do
    close (11,status='keep') 
    
    !  When Woesten Soil Parameters are not available for a particular Soil Class
    !  ,as assumed by tiny triangles in HWSD soil triangle, Woesten Soil
    !  parameters from the nearest available tiny triangle will be substituted.
    	     	  
      do n =1,10
         do k=1,n*2 -1
            
            min_percs%clay_perc = 100. -((n-1)*10 + 5)
            min_percs%sand_perc = 100. -  min_percs%clay_perc -2.-(k-1)*5.
            min_percs%silt_perc = 100. -  min_percs%clay_perc - min_percs%sand_perc
            
            i = soil_class (min_percs)
            
            if(table_map (i,1)== 0) then
               j = GDL_center_pix (a_clay(1:nsoil_pcarbon(1)),a_sand(1:nsoil_pcarbon(1)), &
                    min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)         
               min_percs%clay_perc = a_clay(j)
               min_percs%silt_perc = a_silt(j)
               min_percs%sand_perc = a_sand(j)		   
               table_map (i,1)= table_map (soil_class (min_percs),1)
            endif
            
            min_percs%clay_perc = 100. -((n-1)*10 + 5)
            min_percs%sand_perc = 100. -  min_percs%clay_perc -2.-(k-1)*5.
            min_percs%silt_perc = 100. -  min_percs%clay_perc - min_percs%sand_perc
            
            if(table_map (i,2)== 0) then   
               j = GDL_center_pix(a_clay(nsoil_pcarbon(1)+1 : nsoil_pcarbon(2)),         &
                    a_sand(nsoil_pcarbon(1)+1 : nsoil_pcarbon(2)),         &
                    min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.) 
               min_percs%clay_perc = a_clay(j + nsoil_pcarbon(1))
               min_percs%silt_perc = a_silt(j + nsoil_pcarbon(1))
               min_percs%sand_perc = a_sand(j + nsoil_pcarbon(1))		   
               table_map (i,2)= table_map (soil_class (min_percs),2)	         
            endif
            
            min_percs%clay_perc = 100. -((n-1)*10 + 5)
            min_percs%sand_perc = 100. -  min_percs%clay_perc -2.-(k-1)*5.
            min_percs%silt_perc = 100. -  min_percs%clay_perc - min_percs%sand_perc
            
            if(table_map (i,3)== 0) then   
               j = GDL_center_pix (a_clay(nsoil_pcarbon(2)+1 : nsoil_pcarbon(3)),         &
                    a_sand(nsoil_pcarbon(2)+1 : nsoil_pcarbon(3)),         &
                    min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.) 
               min_percs%clay_perc = a_clay(j + nsoil_pcarbon(2))
               min_percs%silt_perc = a_silt(j + nsoil_pcarbon(2))
               min_percs%sand_perc = a_sand(j + nsoil_pcarbon(2))		   
               table_map (i,3)= table_map (soil_class (min_percs),3)	         	         
            endif
         end do
      end do
      
      ! Now deriving soil types based on NGDC-HWSD-STATSGO merged soil property maps
      
      allocate (soil_class_top (1:maxcat))
      allocate (soil_class_com (1:maxcat))
      soil_class_top =-9999
      soil_class_com =-9999
      
      allocate(low_ind(n_threads))
      allocate(upp_ind(n_threads))
      low_ind(1)         = 1
      upp_ind(n_threads) = maxcat
      
      if (running_omp)  then
         do i=1,n_threads-1  
            upp_ind(i)   = low_ind(i) + (maxcat/n_threads) - 1 
            low_ind(i+1) = upp_ind(i) + 1
         end do
      end if
      
      !$OMP PARALLELDO DEFAULT(NONE)                          &
      !$OMP SHARED( n_threads, low_ind, upp_ind, tileid_vec,  &
      !$OMP         sf,data_vec1,data_vec2,data_vec3,         &
      !$OMP         data_vec4,data_vec5,data_vec6,cF_lim,     &
      !$OMP         table_map,soil_class_top,soil_class_com,  &
      !$OMP         soc_vec,poc_vec,ncells_top,ncells_top_pro,&
      !$OMP         ncells_sub_pro)    &
      !$OMP PRIVATE(n,i,j,k,icount,t_count,i1,i2,ss_clay,     &
      !$OMP         ss_sand,ss_clay_all,ss_sand_all,          &
      !$OMP         ss_oc_all,cFamily,factor,o_cl,o_clp,ktop, &
      !$OMP         min_percs, fac_count, write_file)
      
      DO t_count = 1,n_threads
         DO n = low_ind(t_count),upp_ind(t_count)
            
            write_file = .false.
            
            !	if (n==171010)  write_file = .true.
            
            if(n==low_ind(t_count)) then
               icount = 1
               do k=1,low_ind(t_count) - 1
                  do while (tileid_vec(icount)== k)
                     icount = icount + 1
                  end do
               end do
            endif
            
            i1 = icount 
            
            loop: do while (tileid_vec(icount)== n)
               if(icount <= size(tileid_vec,1)) icount = icount + 1
               if(icount > size(tileid_vec,1)) exit loop
            end do loop
            
            i2 = icount -1
            i = i2 - i1 + 1
            
            allocate(ss_clay    (1:2*i))
            allocate(ss_sand    (1:2*i))
            allocate(ss_clay_all(1:2*i))
            allocate(ss_sand_all(1:2*i))
            allocate(ss_oc_all  (1:2*i))
            
            ss_clay    = 0    
            ss_sand    = 0	
            ss_clay_all= 0
            ss_sand_all= 0
            ss_oc_all  = 0
            
            ss_clay_all (1:i)     = data_vec1(i1:i2)
            ss_sand_all (1:i)     = data_vec2(i1:i2)
            ss_oc_all   (1:i)     = data_vec3(i1:i2)	
            ss_clay_all (1+i:2*i) = data_vec4(i1:i2) 
            ss_sand_all (1+i:2*i) = data_vec5(i1:i2)
            ss_oc_all   (1+i:2*i) = data_vec6(i1:i2)	
            
            cFamily = 0.
            
            do j=1,i
               if(j <= i) factor = 1.
               if((ss_oc_all(j)*sf >=  cF_lim(1)).and. (ss_oc_all(j)*sf < cF_lim(2))) cFamily(1) = cFamily(1) + factor
               if((ss_oc_all(j)*sf >=  cF_lim(2)).and. (ss_oc_all(j)*sf < cF_lim(3))) cFamily(2) = cFamily(2) + factor
               if((ss_oc_all(j)*sf >=  cF_lim(3)).and. (ss_oc_all(j)*sf < cF_lim(4))) cFamily(3) = cFamily(3) + factor
               if((ss_oc_all(j)*sf >=  cF_lim(4) ))                                   cFamily(4) = cFamily(4) + factor
            end do
            
            if (sum(cFamily) == 0.) o_cl  = 1
            if (sum(cFamily)  > 0.) o_cl  = maxloc(cFamily, dim = 1)
            
            cFamily = 0.
            
            do j=1,2*i
               if(j <= i) factor = 1.
               if(j  > i) factor = 2.33
               if((ss_oc_all(j)*sf >=  cF_lim(1)).and. (ss_oc_all(j)*sf < cF_lim(2))) cFamily(1) = cFamily(1) + factor
               if((ss_oc_all(j)*sf >=  cF_lim(2)).and. (ss_oc_all(j)*sf < cF_lim(3))) cFamily(2) = cFamily(2) + factor
               if((ss_oc_all(j)*sf >=  cF_lim(3)).and. (ss_oc_all(j)*sf < cF_lim(4))) cFamily(3) = cFamily(3) + factor
               if((ss_oc_all(j)*sf >=  cF_lim(4) ))                                   cFamily(4) = cFamily(4) + factor
            end do
            
            if (sum(cFamily) == 0.) o_clp = 1
            if (sum(cFamily)  > 0.) o_clp = maxloc(cFamily, dim = 1)
            
            if(o_cl == 4) then 
               soil_class_top(n) = n_SoilClasses
               ktop = 0
               do j=1,i
                  if(ss_oc_all(j)*sf >= cF_lim(4)) then
                     soc_vec (n) = soc_vec(n) + ss_oc_all(j)*sf
                     ktop = ktop + 1
                  endif
               end do
               if(ktop.ne.0) soc_vec (n)   = soc_vec(n)/ktop
               ncells_top(n) = 100.*float(ktop)/float(i)
            else 
               k = 1
               ktop = 1
               
               do j=1,i
                  if((ss_oc_all(j)*sf >= cF_lim(o_cl)).and.(ss_oc_all(j)*sf < cF_lim(o_cl + 1))) then 
                     if((ss_clay_all(j)*sf >= 0.).and.(ss_sand_all(j)*sf >= 0.)) then   
                        ss_clay (k) = ss_clay_all(j)
                        ss_sand (k) = ss_sand_all(j)
                        if((ss_clay (k) + ss_sand (k)) > 9999) then
                           if(ss_clay (k) >= ss_sand (k)) then
                              ss_sand (k) = 10000 - ss_clay (k)
                           else
                              ss_clay (k) = 10000 - ss_sand (k)
                           endif
                        endif
                        soc_vec (n) = soc_vec(n) + ss_oc_all(j)*sf
                        k = k + 1
                        ktop = ktop + 1
                     endif
                  endif
               end do
               
               k = k - 1
               ktop = ktop -1
               if(ktop.ne.0) soc_vec (n) = soc_vec(n)/ktop
               ncells_top(n) = 100.*float(ktop)/float(i)
               if (write_file) write(80+n,*)ktop,o_cl
               if(ktop > 0) then 
                  if (write_file) write (80+n,*)ss_clay(1:ktop)
                  if (write_file) write (80+n,*)ss_sand(1:ktop)
               endif
               j = GDL_center_pix (sf, ktop,ktop, ss_clay(1:ktop),ss_sand(1:ktop))
               if (write_file) write(80+n,*)j
               
               if(j >=1) then 
                  min_percs%clay_perc = ss_clay(j)*sf
                  min_percs%sand_perc = ss_sand(j)*sf
                  min_percs%silt_perc = 100. - ss_clay(j)*sf - ss_sand(j)*sf
                  soil_class_top (n) = table_map(soil_class (min_percs),o_cl)   
               endif
            endif
            if (write_file) write(80+n,*)soil_class_top (n) 
            if(o_clp == 4) then 
               soil_class_com(n) = n_SoilClasses
               fac_count = 0.
               k =0
               ktop =0
               do j=1,2*i
                  if(ss_oc_all(j)*sf >= cF_lim(4)) then
                     if(j <= i) factor = 1.
                     if(j  > i) factor = 2.33
                     if(j  > i) k = k + 1
                     if(j <= i) ktop = ktop + 1
                     poc_vec (n) = poc_vec(n) + ss_oc_all(j)*sf*factor
                     fac_count = fac_count + factor
                  endif
               end do
               if(fac_count.ne.0) poc_vec (n) = poc_vec (n)/fac_count
               ncells_sub_pro(n) = 100.*float(k)/float(i)
               ncells_top_pro(n) = 100.*float(ktop)/float(i)
            else
               k = 1
               ktop = 1
               
               ss_clay=0
               ss_sand=0
               fac_count = 0.
               
               do j=1,2*i
                  if((ss_oc_all(j)*sf >=  cF_lim(o_clp)).and.(ss_oc_all(j)*sf < cF_lim(o_clp + 1))) then 
                     if((ss_clay_all(j)*sf >= 0.).and.(ss_sand_all(j)*sf >= 0.)) then 
                        if(j <= i) factor = 1.
                        if(j  > i) factor = 2.33
                        poc_vec (n) = poc_vec(n) + ss_oc_all(j)*sf*factor
                        fac_count = fac_count + factor
                        if(j <= i) then
                           ss_clay (k) = ss_clay_all(j)
                           ss_sand (k) = ss_sand_all(j)
                           if((ss_clay (k) + ss_sand (k)) > 9999) then
                              if(ss_clay (k) >= ss_sand (k)) then
                                 ss_sand (k) = 10000 - ss_clay (k)
                              else
                                 ss_clay (k) = 10000 - ss_sand (k)
                              endif
                           endif
                           k = k + 1
                           ktop = ktop + 1
                     else
                        ss_clay (k) = ss_clay_all(j)
                        ss_sand (k) = ss_sand_all(j)
                        if((ss_clay (k) + ss_sand (k)) > 9999) then
                           if(ss_clay (k) >= ss_sand (k)) then
                              ss_sand (k) = 10000 - ss_clay (k)
                           else
                              ss_clay (k) = 10000 - ss_sand (k)
                           endif
                        endif
                        k = k + 1                         
		       endif   
                     endif
                  endif
               end do
	    
               k = k - 1
               ktop = ktop -1
               if(fac_count.ne.0) poc_vec (n) = poc_vec(n)/fac_count
               ncells_top_pro(n) = 100.*float(ktop)/float(i)
               ncells_sub_pro(n) = 100.*float(k-ktop)/float(i)
               
               if (write_file) write (80+n,*)ktop,k,o_cl
               if (write_file) write (80+n,*)ss_clay(1:k)
               if (write_file) write (80+n,*)ss_sand(1:k)
               j = GDL_center_pix (sf, ktop,k, ss_clay(1:k),ss_sand(1:k))
               if (write_file) write(80+n,*) j
               if(j >=1) then 
                  min_percs%clay_perc = ss_clay(j)*sf
                  min_percs%sand_perc = ss_sand(j)*sf
                  min_percs%silt_perc = 100. - ss_clay(j)*sf - ss_sand(j)*sf
                  soil_class_com (n) = table_map(soil_class (min_percs),o_clp)  
               endif
               if (write_file) write(80+n,*) soil_class_com (n) 
               if (write_file) close(80+n)          
            endif
            deallocate (ss_clay,ss_sand,ss_clay_all,ss_sand_all,ss_oc_all)
         END DO
      END DO
      !$OMP ENDPARALLELDO

      call process_peatmap (nx, ny, gfiler, pmap)

      inquire(file='clsm/catch_params.nc4', exist=file_exists)

      if(file_exists) then
         status = NF_OPEN ('clsm/catch_params.nc4', NF_WRITE, ncid) ; VERIFY_(STATUS)
         allocate (parms4file (1:maxcat, 1:10))
      endif
    
      fname='clsm/catchment.def'
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      fname ='clsm/soil_param.first'
      open (11,file=trim(fname),form='formatted',status='unknown',action = 'write')

      fname ='clsm/tau_param.dat'
      open (12,file=trim(fname),form='formatted',status='unknown',action = 'write')

      fname ='clsm/mosaic_veg_typs_fracs'
      open (13,file=trim(fname),form='formatted',status='old',action = 'read')

      do n = 1, maxcat

      	 read (10,*) tindex,pfafindex
         read (13,*) tindex,pfafindex,vtype

         ! fill gaps from neighbor for rare missing values came from inconsistent masks
         if ((soil_class_top (n) == -9999).or.(soil_class_com (n) == -9999)) then

            ! if com-layer has data the issues is only with top-layer
            ! -------------------------------------------------------

            if(soil_class_com (n) >= 1) soil_class_top (n) = soil_class_com (n)

            ! if there is nothing look for the neighbor
            ! -----------------------------------------
            
            if (soil_class_com (n) == -9999) then
               do k = 1, maxcat
                  j  = 0
                  i1 = n - k
                  i2 = n + k
                  if((i1 >=     1).and.(soil_class_com (i1) >=1)) j = i1
                  if((i2 <=maxcat).and.(soil_class_com (i2) >=1)) j = i2

                  if (j > 0) then
                     soil_class_com (n) = soil_class_com (j)
                     soil_class_top (n) = soil_class_com (n)
                     grav_vec(n)        = grav_vec(j)
                     soc_vec(n)         = soc_vec (j)
                     poc_vec(n)         = poc_vec (j)
                  endif

                  if (soil_class_com (n) >=1) exit
               end do
            endif

         endif

         fac_surf = soil_class_top(n)
	 fac      = soil_class_com(n)

         wp_wetness = a_wp(fac) /a_poros(fac)

         d_poros = a_poros(fac)
         d_bee   = a_bee(fac)
         d_psis  = a_psis(fac)  
         d_ks    = a_aksat(fac)

         if((pmap (n) > pmap_thresh).or.(fac == 253)) then
            d_poros = p_poros
            d_bee   = p_bee
            d_psis  = p_psis 
            d_ks    = p_ks
            wp_wetness = a_wp(fac)     /d_poros
         endif

!         write (11,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4, f8.4)')tindex,pfafindex,      &
!               soil_class_top(n),soil_class_com(n),d_bee,d_psis,d_poros,               &
!               d_ks/exp(-1.0*zks*gnu),wp_wetness,soildepth(n),                 &
!               grav_vec(n),soc_vec(n),poc_vec(n), &
!               a_sand(fac_surf),a_clay(fac_surf),a_sand(fac),a_clay(fac), &
!	       a_wpsurf(fac_surf)/a_porosurf(fac_surf),a_porosurf(fac_surf), pmap(n)
	       	    
         write (12,'(i8,i8,4f10.7)')tindex,pfafindex, &
	       atau_2cm(fac_surf),btau_2cm(fac_surf),atau(fac_surf),btau(fac_surf)  

!         if (allocated (parms4file)) then
!
!            parms4file (n, 1) = d_bee
!            parms4file (n, 2) = d_ks /exp(-1.0*zks*gnu)
!            parms4file (n, 3) = d_poros
!            parms4file (n, 4) = d_psis
!            parms4file (n, 5) = wp_wetness
!            parms4file (n, 6) = soildepth(n)
!            parms4file (n, 7) = atau_2cm(fac_surf)
!            parms4file (n, 8) = btau_2cm(fac_surf)
!            parms4file (n, 9) = atau(fac_surf)
!            parms4file (n,10) = btau(fac_surf) 
!  
!  	 endif

      end do
      write (11,'(a)')'                    '
      write (11,'(a)')'FMT=i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4'
      write (11,'(a)')'TileIndex PfafID SoilClassTop SoilClassProfile BEE PSIS POROS Ks_at_SURF WPWET SoilDepth %Grav %OCTop %OCProf %Sand_top %Clay_top %Sand_prof %Clay_prof WPWET_SURF POROS_SURF'
      close (10, status = 'keep')	            
      close (11, status = 'keep')	            
      close (12, status = 'keep')	            
      close (13, status = 'keep')

      deallocate (data_vec1, data_vec2,data_vec3, data_vec4,data_vec5, data_vec6)
      deallocate (tileid_vec)
      deallocate (a_sand,a_clay,a_silt,a_oc,a_bee,a_psis,       &
            a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf, &
            atau_2cm,btau_2cm)
      deallocate (soildepth, grav_vec,soc_vec,poc_vec,&
             ncells_top,ncells_top_pro,ncells_sub_pro,soil_class_top,soil_class_com)
      if(file_exists) then
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BEE'  ) ,(/1/),(/maxcat/), parms4file (:, 1)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'COND' ) ,(/1/),(/maxcat/), parms4file (:, 2)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'POROS') ,(/1/),(/maxcat/), parms4file (:, 3)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'PSIS' ) ,(/1/),(/maxcat/), parms4file (:, 4)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'WPWET') ,(/1/),(/maxcat/), parms4file (:, 5)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'DP2BR') ,(/1/),(/maxcat/), parms4file (:, 6)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ATAU2') ,(/1/),(/maxcat/), parms4file (:, 7)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BTAU2') ,(/1/),(/maxcat/), parms4file (:, 8)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ATAU5') ,(/1/),(/maxcat/), parms4file (:, 9)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BTAU5') ,(/1/),(/maxcat/), parms4file (:,10)) ; VERIFY_(STATUS) 
         STATUS   = NF_CLOSE (NCID) ; VERIFY_(STATUS)
         DEALLOCATE (parms4file)
      endif

    END SUBROUTINE soil_para_hwsd


    
    ! --------------------------------------------------------------------------------------

    SUBROUTINE process_peatmap (nc, nr, gfiler, pmap)
      
      implicit none
      integer  , parameter                         :: N_lon_pm = 43200, N_lat_pm = 21600
      integer, intent (in)                         :: nc, nr
      real, pointer, dimension (:), intent (inout) :: pmap
      character(*), intent (in)                    :: gfiler
      integer                                      :: i,j, status, varid, ncid
      integer                                      :: NTILES        
      REAL, ALLOCATABLE, dimension (:)             :: count_pix
      REAL, ALLOCATABLE, dimension (:,:)           :: data_grid, pm_grid
      INTEGER, ALLOCATABLE, dimension (:,:)        :: tile_id
      character*100                                :: fout    
      
      ! Reading number of tiles
      ! -----------------------
      
      open (20, file = 'clsm/catchment.def', form = 'formatted', status = 'old', action =  'read')
      
      read (20, *) NTILES
      
      close (20, status = 'keep')
      
      ! READ PEATMAP source data files and regrid
      ! -----------------------------------------
      
      status  = NF_OPEN ('data/CATCH/PEATMAP_mask.nc4', NF_NOWRITE, ncid)
      
      allocate (pm_grid   (1 : NC      , 1 : NR))
      allocate (data_grid (1 : N_lon_pm, 1 : N_lat_pm)) 
      
      status  = NF_INQ_VARID (ncid,'PEATMAP',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_REAL (ncid,VarID, (/1,1/),(/N_lon_pm, N_lat_pm/), data_grid) ; VERIFY_(STATUS)
      
      call LDT_RegridRaster (data_grid, pm_grid)
      
      status = NF_CLOSE(ncid)
      
      ! Grid to tile
      ! ------------
      
      ! Reading tile-id raster file
      
      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do
      
      close (10,status='keep')     
      
      allocate (pmap      (1:NTILES))
      allocate (count_pix (1:NTILES))
      
      pmap      = 0.
      count_pix = 0.
      
      do j = 1,nr
         do i = 1, nc
            if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then                
               if(pm_grid(i,j) > 0.)  pmap (tile_id(i,j)) = pmap (tile_id(i,j)) + pm_grid(i,j)
               count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
            endif
         end do
      end do
      
      where (count_pix >   0.) pmap = pmap/count_pix
      
      deallocate (count_pix)
      deallocate (pm_grid)
      deallocate (tile_id)
      
    END SUBROUTINE process_peatmap

        
end module mod_HWSD_STATSGO2_texture
