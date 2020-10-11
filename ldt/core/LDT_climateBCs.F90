!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_ClimateBCsReader
!BOP
!
! !MODULE:  LDT_ClimateBCsReader
! USAGE:
! 1) read in from LDT config file  
!    BCS data source:                  SOURCE  
!    BCS climatology interval:         daily, 5day, 8day, monthly
!    BCS map projection:               latlon
!
! 2) In parameter reader
!    use LDT_ClimateBCsReader, ONLY : ClimateBCsReader 
!    type (ClimateBCsReader) :: bcr
!    call bcr%readDataset (nest,SOURCE, INTERVAL, projection, array) to process a single data set OR
!    call bcr%readDataset (nest,SOURCE, INTERVAL, projection, array, array2=array2, array3=array3, array4=array4) to process
!             multiple data for e.g. whitesky/blacksky VISDF, NIRDF can be processed together 
!    where array dimensions are NX (LDT_rc%lnc(nest)),NY (LDT_rc%lnr(nest)),
!    and NT (tsteps per year in clim array)
!
! 3) SUBROUTINE readDataset is the entry point to this module 
!    call ClimateBCsReader%init_bcs - initializes ClimateBCsReader structure
!    call ClimateBCsReader%%run_annual_cycle
!       SUBROUTINE run_annual_cycle
!       initialize data
!       annual_loop
!          call ClimateBCsReader%update_data at BCs data interval
!          update output arrays at INTERVAL
!       SUBROUTINE update_data
!          call read_10deg_geos5tiles/ read_global_nc4 -read in high resolution data
!          call regrid_to_lisgrid
!          call fill_gaps  
!
! 4) Steps to add new data sets
!    follow MCD15A2H climate filename convention and data strcuture  
!    add a new case to SUBROUTINE init_bcs
!    add SOURCE to call read_global_nc4 in SUBROUTINE update_data
!
!  27 August 2020: Sarith Mahanama; Initial implementation

  use ESMF
  use LDT_timeMgrMod, only : LDT_calendar
  use LDT_gridmappingMod
  use LDT_logMod,  only : LDT_logunit, LDT_endrun
  use LDT_coreMod
  use map_utils,   only : ij_to_latlon
  use netcdf
   
  implicit none
  
  PRIVATE
  
  type, public :: ClimateBCsReader
     
     integer                             :: NTIMES, NX, NY
     type(ESMF_TimeInterval)             :: data_dtstep, out_dtstep
     integer, allocatable, dimension (:) :: data_doys
     character*20                        :: var1, var2, source
     real                                :: sf
     integer, dimension(2)               :: range
     character*300                       :: BCS_PATH
     real, dimension(20)                 :: param_gridDesc, subparam_gridDesc
     integer                             :: glpnc, glpnr, subpnc, subpnr   
     integer,allocatable,dimension (:,:) :: lat_line, lon_line, local_mask     
     
   contains 
     
     ! public 
     generic,   public ::  readDataset => readDataset_general, readDataset_clsm
     procedure, public ::  run_annual_cycle
     procedure, public ::  update_data
     procedure, public ::  update_sibinputs

     ! private
     procedure, private :: init_bcs
     procedure, private :: readDataset_general
     procedure, private :: readDataset_clsm
     
  end type ClimateBCsReader

  type, private :: cdata_struc
     real, dimension  (:,:,:), allocatable :: gdata
  end type cdata_struc

  type, private :: sibalb_inputs_struc
     integer, allocatable, dimension (:,:)   :: ITYP
     real,    allocatable, dimension (:,:,:) :: LAI
     real,    allocatable, dimension (:,:,:) :: GREEN
  end type sibalb_inputs_struc

  integer, parameter                     :: ref_year = 2002
  integer                                :: status
  type (sibalb_inputs_struc),save        :: sibinputs
     
contains

  ! --------------------------------------------------------------

  SUBROUTINE readDataset_general (this, nest, SOURCE, out_interval, project, &
       CLIMDATA,CLIMDATA2,CLIMDATA3, CLIMDATA4, maskarray)     

      implicit none

      class (ClimateBCsReader), intent(inout)        :: this
      character(*), intent(in)                       :: SOURCE, out_interval, project
      integer, intent (in)                           :: nest
      real, dimension(:,:,:), intent(inout) :: CLIMDATA
      real, dimension(:,:,:), intent(inout), OPTIONAL :: CLIMDATA2, CLIMDATA3, CLIMDATA4
      real, optional, intent(inout)                  :: maskarray(LDT_rc%lnc(nest),LDT_rc%lnr(nest))
      type (cdata_struc),dimension(:),allocatable    :: cdata
      integer                                        :: NF, NX, NY, NT, i
 
      ! --------------------------------
      ! Initialize parameter description
      ! --------------------------------

      call this%init_bcs (trim(SOURCE),nest)
      
      if (out_interval == "daily"  ) call ESMF_TimeIntervalSet(this%out_dtstep, h=24  , rc=status ) ; VERIFY_(STATUS)
      if (out_interval == "5day"   ) call ESMF_TimeIntervalSet(this%out_dtstep, h=24*5, rc=status ) ; VERIFY_(STATUS)
      if (out_interval == "8day"   ) call ESMF_TimeIntervalSet(this%out_dtstep, h=24*8, rc=status ) ; VERIFY_(STATUS)
      if (out_interval == "monthly") call ESMF_TimeIntervalSet(this%out_dtstep,mm=1   , rc=status ) ; VERIFY_(STATUS)
      if(trim (project) == 'latlon') this%param_gridDesc(1)  = 0.

      ! ------------------------------------------------------------
      !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN BCS DATA
      ! ------------------------------------------------------------

      !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
      this%subparam_gridDesc = 0.
      call LDT_RunDomainPts( nest, project, this%param_gridDesc(:), &
           this%glpnc, this%glpnr, this%subpnc, this%subpnr, this%subparam_gridDesc, &
           this%lat_line, this%lon_line)

      if((this%subpnc /= size(CLIMDATA,1)) .or.(this%subpnr /= size(CLIMDATA,2))) then
          write(LDT_logunit,*)'[ERR] LIS domain dimensions mismatch in READ_CLIM_BCS'
          write(LDT_logunit,*)'  input dimensions : ', size(CLIMDATA,1), size(CLIMDATA,2)
          write(LDT_logunit,*)'  domain dimensions: ', this%subpnc, this%subpnr
          VERIFY_(1)
          call LDT_endrun         
      endif

      ! ------------------------------------------------
      !    READ REGRID INTERPOLATE AND POPULATE CLIMDATA
      ! ------------------------------------------------

      NF = 1
      if (present (CLIMDATA2)) NF = NF + 1
      if (present (CLIMDATA3)) NF = NF + 1
      if (present (CLIMDATA4)) NF = NF + 1

      NX = size(CLIMDATA,1)
      NY = size(CLIMDATA,2)
      NT = size(CLIMDATA,3)

      allocate (cdata (1:NF))
      do i = 1, NF
         allocate (cdata(i)%gdata (1:NX, 1:NY, 1: NT))
      end do

      call this%run_annual_cycle (nest, NF, NX, NY,NT, cdata)

      CLIMDATA  (:,:,:) = cdata(1)%gdata(:,:,:)      
      if (present (CLIMDATA2))  CLIMDATA2 (:,:,:) = cdata(2)%gdata(:,:,:)
      if (present (CLIMDATA3))  CLIMDATA3 (:,:,:) = cdata(3)%gdata(:,:,:)
      if (present (CLIMDATA4))  CLIMDATA4 (:,:,:) = cdata(3)%gdata(:,:,:)
      deallocate (cdata)

    end SUBROUTINE readDataset_general

    ! --------------------------------------------------------------------------------------------
    
    SUBROUTINE readDataset_clsm (this, nest, SOURCE, out_interval, project, CLIMDATA1,CLIMDATA2, maskarray, clsm)     

      implicit none

      class (ClimateBCsReader), intent(inout)        :: this
      character(*), intent(in)                       :: SOURCE, out_interval, project
      integer, intent (in)                           :: nest
      logical, intent (in)                           :: clsm
      real, dimension(:,:,:), intent(inout) :: CLIMDATA1, CLIMDATA2
      real, optional, intent(inout)                  :: maskarray(LDT_rc%lnc(nest),LDT_rc%lnr(nest))
      type (cdata_struc),dimension(:),allocatable    :: cdata
      integer                                        :: NF, NX, NY, NT, i
 
      ! --------------------------------
      ! Initialize parameter description
      ! --------------------------------

      call this%init_bcs (trim(SOURCE), nest)
      
      if (out_interval == "daily"  ) call ESMF_TimeIntervalSet(this%out_dtstep, h=24  , rc=status ) ; VERIFY_(STATUS)
      if (out_interval == "5day"   ) call ESMF_TimeIntervalSet(this%out_dtstep, h=24*5, rc=status ) ; VERIFY_(STATUS)
      if (out_interval == "8day"   ) call ESMF_TimeIntervalSet(this%out_dtstep, h=24*8, rc=status ) ; VERIFY_(STATUS)
      if (out_interval == "monthly") call ESMF_TimeIntervalSet(this%out_dtstep,mm=1   , rc=status ) ; VERIFY_(STATUS)
      if(trim (project) == 'latlon') this%param_gridDesc(1)  = 0.

      ! ------------------------------------------------------------
      !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN BCS DATA
      ! ------------------------------------------------------------

      !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
      this%subparam_gridDesc = 0.
      call LDT_RunDomainPts( nest, project, this%param_gridDesc(:), &
           this%glpnc, this%glpnr, this%subpnc, this%subpnr, this%subparam_gridDesc, &
           this%lat_line, this%lon_line)   

      if((this%subpnc /= size(CLIMDATA1,1)) .or.(this%subpnr /= size(CLIMDATA1,2))) then
          write(LDT_logunit,*)'[ERR] LIS domain dimensions mismatch in READ_CLIM_BCS'
          write(LDT_logunit,*)'  input dimensions : ', size(CLIMDATA1,1), size(CLIMDATA1,2)
          write(LDT_logunit,*)'  domain dimensions: ', this%subpnc, this%subpnr
          VERIFY_(1)
          call LDT_endrun         
      endif

      ! ------------------------------------------------
      !    READ REGRID INTERPOLATE AND POPULATE CLIMDATA
      ! ------------------------------------------------

      NF = 2

      NX = size(CLIMDATA1,1)
      NY = size(CLIMDATA1,2)
      NT = size(CLIMDATA1,3)

      allocate (cdata (1:NF))
      do i = 1, NF
         allocate (cdata(i)%gdata (1:NX, 1:NY, 1: NT))
      end do

      call this%run_annual_cycle (nest, NF, NX, NY,NT, cdata, clsm = .true.)

      CLIMDATA1  (:,:,:) = cdata(1)%gdata(:,:,:)
      CLIMDATA2  (:,:,:) = cdata(2)%gdata(:,:,:)
      deallocate (cdata)

    end SUBROUTINE readDataset_clsm
    
    !###################################################################################

    SUBROUTINE init_bcs (CP, SOURCE, nest)
          
      implicit none
      class (ClimateBCsReader), intent(inout) :: CP
      integer, intent (in)      :: nest
      character(*), intent (in) :: SOURCE
      integer, parameter  :: N_MONTHLY_DATES = 12
      integer, parameter  :: N_MODIS_DATES8  = 46
      integer, parameter  :: N_MODIS_DATES5  = 74

      integer, dimension (N_MODIS_DATES8), target    :: &
           MODIS_DOYS = (/                              &
           1  ,  9, 17, 25, 33, 41, 49, 57, 65,         &
           73 , 81, 89, 97,105,113,121,129,137,         &
           145,153,161,169,177,185,193,201,209,         &
           217,225,233,241,249,257,265,273,281,         &
           289,297,305,313,321,329,337,345,353,361/)

      integer, dimension (N_MONTHLY_DATES), target   :: &
           MONTHLY_DOYS = (/                            &
           1  , 32, 60, 91,121,152,                     &
           182,213,244,274,305,335/) 

      integer, dimension (N_MODIS_DATES5), target ::    &
           MODIS_DOYS5 = (/                             &
           1  ,  6, 11, 16, 21, 26, 31, 36, 41, 46,     &
           51 , 56, 61, 66, 71, 76, 81, 86, 91, 96,     &
           101,106,111,116,121,126,131,136,141,146,     &
           151,156,161,166,171,176,181,186,191,196,     &
           201,206,211,216,221,226,231,236,241,246,     &
           251,256,261,266,271,276,281,286,291,296,     &
           301,306,311,316,321,326,331,336,341,346,     &
           351,356,361,366/)

      integer :: input_cols 
      integer :: input_rows 
      real    :: IN_xres 
      real    :: IN_yres
      integer :: c,r,gr,gc, glpnc, glpnr
      real    :: rlon(LDT_rc%lnc(nest),LDT_rc%lnr(nest)),rlat(LDT_rc%lnc(nest),LDT_rc%lnr(nest)) 
      real    :: param_grid(20)
      character*300         :: c_data = '/discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001/'

      CP%source   = trim(source)

      select case (trim(SOURCE))
         
      case ('GSWPH')
         CP%NTIMES   = N_MONTHLY_DATES
         CP%NX       = 43200
         CP%NY       = 21600
         CP%var1     = 'grnFrac'
         CP%sf       = 0.001
         CP%range(1) = 0
         CP%range(2) = 1000
         CP%BCS_PATH = trim(c_data)//'/GSWP2_30sec_VegParam/GSWP2_VegParam_'
         allocate (CP%data_doys (1: CP%NTIMES))
         CP%data_doys= MONTHLY_DOYS
         call ESMF_TimeIntervalSet(CP%data_dtstep, mm=1, rc=status )   ; VERIFY_(STATUS)
                  
      case ('MCD43GFv5')
         CP%NTIMES = N_MODIS_DATES8
         CP%NX       = 43200
         CP%NY       = 21600
         CP%var1     = 'Alb_0.3_0.7'
         CP%var2     = 'Alb_0.7_5.0'
         CP%sf       = 0.001
         CP%range(1) = 0
         CP%range(2) = 1000
         CP%BCS_PATH = trim(c_data)//'/MODIS-Albedo2/MCD43GF_wsa_'
         allocate (CP%data_doys (1: CP%NTIMES))
         CP%data_doys= MODIS_DOYS
         call ESMF_TimeIntervalSet(CP%data_dtstep, h=24*8, rc=status ) ; VERIFY_(STATUS)    

      case ('MCD43GFv6')
         CP%NTIMES   = N_MODIS_DATES5
         CP%NX       = 43200
         CP%NY       = 21600
         CP%var1     = 'Albedo_Map_0.3-0.7'
         CP%var2     = 'Albedo_Map_0.7-5.0'
         CP%sf       = 0.001
         CP%range(1) = 0
         CP%range(2) = 1000
         CP%BCS_PATH = 'LS_PARAMETERS/MODIS/MCD43GF.006/MCD43GF.006_WSA_YYYY'
         allocate (CP%data_doys (1: CP%NTIMES))
         CP%data_doys= MODIS_DOYS5
         call ESMF_TimeIntervalSet(CP%data_dtstep, h=24*5, rc=status ) ; VERIFY_(STATUS)
         
      case ('GLASSA')
         CP%NTIMES = N_MODIS_DATES8
         CP%NX       = 7200
         CP%NY       = 3600
         CP%var1     = 'LAI'
         CP%sf       = 0.01
         CP%range(1) = 0
         CP%range(2) = 1000
         CP%BCS_PATH = trim(c_data)//'/GLASS-LAI/AVHRR.v4/GLASS01B02.V04.AYYYY'
         allocate (CP%data_doys (1: CP%NTIMES))
         CP%data_doys= MODIS_DOYS
         call ESMF_TimeIntervalSet(CP%data_dtstep, h=24*8, rc=status ) ; VERIFY_(STATUS)         
         
      case ('GLASSM')
         CP%NTIMES = N_MODIS_DATES8
         CP%NX       = 7200
         CP%NY       = 3600
         CP%var1     = 'LAI'
         CP%sf       = 0.01
         CP%range(1) = 0
         CP%range(2) = 1000
         CP%BCS_PATH = trim(c_data)//'/GLASS-LAI/MODIS.v4/GLASS01B01.V04.AYYYY'
         allocate (CP%data_doys (1: CP%NTIMES))
         CP%data_doys = MODIS_DOYS
         call ESMF_TimeIntervalSet(CP%data_dtstep, h=24*8, rc=status ) ; VERIFY_(STATUS)
         
      case ('MCD15A2H')
         CP%NTIMES   = N_MODIS_DATES8
         CP%NX       = 86400
         CP%NY       = 43200
         CP%var1     = 'Lai_500m'
         CP%sf       = 0.1
         CP%range(1) = 0
         CP%range(2) = 100
         CP%BCS_PATH = 'LS_PARAMETERS/MODIS/MCD15A2H.006/MCD15A2H.006_LAI_YYYY'
         allocate (CP%data_doys (1: CP%NTIMES))
         CP%data_doys= MODIS_DOYS
         call ESMF_TimeIntervalSet(CP%data_dtstep, h=24*8, rc=status ) ; VERIFY_(STATUS)

      case default
         
         print *, 'Unknown climatological data label : ', trim (SOURCE)
         VERIFY_(1) 
         
      end select
      
      input_cols = CP%NX
      input_rows = CP%NY
      IN_xres    = 360./REAL(CP%NX)
      IN_yres    = 180./REAL(CP%NY)
      CP%param_gridDesc(1)  = 0.          ! Latlon
      CP%param_gridDesc(2)  = input_cols
      CP%param_gridDesc(3)  = input_rows
      CP%param_gridDesc(4)  = -60.0  + (IN_yres/2) ! LL lat 
      CP%param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon 
      CP%param_gridDesc(6)  = 128
      CP%param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat
      CP%param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon
      CP%param_gridDesc(9)  = IN_yres     
      CP%param_gridDesc(10) = IN_xres     
      CP%param_gridDesc(20) = 64

      param_grid(:) = LDT_rc%mask_gridDesc(nest,:)
      glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
      glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1
      allocate (CP%local_mask (LDT_rc%lnc(nest),LDT_rc%lnr(nest)))

      CP%local_mask = LDT_rc%udef
      
      do r = 1, LDT_rc%lnr(nest)
         do c = 1, LDT_rc%lnc(nest)
            call ij_to_latlon(LDT_domain(nest)%ldtproj,float(c),float(r),&
                 rlat(c,r),rlon(c,r))
            gr = nint((rlat(c,r)-param_grid(4))/param_grid(10)) + 1
            gc = nint((rlon(c,r)-param_grid(5))/param_grid( 9)) + 1
            if(LDT_rc%global_mask(gc,gr) > 0. ) CP%local_mask(c,r) = 1
         end do
      end do
      
    END SUBROUTINE init_bcs
   
    ! ---------------------------------------------------------------------------
    
    subroutine update_data (CP, time_slice, nest, var_subset)
      
      implicit none
      class (ClimateBCsReader), intent(inout)      :: CP
      integer, intent (in)                         :: time_slice, nest
      real, intent (out), dimension(:,:,:)         :: var_subset
      real, allocatable, dimension (:,:,:)         :: datain
      
      allocate (datain (CP%NX, CP%NY,size(var_subset,3)))
      datain     = LDT_rc%udef
      var_subset = LDT_rc%udef

      if((trim(CP%source) == 'GSWPH') .OR. (trim(CP%source) == 'MCD43GFv5')) &
           call read_10deg_geos5tiles (time_slice, datain)
      
      if((trim(CP%source) == 'GLASSA')  .OR. (trim(CP%source) == 'GLASSM') .OR. &
         (trim(CP%source) == 'MCD15A2H').OR. (trim(CP%source) == 'MCD43GFv6'))  &
         call read_global_nc4 (time_slice, datain)
      
      call regrid_to_lisgrid (nest, datain, var_subset)
      call fill_gaps         (var_subset)
      
      deallocate (datain)
      
      contains

        ! ---------------------------------------------------------------------------
    
        subroutine read_10deg_geos5tiles (time_slice,datain)
          
          implicit none
          real, dimension (:,:,:), intent (inout) :: datain
          integer, intent (in)                    :: time_slice
          integer     :: jx, ix, ncid,iLL,jLL,  d_undef, nv,i,j
          character*2 :: hh,vv
          real        :: sf
          character*300 :: fname
          integer, parameter :: nc_10 = 1200
          integer, allocatable, dimension (:,:) :: net_data1,  net_data2

          nv = SIZE (datain,3)
          allocate(net_data1 (1:nc_10,1:nc_10))
          
          if(trim(CP%source) == 'MCD43GFv5') then
             allocate(net_data2 (1:nc_10,1:nc_10))
             if (nv .NE. 2) then
                print *, ' SOURCE : ', trim(CP%source),nv
                print *, '[ERR] in consistent 3rd dimension for albedo .. stopping '
                ASSERT_(.FALSE.)
                call LDT_endrun
             endif
          endif

          do jx = 1,18	
             do ix = 1,36
                write (vv,'(i2.2)')jx
                write (hh,'(i2.2)')ix 
                fname = trim(CP%BCS_PATH)//'H'//hh//'V'//vv//'.nc'
                status = NF90_OPEN(trim(fname),NF90_NOWRITE, ncid)                    ; VERIFY_(STATUS)
                if(status == 0) then
                   status = NF90_GET_att (ncid,NF90_GLOBAL,'i_ind_offset_LL',iLL)     ; VERIFY_(STATUS)
                   status = NF90_GET_att (ncid,NF90_GLOBAL,'j_ind_offset_LL',jLL)     ; VERIFY_(STATUS)
                   status = NF90_GET_att (ncid,VARID(ncid,trim(CP%var1)),'UNDEF',d_undef)  ; VERIFY_(STATUS)
                   status = NF90_GET_att (ncid,VARID(ncid,trim(CP%var1)),'ScaleFactor',sf) ; VERIFY_(STATUS) 
                   if (sf .NE. CP%sf) then
                      ASSERT_ (.FALSE.)
                   endif
                   status = NF90_GET_VAR (ncid, VARID(ncid,trim(CP%var1)),net_data1,start=(/1,1,time_slice/),count=(/nc_10,nc_10,1/)); VERIFY_(STATUS)
                   if (nv == 2) then
                      status = NF90_GET_VAR (ncid, VARID(ncid,trim(CP%var2)),net_data2,start=(/1,1,time_slice/),count=(/nc_10,nc_10,1/))
                      VERIFY_(STATUS)
                   endif
                   status = NF90_CLOSE(ncid)
                   do j = jLL,jLL + nc_10 -1 
                      do i = iLL, iLL + nc_10 -1 
                         if((net_data1(i-iLL +1 ,j - jLL +1) >=  CP%range(1)).and. & 
                              (net_data1(i-iLL +1 ,j - jLL +1) <=  CP%range(2))) then
                            datain (i,j,1) =  sf*net_data1(i-iLL +1 ,j - jLL +1)
                            if (nv == 2)  datain (i,j,2) =  sf*net_data2(i-iLL +1 ,j - jLL +1)
                         endif
                      end do
                   end do
                endif
             end do
          end do

          deallocate (net_data1)
          if(allocated (net_data2)) deallocate (net_data2)

        end subroutine read_10deg_geos5tiles
        
        ! ---------------------------------------------------------------------
    
        subroutine read_global_nc4 (time_slice, datain)
          
          implicit none
          real, dimension (:,:,:), intent (inout) :: datain
          integer, intent (in)                    :: time_slice
          character*200                           :: filename
          character*3                             :: DOY
          integer                                 :: i,j, nv, ncid
          integer*2, allocatable, dimension (:,:) :: int2_arr
          integer, allocatable, dimension (:,:)   :: data_array, data_array2

          nv = SIZE (datain,3)
          if (trim(CP%source) == 'MCD43GFv6') then
             allocate (data_array2 (1:CP%NX,1:CP%NY))
             if (nv .NE. 2) then
                print *, ' SOURCE : ', trim(CP%source),nv
                print *, '[ERR] in consistent 3rd dimension for albedo .. stopping '
                ASSERT_(.FALSE.)
                call LDT_endrun
             endif
          endif

          allocate (data_array (1:CP%NX,1:CP%NY))
          write (DOY,'(i3.3)') CP%data_doys(time_slice)
          filename = trim(CP%BCS_PATH)//DOY//'.nc4'
          status = NF90_OPEN(trim(filename),NF90_NOWRITE, ncid) ; VERIFY_(STATUS)
          
          if (trim(CP%source) == 'MCD15A2H') then
             allocate (int2_arr (1:CP%NX, 1:CP%NY))
             STATUS = NF90_GET_VAR (NCID,VARID(NCID,CP%var1 ), int2_arr)            ; VERIFY_(STATUS)
             data_array  = int2_arr
             if (nv == 2) STATUS = NF90_GET_VAR (NCID,VARID(NCID,CP%var2), int2_arr); VERIFY_(STATUS)
             data_array2 = int2_arr
             deallocate (int2_arr)
          else             
             STATUS = NF90_GET_VAR (NCID,VARID(NCID,CP%var1 ),   data_array)             ; VERIFY_(STATUS)
             if (nv == 2) STATUS = NF90_GET_VAR (NCID,VARID(NCID,CP%var2),   data_array2); VERIFY_(STATUS)
          endif
          STATUS = NF90_CLOSE (NCID)

          do j = 1, CP%NY
             do i = 1, CP%NX
                if((data_array(i,j) >=  CP%range(1)).and.(data_array(i,j) <=  CP%range(2))) then
                   datain (i,j,1) = CP%sf * data_array(i,j) 
                   if (nv == 2)  datain (i,j,2) = CP%sf * data_array2(i,j) 
                endif
             end do
          end do

          deallocate (data_array)
          if(allocated (data_array2)) deallocate (data_array2)         

        end subroutine read_global_nc4

        ! ------------------------------------------------------------
    
        integer function VARID (NCFID, VNAME) 
      
          integer, intent (in)      :: NCFID
          character(*), intent (in) :: VNAME
          
          STATUS = NF90_INQ_VARID (NCFID, trim(VNAME) ,VARID) ; VERIFY_(STATUS)
          
        end function VARID

        ! ---------------------------------------------------------------------

        SUBROUTINE regrid_to_lisgrid (nest, data_in, var_subset)

          ! adapted from params/irrigation/read_GRIPC_irrigfrac.F90 
          
          implicit none
          integer, intent (in)                  :: nest
          real, dimension(:,:,:), intent   (in) :: data_in
          real, dimension(:,:,:), intent(inout) :: var_subset
          
          ! arrays read from LL
          
          integer   :: mi                     ! Total number of input param grid array points
          integer   :: mo                     ! Total number of output LIS grid array points
          integer   :: nc, nr, i, j
          integer, allocatable  :: n11(:)     ! Array that maps the location of each input grid
          !   point in the output grid. 
          real,    allocatable  :: gi(:)      ! Input parameter 1d grid
          logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)
          real, allocatable, dimension (:,:) :: var_in
          real, allocatable, dimension (:,:)      :: go2     ! Output lis 1d grid
          logical*1, allocatable, dimension (:,:) :: lo2  ! Output logical mask (to match go)
          
          mi = CP%subpnc*CP%subpnr
          mo = LDT_rc%lnc(nest)*LDT_rc%lnr(nest)
          
          allocate (var_in(CP%subpnc,CP%subpnr))                
          allocate ( li(mi), gi (mi), n11(mi))
          allocate (go2(mo,1), lo2(mo,1))
          
          !- Create mapping between parameter domain and LIS grid domain:
          call upscaleByAveraging_input( CP%subparam_gridDesc, &
               LDT_rc%gridDesc(nest,:), mi, mo, n11)
          
          Data_Fields: do j = 1, size (data_in, 3)
             
             var_in = LDT_rc%udef
             do nr = 1, CP%subpnr
                do nc = 1, CP%subpnc
                   var_in (nc,nr) = data_in (CP%lon_line(nc,1),CP%lat_line (1,nr),j)
                end do
             end do
             
             ! -------------------------------------------------------------------
             !     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
             ! -------------------------------------------------------------------
             
             gi  = LDT_rc%udef
             li  = .false.
             
             !- Assign 2-D array to 1-D for aggregation routines:
             i = 0
             do nr = 1, CP%subpnr
                do nc = 1, CP%subpnc
                   i = i + 1            
                   if( var_in(nc,nr) .NE. LDT_rc%udef) then  
                      gi(i) = var_in(nc,nr)  
                      li(i) = .true.
                   endif
                end do
             enddo
             
             !- Calculate total counts for valid land pixels in each coarse gridcell:
             call upscaleByCnt( mi, mo, 1, LDT_rc%udef, n11, li, gi, &
                  lo2(:,1), go2 (:,1))
             i = 0
             do nr = 1, LDT_rc%lnr(nest)
                do nc = 1, LDT_rc%lnc(nest)
                   i = i + 1
                   var_subset (nc,nr,j) = go2(i,1)
                enddo
             enddo
          end do Data_Fields
          
          deallocate (var_in, li, gi, n11, go2, lo2)
   
        END SUBROUTINE regrid_to_lisgrid

        ! -------------------------------------------------------------------

        SUBROUTINE fill_gaps (var_subset)

          implicit none
          real, intent (inout), dimension(:,:,:)      :: var_subset
          real, allocatable, target, dimension(:,:,:) :: data_grid
          real,    pointer,     dimension (:,:)       :: subset
          integer :: NX, NY, NF, I,J,F
          INTEGER :: imn,imx,jmn,jmx,mval,l, ix, jx
          
          NX = SIZE (var_subset,1)
          NY = SIZE (var_subset,2)
          NF = SIZE (var_subset,3)
          allocate (data_grid(1:NX, 1:NY, 1:NF))

          data_grid = var_subset

          do F = 1, NF
             do j = 1, NY
                do i = 1, NX
                   if ((CP%local_mask (i,j) == 1).and.(var_subset (i,j,f) == LDT_rc%udef)) then
                      l = 1
                      do
                         imx=i + l
                         imn=i - l
                         jmn=j - l
                         jmx=j + l
                         imn=MAX(imn,1)
                         jmn=MAX(jmn,1)
                         imx=MIN(imx,nx)
                         jmx=MIN(jmx,ny)
                         subset => data_grid(imn: imx,jmn:jmx, f)
                         if(maxval(subset) > 0.) then 
                            var_subset (i,j,f) = sum(subset, subset>0.)/(max(1,count(subset>0.)))
                            exit
                         endif
                         l = l + 1
                         NULLIFY (subset)
                      end do
                   endif
                end do
             end do
          end do
          deallocate (data_grid)
        END SUBROUTINE fill_gaps
         
    end subroutine update_data

    !###################################################################################
      
    subroutine run_annual_cycle (CP, nest, NF, NX, NY,NT, cdata, clsm)
      
      implicit none
      class (ClimateBCsReader), intent(inout)        :: CP
      integer, intent (in)                           :: nest, NF, NX, NY,NT
      type (cdata_struc),dimension(NF),intent(inout) :: cdata
      logical, intent (in), optional                 :: clsm
      type(ESMF_Time)                                :: CURRENT_TIME, OUT_RING, END_TIME, TIME1, TIME2 
      type(ESMF_TimeInterval)                        :: DAY_DT
      integer                                        :: t1, t2, switch,out_time, time_count, i, day_count
      real, allocatable, dimension (:,:,:)           :: var_subset1, var_subset2, out_ave, sib_ave
      real, allocatable, dimension (:,:)             :: daily_array, sibvis, sibnir, sibvis_ave, sibnir_ave
      logical                                        :: last
      real                                           :: time_frac
      logical                                        :: save_lai=.false., save_green = .false.

      if(LDT_rc%lsm.eq."CLSMJ3.2")then
         select case (trim(CP%SOURCE))
         case ('GSWPH')
            save_green = .true.
         case ('GLASSA','GLASSM','MCD15A2H')
            save_lai= .true.      
         case default         
            write(LDT_logunit,*)'[WARNING] Not saving for SiBALB calc though CLSM is the model : '
            write(LDT_logunit,*)'[WARNING] clim data source : ', trim (CP%SOURCE)
            write(LDT_logunit,*)'[WARNING] LDT_rc%lsm :', LDT_rc%lsm         
         end select
      endif
      
      ! initialize time loop
      ! --------------------
      
      call ESMF_TimeSet (CURRENT_TIME, yy=ref_year, mm=1, dd=1, rc=status)     ; VERIFY_(STATUS)
      call ESMF_TimeSet (END_TIME    , yy=ref_year + 1, mm=1, dd=1, rc=status) ; VERIFY_(STATUS)
      call ESMF_TimeIntervalSet (DAY_DT, h=24, rc=status)                      ; VERIFY_(STATUS)
      
      OUT_RING  = CURRENT_TIME
      t1        = CP%NTIMES
      t2        = 1
      out_time  = 1
      time_count= 1
      day_count = 1
      switch    = 0
      last      = .false. 
      TIME1     = MidTime (CP%data_doys(t1),CP%data_doys(1),start=.true.)
      TIME2     = MidTime (CP%data_doys(1) ,CP%data_doys(2))
      
      allocate (var_subset1 (CP%subpnc,CP%subpnr,NF))
      allocate (var_subset2 (CP%subpnc,CP%subpnr,NF))
      allocate (out_ave     (CP%subpnc,CP%subpnr,NF))
      allocate (daily_array (CP%subpnc,CP%subpnr))
      call CP%update_data (t1, nest, var_subset1)
      call CP%update_data (t2, nest, var_subset2)
      out_ave = 0.
      
      if(present (clsm)) then
         allocate (sibvis  (CP%subpnc,CP%subpnr))
         allocate (sibnir  (CP%subpnc,CP%subpnr))
         allocate (sib_ave (CP%subpnc,CP%subpnr,NF))
         sib_ave = 0.
      endif
      
      ! Annual loop to compute climatology at out_interval on the lis grid - this loops through non-leap ref_year 
      
      ANNUAL_LOOP : do while (CURRENT_TIME  <= END_TIME)

         UPDATE_BCS: if (CURRENT_TIME >= TIME2) then

           ! ---------------------------------------------------------------
           ! Update the interpolation limits for climatological data
           ! and get their midmonth times
           ! ---------------------------------------------------------------

            ! update time infor
            if(t2 == 1) t2 = 2
            t1 = t2
            t2 = t2 + 1

            if(t2 > CP%NTIMES) then
               t2 = 1
               t1 = CP%NTIMES 
               switch = 1 
            endif
            
            TIME1 = TIME2
            if(last) then 
               TIME2 = MidTime (CP%data_doys(1),CP%data_doys(2), last = last)
               t2 = 1
            else
               TIME2 = MidTime (CP%data_doys(t1),CP%data_doys(t2))
            endif
            if (switch == 1) last = .true. 
            
            ! read BCs data 
            var_subset1 = var_subset2
            call CP%update_data (t2, nest, var_subset2)

         endif UPDATE_BCS

         ! sum of daily values to compute CLIMDATA data at out_interval
         time_frac = TimeFrac (TIME1, CURRENT_TIME, TIME2)
         do i = 1, NF
            daily_array (:,:) = (var_subset1(:,:,i) * time_frac + & 
                 var_subset2(:,:,i) * (1. - time_frac))
            out_ave(:,:,i) = out_ave(:,:,i) + daily_array
         end do
         
         if(present (clsm)) then
            ! SiB scaling parameters for CLSM
            call sib_scale_params (day_count,sibvis, sibnir)
            sib_ave(:,:,1) = sib_ave(:,:,1) + sibvis
            sib_ave(:,:,2) = sib_ave(:,:,2) + sibnir
         else
            ! NF is 1 for below 2 cases
            if(save_green) call CP%update_sibinputs (nest, GREEN=daily_array, tstep =day_count)
            if(save_lai  ) call CP%update_sibinputs (nest, LAI=daily_array,   tstep =day_count)
         endif

         PROC_OUTPPUT: if (CURRENT_TIME == OUT_RING) then

            if (out_time > NT) then
               write(LDT_logunit,*)'[ERR] out_interval mismatch in READ_CLIM_BCS'
               write(LDT_logunit,*)'  out_intervals : ',NF
               VERIFY_(1)
               call LDT_endrun
            endif

            OUT_RING   = CURRENT_TIME + CP%out_dtstep
            do i = 1, NF
               if(present (clsm)) then
                  out_ave(:,:,i) = out_ave(:,:,i) / real(time_count)
                  sib_ave(:,:,i) = sib_ave(:,:,i) / real(time_count)
                  out_ave(:,:,i) = out_ave(:,:,i) / sib_ave(:,:,i)
                  where (out_ave(:,:,i) <   0.) out_ave(:,:,i) = 1.
                  where (out_ave(:,:,i) > 100.) out_ave(:,:,i) = 1.
                  cdata(i)%gdata(:,:,out_time) = out_ave(:,:,i)
                  sib_ave(:,:,i) = 0.
               else
                  cdata(i)%gdata(:,:,out_time) = out_ave(:,:,i) / real(time_count)
               endif
            end do
            out_time   = out_time + 1
            time_count = 0
            out_ave    = 0.

         end if PROC_OUTPPUT
         
         day_count    = day_count  + 1
         time_count   = time_count + 1
         CURRENT_TIME = CURRENT_TIME + DAY_DT
         
      end do ANNUAL_LOOP

      if(present (clsm)) deallocate (sibinputs%ITYP, sibinputs%LAI, sibinputs%GREEN, sibvis, sibnir)

    end subroutine run_annual_cycle

    ! --------------------------------------------------------------
    
    SUBROUTINE update_sibinputs (this, nest, ITYP, LAI, GREEN, tstep)

      implicit none
      class (ClimateBCsReader), intent(inout)        :: this
      integer, intent (in)                           :: nest
      integer, dimension (:,:), intent(in), optional :: ITYP
      real,    dimension (:,:), intent(in), optional :: LAI, GREEN
      integer, optional                              :: tstep
      integer :: nx,ny

      if (present (ITYP)) then
         NX = SIZE (ITYP,1)
         NY = SIZE (ITYP,2)
         ASSERT_ (NX == LDT_rc%lnc(nest))
         ASSERT_ (NY == LDT_rc%lnr(nest))
         if(.not.allocated (sibinputs%ITYP)) allocate (sibinputs%ITYP(1:NX, 1:NY))
         sibinputs%ITYP = ITYP
      endif

      if (present (LAI)) then
         NX = SIZE (LAI,1)
         NY = SIZE (LAI,2)
         ASSERT_ (NX == LDT_rc%lnc(nest))
         ASSERT_ (NY == LDT_rc%lnr(nest))
         if(.not.allocated (sibinputs%LAI)) allocate (sibinputs%LAI(1:NX, 1:NY,1:366))
         sibinputs%LAI(:,:,tstep) = LAI         
      endif

      if (present (GREEN)) then
         NX = SIZE (GREEN,1)
         NY = SIZE (GREEN,2)
         ASSERT_ (NX == LDT_rc%lnc(nest))
         ASSERT_ (NY == LDT_rc%lnr(nest))
         if(.not.allocated (sibinputs%GREEN)) allocate (sibinputs%GREEN(1:NX, 1:NY,1:366))
         sibinputs%GREEN(:,:,tstep) = GREEN         
      endif

    end SUBROUTINE update_sibinputs
    
    ! ###################################################################

    SUBROUTINE sib_scale_params (day,AVISDF,ANIRDF)
    
      implicit none
      INTEGER, INTENT (IN)  :: DAY
      REAL, DIMENSION (:,:) :: AVISDF,ANIRDF
      INTEGER :: NTILES, NX, NY
      logical, dimension (:,:), allocatable :: landmask
      
      NX = SIZE (sibinputs%ITYP,1)
      NY = SIZE (sibinputs%ITYP,2)      
      NTILES = COUNT (sibinputs%ITYP .NE. NINT(LDT_rc%udef))

      allocate (landmask (1:NX, 1:NY))
      
      landmask = .false.
      AVISDF = LDT_rc%udef
      ANIRDF = LDT_rc%udef

      landmask = sibinputs%ITYP .NE. NINT(LDT_rc%udef)
      
      call SIBALB (NTILES,                                  &
           pack(sibinputs%ITYP, mask = landmask),           &
           pack(sibinputs%LAI(:,:,day),  mask = landmask),  &
           pack(sibinputs%GREEN(:,:,day),mask = landmask),  &
           pack(AVISDF,mask = landmask),                    &
           pack(ANIRDF,mask = landmask))

      deallocate (landmask)
      
    contains
      
      SUBROUTINE SIBALB (NTILES, ITYP, VLAI, VGRN, &
           AVISDF, ANIRDF)
      
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: NTILES 
        REAL,    DIMENSION (:), INTENT (IN)   :: VLAI, VGRN
        INTEGER, DIMENSION (:), INTENT (IN  ) :: ITYP
        REAL,    DIMENSION (:) :: AVISDF, ANIRDF !AVISDR, ANIRDR 
        
        REAL, PARAMETER :: ALVDRS = 0.1, ALIDRS = 0.2
        REAL, PARAMETER :: ALVDRD = 0.3, ALIDRD = 0.350
        REAL, PARAMETER :: ALVDRI = 0.7, ALIDRI = 0.7
        REAL, PARAMETER :: ZTH    = 0.0, SNW    = 0.0
    
        ! ALVDRS:  Albedo of soil for visible   direct  solar radiation.
        ! ALIDRS:  Albedo of soil for infra-red direct  solar radiation.
        ! ALVDFS:  Albedo of soil for visible   diffuse solar radiation.
        ! ALIDFS:  Albedo of soil for infra-red diffuse solar radiation.
        
        INTEGER, PARAMETER :: NLAI = 14, NTYPS = 9
        REAL,    PARAMETER :: EPSLN = 1.E-6, BLAI = 0.5, DLAI = 0.5, ZERO=0., ONE=1.0
        REAL,    PARAMETER :: ALATRM = BLAI + (NLAI - 1) * DLAI - EPSLN
        
        ! OUTPUTS:
        ! AVISDR:   visible, direct albedo.
        ! ANIRDR:   near infra-red, direct albedo.
        ! AVISDF:   visible, diffuse albedo.
        ! ANIRDF:   near infra-red, diffuse albedo.
        
        ! INPUTS:
        ! VLAI:     the leaf area index.
        ! VGRN:     the greenness index.
        ! ZTH:      The cosine of the solar zenith angle => assumed 0. (noon) at daily time step
        ! SNW:      Snow cover in meters water equivalent=> assumed 0. (snow free)
        
        ! ITYP: Vegetation type as follows:
        !                  1:  BROADLEAF EVERGREEN TREES
        !                  2:  BROADLEAF DECIDUOUS TREES
        !                  3:  NEEDLELEAF TREES
        !                  4:  GROUND COVER
        !                  5:  BROADLEAF SHRUBS
        !                  6:  DWARF TREES (TUNDRA)
        !                  7:  BARE SOIL
        !                  8:  DESSERT
        !                  9:  ICE
        
        ! LOCAL VARIABLES
        
        INTEGER :: I, LAI
        REAL    :: FAC, GAMMA, BETA, ALPHA, DX, DY, ALA, GRN (2)
        !, SNWALB (4, NTYPS), SNWMID (NTYPS)
        
        !   Constants used in albedo calculations:
        
        REAL :: ALVDR (NLAI, 2, NTYPS)
        REAL :: BTVDR (NLAI, 2, NTYPS)
        REAL :: GMVDR (NLAI, 2, NTYPS)
        REAL :: ALIDR (NLAI, 2, NTYPS)
        REAL :: BTIDR (NLAI, 2, NTYPS)
        REAL :: GMIDR (NLAI, 2, NTYPS)
        
        ! (Data statements for ALVDR described in full; data statements for
        ! other constants follow same framework.)
        
        ! BROADLEAF EVERGREEN (ITYP=4); GREEN=0.33; LAI: .5-7
        DATA (ALVDR (I, 1, 1), I = 1, 14)                                &
             /0.0808, 0.0796, 0.0792, 0.0790, 10*0.0789/
        
        ! BROADLEAF EVERGREEN (ITYP=4); GREEN=0.67; LAI: .5-7
        DATA (ALVDR (I, 2, 1), I = 1, 14)                                &
             /0.0788, 0.0775, 0.0771, 0.0769, 10*0.0768/
        
        ! BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.33; LAI: .5-7
        DATA (ALVDR (I, 1, 2), I = 1, 14)                                &
             /0.0803, 0.0790, 0.0785, 0.0784, 3*0.0783, 7*0.0782/
        
        ! BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.67; LAI: .5-7
        DATA (ALVDR (I, 2, 2), I = 1, 14)                                &
             /0.0782, 0.0770, 0.0765, 0.0763, 10*0.0762/
        
        ! NEEDLELEAF (ITYP=3); GREEN=0.33; LAI=.5-7
        DATA (ALVDR (I, 1, 3), I = 1, 14)                                &
             /0.0758, 0.0746, 0.0742, 0.0740, 10*0.0739/
        
        ! NEEDLELEAF (ITYP=3); GREEN=0.67; LAI=.5-7
        DATA (ALVDR (I, 2, 3), I = 1, 14)                                &
             /0.0683, 0.0672, 0.0667, 2*0.0665, 9*0.0664/
        
        ! GROUNDCOVER (ITYP=4); GREEN=0.33; LAI=.5-7
        DATA (ALVDR (I, 1, 4), I = 1, 14)                                &
             /0.2436, 0.2470, 0.2486, 0.2494, 0.2498, 0.2500, 2*0.2501,  &
             6*0.2502 /
        
        ! GROUNDCOVER (ITYP=4); GREEN=0.67; LAI=.5-7
        DATA (ALVDR (I, 2, 4), I = 1, 14) /14*0.1637/
        
        ! BROADLEAF SHRUBS (ITYP=5); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 5), I = 1, 14)                                &
             /0.0807, 0.0798, 0.0794, 0.0792, 0.0792, 9*0.0791/
        
        ! BROADLEAF SHRUBS (ITYP=5); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 5), I = 1, 14)                                &
             /0.0787, 0.0777, 0.0772, 0.0771, 10*0.0770/
        
        ! DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 6), I = 1, 14)                                &
             /0.0802, 0.0791, 0.0787, 0.0786, 10*0.0785/
        
        ! DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 6), I = 1, 14)                                &
             /0.0781, 0.0771, 0.0767, 0.0765, 0.0765, 9*0.0764/
        
        ! BARE SOIL
        DATA (ALVDR (I, 1, 7), I = 1, 14) /14*ALVDRS/
        DATA (ALVDR (I, 2, 7), I = 1, 14) /14*ALVDRS/
        
        ! DESERT
        DATA (ALVDR (I, 1, 8), I = 1, 14) /14*ALVDRD/
        DATA (ALVDR (I, 2, 8), I = 1, 14) /14*ALVDRD/
        
        ! ICE
        DATA (ALVDR (I, 1, 9), I = 1, 14) /14*ALVDRI/
        DATA (ALVDR (I, 2, 9), I = 1, 14) /14*ALVDRI/
    
        !**** -----------------------------------------------------------
        
        DATA (BTVDR (I, 1, 1), I = 1, 14)                                      &
             /0.0153, 0.0372, 0.0506, 0.0587, 0.0630, 0.0652, 0.0663,               &
             0.0668, 0.0671, 0.0672, 4*0.0673 /
        DATA (BTVDR (I, 2, 1), I = 1, 14)                                      &
             /0.0135, 0.0354, 0.0487, 0.0568, 0.0611, 0.0633, 0.0644,             &
             0.0650, 0.0652, 0.0654, 0.0654, 3*0.0655 /
        DATA (BTVDR (I, 1, 2), I = 1, 14)                                      &
             /0.0148, 0.0357, 0.0462, 0.0524, 0.0554, 0.0569, 0.0576,             &
             0.0579, 0.0580, 0.0581, 0.0581, 3*0.0582 /
        DATA (BTVDR (I, 2, 2), I = 1, 14)                                      &
             /0.0131, 0.0342, 0.0446, 0.0508, 0.0539, 0.0554, 0.0560,             &
             0.0564, 0.0565, 5*0.0566 /
        DATA (BTVDR (I, 1, 3), I = 1, 14)                                      &
             /0.0108, 0.0334, 0.0478, 0.0571, 0.0624, 0.0652, 0.0666,             &
             0.0673, 0.0677, 0.0679, 4*0.0680 /
        DATA (BTVDR (I, 2, 3), I = 1, 14)                                      &
             /0.0034, 0.0272, 0.0408, 0.0501, 0.0554, 0.0582, 0.0597,             &
             0.0604, 0.0608, 0.0610, 4*0.0611 /
        DATA (BTVDR (I, 1, 4), I = 1, 14)                                      &
             /0.2050, 0.2524, 0.2799, 0.2947, 0.3022, 0.3059, 0.3076,             &
             0.3085, 0.3088, 0.3090, 4*0.3091 /
        DATA (BTVDR (I, 2, 4), I = 1, 14)                                      &
             /0.1084, 0.1404, 0.1617, 0.1754, 0.1837, 0.1887, 0.1915,             &
             0.1931, 0.1940, 0.1946, 0.1948, 0.1950, 2*0.1951  /
        DATA (BTVDR (I, 1, 5), I = 1, 14)                                      &
             /0.0203, 0.0406, 0.0548, 0.0632, 0.0679, 0.0703, 0.0716,             &
             0.0722, 0.0726, 0.0727, 0.0728, 0.0728, 0.0728, 0.0729 /
        DATA (BTVDR (I, 2, 5), I = 1, 14)                                      &
             /0.0184, 0.0385, 0.0526, 0.0611,  0.0658, 0.0683, 0.0696,            &
             0.0702, 0.0705, 0.0707, 4*0.0708 /
        DATA (BTVDR (I, 1, 6), I = 1, 14)                                      &
             /0.0199, 0.0388, 0.0494,  0.0554, 0.0584, 0.0599, 0.0606,            &
             0.0609, 0.0611, 5*0.0612  /
        DATA (BTVDR (I, 2, 6), I = 1, 14)                                      &
             /0.0181, 0.0371, 0.0476, 0.0537,  0.0568, 0.0583, 0.0590,            &
             0.0593, 0.0595, 0.0595, 4*0.0596 /
        DATA (BTVDR (I, 1, 7), I = 1, 14) /14*0./
        DATA (BTVDR (I, 2, 7), I = 1, 14) /14*0./
        DATA (BTVDR (I, 1, 8), I = 1, 14) /14*0./
        DATA (BTVDR (I, 2, 8), I = 1, 14) /14*0./
        DATA (BTVDR (I, 1, 9), I = 1, 14) /14*0./
        DATA (BTVDR (I, 2, 9), I = 1, 14) /14*0./
        
        !**** -----------------------------------------------------------
        
        DATA (GMVDR (I, 1, 1), I = 1, 14)                                      &
             /0.0814, 0.1361, 0.2078, 0.2650, 0.2986, 0.3169,  0.3265,            &
             0.3313, 0.3337, 0.3348, 0.3354, 0.3357, 2*0.3358 /
        DATA (GMVDR (I, 2, 1), I = 1, 14)                                      &
             /0.0760, 0.1336, 0.2034, 0.2622, 0.2969, 0.3159,  0.3259,            &
         	   0.3309, 0.3333, 0.3346, 0.3352, 0.3354, 2*0.3356 /
        DATA (GMVDR (I, 1, 2), I = 1, 14)                                      &
             /0.0834, 0.1252, 0.1558, 0.1927, 0.2131,   0.2237, 0.2290,           &
         	   0.2315, 0.2327, 0.2332, 0.2335, 2*0.2336, 0.2337 /
        DATA (GMVDR (I, 2, 2), I = 1, 14)                                      &
             /0.0789, 0.1235, 0.1531, 0.1912, 0.2122, 0.2232,  0.2286,            &
        	   0.2312, 0.2324, 0.2330, 0.2333, 0.2334, 2*0.2335 /
        DATA (GMVDR (I, 1, 3), I = 1, 14)                                      &
             /0.0647, 0.1342, 0.2215, 0.2968, 0.3432, 0.3696, 0.3838,             &
         	   0.3912, 0.3950, 0.3968, 0.3978, 0.3982, 0.3984, 0.3985 /
        DATA (GMVDR (I, 2, 3), I = 1, 14)                                      &
             /0.0258, 0.1227, 0.1999, 0.2825, 0.3339, 0.3634, 0.3794,             &
         	   0.3877, 0.3919, 0.3940, 0.3950, 0.3956, 0.3958, 0.3959 /
        DATA (GMVDR (I, 1, 4), I = 1, 14)                                      &
             /0.3371, 0.5762, 0.7159, 0.7927, 0.8324, 0.8526,  0.8624,            &
         	   0.8671, 0.8693, 0.8704, 0.8709, 0.8710, 2*0.8712 /
        DATA (GMVDR (I, 2, 4), I = 1, 14)                                      &
             /0.2634, 0.4375, 0.5532, 0.6291, 0.6763, 0.7048, 0.7213,             &
         	   0.7310, 0.7363, 0.7395, 0.7411, 0.7420, 0.7426, 0.7428 /
        DATA (GMVDR (I, 1, 5), I = 1, 14)                                      &
             /0.0971, 0.1544, 0.2511, 0.3157, 0.3548, 0.3768, 0.3886,            &
             0.3948, 0.3978, 0.3994, 0.4001, 0.4006, 0.4007, 0.4008 /
        DATA (GMVDR (I, 2, 5), I = 1, 14)                                      &
             /0.0924, 0.1470, 0.2458, 0.3123, 0.3527, 0.3756, 0.3877,            &
             0.3942, 0.3974, 0.3990, 0.3998, 0.4002, 0.4004, 0.4005 /
        DATA (GMVDR (I, 1, 6), I = 1, 14)                                      &
             /0.0970, 0.1355, 0.1841, 0.2230, 0.2447,  0.2561, 0.2617,           &
             0.2645, 0.2658, 0.2664, 0.2667, 3*0.2669 /
        DATA (GMVDR (I, 2, 6), I = 1, 14)                                      &
             /0.0934, 0.1337, 0.1812, 0.2213, 0.2437, 0.2554, 0.2613,            &
             0.2642, 0.2656, 0.2662, 0.2665, 0.2667, 0.2667, 0.2668 /
        DATA (GMVDR (I, 1, 7), I = 1, 14) /14*1./
        DATA (GMVDR (I, 2, 7), I = 1, 14) /14*1./
        DATA (GMVDR (I, 1, 8), I = 1, 14) /14*1./
        DATA (GMVDR (I, 2, 8), I = 1, 14) /14*1./
        DATA (GMVDR (I, 1, 9), I = 1, 14) /14*1./
        DATA (GMVDR (I, 2, 9), I = 1, 14) /14*1./
        
    !***  -----------------------------------------------------------

        DATA (ALIDR (I, 1, 1), I = 1, 14)                                      &
             /0.2867,  0.2840, 0.2828, 0.2822, 0.2819, 0.2818, 2*0.2817,          &
         	   6*0.2816 /
        DATA (ALIDR (I, 2, 1), I = 1, 14)                                      &
             /0.3564, 0.3573, 0.3577, 0.3580, 2*0.3581, 8*0.3582 /
        DATA (ALIDR (I, 1, 2), I = 1, 14)                                      &
             /0.2848, 0.2819, 0.2804, 0.2798, 0.2795, 2*0.2793, 7*0.2792 /
        DATA (ALIDR (I, 2, 2), I = 1, 14)                                      &
             /0.3544, 0.3550, 0.3553, 2*0.3555, 9*0.3556 /
        DATA (ALIDR (I, 1, 3), I = 1, 14)                                      &
             /0.2350, 0.2311, 0.2293, 0.2285, 0.2281, 0.2280, 8*0.2279 /
        DATA (ALIDR (I, 2, 3), I = 1, 14)                                      &
             /0.2474, 0.2436, 0.2418, 0.2410, 0.2406, 0.2405, 3*0.2404,           &
         	   5*0.2403 /
        DATA (ALIDR (I, 1, 4), I = 1, 14)                                      &
             /0.5816, 0.6157, 0.6391, 0.6556, 0.6673, 0.6758, 0.6820,             &
         	   0.6866, 0.6899, 0.6924, 0.6943, 0.6956, 0.6966, 0.6974 /
        DATA (ALIDR (I, 2, 4), I = 1, 14)                                      &
             /0.5489, 0.5770, 0.5955, 0.6079, 0.6163, 0.6221, 0.6261,             &
         	   0.6288, 0.6308, 0.6321, 0.6330, 0.6337, 0.6341, 0.6344 /
        DATA (ALIDR (I, 1, 5), I = 1, 14)                                      &
             /0.2845, 0.2837, 0.2832, 0.2831, 0.2830, 9*0.2829 /
        DATA (ALIDR (I, 2, 5), I = 1, 14)                                      &
             /0.3532, 0.3562, 0.3578,  0.3586, 0.3590, 0.3592, 0.3594,           &
             0.3594, 0.3594, 5*0.3595 /
        DATA (ALIDR (I, 1, 6), I = 1, 14)                                      &
             /0.2825, 0.2812, 0.2806, 0.2803, 0.2802, 9*0.2801 /
        DATA (ALIDR (I, 2, 6), I = 1, 14)                                      &
             /0.3512, 0.3538,  0.3552, 0.3559, 0.3562, 0.3564, 0.3565,           &
             0.3565, 6*0.3566 /
        DATA (ALIDR (I, 1, 7), I = 1, 14) /14*ALIDRS/
        DATA (ALIDR (I, 2, 7), I = 1, 14) /14*ALIDRS/
        DATA (ALIDR (I, 1, 8), I = 1, 14) /14*ALIDRD/
        DATA (ALIDR (I, 2, 8), I = 1, 14) /14*ALIDRD/
        DATA (ALIDR (I, 1, 9), I = 1, 14) /14*ALIDRI/
        DATA (ALIDR (I, 2, 9), I = 1, 14) /14*ALIDRI/
    
   !*** -----------------------------------------------------------
        DATA (BTIDR (I, 1, 1), I = 1, 14)                                      &
             /0.1291, 0.1707, 0.1969, 0.2125, 0.2216,   0.2267, 0.2295,           &
         	   0.2311, 0.2319, 0.2323, 0.2326, 2*0.2327, 0.2328 /
        DATA (BTIDR (I, 2, 1), I = 1, 14)                                      &
             /0.1939, 0.2357, 0.2598, 0.2735, 0.2810,  0.2851, 0.2874,            &
         	   0.2885, 0.2892, 0.2895, 0.2897, 3*0.2898 /
        DATA (BTIDR (I, 1, 2), I = 1, 14)                                      &
             /0.1217, 0.1522, 0.1713, 0.1820,   0.1879,  0.1910, 0.1926,          &
        	   0.1935, 0.1939, 0.1942, 2*0.1943, 2*0.1944 /
        DATA (BTIDR (I, 2, 2), I = 1, 14)                                      &
             /0.1781, 0.2067, 0.2221, 0.2301,   0.2342,  0.2363, 0.2374,          &
         	   0.2379, 0.2382, 0.2383, 2*0.2384, 2*0.2385 /
        DATA (BTIDR (I, 1, 3), I = 1, 14)                                      &
             /0.0846, 0.1299, 0.1614, 0.1814, 0.1935,   0.2004, 0.2043,           &
             0.2064, 0.2076, 0.2082, 0.2085, 2*0.2087, 0.2088 /
        DATA (BTIDR (I, 2, 3), I = 1, 14)                                      &
             /0.0950, 0.1410, 0.1722, 0.1921, 0.2042, 0.2111,  0.2151,            &
         	   0.2172, 0.2184, 0.2191, 0.2194, 0.2196, 2*0.2197 /
        DATA (BTIDR (I, 1, 4), I = 1, 14)                                      &
             /0.5256, 0.7444, 0.9908, 1.2700, 1.5680, 1.8505, 2.0767,             &
         	   2.2211, 2.2808, 2.2774, 2.2362, 2.1779, 2.1160, 2.0564 /
        DATA (BTIDR (I, 2, 4), I = 1, 14)                                      &
             /0.4843, 0.6714, 0.8577, 1.0335, 1.1812, 1.2858, 1.3458,             &
         	   1.3688, 1.3685, 1.3546, 1.3360, 1.3168, 1.2989, 1.2838 /
        DATA (BTIDR (I, 1, 5), I = 1, 14)                                      &
             /0.1498, 0.1930, 0.2201, 0.2364, 0.2460, 0.2514, 0.2544,            &
             0.2560, 0.2569, 0.2574, 0.2577, 0.2578, 0.2579, 0.2579 /
        DATA (BTIDR (I, 2, 5), I = 1, 14)                                      &
             /0.2184, 0.2656, 0.2927, 0.3078, 0.3159,  0.3202, 0.3224,           &
             0.3235, 0.3241, 0.3244, 0.3245, 3*0.3246 /
        DATA (BTIDR (I, 1, 6), I = 1, 14)                                      &
             /0.1369, 0.1681, 0.1860, 0.1958, 0.2010,  0.2038, 0.2053,           &
             0.2060, 0.2064, 0.2066, 0.2067, 3*0.2068 /
        DATA (BTIDR (I, 2, 6), I = 1, 14)                                      &
             /0.1969, 0.2268, 0.2416,  0.2488, 0.2521, 0.2537, 0.2544,           &
             0.2547, 0.2548, 5*0.2549 /
        DATA (BTIDR (I, 1, 7), I = 1, 14) /14*0./
        DATA (BTIDR (I, 2, 7), I = 1, 14) /14*0./
        DATA (BTIDR (I, 1, 8), I = 1, 14) /14*0./
        DATA (BTIDR (I, 2, 8), I = 1, 14) /14*0./
        DATA (BTIDR (I, 1, 9), I = 1, 14) /14*0./
        DATA (BTIDR (I, 2, 9), I = 1, 14) /14*0./
        
    !*** --------------------------------------------------------------
        DATA (GMIDR (I, 1, 1), I = 1, 14)                                      &
             /0.1582, 0.2581, 0.3227, 0.3635, 0.3882, 0.4026, 0.4108,             &
         	   0.4154, 0.4179, 0.4193, 0.4200, 0.4204, 0.4206, 0.4207 /
        DATA (GMIDR (I, 2, 1), I = 1, 14)                                      &
             /0.1934, 0.3141, 0.3818, 0.4200, 0.4415, 0.4533, 0.4598,             &
         	   0.4633, 0.4651, 0.4662, 0.4667, 0.4671, 2*0.4672 /
        DATA (GMIDR (I, 1, 2), I = 1, 14)                                      &
             /0.1347, 0.1871, 0.2277, 0.2515, 0.2651, 0.2727, 0.2768,             &
         	   0.2790, 0.2801, 0.2808, 0.2811, 0.2812, 0.2813, 0.2814 /
        DATA (GMIDR (I, 2, 2), I = 1, 14)                                      &
             /0.1440, 0.2217, 0.2629, 0.2839, 0.2947, 0.3003, 0.3031,             &
         	   0.3046, 0.3054, 0.3058, 0.3060, 2*0.3061, 0.3062 /
        DATA (GMIDR (I, 1, 3), I = 1, 14)                                      &
             /0.1372, 0.2368, 0.3235, 0.3839, 0.4229, 0.4465, 0.4602,             &
         	   0.4679, 0.4722, 0.4745, 0.4758, 0.4764, 0.4768, 0.4770 /
        DATA (GMIDR (I, 2, 3), I = 1, 14)                                      &
             /0.1435, 0.2524, 0.3370, 0.3955, 0.4332, 0.4563, 0.4697,             &
         	   0.4773, 0.4815, 0.4839, 0.4851, 0.4858, 0.4861, 0.4863 /
        DATA (GMIDR (I, 1, 4), I = 1, 14)                                      &
             /0.4298, 0.9651, 1.6189, 2.4084, 3.2992, 4.1928, 4.9611,             &
         	   5.5095, 5.8085, 5.9069, 5.8726, 5.7674, 5.6346, 5.4944 /
        DATA (GMIDR (I, 2, 4), I = 1, 14)                                      &
             /0.4167, 0.8974, 1.4160, 1.9414, 2.4147, 2.7803, 3.0202,             &
        	   3.1468, 3.1954, 3.1932, 3.1676, 3.1328, 3.0958, 3.0625 /
        DATA (GMIDR (I, 1, 5), I = 1, 14)                                      &
             /0.1959, 0.3203, 0.3985, 0.4472, 0.4766, 0.4937, 0.5034,            &
             0.5088, 0.5117, 0.5134, 0.5143, 0.5147, 0.5150, 0.5152 /
        DATA (GMIDR (I, 2, 5), I = 1, 14)                                      &
             /0.2328, 0.3859, 0.4734, 0.5227, 0.5498, 0.5644, 0.5720,            &
             0.5761, 0.5781, 0.5792, 0.5797, 0.5800, 0.5802, 0.5802 /
        DATA (GMIDR (I, 1, 6), I = 1, 14)                                      &
             /0.1447, 0.2244, 0.2698, 0.2953, 0.3094, 0.3170, 0.3211,            &
             0.3233, 0.3244, 0.3250, 0.3253, 0.3255, 0.3256, 0.3256 /
        DATA (GMIDR (I, 2, 6), I = 1, 14)                                      &
             /0.1643, 0.2624, 0.3110, 0.3347, 0.3461, 0.3517, 0.3543,            &
             0.3556, 0.3562, 0.3564, 0.3565, 0.3566, 0.3566, 0.3566 /
        DATA (GMIDR (I, 1, 7), I = 1, 14) /14*1./
        DATA (GMIDR (I, 2, 7), I = 1, 14) /14*1./
        DATA (GMIDR (I, 1, 8), I = 1, 14) /14*1./
        DATA (GMIDR (I, 2, 8), I = 1, 14) /14*1./
        DATA (GMIDR (I, 1, 9), I = 1, 14) /14*1./
        DATA (GMIDR (I, 2, 9), I = 1, 14) /14*1./
        
        DATA GRN /0.33, 0.67/
        
        DO I=1,NTILES
           
           ALA = AMIN1 (AMAX1 (ZERO, VLAI(I)), ALATRM)
           LAI = 1 + MAX(0, INT((ALA-BLAI)/DLAI) )
           DX = (ALA - (BLAI+(LAI-1)*DLAI)) * (ONE/DLAI)
           DY = (VGRN(I)- GRN(1)) * (ONE/(GRN(2) - GRN(1)))
           
           ALPHA = COEFF (ALVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
           BETA  = COEFF (BTVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
           GAMMA = COEFF (GMVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
           
           GAMMA = MAX(GAMMA,0.01)
           
           !	  AVISDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
           AVISDF(I) = ALPHA-BETA                                  &
                + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))
           
           ALPHA = COEFF (ALIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
           BETA  = COEFF (BTIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
           GAMMA = COEFF (GMIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
           
           GAMMA = MAX(GAMMA,0.01)
           
           !	  ANIRDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
           ANIRDF(I) = ALPHA-BETA                                  &
                + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))
           
        END DO
    
      END SUBROUTINE SIBALB
        ! -------------------------------------------------------------------------
      
      REAL FUNCTION COEFF (TABLE, NTABL, LAI ,DX, DY)
    
        implicit none
        INTEGER NTABL, LAI       
        REAL TABLE (NTABL, 2), DX, DY
        
        COEFF = (TABLE(LAI,  1)                                      &
             + (TABLE(LAI  ,2) - TABLE(LAI  ,1)) * DY ) * (1.0-DX) &
             + (TABLE(LAI+1,1)                                     &
             + (TABLE(LAI+1,2) - TABLE(LAI+1,1)) * DY ) * DX
        
      END FUNCTION COEFF
    
  END SUBROUTINE sib_scale_params
        
    ! -------------------------------------------------------------------------

    function MidTime (DOY1, DOY2, start, last)  result (MT)
      
      implicit none
      integer, intent (in)    :: DOY1, DOY2
      type(ESMF_Time)         :: TIME1, MT, CYB
      type(ESMF_TimeInterval) :: DAY_SHIFT
      integer                 :: yyyy
      logical, intent(in), optional :: start, last
      
      if (present(start)) then
         yyyy = ref_year -1
      else
         yyyy = ref_year
      endif
      
      call ESMF_TimeSet (CYB, yy=yyyy, mm=1, dd=1, rc=status)                  ; VERIFY_(STATUS)
      
      if(present(last)) then
         yyyy = ref_year + 1
         call ESMF_TimeSet (CYB, yy=yyyy, mm=1, dd=1, rc=status)               ; VERIFY_(STATUS)
      endif
      
      call ESMF_TimeIntervalSet (DAY_SHIFT, h=24*(DOY1-1), rc=status )         ; VERIFY_(STATUS)
      TIME1 = CYB + DAY_SHIFT
      
      if (DOY2 > DOY1) then
         call ESMF_TimeIntervalSet (DAY_SHIFT, h=24*(DOY2-DOY1)/2, rc=status ) ; VERIFY_(STATUS)
         MT = TIME1 + DAY_SHIFT
      else
         ! climatologies don't account for leap years)
         call ESMF_TimeIntervalSet (DAY_SHIFT, h=24*(366-DOY1)/2, rc=status )  ; VERIFY_(STATUS)
         MT = TIME1 + DAY_SHIFT
         if (DOY2 > 1) then
            call ESMF_TimeIntervalSet (DAY_SHIFT, h=24*(DOY2-1)/2, rc=status ) ; VERIFY_(STATUS)
            MT = MT + DAY_SHIFT
         endif
      endif
          
    end function MidTime
  
    ! ---------------------------------------------------------------------------
    
    real function TimeFrac (PTIME, CTIME, NTIME) result (tf)
      
      implicit none
      type(ESMF_Time),INTENT(IN) :: PTIME, CTIME, NTIME 
      real(ESMF_KIND_R8)         :: tdiff_21, tdiff_20, tfrac      
      type(ESMF_TimeInterval)    :: T21, T20
      
      T21 = NTIME - PTIME
      T20 = NTIME - CTIME
      
      call ESMF_TimeIntervalGet(T21, s_r8=tdiff_21,rc=status) ; VERIFY_(STATUS)
      call ESMF_TimeIntervalGet(T20,s_r8=tdiff_20,rc=status)  ; VERIFY_(STATUS)
      tfrac = tdiff_20/tdiff_21
      tf    = tfrac
      
    end function TimeFrac
               
    ! -------------------------------------------------------------------------
        
  end module LDT_ClimateBCsReader

! /discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001//GSWP2_30sec_VegParam/GSWP2_VegParam_H26V10.nc
! /discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001//MODIS-Albedo2/MCD43GF_wsa_H26V10.nc
! CP%var1     = 'Alb_0.3_0.7'
! CP%var2     = 'Alb_0.7_5.0'
! /discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001//GLASS-LAI/AVHRR.v4/GLASS01B02.V04.AYYYY001.nc4
! /discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001//GLASS-LAI/MODIS.v4/GLASS01B01.V04.AYYYY001.nc4
! /discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/MCD15A2H.006/MCD15A2H.006_LAI_YYYY001.nc4
! /discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/MCD43GF.006/MCD43GF_wsa_001_YYYY_V006.nc4
! /discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/MCD43GF.006/MCD43GF.006_WSA_YYYY001.nc4
! CP%var1     = 'Albedo_Map_0.3-0.7'
! CP%var2     = 'Albedo_Map_0.7-5.0'
