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
! 27 August 2020: Sarith Mahanama; Initial implementation

  use ESMF
  use LDT_timeMgrMod, only : LDT_calendar
  use LDT_gridmappingMod
  use LDT_logMod,  only : LDT_logunit, LDT_endrun
  use LDT_coreMod
  use map_utils,   only : ij_to_latlon
  use CLSM_param_routines, only : sibalb
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
      real, dimension(:,:,:), intent(inout)          :: CLIMDATA
      real, dimension(:,:,:), intent(inout),OPTIONAL :: CLIMDATA2, CLIMDATA3, CLIMDATA4
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
      if (trim(project)== 'latlon' ) this%param_gridDesc(1)  = 0.

      ! ------------------------------------------------------------
      !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN BCS DATA
      ! ------------------------------------------------------------

      !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
      this%subparam_gridDesc = 0.
      call LDT_RunDomainPts( nest, project, this%param_gridDesc(:), &
           this%glpnc, this%glpnr, this%subpnc, this%subpnr, this%subparam_gridDesc, &
           this%lat_line, this%lon_line)

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
    
    SUBROUTINE readDataset_clsm (this, nest, SOURCE, out_interval, project, CLIMDATA1,CLIMDATA2, clsm, maskarray)     

      implicit none

      class (ClimateBCsReader), intent(inout)        :: this
      character(*), intent(in)                       :: SOURCE, out_interval, project
      integer, intent (in)                           :: nest
      logical, intent (in)                           :: clsm
      real, dimension(:,:,:), intent(inout)          :: CLIMDATA1, CLIMDATA2
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
      if (trim (project)== 'latlon') this%param_gridDesc(1)  = 0.

      ! ------------------------------------------------------------
      !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN BCS DATA
      ! ------------------------------------------------------------

      !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
      this%subparam_gridDesc = 0.
      call LDT_RunDomainPts( nest, project, this%param_gridDesc(:), &
           this%glpnc, this%glpnr, this%subpnc, this%subpnr, this%subparam_gridDesc, &
           this%lat_line, this%lon_line)   

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
      integer, parameter  :: N_MODIS_DATES5  = 73

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
           351,356,361/)

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
                  
      case ('MCD43GFv5','MCD43GFv5-CLSM')
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

      case ('MCD43GFv6','MCD43GFv6-CLSM')
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
      CP%param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat 
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
      
      if((trim(CP%source) == 'GLASSA')  .OR. (trim(CP%source) == 'GLASSM') .OR.  &
           (trim(CP%source) == 'MCD15A2H').OR. (trim(CP%source) == 'MCD43GFv6').OR.&
           (trim(CP%source) == 'MCD43GFv6-CLSM'))  &
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
                status = NF90_OPEN(trim(fname),NF90_NOWRITE, ncid)
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
          integer,   allocatable, dimension (:,:) :: data_array, data_array2

          nv = SIZE (datain,3)
          if ((trim(CP%source) == 'MCD43GFv6').OR.(trim(CP%source) == 'MCD43GFv6-CLSM')) then
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
             deallocate (int2_arr)
          else             
             STATUS = NF90_GET_VAR (NCID,VARID(NCID,trim(CP%var1 )),   data_array)             ; VERIFY_(STATUS)
             if (nv == 2) STATUS = NF90_GET_VAR (NCID,VARID(NCID,trim(CP%var2)),   data_array2); VERIFY_(STATUS)
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
          integer,parameter :: k10 = selected_int_kind(10)
          integer(kind=8) :: mi               ! Total number of input param grid array points
          integer   :: mo                     ! Total number of output LIS grid array points
          integer   :: nc, nr, i, j
          integer, allocatable  :: n11(:)     ! Array that maps the location of each input grid
          !   point in the output grid. 
          real,    allocatable  :: gi(:)      ! Input parameter 1d grid
          logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)
          !integer, dimension (3732480000)  :: n11
          !real,    dimension (3732480000)  :: gi
          !logical*1,dimension (3732480000) :: li
          real, allocatable, dimension (:,:)    :: var_in
          real, allocatable, dimension (:)      :: go2     ! Output lis 1d grid
          logical*1, allocatable, dimension (:) :: lo2  ! Output logical mask (to match go)

          mi = INT(dble(CP%NX)*dble(CP%NY), 8)
          mo = LDT_rc%lnc(nest)*LDT_rc%lnr(nest)

          allocate (var_in(CP%subpnc,CP%subpnr))                
          allocate ( li(mi), gi (mi), n11(mi))
          allocate (go2(mo), lo2(mo))

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

             !- Spatial average within each coarse gridcell:
             call upscaleByAveraging ( mi, mo, LDT_rc%udef, n11, li, gi, &
                  lo2(:), go2 (:))
             i = 0
             do nr = 1, LDT_rc%lnr(nest)
                do nc = 1, LDT_rc%lnc(nest)
                   i = i + 1
                   var_subset (nc,nr,j) = go2(i)
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
      logical                                        :: save_lai, save_green

      save_lai   = .false.
      save_green = .false.

      if(LDT_rc%lsm.eq."CLSMJ3.2")then
         select case (trim(CP%SOURCE))
         case ('GSWPH')
            save_green = .true.
            if(allocated (sibinputs%GREEN )) deallocate (sibinputs%GREEN)
            allocate (sibinputs%GREEN(1:LDT_rc%lnc(nest), 1:LDT_rc%lnr(nest),1:366))
         case ('GLASSA','GLASSM','MCD15A2H')
            save_lai= .true.
            if(allocated (sibinputs%LAI )) deallocate (sibinputs%LAI)
            allocate (sibinputs%LAI(1:LDT_rc%lnc(nest), 1:LDT_rc%lnr(nest),1:366))            
         case default         
            write(LDT_logunit,*)'[WARNING] Not saving for SiBALB calc though CLSM is the model : '
            write(LDT_logunit,*)'[WARNING] clim data source : ', trim (CP%SOURCE)
            write(LDT_logunit,*)'[WARNING] LDT_rc%lsm :', LDT_rc%lsm         
         end select
      endif
      
      ! initialize time loop
      ! --------------------
      
      call ESMF_TimeSet (CURRENT_TIME, yy=ref_year, mm=1, dd= 1, rc=status) ; VERIFY_(STATUS)
      call ESMF_TimeSet (END_TIME    , yy=ref_year, mm=12,dd=31, rc=status) ; VERIFY_(STATUS)
      call ESMF_TimeIntervalSet (DAY_DT, h=24, rc=status)                   ; VERIFY_(STATUS)
      
      OUT_RING  = CURRENT_TIME  + CP%out_dtstep
      t1        = CP%NTIMES
      t2        = 1
      out_time  = 1
      time_count= 1
      day_count = 1
      switch    = 0
      last      = .false. 
      TIME1     = MidTime (CP%data_doys(t1),CP%data_doys(1),start=.true.)
      TIME2     = MidTime (CP%data_doys(1) ,CP%data_doys(2))
      
      allocate (var_subset1 (NX,NY,NF))
      allocate (var_subset2 (NX,NY,NF))
      allocate (out_ave     (NX,NY,NF))
      allocate (daily_array (NX,NY))
      call CP%update_data (t1, nest, var_subset1)
      call CP%update_data (t2, nest, var_subset2)
      out_ave = 0.
      
      if(present (clsm)) then
         allocate (sibvis  (NX,NY))
         allocate (sibnir  (NX,NY))
         allocate (sib_ave (NX,NY,NF))
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
            !print *, 'UPDATE_BCS', t1,t2
            !CALL ESMF_TimePrint (TIME1, OPTIONS="string", RC=STATUS )
            !CALL ESMF_TimePrint (TIME2, OPTIONS="string", RC=STATUS )
            ! read BCs data 
            var_subset1 = var_subset2
            call CP%update_data (t2, nest, var_subset2)

         endif UPDATE_BCS

         ! sum of daily values to compute CLIMDATA data at out_interval
         time_frac = TimeFrac (TIME1, CURRENT_TIME, TIME2)
         !print *, time_frac, 1. - time_frac
         do i = 1, NF
            daily_array (:,:) = (var_subset1(:,:,i) * time_frac + & 
                 var_subset2(:,:,i) * (1. - time_frac))            
            out_ave(:,:,i) = out_ave(:,:,i) + daily_array
         end do
         
         if(present (clsm)) then
            ! SiB scaling parameters for CLSM
            call sib_scale_params (nest,day_count,sibvis, sibnir)
            ! print *, day_count, maxval(sibvis), maxval (sibnir)
            sib_ave(:,:,1) = sib_ave(:,:,1) + sibvis
            sib_ave(:,:,2) = sib_ave(:,:,2) + sibnir
         else
            ! NF is 1 for below 2 cases
            if(save_green) call CP%update_sibinputs (nest, GREEN=daily_array, tstep =day_count)
            if(save_lai  ) call CP%update_sibinputs (nest, LAI=daily_array,   tstep =day_count)
         endif

         PROC_OUTPPUT: if ((CURRENT_TIME == OUT_RING).OR.(CURRENT_TIME == END_TIME)) then
            !print *, 'OUT_TIME : ', out_time
            !CALL ESMF_TimePrint (CURRENT_TIME, OPTIONS="string", RC=STATUS )
            if (out_time > NT) then
               write(LDT_logunit,*)'[ERR] out_interval mismatch in READ_CLIM_BCS'
               write(LDT_logunit,*)'  out_intervals : ',NT, out_time
               VERIFY_(1)
               call LDT_endrun
            endif

            OUT_RING   = CURRENT_TIME + CP%out_dtstep
            do i = 1, NF
               if(present (clsm)) then
                  out_ave(:,:,i) = out_ave(:,:,i) / real(time_count)
                  sib_ave(:,:,i) = sib_ave(:,:,i) / real(time_count)
                  out_ave(:,:,i) = out_ave(:,:,i) / (sib_ave(:,:,i) + 1.e-15)
                  where (out_ave(:,:,i) <   0.) out_ave(:,:,i) = 1.
                  where (out_ave(:,:,i) > 100.) out_ave(:,:,i) = 1.
                  cdata(i)%gdata(:,:,out_time) = out_ave(:,:,i)
                  sib_ave(:,:,i) = 0.
 !                 print *, 'SCALE PARAMS :', i, maxval(cdata(i)%gdata(:,:,out_time))
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
 
      if (present (ITYP)) then
         if(allocated (sibinputs%ITYP )) deallocate (sibinputs%ITYP)
         allocate (sibinputs%ITYP(1:LDT_rc%lnc(nest), 1:LDT_rc%lnr(nest)))
         sibinputs%ITYP = ITYP
      endif

      if (present (LAI)  ) sibinputs%LAI  (:,:,tstep) = LAI
      if (present (GREEN)) sibinputs%GREEN(:,:,tstep) = GREEN

    end SUBROUTINE update_sibinputs
    
    ! ###################################################################

    SUBROUTINE sib_scale_params (nest, day,AVISDF,ANIRDF)
    
      implicit none
      INTEGER, INTENT (IN)              :: DAY, nest
      REAL, DIMENSION (:,:),INTENT(OUT) :: AVISDF,ANIRDF
      INTEGER                           :: NTILES, NX, NY
      logical, dimension (:,:), allocatable :: landmask
      real, dimension (:), allocatable      :: VISDF, NIRDF
      real, dimension (:,:), allocatable    :: field
      
      NX = SIZE (sibinputs%ITYP,1)
      NY = SIZE (sibinputs%ITYP,2)      
      NTILES = COUNT (sibinputs%ITYP .NE. NINT(LDT_rc%udef))

      allocate (landmask (1:NX, 1:NY))
      allocate (field    (1:NX, 1:NY))
      allocate (VISDF (1:NTILES))
      allocate (NIRDF (1:NTILES))
      
      landmask = .false.
      AVISDF = LDT_rc%udef
      ANIRDF = LDT_rc%udef

      landmask = sibinputs%ITYP .NE. NINT(LDT_rc%udef)

      call SIBALB (NTILES,                                  &
           pack(sibinputs%ITYP, mask = landmask),           &
           pack(sibinputs%LAI  (:,:,day),mask = landmask),  &
           pack(sibinputs%GREEN(:,:,day),mask = landmask),  &
           VISDF, NIRDF)

      field = 0.
      AVISDF = unpack(VISDF,landmask, field)
      field = 0.
      ANIRDF = unpack(NIRDF,landmask, field)
      deallocate (landmask, visdf, nirdf, field)
     
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
