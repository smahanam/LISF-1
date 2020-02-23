#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"

module CLSM_util

  use LDT_coreMod,  only            : LDT_rc, LDT_domain
  use LDT_numericalMethodsMod, only : LDT_quicksort
  use LDT_constantsMod, ONLY:      &
       RADIUS => LDT_CONST_REARTH, &
       PI => LDT_CONST_PI
  use LDT_paramDataMod

  implicit none
  include 'netcdf.inc'	

  private
  public  LDT_RegridRaster, NC_VarID, GEOS2LIS, LDT_g5map, init_geos2lis_mapping, SRTM_maxcat, & 
       nc_g5_rst, nr_g5_rst, G5_BCSDIR, G52LIS, LISv2g, write_clsm_files, histogram
  
  interface LDT_RegridRaster
     module procedure RegridRaster
     module procedure RegridRaster1
     module procedure RegridRaster2
     module procedure RegridRasterReal
  end interface

  type :: GEOS2LIS

     integer :: NX
     integer :: NY
     integer :: NT_LIS
     integer :: NT_GEOS
     real,    allocatable, dimension (:)   :: lat, lon
     integer, allocatable, dimension (:)   :: ID_LOC, catid_index
     integer, allocatable, dimension (:,:) :: rst

  end type GEOS2LIS

  integer  , parameter  :: nc_esa = 129600, nr_esa = 64800, SRTM_maxcat = 291284, nc_g5_rst = 43200, nr_g5_rst = 21600
  character*300         :: G5_BCSDIR = '/discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001/'
  type (GEOS2LIS), save :: LDT_g5map
  logical               :: write_clsm_files = .true.

  contains

    ! -----------------------------------------------------------

    subroutine init_geos2lis_mapping 

      implicit none

      real                             :: dx, dy, x0, y0
      integer                          :: i, j, n, r, c, status, ncid, dx_esa, dy_esa, NBINS, &
           NPLUS, nc_global,nr_global,msk2rst, catNo,ix1,ix2,iy1,iy2, ii, jj, catCount,ncells,  glpnr, glpnc
      real   , dimension (:), allocatable           :: lat, lon, lat_g5, lon_g5
      integer, allocatable, target, dimension (:,:) :: geos_msk, high_msk 
      real (kind =8)                                :: dxh, dyh
      logical                                       :: regrid, counted
      REAL,    allocatable, DIMENSION (:)           :: loc_val
      INTEGER, ALLOCATABLE, DIMENSION (:)           :: density, loc_int
      logical, dimension (:), allocatable           :: unq_mask 
      logical                             :: tile_found
      logical, allocatable, dimension(:)  :: mask
      integer, allocatable, dimension (:) :: sub_tid, id_loc, tid_geos 
      real   , allocatable, dimension (:) :: sub_lon, sub_lat, rev_dist
      character(10)                       :: tmpstr
      real                                :: dw, dx2, dy2, min_lon, max_lon, min_lat, max_lat     
      real                                :: param_grid(20)

      if (write_clsm_files) call system('mkdir -p LDT_clsm')

      ! Talk to Kristi about the right full domain
      ! ------------------------------------------

      i = 1 ! for now 1 nest
      param_grid(:) = LDT_rc%mask_gridDesc(i,:)
      glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
      glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1      
      
      dx = LDT_rc%mask_gridDesc(i, 9)  ! LDT_fileIOMod.F90 L702
      dy = LDT_rc%mask_gridDesc(i,10)  ! LDT_fileIOMod.F90 L708     
            
      LDT_g5map%NX = nc_g5_rst
      LDT_g5map%NY = nr_g5_rst
      LDT_g5map%NT_LIS = NINT(SUM (LDT_LSMparam_struc(i)%landmask%value(:, :,i)))
      
      allocate (lat    (1: LDT_g5map%NT_LIS))
      allocate (lon    (1: LDT_g5map%NT_LIS))

      n = 0

      do r = 1, glpnr
         do c = 1, glpnc
            if( LDT_rc%global_mask(c,r) > 0. ) then 
               n = n + 1
               lat(n) = LDT_LSMparam_struc(i)%xlat%value(c,r,i)
               lon(n) = LDT_LSMparam_struc(i)%xlon%value(c,r,i)               
            endif
         end do
      end do
      
      ! (1) rst/til arrays for GEOS5
      ! ----------------------------
      
      status    = NF_OPEN (trim(G5_BCSDIR)//'GEOS5_10arcsec_mask.nc', NF_NOWRITE, ncid) ; VERIFY_(STATUS)

      
      x0 = -180. + dx/2.
      y0 = -90.  + dy/2.
      nc_global = nint(360./dx)
      nr_global = nint(180./dy)
      dx_esa = ceiling(real(nc_esa) / real(nc_global)) ! x-dimension (or # of ESA columns within the raster grid cell)
      dy_esa = ceiling(real(nr_esa) / real(nr_global)) ! y-dimension (or # of ESA rows within the raster grid cell)
      msk2rst = nc_esa/nc_g5_rst
      
      allocate(geos_msk              (1: dx_esa,1: dy_esa))
      allocate(high_msk              (1:msk2rst,1:msk2rst))
      allocate(LDT_g5map%rst         (1:nc_g5_rst, 1: nr_g5_rst))
      allocate(LDT_g5map%catid_index (1:LDT_g5map%NT_LIS)) 
      allocate(LDT_g5map%ID_LOC      (1:LDT_g5map%NT_LIS))
      allocate(LDT_g5map%lat         (1:LDT_g5map%NT_LIS))
      allocate(LDT_g5map%lon         (1:LDT_g5map%NT_LIS))
      
      LDT_g5map%rst(:,:)       = 0
      LDT_g5map%catid_index(:) = 0
      
      dxh = 360.d0/nc_g5_rst
      dyh = 180.d0/nr_g5_rst
      
      catCount = 0
      NCELLS   = 1
      
      LDT_g5map%NT_GEOS = LDT_g5map%NT_LIS
      allocate (lat_g5 (1:LDT_g5map%NT_LIS))
      allocate (lon_g5 (1:LDT_g5map%NT_LIS))
      
      do catNo = 1,  LDT_g5map%NT_LIS

         i = nint((lon(catNo) - x0)/dx) + 1
         j = nint((lat(catNo) - y0)/dy) + 1   
         
         ix1 = NINT ((lon(catNo) -dx/2. + 180.)/dxh) + 1 
         ix2 = NINT ((lon(catNo) +dx/2. + 180.)/dxh) 
         
         if(ix1 < 0) ix1 = 1 ! a vertical strip of 0.5 dx that lies on the Eastern hemesphere in DC grid is discarded
         
         iy1 = NINT ((lat(catNo) - dy/2. + 90.)/dyh) + 1
         iy2 = NINT ((lat(catNo) + dy/2. + 90.)/dyh) 
         
         counted = .false.
         
         do jj = iy1,iy2
            do ii = ix1,ix2
               status  = NF_GET_VARA_INT (ncid,4,(/(ii-1)*msk2rst+1,(jj-1)*msk2rst +1/),  &
                    (/msk2rst,msk2rst/),high_msk) 
               if(count(high_msk >= 1 .and. high_msk <= SRTM_maxcat) > 0) then
                  if (.not. counted) catCount       = catCount + 1
                  LDT_g5map%rst (ii,jj) = catCount
                  counted = .true.
               endif
            end do
         end do
         
         if (counted) then
            status  = NF_GET_VARA_INT (ncid,4,(/(ix1-1)*msk2rst+1,(iy1-1)*msk2rst +1/), &
                 (/dx_esa,dy_esa/),geos_msk)  ; VERIFY_(STATUS)
            
            if(maxval (geos_msk) > SRTM_maxcat) then
               if(maxval (geos_msk) == 200000000) where (geos_msk == 200000000) geos_msk = SRTM_maxcat + 2 
               if(maxval (geos_msk) == 190000000) where (geos_msk == 190000000) geos_msk = SRTM_maxcat + 1
            endif
            
            NPLUS = count(geos_msk >= 1 .and. geos_msk <= SRTM_maxcat) ! Count non-ocean ESA pixels within 
            
            if (NPLUS > 0) then ! check whether there are Non-ocean ESA pixels 
               ! catID of the largest Pfafstetter catchment within the grid cell                
               allocate (loc_int (1:NPLUS))
               allocate (unq_mask(1:NPLUS))
               loc_int = pack(geos_msk,mask = (geos_msk >= 1 .and. geos_msk <= SRTM_maxcat)) ! loc_int contains catch_indices of non-ocean ESA pixels 
               call LDT_quicksort (loc_int)
               unq_mask = .true.
               do n = 2,NPLUS 
                  unq_mask(n) = .not.(loc_int(n) == loc_int(n-1)) ! count number of unique numbers in loc_int for binning
               end do
               NBINS = count(unq_mask)
               
               if (NBINS >= 1) then
                  allocate(loc_val (1:NBINS))
                  allocate(density (1:NBINS))
                  loc_val = 1.*pack(loc_int,mask =unq_mask)                               ! loc_val contains available non-ocean catch_indices within the i,j grid cell,
                  ! Those numbers will be used as bin values
                  call histogram (dx_esa*dy_esa, NBINS, density, loc_val, real(geos_msk)) ! density is the pixel count for each bin value
                  LDT_g5map%catid_index(ncells) = loc_val (maxloc(density,1))                               ! picks maximum density as the dominant catchment_index at i,j
                  
                  deallocate (loc_val, density)
               else
                  print *,'Check suspicous'
                  LDT_g5map%catid_index(ncells) = loc_int (1)
               endif
               deallocate (loc_int, unq_mask)
               
            endif
            lon_g5 (ncells) = lon (CatNo)
            lat_g5 (ncells) = lat (CatNo)
            LDT_g5map%lat(ncells) = lat (CatNo)
            LDT_g5map%lon(ncells) = lon (CatNo)
            ncells = ncells + 1
!        write(20,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)') 100,catid_index,lon,lat,i-i_offset,j-j_offset,    &
!             real(count(tileid (ix1:ix2,iy1:iy2) > 0))/real((ix2-ix1+1)*(iy2-iy1+1)),catid_index,da*sum(pfaf_area (1:nbins)),  &
!             SRTM_catid(catid_index), catNo 
         else
            LDT_g5map%NT_GEOS = LDT_g5map%NT_GEOS -1
         endif
      end do

      if (write_clsm_files) then
         open (10, file = 'LDT/clsm/CLSM.til', form = 'formatted', action = 'write')
         do n = 1,  LDT_g5map%NT_GEOS
            write (10, '(2i10, 2f10.4)') n, LDT_g5map%catid_index(n), LDT_g5map%lon(n), LDT_g5map%lat(n)
         end do
         close (10, status = 'keep')
      endif

      ncells = ncells - 1

      if (LDT_g5map%NT_GEOS == LDT_g5map%NT_LIS) then 

         ! NO Gaps, GEOS5 and LIS masks are consistent.
         do n = 1,  LDT_g5map%NT_LIS
            LDT_g5map%ID_LOC(n) = n
         end do

      else

         ! (2) gap filling
         ! ---------------
         
         LDT_g5map%ID_LOC(:) = 9999
         
         LIS_TILES : do n = 1,  LDT_g5map%NT_LIS
            
            dw = 0.5
            
            ZOOMOUT : do  
               
               tile_found = .false. 
               
               ! Min/Max lon/lat of the working window
               ! -------------------------------------
               
               min_lon = MAX(lon (n) - dw, -180.)
               max_lon = MIN(lon (n) + dw,  180.)
               min_lat = MAX(lat (n) - dw,  -90.)
               max_lat = MIN(lat (n) + dw,   90.) 
               mask = .false.
               mask =  ((lat_g5 >= min_lat .and. lat_g5 <= max_lat).and.(lon_g5 >= min_lon .and. lon_g5 <= max_lon))
               nplus =  count(mask = mask)
               
               if(nplus < 0) then
                  dw = dw + 0.5
                  CYCLE
               endif
               
               allocate (sub_tid (1:nplus))
               allocate (sub_lon (1:nplus))
               allocate (sub_lat (1:nplus))
               allocate (rev_dist  (1:nplus))
               
               sub_tid = PACK (tid_geos    , mask= mask) 
               sub_lon = PACK (lon_g5    , mask= mask)
               sub_lat = PACK (lat_g5    , mask= mask)
               
               ! compute distance from the tile
               
               sub_lat = sub_lat * PI/180.
               sub_lon = sub_lon * PI/180.
               
               SEEK : if(LDT_g5map%ID_LOC(n) < 0) then
                  
                  rev_dist  = 1.e20
                  
                  do i = 1,nplus
                     
                     rev_dist(i) = haversine(to_radian(lat(n)), to_radian(lon(n)), &
                          sub_lat(i), sub_lon(i))
                     
                  end do
                  
                  FOUND : if(minval (rev_dist) < 1.e19) then               
                     LDT_g5map%ID_LOC(n) = sub_tid(minloc(rev_dist,1)) 
                     tile_found = .true.                  
                  endif FOUND
                  
               endif SEEK
               
               deallocate (sub_tid, sub_lon, sub_lat, rev_dist)
               
               if(tile_found) GO TO 100
               
               ! if not increase the window size
               dw = dw + 0.5
               
            end do ZOOMOUT
            
100         continue
            
            !     write (10,'(2I6, 4f8.2)') n, Id_loc(n), lis_lat(n), geos_lat( Id_loc(n)),  lis_lon(n), geos_lon( Id_loc(n))
            
         END do LIS_TILES
      endif

      deallocate (lat, lon)
            
    end subroutine init_geos2lis_mapping

    ! -----------------------------------------------------------

    subroutine RegridRaster(Rin,Rout)

      implicit none

      integer, intent(IN)  :: Rin(:,:)
      integer, intent(OUT) :: Rout(:,:)
      
      REAL(KIND=8)  :: xx, yy
      integer :: i,j,ii,jj
      
      xx = size(Rin ,1)/float(size(Rout,1))
      yy = size(Rin ,2)/float(size(Rout,2))
      
      do j=1,size(Rout,2)
         jj = (j-1)*yy + 1
         do i=1,size(Rout,1)
            ii = (i-1)*xx + 1
            Rout(i,j) = Rin(ii,jj)
         end do
      end do
      
    end subroutine RegridRaster
    
    ! --------------------------------------------------------------
    
    subroutine RegridRaster1(Rin,Rout)
      
      integer*1, intent(IN)  :: Rin(:,:)
      integer*1, intent(OUT) :: Rout(:,:)
      
      REAL(KIND=8)  :: xx, yy
      integer :: i,j,ii,jj
      
      xx = size(Rin ,1)/float(size(Rout,1))
      yy = size(Rin ,2)/float(size(Rout,2))
      
      do j=1,size(Rout,2)
         jj = (j-1)*yy + 1
         do i=1,size(Rout,1)
            ii = (i-1)*xx + 1
            Rout(i,j) = Rin(ii,jj)
         end do
      end do
      
    end subroutine RegridRaster1
        
    ! --------------------------------------------------------------
    
    subroutine RegridRaster2(Rin,Rout)

      integer(kind=2), intent(IN)  :: Rin(:,:)
      integer(kind=2), intent(OUT) :: Rout(:,:)
      
      REAL(KIND=8)  :: xx, yy
      integer :: i,j,ii,jj
      
      xx = size(Rin ,1)/float(size(Rout,1))
      yy = size(Rin ,2)/float(size(Rout,2))
      
      do j=1,size(Rout,2)
         jj = (j-1)*yy + 1
         do i=1,size(Rout,1)
            ii = (i-1)*xx + 1
            Rout(i,j) = Rin(ii,jj)
         end do
      end do
    end subroutine RegridRaster2
    
    ! --------------------------------------------------------------
    
    subroutine RegridRasterReal(Rin,Rout)
      
      real, intent(IN)  :: Rin(:,:)
      real, intent(OUT) :: Rout(:,:)
      
      REAL(KIND=8) :: xx, yy
      integer :: i,j,ii,jj
      
      xx = size(Rin ,1)/float(size(Rout,1))
      yy = size(Rin ,2)/float(size(Rout,2))
      
      do j=1,size(Rout,2)
         jj = (j-1)*yy + 1
         do i=1,size(Rout,1)
            ii = (i-1)*xx + 1
            Rout(i,j) = Rin(ii,jj)
         end do
      end do
      
    end subroutine RegridRasterReal
            
    ! ------------------------------------------------------------
    
    integer function NC_VarID (NCFID, VNAME) 
      
      integer, intent (in)      :: NCFID
      character(*), intent (in) :: VNAME
      integer                   :: status
      
      STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,NC_VarID)
      IF (STATUS .NE. NF_NOERR) &
           CALL HANDLE_ERR(STATUS, trim(VNAME))  
      
    end function NC_VarID
    
    ! ------------------------------------------------------------
    
    SUBROUTINE HANDLE_ERR(STATUS, Line)
      
      INTEGER,      INTENT (IN) :: STATUS
      CHARACTER(*), INTENT (IN) :: Line
      
      IF (STATUS .NE. NF_NOERR) THEN
         PRINT *, trim(Line),': ',NF_STRERROR(STATUS)
         STOP 'Stopped'
      ENDIF
      
    END SUBROUTINE HANDLE_ERR

    !----------------------------------------------------------------------
    
    SUBROUTINE HISTOGRAM (NLENS, NBINS, density, loc_val, x, BIN)
      
      implicit none
      
      integer, intent (in) :: NBINS, NLENS
      real,    intent (in) :: x (NLENS)
      integer, intent (out):: density (NBINS)
      real,    intent (inout) :: loc_val (NBINS)
      real,    intent (in), optional :: bin
      real :: xdum(NLENS), xl, xu, min_value
      integer :: n
      
      if(present (bin)) min_value =  real(floor(minval(x)))
      
      DO N = 1, NBINS
         if(present (bin)) then
            xl = (N - 1)*BIN + min_value
            loc_val (n) = xl
            xu = xl + bin
            XDUM = 0.
            where((x >= xl).and.(x < xu))XDUM = 1
         else
            XDUM = 0.
            where(x == loc_val (n)) XDUM = 1
         endif
         density(n) = int(sum(XDUM))       
    END DO
    
  END SUBROUTINE HISTOGRAM

   ! *****************************************************************************

   function to_radian(degree) result(rad)

     ! degrees to radians
     real,intent(in) :: degree
     real :: rad

     rad = degree*PI/180.

   end function to_radian

   ! *****************************************************************************
   
   real function haversine(deglat1,deglon1,deglat2,deglon2)
     ! great circle distance -- adapted from Matlab 
     real,intent(in) :: deglat1,deglon1,deglat2,deglon2
     real :: a,c, dlat,dlon,lat1,lat2
     
!     dlat = to_radian(deglat2-deglat1)
!     dlon = to_radian(deglon2-deglon1)
     !     lat1 = to_radian(deglat1)
!     lat2 = to_radian(deglat2)
     dlat = deglat2-deglat1
     dlon = deglon2-deglon1
     lat1 = deglat1
     lat2 = deglat2     
     a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
     if(a>=0. .and. a<=1.) then
        c = 2*atan2(sqrt(a),sqrt(1-a))
        haversine = radius*c / 1000.
     else
        haversine = 1.e20
     endif
   end function

!
! ====================================================================
!
   function G52LIS (g5_array) result (lis_array)

     implicit none

     real, dimension (:), intent (in)   :: g5_array
     real, dimension (LDT_g5map%NT_LIS) :: lis_array
     integer :: n

     do n = 1, LDT_g5map%NT_LIS
        lis_array (n) = g5_array(LDT_g5map%ID_LOC(n))
     end do
     
   end function G52LIS
   
!
! ====================================================================
!
   
   function LISv2g (nc,nr, lis_array) result (lis_2D)

     implicit none 
     integer, intent (in)             :: nc, nr 
     real, dimension (:), intent (in) :: lis_array
     real, dimension (nc,nr)          :: lis_2D
     integer                          :: c, r, i
     
     i = 1
     do r = 1, nr
        do c = 1, nc
           if( LDT_rc%global_mask(c,r) > 0. ) then
              lis_2D (c,r) = lis_array(i)
              i = i + 1
           endif
        end do
     end do

   end function LISv2g
 
end module CLSM_util
