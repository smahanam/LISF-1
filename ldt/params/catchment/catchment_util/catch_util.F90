module LDT_catch_util

  use LDT_coreMod,  only : LDT_rc, LDT_domain

  implicit none
  include 'netcdf.inc'	

  private
  public  LDT_RegridRaster, NC_VarID, GEOS2LIS, LDT_g5map, init_geos2lis_map 

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
     integer, allocatable, dimension (:)   :: ID_LOC
     integer, allocatable, dimension (:)   :: til
     integer, allocatable, dimension (:,:) :: rst

  end type GEOS2LIS

  type (GEOS2LIS), save :: LDT_g5map

  contains

    ! -----------------------------------------------------------

    subroutine init_geos2lis_map 
      
      implicit none
      
      LDT_g5map%NX = 43200
      LDT_g5map%NY = 21600
         
    end subroutine init_geos2lis_map

    ! -----------------------------------------------------------

    subroutine RegridRaster(Rin,Rout)
      
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
 
end module LDT_catch_util
