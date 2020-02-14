module get_DeLannoy_SoilClass

  implicit none

  private

  public mineral_perc,  n_DeLannoy_classes, DeLannoy_class, GDL_center_pix
  
  type :: mineral_perc
     real :: clay_perc
     real :: silt_perc
     real :: sand_perc
  end type mineral_perc

  interface GDL_center_pix
     module procedure center_pix
     module procedure center_pix_int0
  end interface

  integer, parameter :: n_DeLannoy_classes = 253
  contains
  
    !-------------------------------------------------------------
    
    
    INTEGER FUNCTION DeLannoy_class (min_perc)
      
      ! Function returns a unique soil class [1-100], 
      
      IMPLICIT NONE
      type(mineral_perc), intent (in)  :: min_perc
      !real, intent (in) :: clay_perc,silt_perc,sand_perc
      integer :: clay_row, sand_row, silt_row
      
      clay_row = ceiling((100.- min_perc%clay_perc)/10.)
      if(clay_row == 0 ) clay_row = 1
      if(clay_row == 11) clay_row = 10
      
      sand_row = ceiling((min_perc%sand_perc)/10.)
      if(sand_row == 0 ) sand_row = 1
      if(sand_row == 11) sand_row = 10
      
      silt_row = ceiling((min_perc%silt_perc)/10.)
      if(silt_row == 0 ) silt_row = 1
      if(silt_row == 11) silt_row = 10
      
      if(clay_row == 1) DeLannoy_class=1
      
      if(clay_row > 1) DeLannoy_class=   &
           (clay_row - 1)*(clay_row - 1) + (clay_row - sand_row) + silt_row
      
    end FUNCTION DeLannoy_class
    ! --------------------------------------------------------------------------------------------------------

    INTEGER FUNCTION center_pix_int (sf,ktop, ktot, x,y,x0,y0,z0,ext_point)
      
      implicit none
      
      integer (kind =2), dimension (:), intent (in) :: x,y
      integer, intent (in) :: ktop,ktot
      real, intent (in) :: sf
      real :: xi,xj,yi,yj,xx0,yy0,zz0
      real, allocatable, dimension (:,:) :: length_m
      real, allocatable, dimension (:) :: length
      real, intent (inout) :: x0,y0,z0
      integer :: i,j,npix
      logical, intent(in) :: ext_point
      real :: zi, zj
      
      allocate (length_m (1:ktot,1:ktot))
      allocate (length   (1:ktot))
      length_m =0.
      length   =0.
      
      center_pix_int = -9999
      if(ktot /= 0) then
         do i = 1,ktot
            xi = sf*x(i)
            yi = sf*y(i)
            zi = 100. - xi - yi
            if (.not. ext_point) then
               x0 = xi
               y0 = yi
               z0 = zi
            endif
            
            do j = 1,ktot
               xj = sf*x(j)
               yj = sf*y(j)
               zj = 100. - xj - yj
               xx0= xj - x0
               yy0= yj - y0
               zz0= zj - z0
               
               if(ktot > ktop) then 
                  if(j <= ktop) then
                     length_m (i,j) = (xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5
                  else
                     length_m (i,j) = 2.33*((xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5)
                  endif
               else
                  length_m (i,j) = (xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5
               endif
            end do
            length (i) = sum(length_m (i,:))
         end do
         
         center_pix_int = minloc(length,dim=1)
      endif
      
    END FUNCTION center_pix_int
      
    !
    ! ====================================================================
    !
    
    INTEGER FUNCTION center_pix_int0 (sf,ktop, ktot, x,y)
      
      implicit none
      ! sf = 0.01 (integer to real scale factor), ktop = # of pixels in top layer
      ! ktot = total # of pixels, top + subsurface combined
      ! x (clay), y (sand_
      integer (kind =2), dimension (:), intent (in) :: x,y
      integer, intent (in) :: ktop,ktot
      real, intent (in) :: sf
      real :: xi,xj,yi,yj
      real :: length
      
      integer :: i,j,npix
      real :: zi, zj, mindist,xc,yc,zc
      
      length   =0.
      
      center_pix_int0 = -9999
      
      if(ktot /= 0) then
         ! There should be some data pixels
         if(ktot > ktop) then 
            ! Have both layers
            if(ktop > 0) then
               ! There are data in top layer
               xc = sf*0.3*sum(real(x(1:ktop)))/real(ktop) + sf*0.7*sum(real(x(ktop + 1 : ktot)))/real(ktot - ktop)  
               yc = sf*0.3*sum(real(y(1:ktop)))/real(ktop) + sf*0.7*sum(real(y(ktop + 1 : ktot)))/real(ktot - ktop)
            else
               ! There are no data in top layer
               xc = sf*sum(real(x(1:ktot)))/real(ktot)  
               yc = sf*sum(real(y(1:ktot)))/real(ktot)         
            endif
         else
            ! working on Top layer alone
            xc = sf*sum(real(x(1:ktot)))/real(ktot)  
            yc = sf*sum(real(y(1:ktot)))/real(ktot)
         endif
         zc = 100. - xc - yc
      endif
      
      mindist=100000.*100000.
      
      do i = 1,ktot
         xi = sf*x(i)
         yi = sf*y(i)
         zi = 100. - xi - yi
         length = (xi-xc)**2+(yi-yc)**2+(zi-zc)**2
         if(mindist>length)then
            mindist=length
            center_pix_int0=i
         end if
      end do
      !print *,ktop,ktot,center_pix_int0
      
    END FUNCTION center_pix_int0
    !
    ! ====================================================================
    !
    
    INTEGER FUNCTION center_pix (x,y,x0,y0,z0,ext_point)
      
      implicit none

      real, dimension (:), intent (in) :: x,y
      real, allocatable, dimension (:,:) :: length_m
      real, allocatable, dimension (:) :: length
      real, intent (inout) :: x0,y0,z0
      integer :: i,j,npix,ii
      logical, intent(in) :: ext_point
      real :: zi, zj
      
      npix = size (x)
      allocate (length_m (1:npix,1:npix))
      allocate (length   (1:npix))
      length_m =0.
      length   =0.
      
      do i = 1,npix
         zi = 100. - x(i) - y(i)
         if (.not. ext_point) then
            x0 = x(i)
            y0 = y(i)
            z0 = zi
         endif
         
         do j = i,npix
            zj = 100. - x(j) - y(j)
            !      length_m (i,j) = abs (x(j) - x0) + &
            !            abs (y(j) - y0) +  abs (zj - z0)
            !
            length_m (i,j) = ((x(j) - x0)*(x(j) - x0) &
                 +  (y(j) - y0)*(y(j) - y0) &
                 +  (zj - z0)*(zj - z0))**0.5
            length_m (j,i) = length_m (i,j)
         end do
         length (i) = sum(length_m (i,:))
      end do
      
      center_pix = minloc(length,dim=1)
      
    END FUNCTION center_pix

end module get_DeLannoy_SoilClass
