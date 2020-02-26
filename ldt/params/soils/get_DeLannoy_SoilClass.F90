module get_DeLannoy_SoilClass

  use CLSM_util, only : c_data => G5_BCSDIR
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber

  implicit none

  private

  public mineral_perc,  n_DeLannoy_classes, DeLannoy_class, GDL_center_pix, GDL_TABLE, OC_LIMITS, &
       get_GDL_soil_table

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
  character*300, parameter       :: GDL_TABLE = 'SoilClasses-SoilHyd-TauParam.dat' 
  real, parameter, dimension (5) :: OC_LIMITS = (/0., 0.4, 0.64, 15./1.72, 100.0/)
  integer, dimension (3)         :: nsoil_pcarbon = (/84, 2*84, 3*84/)

  contains
  
    !-------------------------------------------------------------

    subroutine get_GDL_soil_table (table_map,a_sand, a_clay,a_silt,a_oc,a_bee,a_psis, &
         a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf, atau_2cm,btau_2cm)

      implicit none
      integer, dimension (100,3), intent (inout) :: table_map
      real, dimension(:),pointer, intent (inout) :: a_sand,a_clay,a_silt,a_oc,a_bee,a_psis, &
           a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf, atau_2cm,btau_2cm 

      character*100 :: fout
      integer       :: i,j,k,n,ftbl, n_SoilClasses = n_DeLannoy_classes
      type (mineral_perc)        :: min_percs
     
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
      table_map = 0

      ftbl = LDT_getNextUnitNumber()
      open (ftbl, file=trim(c_data)//trim(GDL_TABLE), form='formatted',status='old', action = 'read')      
      read (ftbl,'(a)')fout

      do n =1,n_SoilClasses 
         read (ftbl,'(4f7.3,4f8.4,e13.5,2f12.7,2f8.4,4f12.7)')a_sand(n),a_clay(n),a_silt(n),a_oc(n),a_bee(n),a_psis(n), &
              a_poros(n),a_wp(n),a_aksat(n),atau(n),btau(n),a_wpsurf(n),a_porosurf(n),atau_2cm(n),btau_2cm(n)
         
         min_percs%clay_perc = a_clay(n)
         min_percs%silt_perc = a_silt(n)
         min_percs%sand_perc = a_sand(n)
         if(n <= nsoil_pcarbon(1))                              table_map(DeLannoy_class (min_percs),1) = n  
         if((n > nsoil_pcarbon(1)).and.(n <= nsoil_pcarbon(2))) table_map(DeLannoy_class (min_percs),2) = n  
         if((n > nsoil_pcarbon(2)).and.(n <= nsoil_pcarbon(3))) table_map(DeLannoy_class (min_percs),3) = n 
         
      end do
      close (ftbl,status='keep') 
      call LDT_releaseUnitNumber(ftbl)

      !  When Woesten Soil Parameters are not available for a particular Soil Class
      !  ,as assumed by tiny triangles in HWSD soil triangle, Woesten Soil
      !  parameters from the nearest available tiny triangle will be substituted.
      
      do n =1,10
         do k=1,n*2 -1
            
            min_percs%clay_perc = 100. -((n-1)*10 + 5)
            min_percs%sand_perc = 100. -  min_percs%clay_perc -2.-(k-1)*5.
            min_percs%silt_perc = 100. -  min_percs%clay_perc - min_percs%sand_perc
            
            i = DeLannoy_class (min_percs)
            
            if(table_map (i,1)== 0) then
               j = GDL_center_pix (a_clay(1:nsoil_pcarbon(1)),a_sand(1:nsoil_pcarbon(1)), &
                    min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)         
               min_percs%clay_perc = a_clay(j)
               min_percs%silt_perc = a_silt(j)
               min_percs%sand_perc = a_sand(j)		   
               table_map (i,1)= table_map (DeLannoy_class (min_percs),1)
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
               table_map (i,2)= table_map (DeLannoy_class (min_percs),2)	         
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
               table_map (i,3)= table_map (DeLannoy_class (min_percs),3)	         	         
            endif
         end do
      end do

    END subroutine get_GDL_soil_table

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
