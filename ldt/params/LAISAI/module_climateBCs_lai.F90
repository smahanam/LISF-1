
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: module_climateBCs_lai
! \label{module_climateBCs_lai}
!
! !REVISION HISTORY:
!  20 Apr 2020: Sarith Mahanama Adapted from read_CLSMF25_lai.F90, 

! !INTERFACE:

module module_climateBCs_lai

  use LDT_laisaiMod,only :  LDT_laisai_struc
  use LDT_coreMod, only :  LDT_rc
  use LDT_ClimateBCsReader, ONLY : ClimateBCsReader

  implicit none
  
  private

  public read_ClimateBCs_lai, read_ClimateBCs_laiMax, read_ClimateBCs_laiMin
  real, allocatable, dimension (:,:,:), save :: clim_data
  
  contains

    subroutine read_ClimateBCs_lai(nest, array, maskarray)
  
      implicit none
      ! !ARGUMENTS: 
      integer, intent(in)    :: nest    
      real,    intent(inout) :: array (:,:,:)   
      type (ClimateBCsReader)       :: bcr
      real, optional, intent(inout) :: maskarray(LDT_rc%lnc(nest),LDT_rc%lnr(nest))

      ! !DESCRIPTION:
      !  This subroutine retrieves the greenness fraction climatology for the 
      !  specified month and returns the values in the latlon projection
      !  
      !  The arguments are:
      !  \begin{description}
      !  \item[nest]
      !   index of the nest
      !  \item[array]
      !   output field with the retrieved greenness fraction
      !  \item[maskarray]
      !   optional input field of reading in the mask array
      !  \end{description}
      !
      !EOP      
      ! __________________________________________________________
      
      array = LDT_rc%udef
      if (allocated (clim_data)) deallocate (clim_data)
      call bcr%readDataset (nest,trim(LDT_laisai_struc(nest)%lai%source), &
           LDT_laisai_struc(nest)%laisaiInterval, LDT_laisai_struc(nest)%laisai_proj, array)

      allocate (clim_data (size (array,1),size (array,2), size (array,3)))
      clim_data = array

    end subroutine read_ClimateBCs_lai

    ! --------------------------------------------------------------

    subroutine read_ClimateBCs_laiMax (nest, array, maskarray)

      implicit none
      ! !ARGUMENTS: 
      integer, intent(in)    :: nest    
      real,    intent(inout) :: array (:,:)
      real, optional, intent(inout) :: maskarray(LDT_rc%lnc(nest),LDT_rc%lnr(nest))
      integer :: i,j
      
      array = LDT_rc%udef
      do j = 1, LDT_rc%lnr(nest)
         do i = 1, LDT_rc%lnc(nest)
            if (maxval (clim_data (i,j,:)) .ne. LDT_rc%udef) &
                 array (i,j) = maxval (clim_data (i,j,:))
         end do
      end do

    end subroutine read_ClimateBCs_laiMax

    ! --------------------------------------------------------------

    subroutine read_ClimateBCs_laiMin (nest, array, maskarray)

      implicit none
      ! !ARGUMENTS: 
      integer, intent(in)    :: nest    
      real,    intent(inout) :: array (:,:)
      real, optional, intent(inout) :: maskarray(LDT_rc%lnc(nest),LDT_rc%lnr(nest))
      integer :: i,j
      
      array = LDT_rc%udef
      do j = 1, LDT_rc%lnr(nest)
         do i = 1, LDT_rc%lnc(nest)
            if (minval (clim_data (i,j,:)) .ne. LDT_rc%udef) &
                 array (i,j) = minval (clim_data (i,j,:))
         end do
      end do

    end subroutine read_ClimateBCs_laiMin
    
  end module module_climateBCs_lai
