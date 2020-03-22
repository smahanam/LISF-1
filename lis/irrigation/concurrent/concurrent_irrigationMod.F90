!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module concurrent_irrigationMod
!BOP
!
! !MODULE: concurrent_irrigationMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  11 Nov 2012: Sujay Kumar; Initial implementation
!  18 Jun 2014: Ben Zaitchik; Modified for flood
!
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod

  implicit none

  PRIVATE
  
  PUBLIC  :: concurrent_irrigation_init
  PUBLIC  :: concurrent_irrigation_updates  
  integer, parameter, PUBLIC :: N_IRRIG_TYPES = 3 

contains
  
  subroutine concurrent_irrigation_init(irrigState)


    type(ESMF_State) :: irrigState(LIS_rc%nnest)

    integer              :: n 
    integer              :: rc, status
    type(ESMF_ArraySpec) :: arrspec1, arrspec2
    type(ESMF_Field)     :: irrigRateField, irrigFracField
    type(ESMF_Field)     :: irrigRootDepthField, irrigScaleField
    real,  pointer       :: irrigrate(:)
    real,  pointer       :: irrigFrac(:,:), frac(:,:)
    real,  pointer       :: irrigRootdepth(:), rootdepth(:)
    real,  pointer       :: irrigScale(:,:),scale(:,:)
    character*100        :: maxrootdepthfile

    do n=1,LIS_rc%nnest
       allocate(irrigFrac(LIS_rc%npatch(n,LIS_rc%lsm_index),N_IRRIG_TYPES))
       allocate(irrigRootDepth(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irrigScale(LIS_rc%npatch(n,LIS_rc%lsm_index),N_IRRIG_TYPES))

       write(LIS_logunit,*) " Running the 'Concurrent' irrigation method ... "

       call ESMF_ConfigGetAttribute(LIS_config,maxrootdepthfile,&
            label="Concurrent irrigation max root depth file:",&
            rc=rc)
       call LIS_verify(rc,&
            'Concurrent irrigation max root depth file: option not specified in the config file')

       call read_irrigfrac(n, irrigFrac)
       call read_irrigrootdepth(n, maxrootdepthfile, irrigRootDepth)

       call compute_irrigScale(n,irrigFrac, irrigScale)

       call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_ArraySpecSet failed in concurrent_irrigation_init")
       call ESMF_ArraySpecSet(arrspec2,rank=2,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_ArraySpecSet arrspec2 failed in concurrent_irrigation_init")

       irrigRateField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation rate", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in concurrent_irrigation_init")
       
       call ESMF_StateAdd(irrigState(n),(/irrigRateField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irrigRate failed in concurrent_irrigation_init")
       
       irrigFracField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec2,&
            name="Irrigation frac", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in concurrent_irrigation_init")

       call ESMF_FieldGet(irrigFracField,localDE=0,&
            farrayPtr=frac,rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for IrrigFrac')
       
       frac = irrigFrac

       call ESMF_StateAdd(irrigState(n),(/irrigFracField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irrigFrac failed in concurrent_irrigation_init")
       deallocate(irrigFrac)

       irrigRootdepthField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation max root depth", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in concurrent_irrigation_init")

       call ESMF_FieldGet(irrigRootdepthField,localDE=0,&
            farrayPtr=rootdepth,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for root depth')
       rootdepth=irrigRootDepth
       
       call ESMF_StateAdd(irrigState(n),(/irrigRootdepthField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for max root depth failed in concurrent_irrigation_init")
       deallocate(irrigRootDepth)

       irrigScaleField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec2,&
            name="Irrigation scale",rc=status)
       call LIS_verify(status, &
            'ESMF_FieldCreate failed in concurrent_irrigation_init')
       
       call ESMF_FieldGet(irrigScaleField,localDE=0,&
            farrayPtr = scale,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for irrigation scale')

       scale = irrigScale
       call ESMF_StateAdd(irrigState(n),(/irrigScaleField/),rc=status)
       call LIS_verify(status,&
            'ESMF_StateAdd for irrigation scale failed in concurrent_irrigation_init')
       deallocate(irrigScale)
    enddo

  end subroutine concurrent_irrigation_init

  ! -----------------------------------------------------------------------------------

  subroutine concurrent_irrigation_updates(n, irrigState)

    use LIS_FORC_AttributesMod 
    use LIS_histDataMod
    use LIS_metforcingMod, only : LIS_FORC_State    

    implicit none

    integer, intent(in) :: n 
    type(ESMF_State)    :: irrigState
    
    ! local trigger time for Drip, Sprinkler, and Flood is 6 AM
    ! Duration is 0.5 for flood and 10 for other 2 methods, thus we use the maximum of the 2 durations
    ! since irrigrate get updated from drip and sprinkler during that time.
    real, parameter     :: otimes = 6.0 ! local trigger check start time [hour]
    real, parameter     :: irrhrf = 10. !duration of sprinkler and drip irrigation [hour]
    
    integer             :: t,gid
    real                :: ltime,otimee
    integer             :: chhr, lhr
    integer             :: status

    type(ESMF_Field)    :: irrigRateField,prcpField
    real,    pointer    :: prcp(:)
    real                :: irrigAmt(LIS_rc%npatch(n,LIS_rc%lsm_index))
    real,    pointer    :: irrigRate(:)

    call ESMF_StateGet(irrigState,&
         "Irrigation rate",&
         irrigRateField,rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed for Irrigation rate')
    
    call ESMF_FieldGet(irrigRateField,localDE=0,&
         farrayPtr=irrigrate,rc=status)
    call LIS_verify(status,'ESMF_FieldGet failed for irrigrate ')

!    call ESMF_StateGet(LIS_FORC_State(n),&
!         trim(LIS_FORC_Rainf%varname(1)),prcpField,&
!         rc=status)
!    call LIS_verify(status,&
!         'ESMF_StateGet failed for rainf in concurrent_irrigation')

    
!    call ESMF_FieldGet(prcpField,localDE=0, farrayPtr=prcp,rc=status)
!    call LIS_verify(status,&
!         'ESMF_FieldGet failed for rainf in concurrent_irrigation')

    irrigAmt = 0.0
    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

       gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
       chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
       if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
            (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
          chhr = chhr -1
       lhr = LIS_rc%hr +chhr
       if(lhr.ge.24) lhr = lhr-24
       if(lhr.lt.0) lhr = lhr+24
       
       ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0
       otimee = otimes + irrhrf
       if((ltime.ge.otimes).and.(ltime.lt.otimee)) then           
          irrigAmt(t) = irrigRate(t)
       endif
       call LIS_diagnoseIrrigationOutputVar(n,t,LIS_MOC_IRRIGATEDWATER,&
            value=irrigAmt(t),unit="kg m-2 s-1",direction="-",vlevel=1)
    enddo
    
  end subroutine concurrent_irrigation_updates


  subroutine read_irrigFrac(n,frac)
    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,      intent(in) :: n 
    real                     :: frac(LIS_rc%npatch(n,LIS_rc%lsm_index), N_IRRIG_TYPES)

    integer                  :: t,col,row
    integer                  :: nid,ios,status,fracId,itid, irrigtypes
    logical                  :: file_exists    
    real, allocatable        :: l_frac(:,:,:)
    real,         pointer    :: glb_frac(:,:,:)
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       allocate (l_frac(LIS_rc%lnc(n),LIS_rc%lnr(n),N_IRRIG_TYPES))
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

       write(LIS_logunit,*) " Reading in the irrigation fraction field ... "
       
       allocate(glb_frac(LIS_rc%gnc(n),LIS_rc%gnr(n),N_IRRIG_TYPES))
       
       ios = nf90_inq_varid(nid,'IRRIGTYPE',fracId)
       call LIS_verify(ios,'nf90_inq_varid failed for IRRIGFRAC')
       ios = nf90_inq_dimid  (nid, 'irrigtypes', itid)
       call LIS_verify(ios,'nf90_inq_varid failed for irrigtypes dimension')
!       ios = nf90_inq_dimlen (nid, itid, irrigtypes)  
!       call LIS_verify(ios,'nf90_inq_varid failed for irrigtypes dimension value')

       if(irrigtypes < N_IRRIG_TYPES) then
          write(LIS_logunit,*) 'irrigtypes < N_IRRIG_TYPES '
          write(LIS_logunit,*) 'program stopping ...'
          call LIS_endrun         
       endif

       ios = nf90_get_var(nid,fracId, glb_frac, start =  (/1,1,1/), count = (/LIS_rc%gnc(n),LIS_rc%gnr(n),N_IRRIG_TYPES/))
       call LIS_verify(ios,'nf90_get_var failed for in concurrent_irrigationMod')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in concurrent_irrigationMod')
       
       l_frac(:,:,:) = glb_frac(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1),:)
       deallocate(glb_frac)

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          
          frac(t,:) = l_frac(col,row,:)
       enddo

    else
       write(LIS_logunit,*) 'irrigation frac map: ',&
            LIS_rc%paramfile(n), ' does not exist'
       write(LIS_logunit,*) 'program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine read_irrigFrac

  subroutine read_irrigRootdepth(n, rdfile, rootdepth)

    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,    intent(in) :: n 
    character(len=*)       :: rdfile
    real                   :: rootdepth(LIS_rc%npatch(n,LIS_rc%lsm_index))
    integer                :: total_vegtypes
!    integer, parameter     :: nt = 32 ! used to be hardcoded 
    real, allocatable      :: rootd(:)
    integer                :: ftn
    integer                :: t,j,col,row
    integer                :: nid, ios, status, croptypeId
    logical                :: file_exists    
    real,   allocatable    :: l_croptype(:,:)
    real,   allocatable    :: glb_croptype(:,:)
    real,   allocatable    :: glb_croptype1(:,:)
! __________________________________________________________________________
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

 !- Read in LDT input crop classification information (done here for now):
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
     ! Read in LDT-generated netcdf file information:
       write(LIS_logunit,*)"[INFO] Reading crop classification information ..."
       ios = nf90_open(path=LIS_rc%paramfile(n),&
                       mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_irrigRootdepth (concurrent)')
 
       ios = nf90_get_att(nid, NF90_GLOBAL, 'CROPCLASS_SCHEME', LIS_rc%cropscheme)
       call LIS_verify(ios,'Error in nf90_get_att in read_irrigRootdepth (concurrent)')
 
       ios = nf90_get_att(nid, NF90_GLOBAL, 'CROPCLASS_NUMBER', LIS_rc%numbercrops)
       call LIS_verify(ios,'Error in nf90_get_att in LIS_irrigation_init')
       write(LIS_logunit,*)"[INFO] Read in crop classfication: ",trim(LIS_rc%cropscheme),&
                          ", with the number of crop types:",LIS_rc%numbercrops
       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in read_irrigRootdepth (concurrent)')
    endif

  ! Estimate total crop types, added to landcover scheme class total:
    select case ( LIS_rc%lcscheme )
     case( "UMD" )
       total_vegtypes = 13 + LIS_rc%numbercrops
     case( "UMD+MIRCA" )
       total_vegtypes = 14 + LIS_rc%numbercrops
!     case( "IGBP", "IGBPNCEP" )
     case( "IGBP", "IGBPNCEP", "IGBP+MIRCA", "IGBPNCEP+MIRCA" )
       total_vegtypes = 20 + LIS_rc%numbercrops
     case( "USGS" )
       total_vegtypes = 24 + LIS_rc%numbercrops
     case default
       write(LIS_logunit,*) "[ERR] The landcover scheme, ",trim(LIS_rc%lcscheme),","
       write(LIS_logunit,*) "[ERR] is not supported for sprinkler irrigation Stopping program ... "
       call LIS_endrun()
    end select
  ! Assign default 32 UMD+CROPMAP for now, due to indexing for max root depth input files:
!    total_vegtypes = 13 + LIS_rc%numbercrops ! KRA: Commenting out to account for crops

    allocate(l_croptype(LIS_rc%lnc(n),LIS_rc%lnr(n)))

 !- Read the max root depth table file:
    inquire(file=rdfile,exist=file_exists)
    if(file_exists) then 
       write(LIS_logunit,*) "[INFO] Reading in the max root depth file: ",trim(rdfile)
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=rdfile,status='old')
       allocate( rootd(total_vegtypes) )
       read(ftn,*) (rootd(j),j=1,total_vegtypes)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) "[ERR] Max root depth file, ",trim(rdfile),", not found."
       write(LIS_logunit,*) "[ERR] Stopping program ..."
       call LIS_endrun()
    endif

 !- Read in crop type map file (specified in LIS parameter input file)
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       ios = nf90_open(path=LIS_rc%paramfile(n),&
                       mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_irrigRootdepth (sprinkler)')

       write(LIS_logunit,*) "[INFO] Reading in the crop type field ... "
       
       allocate(glb_croptype(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       
       ios = nf90_inq_varid(nid,'CROPTYPE',croptypeId)
       call LIS_verify(ios,'nf90_inq_varid failed for CROPTYPE (sprinkler)')
       
       ios = nf90_get_var(nid, croptypeId, glb_croptype)
       call LIS_verify(ios,'nf90_get_var failed for CROPTYPE  (sprinkler)')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in read_irrigRootdepth (sprinkler)')
       
       l_croptype(:,:) = glb_croptype(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1))
       deallocate(glb_croptype)

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          
          if(l_croptype(col,row).gt.0) then 
             rootdepth(t) = rootd(nint(l_croptype(col,row)))
          else
             rootdepth(t) = 0 
          endif
       enddo
       deallocate( rootd )

    else
       write(LIS_logunit,*) "[ERR] The irrigation croptype map: ",&
             LIS_rc%paramfile(n)," does not exist."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun
    endif
    deallocate(l_croptype)
#endif

  end subroutine read_irrigRootdepth

  subroutine compute_irrigScale(n,irrigFrac, irrigScale)

    integer                :: n, i 
    real                   :: irrigFrac (LIS_rc%npatch(n,LIS_rc%lsm_index),N_IRRIG_TYPES)
    real                   :: irrigScale(LIS_rc%npatch(n,LIS_rc%lsm_index),N_IRRIG_TYPES)
    
    integer                :: t,gid,vegt
    real                   :: crppix(LIS_rc%ngrid(n))
    real                   :: grasspix(LIS_rc%ngrid(n))
    real                   :: restpix(LIS_rc%ngrid(n))
    real                   :: irrpix,excess
    integer                :: crop1,crop2,grass,shrub1,shrub2
    crppix    = 0.0
    grasspix  = 0.0
    restpix   = 0.0

!------------------------------------------------------------------------
! WARNING: The following code is valid only for the no-tiling or the 
! vegetation only tiling. The fgrd values are not valid when multiple
! modes of tiling are turned on. 
!------------------------------------------------------------------------
    if(LIS_rc%lcscheme.eq."UMD") then !UMD
      crop1 = 11
      crop2 = 11
      grass = 10 
      shrub1 = 6
      shrub2 = 9
   elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP" .or. &
       LIS_rc%lcscheme.eq."IGBP+MIRCA" .or.  LIS_rc%lcscheme.eq."IGBPNCEP+MIRCA" ) then 
      crop1 = 12
      crop2 = 14
      grass = 10 
      shrub1 = 6
      shrub2 = 9
   elseif(LIS_rc%lcscheme.eq."USGS") then !UMD
      crop1 = 2
      crop2 = 6
      grass = 7 
      shrub1 = 8
      shrub2 = 10
!   elseif(LIS_rc%lcscheme.eq."UMD+MIRCAIrrig")  ! UMD + MIRCA crop types
!      crop1 = 15
!      crop2 = 16
!      grass = 10
!      shrub1 = 6
!      shrub2 = 9
   else
      write(LIS_logunit,*) "The landcover scheme, ",trim(LIS_rc%lcscheme)
      write(LIS_logunit,*) "is not supported for 'Concurrent' irrigation at this time. "
      write(LIS_logunit,*) " Stopping ..."
      call LIS_endrun()
   endif
   
   do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
      if(LIS_domain(n)%tile(t)%vegt.ne.LIS_rc%waterclass.and.&
           LIS_domain(n)%tile(t)%vegt.ne.LIS_rc%urbanclass) then 
         if(LIS_domain(n)%tile(t)%vegt.ge.crop1.and.&
              LIS_domain(n)%tile(t)%vegt.le.crop2) then !crop tiles
            crppix(gid) = crppix(gid) + LIS_domain(n)%tile(t)%fgrd* & 
                 LIS_domain(n)%tile(t)%pens
         elseif(LIS_domain(n)%tile(t)%vegt.eq.grass) then !grassland
            grasspix(gid) = LIS_domain(n)%tile(t)%fgrd*& 
                 LIS_domain(n)%tile(t)%pens
         elseif(LIS_domain(n)%tile(t)%vegt.ge.shrub1.and.&
              LIS_domain(n)%tile(t)%vegt.le.shrub2) then !shrubs
            restpix(gid) = restpix(gid) + LIS_domain(n)%tile(t)%fgrd* & 
                 LIS_domain(n)%tile(t)%pens
         endif
! logic for the more detailed (32 category map)
!             if(LIS_domain(n)%tile(t)%vegt.ge.14) then !crop tiles
!                crppix(gid) = crppix(gid) + LIS_domain(n)%tile(t)%fgrd
!             elseif(LIS_domain(n)%tile(t)%vegt.eq.10) then !grassland
!                grasspix(gid) = LIS_domain(n)%tile(t)%fgrd
!             elseif(LIS_domain(n)%tile(t)%vegt.gt.5.and.&
!                  LIS_domain(n)%tile(t)%vegt.lt.10) then !shrubs
!                restpix(gid) = restpix(gid) + LIS_domain(n)%tile(t)%fgrd
!             endif
      endif
   enddo
   
   irrigScale =1.0
   do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
      vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt

      IRR_TYPE : do i = 1, N_IRRIG_TYPES
         irrpix = irrigFrac(t,i)*0.01
         if(vegt.ne.LIS_rc%waterclass) then 
            if(irrpix < crppix(gid)) then 
               if(vegt.ge.crop1.and.vegt.le.crop2) then 
                  irrigScale(t,i) = irrpix/crppix(gid)
               else
                  irrigScale(t,i) = 0.0
               endif
            else
               excess = irrpix -crppix(gid)
               if(excess.gt.grasspix(gid)) then 
                  if(vegt.ge.shrub1.and.vegt.le.shrub2) then 
                     irrigScale(t,i) = (excess - grasspix(gid)) /restpix(gid)
                  elseif(vegt.eq.grass) then !grass
                     irrigScale(t,i) = 1.0
                  elseif(vegt.ge.crop1.and.vegt.le.crop2) then !crop
                     irrigScale(t,i) = 1.0
                  else
                     irrigScale(t,i) = 0.0
                  endif
               elseif(excess.lt.grasspix(gid)) then 
                  if(vegt.eq.grass) then 
                     irrigScale(t,i) = excess/grasspix(gid)
                  elseif(vegt.ge.crop1.and.vegt.le.crop2) then 
                     irrigScale(t,i) = 1.0
                  else
                     irrigScale(t,i) = 0.0
                  endif
               else
                  if(vegt.eq.grass) then 
                     irrigScale(t,i) =1.0
                  elseif(vegt.ge.crop1.and.vegt.le.crop2) then 
                     irrigScale(t,i) =1.0
                  else
                     irrigScale(t,i) =0.0
                  endif
               endif
            endif
         endif
      end do IRR_TYPE
   enddo
#if 0        
       irrigScale =1.0
       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
          vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
          
          IRR_TYPE2 : do i = 1, N_IRRIG_TYPES          
             irrpix = irrigFrac(t,i)*0.01
             if(vegt.ne.LIS_rc%waterclass) then 
                if(irrpix < crppix(gid)) then 
                   if(vegt.ge.14) then 
                      irrigScale(t,i) = irrpix/crppix(gid)
                   else
                      irrigScale(t,i) = 0.0
                   endif
                else
                   excess = irrpix -crppix(gid)
                   if(excess.gt.grasspix(gid)) then 
                      if(vegt.gt.5.and.vegt.lt.10) then 
                         irrigScale(t,i) = (excess - grasspix(gid)) /restpix(gid)
                      elseif(vegt.eq.10) then !grass
                         irrigScale(t,i) = 1.0
                      elseif(vegt.ge.14) then !crop
                         irrigScale(t,i) = 1.0
                      else
                         irrigScale(t,i) = 0.0
                      endif
                   elseif(excess.lt.grasspix(gid)) then 
                      if(vegt.eq.10) then 
                         irrigScale(t,i) = excess/grasspix(gid)
                      elseif(vegt.ge.14) then 
                         irrigScale(t,i) = 1.0
                      else
                         irrigScale(t,i) = 0.0
                      endif
                   else
                      if(vegt.eq.10) then 
                         irrigScale(t,i) =1.0
                      elseif(vegt.ge.14) then 
                         irrigScale(t,i) =1.0
                      else
                         irrigScale(t,i) =0.0
                      endif
                   endif
                endif
             endif
          end do IRR_TYPE2
       enddo
#endif

  end subroutine compute_irrigScale

end module concurrent_irrigationMod
