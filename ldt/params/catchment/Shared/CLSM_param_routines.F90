#include "LDT_misc.h"
module CLSM_param_routines

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_gfracMod
  use LDT_albedoMod
  use LDT_paramDataMod

  use CLSM_util, ONLY: NC_VarID, &
       c_data => G5_BCSDIR,      &
       LDT_g5map, write_clsm_files, &
       NC => nc_g5_rst,             &
       NR => nr_g5_rst

  use get_DeLannoy_SoilClass, ONLY:          & 
       mineral_perc, GDL_center_pix,         &
       n_SoilClasses => n_DeLannoy_classes,  &
       soil_class    => DeLannoy_class,      &
       GDL_TABLE
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber,LDT_verify 
  
  implicit none

  private

  public create_CLSM_parameters, catchmentParms_writeHeader, catchmentParms_writeData,  clsm_type_dec, &
       set_CLSM_param_attribs, CLSMF25_struc, CLSMJ32_struc, sibalb, jpl_canoph

  type :: clsm_type_dec

     character*100 :: albnirfile
     character*100 :: albvisfile
     real          :: catchparms_gridDesc(20)
     character*50  :: catchparms_proj
     character*50  :: catchparms_gridtransform
!  - Catchment (F2.5):
!     character*100 :: modisdir        ! sub-directory for MODIS files
!     character*140 :: tile_coord_file ! tile coordinate file
!     character*140 :: tile_veg_file   ! tile vegetation file
     character*140 :: soilparamfile   ! soil parameters file
     character*140 :: sltsfile        ! surface layer timescales file
     character*140 :: topo_ar_file    ! topography parameters file
     character*140 :: topo_bf_file    ! topography parameters file
     character*140 :: topo_ts_file    ! topography parameters file
     character*140 :: catchgreenfile  ! greenness climatology file
     character*140 :: catchlaifile    ! LAI climatology file
     real          :: dzsfcrd         ! CLSM top soil layer depth from ldt.config
     real          :: addbdrckcrd     ! CLSM add to bedrock depth from ldt.config
     type(LDT_paramEntry) :: porosity    ! soil porosity
     type(LDT_paramEntry) :: psisat      ! saturated soil moisture potential
     type(LDT_paramEntry) :: bexp        ! Clapp-Hornberger parameter
     type(LDT_paramEntry) :: wpwet       ! wilting point wetness
     type(LDT_paramEntry) :: bdrckdpth   ! depth to bedrock 
     type(LDT_paramEntry) :: ksat        ! saturated hydraulic conductivity (SATDK; m s-1)
     type(LDT_paramEntry) :: gnu         ! vertical decay factor for transmissivity
     type(LDT_paramEntry) :: ars1        ! Wetness parameters
     type(LDT_paramEntry) :: ars2         
     type(LDT_paramEntry) :: ars3
     type(LDT_paramEntry) :: ara1        ! Shape parameters
     type(LDT_paramEntry) :: ara2
     type(LDT_paramEntry) :: ara3
     type(LDT_paramEntry) :: ara4
     type(LDT_paramEntry) :: arw1        ! Minimum Theta parameters:
     type(LDT_paramEntry) :: arw2
     type(LDT_paramEntry) :: arw3
     type(LDT_paramEntry) :: arw4
     type(LDT_paramEntry) :: bf1         ! Baseflow topographic params
     type(LDT_paramEntry) :: bf2
     type(LDT_paramEntry) :: bf3
     type(LDT_paramEntry) :: tsa1        ! Water transfer parameters
     type(LDT_paramEntry) :: tsa2
     type(LDT_paramEntry) :: tsb1        
     type(LDT_paramEntry) :: tsb2
     type(LDT_paramEntry) :: atau
     type(LDT_paramEntry) :: btau
     type(LDT_paramEntry) :: z2
     type(LDT_paramEntry) :: albnirdir   ! Albedo NIR direct scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albnirdif   ! Albedo NIR diffuse scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albvisdir   ! Albedo VIS direct scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albvisdif   ! Albedo VIS diffuse scale factor (CLSM F2.5)

  end type clsm_type_dec
  type(clsm_type_dec), save, allocatable :: CLSMF25_struc(:), CLSMJ32_struc(:)
  real, parameter :: gnu = 1.0, zks = 2.0
  integer, PARAMETER :: nbdep=150, NAR=1000,nwt=81,nrz=41
  logical :: preserve_soiltype = .false.
  logical, parameter ::  bug =.false.
  logical, parameter :: error_file=.false. 
  real, parameter :: slice=0.1, lim =5.,grzdep =1.1

contains

  !--------------------------------------------------------------------

  SUBROUTINE create_CLSM_parameters (nbcatch, BEE, POROS, WPWET, PSIS, COND,  &
       SOILDEPTH, SOIL_CLASS_TOP, SOIL_CLASS_COM,                             &
       AGNU, ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4,                             &
       ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3, tsa1, tsa2,tsb1, tsb2)

    implicit none
    integer, intent (in) :: nbcatch
    real, dimension (:), intent (inout) :: BEE, POROS, WPWET, PSIS, COND, SOILDEPTH, &
         AGNU, ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4, &
         ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3, tsa1, tsa2,tsb1, tsb2 
    integer, dimension (:), intent (inout) :: SOIL_CLASS_TOP, SOIL_CLASS_COM  

    REAL, allocatable, dimension(:) :: TOPMEAN, TOPVAR, TOPSKEW
    real, allocatable, dimension(:) :: a_sand,a_clay,a_silt,a_oc,  &
         tile_lon, tile_lat
    
    real, allocatable, dimension (:,:)    :: good_clay, good_sand
    integer, allocatable, dimension (:,:) :: tile_add, tile_pick
    type (mineral_perc) :: min_percs
    integer :: CF1, CF2, CF3, CF4
    integer i,j,n,k,ftbl, far, fbf, fts, fso, ftcl
    integer soil_gswp
    real rzdep, atile_sand,atile_clay 
     
    REAL ST(NAR), AC(NAR),COESKEW
    REAL, allocatable, dimension (:) ::       &
         taberr1,taberr2,normerr1,normerr2,   &
         taberr3,taberr4,normerr3,normerr4

    real watdep(nwt,nrz),wan(nwt,nrz),rzexcn(nwt,nrz),frc(nwt,nrz)
    real, allocatable, dimension  (:,:,:) :: &
         gwatdep,gwan,grzexcn,gfrc
    real :: wtdep,wanom,rzaact,fracl,profdep,dist_save,     &
         ncells_top, ncells_top_pro,ncells_sub_pro,tile_distance
    character*400 :: fname,fout,losfile
    character*6 rdep,ext
    integer :: iwt,irz,group
    logical :: picked
    logical :: preserve_soiltype = .false.
 
! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
!
! NOTE: "!$" is for conditional compilation
!
logical :: running_omp = .false.
!
!$ integer :: omp_get_thread_num, omp_get_num_threads
!
integer :: n_threads=1, li, ui
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
      
!c-------------------------------------------------------------------------

      write(LDT_logunit,'(A60)')'[CLSM create_CLSM_parameters] creating model params .....'

      ! Scale saturated hydraulic conductivity  to surface 
      ! --------------------------------------------------

      COND = COND/exp(-1.0*zks*gnu)
      AGNU = GNU

      allocate (TOPMEAN (1: nbcatch))
      allocate (TOPVAR  (1: nbcatch))
      allocate (TOPSKEW (1: nbcatch))
      call read_cti_stats (TOPMEAN, TOPVAR, TOPSKEW)
      ftbl = LDT_getNextUnitNumber()
      open (ftbl, file=trim(c_data)//trim(GDL_TABLE), form='formatted',status='old', &
           action = 'read')
      read (ftbl,'(a)')fout           
      losfile =trim(c_data)//'/Woesten_SoilParam/loss_pd_top/loss_perday_rz1m_'

      allocate (a_sand (1:n_SoilClasses))
      allocate (a_silt (1:n_SoilClasses))
      allocate (a_clay (1:n_SoilClasses))
      allocate (a_oc   (1:n_SoilClasses))
      allocate (gwatdep(1:nwt,1:nrz,1:n_SoilClasses))
      allocate (gwan   (1:nwt,1:nrz,1:n_SoilClasses))
      allocate (grzexcn(1:nwt,1:nrz,1:n_SoilClasses))
      allocate (gfrc   (1:nwt,1:nrz,1:n_SoilClasses))
      
      do n =1,n_SoilClasses
         read (ftbl,'(4f7.3)')a_sand(n),a_clay(n),a_silt(n),a_oc(n)
         write (fout,'(i2.2,i2.2,i4.4)')nint(a_sand(n)),nint(a_clay(n)),nint(100*a_oc(n))
         ftcl = LDT_getNextUnitNumber()
         open (ftcl,file=trim(losfile)//trim(fout),  &
              form='formatted',status='old')
         
         do iwt=1,nwt
            do irz=1,nrz
               read(ftcl,2000) wtdep,wanom,rzaact,fracl
2000           format(1x,4e16.8)
               gwatdep(iwt,irz,n)= wtdep
               gwan(iwt,irz,n)   = wanom
               grzexcn(iwt,irz,n)= rzaact
               gfrc(iwt,irz,n)   = amin1(fracl,1.)
            enddo
         enddo
         close (ftcl,status='keep')	   
         call LDT_releaseUnitNumber(ftcl)
      end do
      close (ftbl,status='keep')  
      call LDT_releaseUnitNumber(ftbl)

      allocate (tile_lon  (1:nbcatch))
      allocate (tile_lat  (1:nbcatch))
      allocate (TABERR1   (1:nbcatch))
      allocate (TABERR2   (1:nbcatch))
      allocate (TABERR3   (1:nbcatch))
      allocate (TABERR4   (1:nbcatch))
      allocate (NORMERR1  (1:nbcatch))
      allocate (NORMERR2  (1:nbcatch))
      allocate (NORMERR3  (1:nbcatch))
      allocate (NORMERR4  (1:nbcatch))
      allocate (good_clay (1:100,4))
      allocate (good_sand (1:100,4))
      allocate (tile_add  (1:100,4))   
      allocate (tile_pick (1:100,4)) 
      tile_add = 0
      tile_pick= 0
      good_clay =0.
      good_sand =0.

      tile_lon = LDT_g5map%lon
      tile_lat = LDT_g5map%lat

      allocate(low_ind(n_threads))
      allocate(upp_ind(n_threads))
      low_ind(1)         = 1
      upp_ind(n_threads) = nbcatch
      
      if (running_omp)  then
         do i=1,n_threads-1
            
            upp_ind(i)   = low_ind(i) + (nbcatch/n_threads) - 1 
            low_ind(i+1) = upp_ind(i) + 1
            
         end do
      end if

!$OMP PARALLELDO DEFAULT(NONE)                          &
!$OMP SHARED( BEE, PSIS,POROS,COND,WPWET,soildepth,     &
!$OMP        TOPMEAN, TOPVAR, TOPSKEW,                  &
!$OMP        ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4,       &
!$OMP        ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,         &
!$OMP        tsa1, tsa2,tsb1, tsb2,                     &
!$OMP        taberr1,taberr2,normerr1,normerr2,         &
!$OMP        taberr3,taberr4,normerr3,normerr4,         &
!$OMP        gwatdep,gwan,grzexcn,gfrc,soil_class_com,  &
!$OMP        n_threads, low_ind, upp_ind )              &
!$OMP PRIVATE(k,li,ui,n,i,watdep,wan,rzexcn,frc,ST,AC,  &
!$OMP COESKEW,profdep)

     do k=1,n_threads

        li = low_ind(k)
        ui = upp_ind(k)

      do n=li,ui

      CALL TGEN (                              &
          TOPMEAN(n),TOPVAR(n),TOPSKEW(n),     &
          ST,AC,COESKEW)
      
!c Areal fractioning parameters

      CALL SAT_PARAM(                                              &
                     BEE(n),PSIS(n),POROS(n),COND(n),              & 
                     WPWET(n), ST, AC, COESKEW,n,                  &
                     soildepth(n),                                 &
                     ars1(n),ars2(n),ars3(n),                      &
                     ara1(n),ara2(n),ara3(n),ara4(n),              &
                     arw1(n),arw2(n),arw3(n),arw4(n),              &
                     taberr1(n),taberr2(n),taberr3(n),taberr4(n),  &
                     normerr1(n),normerr2(n),normerr3(n),normerr4(n))


      CALL BASE_PARAM(                                  &
          BEE(n),PSIS(n),POROS(n),COND(n),              &
          ST, AC,                                       &
          bf1(n),bf2(n),bf3(n),                         &
          taberr1(n),taberr2(n),normerr1(n),normerr2(n) &
          )
            

	  watdep (:,:) =  gwatdep (:,:,soil_class_com(n))   
	  wan    (:,:) =  gwan    (:,:,soil_class_com(n))
	  rzexcn (:,:) =  grzexcn (:,:,soil_class_com(n))
	  frc    (:,:) =  gfrc    (:,:,soil_class_com(n))

      CALL TS_PARAM(                       &
          BEE(n),PSIS(n),POROS(n),         &
          ST, AC,                          &
          watdep,wan,rzexcn,frc,           &
          tsa1(n),tsa2(n),tsb1(n),tsb2(n)  &
          )

   END DO
   END DO
          !$OMP ENDPARALLELDO

     CF1 =0
     CF2 =0
     CF3 =0
     CF4 =0

     if (write_clsm_files) then
        far = LDT_getNextUnitNumber()
        fts = LDT_getNextUnitNumber()
        fbf = LDT_getNextUnitNumber()
        fso = LDT_getNextUnitNumber()

        open (far, file = 'LDT_clsm/ar.new'        , form = 'formatted', action = 'write')
        open (fts, file = 'LDT_clsm/ts.dat'        , form = 'formatted', action = 'write')
        open (fbf, file = 'LDT_clsm/bf.dat'        , form = 'formatted', action = 'write')
        open (fso, file = 'LDT_clsm/soil_param.dat', form = 'formatted', action = 'write')
     endif

     DO n=1,nbcatch
        atile_clay = a_clay(SOIL_CLASS_COM(n))
        atile_sand = a_sand(SOIL_CLASS_COM(n))
     	if((ars1(n).ne.9999.).and.(arw1(n).ne.9999.))then

         if ((soil_class_com(n)>=1).and.(soil_class_com(n)<=84)) then	
            group=1
         else if ((soil_class_com(n) > 84).and.(soil_class_com(n)<=168)) then
            group=2
         else if ((soil_class_com(n) >168).and.(soil_class_com(n)< N_SoilClasses)) then
            group=3
         else
            group=4
         endif
         
         min_percs%clay_perc = atile_clay
         min_percs%sand_perc = atile_sand
         min_percs%silt_perc = 100. - min_percs%clay_perc - min_percs%sand_perc
         
         if(tile_pick(soil_class (min_percs),group) == 0) then
            tile_pick(soil_class (min_percs),group) = n 
            
            select case (group)
               
            case (1)
               
               CF1 = CF1 + 1
               good_clay (CF1,group) = atile_clay
               good_sand (CF1,group) = atile_sand
               tile_add  (CF1,group) = n
               
            case (2)
               CF2 = CF2 + 1
               good_clay (CF2,group) = atile_clay
               good_sand (CF2,group) = atile_sand
               tile_add  (CF2,group) = n
               
            case (3)
               CF3 = CF3 + 1
               good_clay (CF3,group) = atile_clay
               good_sand (CF3,group) = atile_sand
               tile_add  (CF3,group) = n
               
            case (4)
               CF4 = CF4 + 1
               good_clay (CF4,group) = atile_clay
               good_sand (CF4,group) = atile_sand
               tile_add  (CF4,group) = n
               
            end select
         endif
      endif
   END DO
   
   DO n=1,nbcatch
      atile_clay = a_clay(SOIL_CLASS_COM(n))
      atile_sand = a_sand(SOIL_CLASS_COM(n))
      if((ars1(n).ne.9999.).and.(arw1(n).ne.9999.))then 
         k = n
      else
         if(preserve_soiltype) then 
            if ((soil_class_com(n)>=1).and.(soil_class_com(n)<=84)) then	
               group=1
            else if ((soil_class_com(n)>  84).and.(soil_class_com(n)<=168)) then
               group=2
            else if ((soil_class_com(n)> 168).and.(soil_class_com(n)< N_SoilClasses)) then
               group=3
            else
               group=4
            endif
            
            min_percs%clay_perc = atile_clay
            min_percs%sand_perc = atile_sand
            min_percs%silt_perc = 100. - min_percs%clay_perc - min_percs%sand_perc
            
            if(tile_pick(soil_class (min_percs),group) > 0) then
               k = tile_pick(soil_class (min_percs),group)               
            else
               select case (group)
                  
               case (1)
                  j = GDL_center_pix (good_clay(1:CF1,group),good_sand(1:CF1,group),       &
                       min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)
                  k = tile_add  (j,group) 
               case (2)
                  j = GDL_center_pix (good_clay(1:CF2,group),good_sand(1:CF2,group),       &
                       min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)
                  k = tile_add  (j,group)   
               case (3)
                  j = GDL_center_pix (good_clay(1:CF3,group),good_sand(1:CF3,group),       &
                       min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)
                  k = tile_add  (j,group) 
               case (4)
                  j = GDL_center_pix (good_clay(1:CF4,group),good_sand(1:CF4,group),       &
                       min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)
                  k = tile_add  (j,group)   
               end select
               print *,'NO Similar SoilClass :',soil_class (min_percs),group,n,k            
            endif

         else 
            dist_save = 1000000.
            k = 0
            do i = 1,nbcatch
               if(i /= n) then
                  if((ars1(i).ne.9999.).and.(arw1(i).ne.9999.)) then
                  
                     tile_distance = (tile_lon(i) - tile_lon(n)) * (tile_lon(i) - tile_lon(n)) + &
                          (tile_lat(i) - tile_lat(n)) * (tile_lat(i) - tile_lat(n))
                     if(tile_distance < dist_save) then
                        k = i
                        dist_save = tile_distance
                     endif
                  endif
               endif
            enddo
         endif

         BEE(n)      =  BEE(k)       
         PSIS(n)     =  PSIS(k)      
         POROS(n)    =  POROS(k)     
         COND(n)     =  COND(k)      
         WPWET(n)    =  WPWET(k)     
         soildepth(n)=  soildepth(k) 
         ars1(n)     =  ars1(k)      
         ars2(n)     =  ars2(k)      
         ars3(n)     =  ars3(k)                  
         ara1(n)     =  ara1(k)      
         ara2(n)     =  ara2(k)      
         ara3(n)     =  ara3(k)      
         ara4(n)     =  ara4(k)      
         arw1(n)     =  arw1(k)      
         arw2(n)     =  arw2(k)      
         arw3(n)     =  arw3(k)      
         arw4(n)     =  arw4(k)      
         bf1(n)      =  bf1(k)       
         bf2(n)      =  bf2(k)       
         bf3(n)      =  bf3(k)       
         tsa1(n)     =  tsa1(k)      
         tsa2(n)     =  tsa2(k)      
         tsb1(n)     =  tsb1(k)      
         tsb2(n)     =  tsb2(k)      
         
      endif

      if (write_clsm_files) then
         write(far,'(i8,i8,f5.2,11(2x,e14.7))')           &
              n,LDT_g5map%catid_index(n), gnu,           &
              ars1(n),ars2(n),ars3(n),                   &
              ara1(n),ara2(n),ara3(n),ara4(n),           &
              arw1(n),arw2(n),arw3(n),arw4(n) 
         write(fts,'(i8,i8,f5.2,4(2x,e13.7))')            &
              n,LDT_g5map%catid_index(n), gnu,           &
              tsa1(n),tsa2(n),tsb1(n),tsb2(n)
         write(fbf,'(i8,i8,f5.2,3(2x,e13.7))')            &
              n,LDT_g5map%catid_index(n), gnu,           &
              bf1(n),bf2(n),bf3(n)
         write(fso,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,2f7.3)')        &
              n,LDT_g5map%catid_index(n),                              &
              soil_class_top(k),soil_class_com(k),                     &
              BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n),  &
              atile_sand, atile_clay
      endif
   END DO

   if (write_clsm_files) then
      close (far, status = 'keep')
      close (fts, status = 'keep')
      close (fbf, status = 'keep')
      close (fso, status = 'keep')
      call LDT_releaseUnitNumber(far)
      call LDT_releaseUnitNumber(fts)
      call LDT_releaseUnitNumber(fbf)
      call LDT_releaseUnitNumber(fso)
   endif
      
END SUBROUTINE create_CLSM_parameters

!----------------------------------------------------------------------

  SUBROUTINE read_cti_stats (TOPMEAN, TOPVAR, TOPSKEW)

    IMPLICIT NONE
    real, dimension (:), intent (inout) :: TOPMEAN, TOPVAR, TOPSKEW
    INTEGER, PARAMETER :: nofvar= 3
    INTEGER :: n, i, ncat,idum2
    REAL, allocatable, dimension (:,:) :: var
    real :: dum1,dum2,dum3,dum4,dum5,dum6
    INTEGER*8 :: idum8

    open (1010,file=trim(c_data)//'/SRTM-TopoData/SRTM_cti_stats.dat',       &
         form='formatted', status='old',action='read')

    read (1010,*)ncat
    allocate(var(1:ncat,1:nofvar))

    do i = 1, ncat
       read(1010,'(i8,i15,5(1x,f8.4),i5,e18.3)')n,idum8,dum1,dum2,dum3,dum4,dum5,idum2,dum6
          var(i,1)=dum1
          var(i,2)=dum2
          var(i,3)=dum5
    end do

    close (1010,status='keep')

    do i = 1,LDT_g5map%NT_GEOS 
       TOPMEAN(i) = var(LDT_g5map%catid_index(i),1)
       TOPVAR (i) = var(LDT_g5map%catid_index(i),2) * var(LDT_g5map%catid_index(i),2)
       TOPSKEW(i) = var(LDT_g5map%catid_index(i),3) * var(LDT_g5map%catid_index(i),2) * TOPVAR (i)
    end do

  END SUBROUTINE read_cti_stats
 
! ----------------------------------------------------------------

!BOP
! !ROUTINE:  set_CLSM_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_CLSM_param_attribs(paramEntry, short_name, source, vlevels, &
                 units, full_name )

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name, source
   integer, optional                  :: vlevels
   character(len=*),     optional     :: units 
   character(len=*),     optional     :: full_name

   integer   :: v_temp
   character(20) :: unit_temp
   character(100):: name_temp

! ____________________________________________________
    
   if(present(vlevels)) then 
      v_temp = vlevels
   else
      v_temp = 1
   endif

   if(present(units)) then
      unit_temp = units
   else
      unit_temp = "none"
   endif

   if(present(full_name)) then
      name_temp = full_name
   else
      name_temp = trim(short_name)
   endif

   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = v_temp
   paramEntry%selectOpt = 1
   paramEntry%source = trim(source)
   paramEntry%units = trim(unit_temp)
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

 end subroutine set_CLSM_param_attribs

! -----------------------------------------------------------------

  subroutine catchmentParms_writeHeader (n,ftn,dimID,monthID)

    integer     :: n
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: monthID

    if (LDT_rc%lsm.eq."CLSMF2.5") call writeHeader (n,CLSMF25_struc(:),ftn,dimID,monthID)
    if (LDT_rc%lsm.eq."CLSMJ3.2") call writeHeader (n,CLSMJ32_struc(:),ftn,dimID,monthID)

    contains 

      subroutine writeHeader (n, CLSM_struc, ftn,dimID,monthID)

        type (clsm_type_dec), intent (in) ::  CLSM_struc(:)
        integer     :: n
        integer     :: ftn
        integer     :: dimID(3)
        integer     :: monthID
        
        integer     :: t_dimID(3)
        integer     :: tdimID(3)
        
        tdimID(1) = dimID(1)
        tdimID(2) = dimID(2)
        
        t_dimID(1) = dimID(1)
        t_dimID(2) = dimID(2)

        if(LDT_gfrac_struc(n)%gfrac%selectOpt.gt.0) then
           if(have_str (trim(LDT_gfrac_struc(n)%gfracInterval))) t_dimID(3) = monthID
        end if
        
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%gnu)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%ars1)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%ars2)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%ars3)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%ara1)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%ara2)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%ara3)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%ara4)
        
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%arw1)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             clsm_struc(n)%arw2)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             cLSM_struc(n)%arw3)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%arw4)
        
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%bf1)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%bf2)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%bf3)
        
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%tsa1)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%tsa2)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%tsb1)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%tsb2)
        
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%atau)
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Clsm_struc(N)%btau)
        if (LDT_rc%lsm.eq."CLSMJ3.2") call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
             Clsm_struc(N)%porosity) 
        call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
             Clsm_struc(N)%psisat)
        call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
             Clsm_struc(N)%ksat)
        call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
             Clsm_struc(N)%bexp)
        call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
             Clsm_struc(N)%wpwet)
        !    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
        !         LDT_LSMparam_struc(n)%quartz)
        !    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
        !         LDT_LSMparam_struc(n)%soildepth)
        call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
             Clsm_struc(N)%bdrckdpth)
        
        !- Albedo NIR scale factors:
        if(have_str (trim(LDT_albedo_struc(n)%albInterval))) t_dimID(3) = monthID
        call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
             Clsm_struc(N)%albnirdir)
        call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
             Clsm_struc(N)%albnirdif)
        
        !- Albedo VIS scale factors:
        if(have_str (trim(LDT_albedo_struc(n)%albInterval))) t_dimID(3) = monthID
        call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
             Clsm_struc(N)%albvisdir)
        call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
             Clsm_struc(N)%albvisdif)
        
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ALBEDO_DATA_INTERVAL", &
             LDT_albedo_struc(n)%albInterval))
      end subroutine writeHeader

      ! ---------------------------------------------------------
      
      logical function have_str (str)
        character(*), intent (in) :: str
        have_str = str == ''
        if(have_str) then
           have_str = .FALSE.
        else
           have_str = .TRUE.
        endif
      end function have_str
      
  end subroutine catchmentParms_writeHeader

 ! ----------------------------------------------------------------------------------

  subroutine catchmentParms_writeData (n,ftn)

    integer   :: ftn
    integer   :: n
    
    if (LDT_rc%lsm.eq."CLSMF2.5") call writeData (n,ftn,  CLSMF25_struc(:))
    if (LDT_rc%lsm.eq."CLSMJ3.2") call writeData (n,ftn,  CLSMJ32_struc(:))

    CONTAINS 

      subroutine writeData (n,ftn, CLSM_struc)

        integer :: ftn, n
        type(clsm_type_dec), intent (in) :: CLSM_struc(:)
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%gnu)
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%ars1)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%ars2)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%ars3)
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%ara1)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%ara2)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%ara3)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%ara4)
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%arw1)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%arw2)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%arw3)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%arw4)
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%bf1)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%bf2)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%bf3)
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%tsa1)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%tsa2)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%tsb1)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%tsb2)
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%atau)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%btau)
        
        if (LDT_rc%lsm.eq."CLSMJ3.2") call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%porosity)
        if (LDT_rc%lsm.eq."CLSMJ3.2") call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%z2)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%psisat)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%ksat)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%bexp)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%wpwet)        
        !    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%quartz)
        !    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%soildepth)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%bdrckdpth)
        
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%albnirdir)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%albnirdif)
        
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%albvisdir)
        call LDT_writeNETCDFdata(n,ftn,Clsm_struc(N)%albvisdif)
      END subroutine writeData

  end subroutine catchmentParms_writeData
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
    
      contains

        REAL FUNCTION COEFF (TABLE, NTABL, LAI ,DX, DY)
    
          implicit none
          INTEGER NTABL, LAI       
          REAL TABLE (NTABL, 2), DX, DY
          
          COEFF = (TABLE(LAI,  1)                                      &
               + (TABLE(LAI  ,2) - TABLE(LAI  ,1)) * DY ) * (1.0-DX) &
               + (TABLE(LAI+1,1)                                     &
               + (TABLE(LAI+1,2) - TABLE(LAI+1,1)) * DY ) * DX
          
        END FUNCTION COEFF
        
      END SUBROUTINE SIBALB
      
 ! ------------------------------------------------------------------------------------------

    SUBROUTINE jpl_canoph (ntiles, z2)

      implicit none
      include 'netcdf.inc'
      ! 1) JPL Canopy Height 
      ! /discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001//Simard_Pinto_3DGlobalVeg_JGR.nc4
     
      integer, intent (in)               :: ntiles
      real, pointer, dimension (:), intent (inout) :: z2
      integer  , parameter               :: N_lon_jpl = 43200, N_lat_jpl = 21600
      integer                            :: i,j, status, ncid,varid
      REAL, ALLOCATABLE, dimension (:)   :: count_pix
      INTEGER, ALLOCATABLE, dimension (:,:) :: z2_grid
      INTEGER, ALLOCATABLE, dimension (:,:) :: tile_id
      character*100                      :: fout

      ! READ JPL source data files and regrid
      ! -------------------------------------

      status  = NF_OPEN (trim(c_data)//'/Simard_Pinto_3DGlobalVeg_JGR.nc4', NF_NOWRITE, ncid)
            
      allocate (z2_grid   (1 : NC         , 1 : NR))

      status  = NF_INQ_VARID (ncid,'CanopyHeight',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_INT (ncid,VarID, (/1,1/),(/N_lon_jpl, N_lat_jpl/), z2_grid) ; VERIFY_(STATUS)
      
      status = NF_CLOSE(ncid)

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file

      allocate(tile_id(1:nc,1:nr))
      tile_id = LDT_g5map%rst
      
      allocate (z2        (1:NTILES))
      allocate (count_pix (1:NTILES))
      
      z2        = 0.
      count_pix = 0.

      do j = 1,nr
         do i = 1, nc
            if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then

               if(z2_grid(i,j) >= 0.) then 
                 z2 (tile_id(i,j)) = z2 (tile_id(i,j)) + real (z2_grid(i,j))
                 count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
               endif

            endif
         end do
      end do
      
      where (count_pix >   0.) z2 = z2/count_pix
      where (z2        <  0.5) z2 = 0.5            ! to ensure Z2 >= MIN_VEG_HEIGHT

      deallocate (count_pix)
      deallocate (z2_grid)
      deallocate (tile_id)
      
    END SUBROUTINE jpl_canoph
    
!---------------------------------------------------------------------

      SUBROUTINE TS_PARAM(                                &
                           BEE,PSIS,POROS,                &
                           VALX, PX,                      &
                           watdep,wan,rzexcn,frc,         &
                           tsa1,tsa2,tsb1,tsb2            &
                          )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                         c
!c Given pre-computed 1-D relationships between a "local" root zone excess c 
!c and a "local" catchment deficit, the timescale of the bulk vertical     c
!c transfer between the two bulk prognostic variables is computed using    c
!c the distribution of the local deficit established from the distribution c
!c of the topographic index, then an approximated function of catdef and   c
!c rzex is derived.                                                        c
!c                                                                         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER NAR0
      REAL, intent (in) :: BEE, PSIS, POROS
      REAL, intent (in) :: VALX(NAR), PX(NAR)
      real, intent (inout) :: watdep(nwt,nrz),wan(nwt,nrz),  &
                      rzexcn(nwt,nrz),frc(nwt,nrz)
      real, intent (out) ::  tsa1, tsa2 ,tsb1, tsb2

      integer :: tex,iwt,irz,n,idep,k, index1,i0
      REAL VALX0(NAR), PX0(NAR),sumta,sumta2,timean,zbar, rzw
      REAL :: term1, term2, sumdef, suma, frcsat,rzexc, rzact
      real zdep(nar),def(nar),wrz(nar),wbin(500),rze(nar)
      real catd(2,2),tsc(2,2), satfrc,sumfrac,sumz,frac
      real, parameter ::  frcmax = .041
      real  wtdep,wanom,rzaact,fracl,profdep,rzdep

!      logical bug

!c----------------------------------------------------------------
!c Is loss.dat compatible with rzdep = 0.49 ???

      rzdep = grzdep

!c Convert fractions to "per-hour" values
      do iwt=1,nwt
         do irz=1,nrz
            frc(iwt,irz)=1.-((1.-frc(iwt,irz))**(1./24.))
         enddo
      enddo

         nar0=0
         do n=1,nar
            if (px(n) .ne. 0.) then
               nar0=nar0+1
               valx0(nar0)=valx(n)
               px0(nar0)=px(n)
            endif
         enddo

         sumta=0.
         sumta2=0.
         suma=0.
         do n=1,nar0
            sumta=sumta+px0(n)*valx0(n)
            sumta2=sumta2+px0(n)*valx0(n)*valx0(n)
            suma=suma+px0(n)
         enddo

         timean=sumta/suma

!c**** Loop over two water table depths
         do idep=1,2
            if(idep.eq.1) zbar=1.5 ! zbar in meters
            if(idep.eq.2) zbar=2.0

!c**** Compute array of water table depths:
            do k=1,nar0
               term1=(1/gnu)*(valx0(k)-timean)
               zdep(k)=zbar-term1
               if(zdep(k) .lt. 0.) zdep(k)=0.
            enddo
!c            write(*,*)"  End water table depth"
!c**** Compute array of moisture deficits:
            do k=1,nar0
               term1=(psis-zdep(k))/psis
               term1=term1**(1.-1./bee)
               term2=-psis*(bee/(bee-1.))*(term1-1.)
               def(k)=poros*(zdep(k)-term2)
            enddo

!c**** Add deficits to produce catdef:
            sumdef=0.
            do k=1,nar0
               sumdef=sumdef+def(k)*px0(k)*1000.
            enddo
!c            write(*,*)"  End catchment deficit"
!c**** Compute array of root zone moisture (degree of wetness in root zone):
            do k=1,nar0

               if(zdep(k).eq.0.) then
                  wrz(k)=1.
               elseif(zdep(k)-rzdep.lt.0.) then
                  term1=((psis-zdep(k))/psis)**(1.-1./bee)
                  wrz(k)=(-psis/zdep(k))*(bee/(bee-1.))   &
                      *(term1-1.)
                  frcsat=1.-zdep(k)/rzdep
                  wrz(k)=(1.-frcsat)*wrz(k)+frcsat*1.
               else
                  term1=((psis-zdep(k))/psis)**(1.-1./bee)
                  term2=((psis-zdep(k)+rzdep)/psis)    &
                      **(1.-1./bee)
                  wrz(k)=(-psis/rzdep)*(bee/(bee-1.))  &
                      *(term1-term2)
               endif
            enddo

!c       Loop over two root zone excess values:
            do irz=1,2
               if(irz.eq.1) rzexc=-0.1*poros
               if(irz.eq.2) rzexc=0.1*poros

!c       Determine actual root zone excess
               rzact=0.
               do k=1,nar0
                  rze(k)=rzexc
                  rzw=wrz(k)*poros
                  if(rzw+rze(k) .gt. poros) rze(k)=poros-rzw
                  if(rzw+rze(k) .lt. 0.) rze(k)=rzw
                  rzact=rzact+rze(k)*px0(k)
               enddo
!c            write(*,*)"  End root zone excess"
!c       Compute the average timescale

               satfrc=0.
               do k=1,nar0
                  if(zdep(k).lt.0.) satfrc=satfrc+px0(k)
               enddo

               sumfrac=0.
               sumz=0.
               do k=1,nar0
                  sumz=sumz+zdep(k)*px0(k)
                  if(zdep(k) .lt. 1.) frac=frcmax
                  if(zdep(k) .ge. 1.) then
                     index1=1+int(((zdep(k)*100.)-99)/5.)
                     if(index1.gt.nwt) index1 = nwt
                     frac=amin1(frc(index1,1),frcmax)
                     do i0=2,nrz
                        if(rze(k) .ge. rzexcn(index1,i0))  &
                            frac=amin1(frc(index1,i0),frcmax)
                     enddo
                  endif
                  sumfrac=sumfrac+frac*px0(k)
               enddo
!c            write(*,*)"  End average time scale"
               catd(idep,irz)=sumdef
               tsc(idep,irz)=sumfrac

            enddo
         enddo

         tsb1=(alog(tsc(2,2))-alog(tsc(1,2)))/(catd(2,2)-catd(1,2))
         tsb2=(alog(tsc(2,1))-alog(tsc(1,1)))/(catd(2,1)-catd(1,1))
         tsa1=alog(tsc(2,2))-tsb1*catd(2,2)
         tsa2=alog(tsc(2,1))-tsb2*catd(2,1)

       END SUBROUTINE TS_PARAM

!*********************************************************************

      SUBROUTINE BASE_PARAM(                                  &
                           BEE,PSIS,POROS,COND,               &
                           VALX, PX,                          &
                           bf1,bf2,bf3,                       &
                           taberr1,taberr2,normerr1,normerr2  &
                          )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                   c
!c New way to get baseflow: we parametrize the relationship between  c 
!c catdef and zbar (two parameters bf1 and bf2).                     c
!c Then, in the LSM/catchment.f/base.f, we use the original relation c 
!c from TOPMODEL to infer baseflow from catdef and the mean of the   c
!c topographic index (topmean=bf3, a third parameter).               c
!c                                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER IDMAX,i1,i2,i,icount

      REAL, intent (in) :: BEE, PSIS,POROS,COND,VALX(NAR),PX(NAR)
      real zbar(nbdep),catdef(nbdep),bflow(nbdep)
      real, intent (out) :: bf1,bf2,bf3,taberr1,taberr2,normerr1,normerr2
      integer :: n,idep
      real suma,sumta,timean

      real catfit(nbdep),bfit(nbdep),dfit(nbdep),catmean,bfmean
      real catref(nbdep),bref(nbdep)
      real err1, err2
!      logical, intent (in) :: bug

         sumta=0.
         suma=0.
         do n=1,nar
            sumta=sumta+px(n)*valx(n)
            suma=suma+px(n)
         enddo
         timean=sumta/suma
         bf3 = timean

!c**** Loop over water table depths

         do idep=1,nbdep

!c           write(*,*) 'idep=',idep

            CALL BASIDEP(                  &
                IDEP,                      &
                BEE,PSIS,POROS,COND,       &
                VALX,PX,TIMEAN,SUMA,       &
                ZBAR,CATDEF,BFLOW)

         enddo


         i1=10   ! zbar= 0 m
         i2=35   ! zbar= 2.5 m

         bf2=zbar(i2)*SQRT(catdef(i1))               &
             /(SQRT(catdef(i2))-SQRT(catdef(i1)))
         bf1=catdef(i1)/(bf2*bf2)

         if (bf1 .le. 0) write(*,*) 'bf1 le 0 for i=',i
         if (bf2 .le. 0) write(*,*) 'bf2 le 0 for i=',i

!c Errors: Root mean square errors: only for points where catdef GT 0.5mm

         do idep=1,nbdep
            catref(idep)=0.
            bref(idep)=0.
         enddo
         catmean=0.
         bfmean=0.
         icount=0
         do idep=1,nbdep
            if (catdef(idep) .gt. lim) then
               icount=icount+1
               catref(icount)=catdef(idep)
               bref(icount)=bflow(idep)
               catfit(icount)=bf1*(zbar(idep)+bf2)             &
                   *(zbar(idep)+bf2)
               dfit(icount)=SQRT(catdef(idep)/bf1)-bf2
               bfit(icount)=cond*exp(-timean-gnu*dfit(icount)) &
                   /gnu
               catmean=catmean+catdef(idep)
               bfmean=bfmean+bflow(idep)
            endif
         enddo
         catmean=catmean/icount
         bfmean=bfmean/icount
	 if (icount.gt.1) then
         call RMSE(catref,catfit,icount,err1)
         call RMSE(bref,bfit,icount,err2)

         taberr1=err1
         taberr2=err2
         normerr1=err1/catmean
         normerr2=err2/bfmean
	 endif

contains

  ! ************************************************************************

      SUBROUTINE BASIDEP(                          &
                         IDEP,                     &
                         BEE,PSIS,POROS,COND,      &
                         VALX,PX,TIMEAN,SUMA,      &
                         ZBAR,CATDEF,BFLOW)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c This program returns the eight parameters for the areal fractioning  c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
     implicit none
      INTEGER, intent (in) :: idep
      integer nref, nind,nmax,indmin,locmax,shift,ord,locmin,ordref,width,k
      REAL, intent (in) :: BEE, PSIS, POROS, COND,VALX(NAR), PX(NAR), &
           suma,timean
      real :: dx,sumdef,dz
      real, intent (out) ::  catdef(nbdep),bflow(nbdep),zbar(idep)
      real term1,term2,sum
      real zdep(nar),locdef(nar)
!      logical bug

!c-------------------------------------------------------------------------
!c integral(f(x)dx)=1. for a pdf
!c here px=f(x)dx

      dx=valx(1)-valx(2)

      if (bug) write(*,*) 'IDEP=',IDEP,' dx=',dx, 'gnu=',gnu

!c the loops over idmax and nbdep are initiated in sta_params4.f

      zbar(idep)=float(idep-10)*slice ! zdep in meters
            
!c**** Compute array of water table depths:
      do k=1,nar
         term1=(1/gnu)*(valx(k)-timean)
         zdep(k)=AMAX1(0.,zbar(idep)-term1)
      enddo
            
!c variable change must be reflected in dx
      dz=dx/gnu

      if (bug) write(*,*) 'basidep: ok1'
      
!c**** Compute array of moisture deficits:
      do k=1,nar
         term1=(psis-zdep(k))/psis
         term1=term1**(1.-1./bee)
         term2=-psis*(bee/(bee-1.))*(term1-1.)
         locdef(k)=zdep(k)-term2
      enddo

!c**** Add deficits to produce catdef:
      sumdef=0.
      do k=1,nar
         sumdef=sumdef+locdef(k)*px(k)
      enddo
      catdef(idep)=poros*1000.*sumdef/suma

      if (bug) write(*,*) 'basidep: ok2'

      bflow(idep)=cond*exp(-timean-gnu*zbar(idep))/gnu

      if (bug) write(*,*) 'basidep: ok3'

    END SUBROUTINE BASIDEP
  
END SUBROUTINE BASE_PARAM


!*****************************************************************************

      SUBROUTINE SAT_PARAM(                                                   &
                           BEE,PSIS,POROS,COND,                               &
                           WPWET,VALX, PX, COESKEW,PFC,                       &
                           soildepth,                                         &
                           ARS1,ARS2,ARS3,                                    &
                           ARA1,ARA2,ARA3,ARA4,                               &
                           ARW1,ARW2,ARW3,ARW4,                               &
                           taberr1,taberr2,taberr3,taberr4,                   &
                           normerr1,normerr2,normerr3,normerr4,               &
                           DBG_UNIT)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c This program returns the eleven parameters for the areal fractioning c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IMPLICIT NONE

      INTEGER, intent (in) :: pfc
      REAL, intent (in) :: BEE,PSIS,POROS,COND,WPWET, &
           VALX(NAR), PX(NAR)
      REAL, intent (in) :: soildepth, COESKEW
      REAL, intent (inout) :: ARS1,ARS2,ARS3,                                 &
                           ARA1,ARA2,ARA3,ARA4,                               &
                           ARW1,ARW2,ARW3,ARW4,                               &
                           taberr1,taberr2,taberr3,taberr4,                   &
                           normerr1,normerr2,normerr3,normerr4
      INTEGER idep,n,k,i,icount,iref
      integer nar0
      integer nref, nind,nmax,indmin,locmax,shift,ord,locmin
      integer loc1,loc2,loc3,loc0,flag
      REAL VALX0(NAR), PX0(NAR)
      integer :: adjust,loc2save,inc,dec
      real sumta,suma,timean,upval,loval,profdep
      real rjunk,rjunk2
      integer, intent (in), optional :: DBG_UNIT
      real catdef(nbdep),wmin(nbdep),ar1(nbdep),aa(nbdep),aabis(nbdep)
      real ar2(nbdep),ar3(nbdep),swsrf2(nbdep),swsrf3(nbdep),rzeq(nbdep)
      real zbar0,catdef0,wmin0,RZDEP,wminsave(nbdep)

      real x1,x2,x3,x4,w1,w1_0,w2,w3,w4,ref1
      real y0,f1,f2,f3,g1,g2,g3,df,dg,dx,bf,bg,delta,z1,z2

      real nar1(nbdep),nar2(nbdep),nmean2(nbdep),neq(nbdep)
      real shape, nwm, area1,cdi,nar3(nbdep),nmean3
      real err1,err2,err3,err4,sum
      real tabact(nbdep),tabfit(nbdep)

      integer :: mp,isvd,j,first_loop
!      REAL*8, allocatable :: A(:,:),AP(:,:)
!      REAL*8, allocatable :: B(:)
      REAL*8, allocatable, target :: A(:,:)
      REAL*8, allocatable, target :: B(:)
      REAL*8, pointer             :: AP(:,:)
      REAL*8, pointer             :: BP(:)
      REAL*8 V(3,3),W(3),ANS(3),sdmax,sdmin,wbrac

      real :: cdcr1,cdcr2,term1,term2,zmet
      logical :: smooth,ars_svd_loop
      logical, parameter ::  bug=.false.
      logical, parameter :: SingValDecomp = .true.
      integer, parameter :: nl=4, nr=4, m=4, NP=50
      real :: savgol_coeff(NP)  
      integer :: savgol_ind(NP)
      integer :: nbdepl,istart

      ref1 = 100.
!      print *,'PFC', pfc   
      if (bug) write(*,*) 'starting sat_param'

      if(SingValDecomp) then
           savgol_ind(1)=0 
	   j=3
           do i=2, nl+1
              savgol_ind(i)=i-j
	      j=j+2
           end do

           j=2
           do i=nl+2, nl+nr+1
              savgol_ind(i)=i-j
	      j=j+2
           end do   
         call savgol(savgol_coeff,nl+nr+1,nl,nr,0,m)
      endif

      profdep = soildepth
      rzdep =grzdep
      profdep=profdep/1000.
      profdep=amax1(1.,profdep)
      if (rzdep .gt. .75*profdep) then
        rzdep=0.75*profdep
      end if
      
      zmet=profdep
      term1=-1.+((psis-zmet)/psis)**  &
           ((bee-1.)/bee)
      term2=psis*bee/(bee-1)
      cdcr1=1000.*poros*(zmet-(-term2*term1))
      cdcr2=(1-wpwet)*poros*1000.*zmet
!c mean of the topographic index distribution

      nar0=0
      do n=1,nar
         if (px(n) .ne. 0.) then
            nar0=nar0+1
            valx0(nar0)=valx(n)
            px0(nar0)=px(n)
         endif
      enddo

      sumta=0.
      suma=0.
      do n=1,nar0
         sumta=sumta+px0(n)*valx0(n)
         suma=suma+px0(n)
      enddo
      timean=sumta/suma

      if (bug) write(*,*) 'ok 0: sumta,suma,nar0=',sumta,suma,nar0

!c**** Loop over water table depths

         do idep=1,nbdep

            CALL FUNCIDEP(                                       &
                         NAR0,IDEP,                              &
                         BEE,PSIS,POROS,COND,RZDEP,WPWET,        &
                         VALX0,PX0,COESKEW,TIMEAN,SUMA,          &
                         CATDEF,AR1,WMIN,AA,AABIS,               &
                         AR2,AR3,SWSRF2,SWSRF3,RZEQ)             
         enddo

         nbdepl = 100
         if(catdef(50) > cdcr1 + 20.) nbdepl = 50
         if(soildepth > 6500.)  nbdepl = nbdep

         if (bug) write(*,*) 'funcidep loop ok'

!c**** for wmin's adjustment, we need an estimate of its limit toward INF
         adjust =0
         ZBAR0=10.
         CALL FUNCZBAR(                                        &
                         NAR0,ZBAR0,                           &
                         BEE,PSIS,POROS,COND,RZDEP,WPWET,      &
                         VALX0,PX0,COESKEW,TIMEAN,SUMA,        &
                         CATDEF0,WMIN0)

         if (bug) write(*,*) 'funczbar ok'

	 if (wmin0 == 0.9999900) then
	       do idep=1,nbdep-1
	         if(catdef(idep).le.cdcr1+10.) then
		 if((wmin(idep) - wmin(idep +1)) > -0.01) then 
	           wmin0=wmin(idep)
                 endif
		 endif  
	       enddo
              wmin0 = 0.1*(nint(wmin0*100000.)/10000) -0.02		    
	 endif
        
       if(present(dbg_unit)) then
          write (dbg_unit,*) nbdep,nbdepl,wmin0,cdcr1,cdcr2
          write (dbg_unit,*) catdef
          write (dbg_unit,*) ar1
          write (dbg_unit,*) wmin
       endif

!c**** AR1 adjustment: 3 points + limit in INF = 0.

         if (bug) write(*,*) 'STARTING AR1'

         ! Singular value decomposition               
         loc1=1
         loc3=nbdepl
	 loc2=loc3	 

	 do idep = 1,loc2
	    if(ar1(idep) < 1.e-10) then	
	       loc3 = idep - 1
	       exit
	    endif
	 end do	 
  
         first_loop = 0      
         ars_svd_loop = .TRUE.
         DO while (ars_svd_loop)
         
         first_loop = first_loop + 1
         mp = loc3-loc1+1
           
         allocate(A(mp,3))
         allocate(AP(mp,3))
         allocate(B(mp))
              
         a=0.
         ap=0.
         b=0.
         v=0.
         w=0.
         ans=0.
              
         do isvd=loc1,loc3
            A(isvd-loc1+1,1)=catdef(isvd)
            A(isvd-loc1+1,2)=-catdef(isvd)*ar1(isvd)
            A(isvd-loc1+1,3)=-ar1(isvd)*((catdef(isvd))**2.)
            B(isvd-loc1+1)=ar1(isvd)-1.
         end do
              
         ap = a
         call svdcmp(ap,mp,3,w,v)
         sdmax=0.
         do j=1,3
            if(w(j).gt.sdmax)sdmax=w(j)
         end do
         sdmin=sdmax*1.0e-6
         do j=1,3
            if(w(j).lt.sdmin)w(j)=0.
         end do

         call svbksb(ap,w,v,mp,3,b,ans)
              
         ars1 = real(ans(1))
         ars2 = real(ans(2))
         ars3 = real(ans(3))  

         flag=0
         call curve1(ars1,ars2,ars3,cdcr2,flag)
         deallocate (A, AP, B)

         IF(FLAG == 1) THEN
            LOC3 = NBDEP
            LOC1 =1
            IF(first_loop > 1) ars_svd_loop=.FALSE.
            ELSE
            ars_svd_loop=.FALSE.  
         ENDIF
         END DO

         IF (FLAG.EQ.1) then
  
         flag=0
	 loc1=1
         do idep=1,nbdepl
            if (catdef(idep) .le. 20.) loc1=idep
         enddo

         loc3=1
         do idep=1,nbdepl -1
            if ((ar1(idep) >= 0.0001).and.(catdef(idep) <= cdcr1)) loc3=idep + 1
         enddo

         if (loc3.le.loc1+1) then
            loc1=MIN(loc3-4,loc1-4)
            loc1=MAX(1,loc1)
         endif

!c below is what was used for no regression, but it's not equivalent to the 
!c IDL program
         loc2=AINT(loc1-1+(loc3-loc1)*3./5.)+1
         
         w1=ar1(loc1)
         w2=ar1(loc2)
         w3=ar1(loc3)

         if(w3.eq.0.)then
 95         loc3=loc3-1
            if(loc3.eq.loc2)loc2=loc2-1
               w3=ar1(loc3)
               w2=ar1(loc2)
               if(w3.eq.0.)goto 95
            endif
            w4=0.

            if((loc1.ge.loc2).or.(loc2.ge.loc3))then
               loc1=10
               loc2=14
               loc3=18
            endif

 115        x1=catdef(loc1)
            x2=catdef(loc2)
            x3=catdef(loc3)
            w1=ar1(loc1)
            w2=ar1(loc2)
            w3=ar1(loc3)
               
            if (bug) then
               write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
               write(*,*) 'x1,x2,x3=',x1,x2,x3
               write(*,*) 'w1,w2,w3=',w1,w2,w3
            endif

            y0=w4
            f1=(1.-w1)/(w1-y0)/x1
            f2=(1.-w2)/(w2-y0)/x2
            f3=(1.-w3)/(w3-y0)/x3
            g1=(1.-y0)/(w1-y0)
            g2=(1.-y0)/(w2-y0)
            g3=(1.-y0)/(w3-y0)
            df=f2-f1
            dg=g2-g1
            dx=x2-x1
            bf=f1-x1*df/dx
            bg=g1-x1*dg/dx

            ars1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
            ars2 = bf+ars1*bg
            ars3 = (df+ars1*dg)/dx

            delta=ars2*ars2-4*ars3
            upval=1.+200.*ars1
            loval=1.+200.*ars2+40000.*ars3
	    z1=0.
            z2=0.

            if (delta .ge. 0.) then !if 8
               z1=(-ars2-SQRT(delta))/2./ars3
               z2=(-ars2+SQRT(delta))/2./ars3 
            endif

            if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
               (z2 .gt. 0. .and. z1 .lt. cdcr1) .or.  &
               ((upval/loval).lt.-.01)) then   !if 7
               z1=0.
               z2=0.
               if (loc1 .eq. 10) then 
                  loc1=1
1              else  
                  loc1=1
                  do idep=1,nbdepl
                     if (catdef(idep) .gt. 60.) then
                        loc1=idep
                        if(loc1.ge.loc3-1)then
                        !                           write(*,*)'Loc1 exceeded loc3 in 2nd attempt'
                           loc1=loc3-5
                        endif
                        goto 46
                      endif
                  enddo
            endif
46          loc2=loc1+AINT(float(loc3-loc1)*3./5.)+1
            if(loc2.ge.loc3)loc2=loc3-1
            loc2save=loc2
            INC=1
            DEC=0

47          w1=ar1(loc1)
            w2=ar1(loc2)
            x1=catdef(loc1)
            x2=catdef(loc2)
            
            if (bug) then
               write(*,*) 'z1,z2=',z1,z2,' -> ar1, 2nd try'
               write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
               write(*,*) 'x1,x2,x3=',x1,x2,x3
               write(*,*) 'w1,w2,w3=',w1,w2,w3
            endif
            
            f1=(1.-w1)/(w1-y0)/(x1 + 1.e-20)
            f2=(1.-w2)/(w2-y0)/(x2 + 1.e-20)
            g1=(1.-y0)/(w1-y0 + 1.e-20 )
            g2=(1.-y0)/(w2-y0 + 1.e-20)
            df=f2-f1
            dg=g2-g1
            dx=x2-x1
            bf=f1-x1*df/dx
            bg=g1-x1*dg/dx
            
            ars1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx  + 1.e-10)
            ars2 = bf+ars1*bg
            ars3 = (df+ars1*dg)/dx
            delta=ars2*ars2-4*ars3
            upval=1.+200.*ars1
            loval=1.+200.*ars2+40000.*ars3

            if (delta .ge. 0.) then   !if 6
               z1=(-ars2-SQRT(delta))/2./ars3
               z2=(-ars2+SQRT(delta))/2./ars3
            end if

            if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
                 (z2 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
                 ((upval/loval).lt.-.01)) then  !if 5
               !c Sarith ---
               z1=0.
               z2=0.
               IF(INC.EQ.1)loc2=loc2+1
               IF(DEC.EQ.1)LOC2=LOC2-1
               if(inc.eq.1)then   !if 4
                  if(loc2.ge.loc3)then   !if 3
                     !                     WRITE(*,*)'INCREASING LOC2 FAILED'
                     INC=0
                     DEC=1
                     loc2=loc2save
                  else
                     adjust=ADJUST+1
                     goto 47
                  end if    !if 3
               endif    !if 4

               if(dec.eq.1)then   !if 2  
                  if(loc2.eq.loc1)then  !if 1
                     !                     WRITE(*,*)'Decreasing too failed'
                     INC=1
                     DEC=0
                     ars1=9999. !ars1old
                     ars2=9999. !ars2old
                     ars3=9999. !ars3old
                     !                     write(*,*) 'AR1: PROBLEM for pfc=',pfc
                  else
                     adjust=ADJUST+1
                     !c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                     goto 47
                  end if   !if 1
               endif  !if 2 
            endif     !if 5
            !c               endif    !if 6
         endif            !if 7
         
         !c         endif  !if 8
         flag=0
         call curve1(ars1,ars2,ars3,cdcr2,flag)

         IF (FLAG.EQ.1)then
            !            WRITE(*,*)'Curve problem in the catchment pfc=',pfc
            ars1=9999. 
            ars2=9999. 
            ars3=9999. 
            !                     write(*,*) 'Pick values from icatch-1'
            flag=0
         end if
      endif

         adjust=0

         if (bug) write(*,*) 'ar1 adjustment ok'

!c**** WMIN adjustment: 3 points + limit in INF = wmin0

         if (bug) write(*,*) 'STARTING WMIN'

         w4=wmin0
         y0=w4

!         write(*,*) 'wmin=',(wmin(idep),idep=1,50)

            loc1=1
            do idep=1,nbdepl
               if (catdef(idep) <= 10.) loc1=idep
             enddo
                  
            loc3=1
            do idep=1,nbdepl - 2
               if ((wmin(idep) >= wmin0).and.(catdef(idep) <= cdcr1)) loc3=idep + 2
            enddo
 
            loc2=loc1 + 2
            do idep=1,nbdepl -1 
               if ((wmin(idep) >= wmin0).and.(catdef(idep) <= cdcr1/2.))loc2=idep + 1
            enddo 

!c For global catch         
            INC=1
            DEC=0

          if(loc3.eq.loc2)loc2=loc2-2
          if(loc2 <= loc1) loc1= loc1-2
 44         loc2save=loc2
            if(loc1 < 1) then
               loc1 =1
               loc2 =2 
               loc3 =3
            endif

            w1=wmin(loc1)
            w2=wmin(loc2)
            w3=wmin(loc3)
            x1=catdef(loc1)
            x2=catdef(loc2)
            x3=catdef(loc3)
         
            if (bug) then
               write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
               write(*,*) 'x1,x2,x3=',x1,x2,x3
               write(*,*) 'w1,w2,w3,w4=',w1,w2,w3,w4
            endif

            f1=(1.-w1)/(w1-y0)/x1
            f2=(1.-w2)/(w2-y0)/x2
            f3=(1.-w3)/(w3-y0)/x3
            g1=(1.-y0)/(w1-y0)
            g2=(1.-y0)/(w2-y0)
            g3=(1.-y0)/(w3-y0)
            df=f2-f1
            dg=g2-g1
            dx=x2-x1
            bf=f1-x1*df/dx
            bg=g1-x1*dg/dx
            
            arw1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
            arw2 = bf+arw1*bg
            arw3 = (df+arw1*dg)/dx
            arw4 = y0

!c wmin=arw4+(1.-arw4)*(1.+arw1*catdef(idep))
!c     /(1.+arw2*catdef(idep)+arw3*catdef(idep)*catdef(idep))
!c we want to check the roots of the denominator

            delta=arw2*arw2-4*arw3

            if (delta .ge. 0.) then !if 8

               z1=(-arw2-SQRT(delta))/2./arw3
               z2=(-arw2+SQRT(delta))/2./arw3

               if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.           &
                   (z2 .gt. 0. .and. z1 .lt. cdcr1)) then !if 7

                  w1_0=w1
                  w1=(1.+w1_0)/2.
                  x1=x1/4.

!                  if (gnu .eq. 3.26/1.5) then 
!                        w1=(1.+w1_0)/3.               ! already difficult
!                        w3=wmin(nint(cdcr1))                ! with gnu=3.26
!                        x3=catdef(nint(cdcr1))
!                        f3=(1.-w3)/(w3-y0)/x3
!                        g3=(1.-y0)/(w3-y0)
!                  endif
                 
                  f1=(1.-w1)/(w1-y0)/x1
                  g1=(1.-y0)/(w1-y0)
                  df=f2-f1
                  dg=g2-g1
                  dx=x2-x1
                  bf=f1-x1*df/dx
                  bg=g1-x1*dg/dx

                  if (bug) then
                     write(*,*) 'z1,z2=',z1,z2,' -> wmin, 2nd try'
                     write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
                     write(*,*) 'x1,x2,x3=',x1,x2,x3
                     write(*,*) 'w1,w2,w3=',w1,w2,w3
                     write(*,*) 'wmin0=',wmin0
                  endif
               
                  arw1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
                  arw2 = bf+arw1*bg
                  arw3 = (df+arw1*dg)/dx
                  arw4 = y0
                  
                  delta=arw2*arw2-4*arw3
                  
                  if (delta .ge. 0.) then  !if 6 
                     z1=(-arw2-SQRT(delta))/2./arw3
                     z2=(-arw2+SQRT(delta))/2./arw3

                     if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.         &
                         (z2 .gt. 0. .and. z1 .lt. cdcr1)) then    !if 5  
!c Sarith ---
                     IF(INC.EQ.1)loc2=loc2+1
                     IF(DEC.EQ.1)LOC2=LOC2-1
                     if(inc.eq.1)then   !if 4
                     if(loc2.eq.loc3)then   !if 3
!                     WRITE(*,*)'INCREASING LOC2 FAILED: WMIN'
                     INC=0
                     DEC=1
                     loc2=loc2save
                     else
                        adjust=ADJUST+1
!c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                        goto 44
                        end if    !if 3
                        endif    !if 4
                     if(dec.eq.1)then   !if 2  
                     if(loc2.eq.loc1)then  !if 1
!                     WRITE(*,*)'Decreasing too failed: WMIN'
                     INC=1
                     DEC=0

                     arw1=9999. 
                     arw2=9999. 
                     arw3=9999. 
                     arw4=9999. 

                     else
                        adjust=ADJUST+1
!c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                        goto 44
                        end if   !if 1
                     endif  !if 2                         
                     endif !if 5     
                  endif   !if 6
                  
               endif !if 7
            endif !if 8
         adjust=0
!         endif ! pfc=12821
         flag=0

         call curve2(arw1,arw2,arw3,arw4,cdcr1,WPWET,flag)

         IF (FLAG.EQ.1) THEN
	     arw1=9999. !arw1old
             arw2=9999. !arw2old
             arw3=9999. !arw3old
             arw4=9999. !arw4old
             flag=0
         endif

         if(arw1==9999.) then 
! Singular Value Decomposition

         w4=wmin0
         y0=w4

            loc1=1
            loc3=nbdepl

         mp = loc3-loc1+1   

         if(mp.lt.3)then

            write(*,*)'WMIN Note: not sufficient points MP = ',mp
	    print *,w4,cdcr1,catdef(loc3),wmin(loc3)
                     arw1 = 9999.
                     arw2 = 9999.
                     arw3 = 9999.
                     arw4 = 9999.
            else
	    
            mp = 1
            istart =1
            w4 = wmin(istart)

            if(w4 <=0) then
               do idep=2,nbdepl
                  if(wmin(idep) > 0.) istart = idep
                  if(wmin(idep) > 0.) exit
               enddo
            endif

            w4 = wmin(istart)

  	    do idep=istart+1,nbdepl
!	       if(wmin(idep).lt.w4) then
               if((wmin(idep) - w4).lt.0.0005) then
	       	w4 = wmin(idep)
		mp = mp +1
	       endif		    		    	          
	    enddo             
	    loc3 = mp   
 	    allocate(A(mp,3))
            allocate(AP(mp,3))
            allocate(B(mp))
            allocate(BP(mp))             
	       smooth = .false.
	       do idep=istart,nbdepl-1
	         if(catdef(idep).le.cdcr1+10.) then
		 if((wmin(idep) - wmin(idep +1)) < -0.01) smooth = .true.   
		 endif  
	       enddo
	       if(smooth) then
	       wminsave = wmin
               ! Apply filter to input data
               do i=istart, nbdepl-nr
    	            wmin(i)=0.
                    do j=1, nl+nr+1
	                if (i+savgol_ind(j).gt.0) then  !skip left points that do not exist
	    	           wmin(i)=wmin(i)+savgol_coeff(j)*wminsave(i+savgol_ind(j))
                        endif
                    end do
               enddo	         
	       wmin (istart:istart+4) = wminsave (istart:istart+4)

	       endif

	       j = 1
	       w4 = wmin(istart)
               do isvd=1,size(wmin)
                  if (j <= mp) then 
                     if(isvd == 1) then
                        wbrac=(wmin(isvd + istart -1)-y0)/(1.-y0 + 1.e-20)
                        A(j,1)=catdef(isvd + istart -1)
                        A(j,2)=-catdef(isvd + istart -1)*wbrac
                        A(j,3)=-wbrac*((catdef(isvd + istart -1))**2.)
                        B(j)=wbrac-1.
                        j = j + 1
                     else
                        if((wmin(isvd + istart -1).lt.w4).and.(wmin(isvd + istart -1).gt.y0)) then
                           wbrac=(wmin(isvd + istart -1)-y0)/(1.-y0 + 1.e-20)
                           A(j,1)=catdef(isvd + istart -1)
                           A(j,2)=-catdef(isvd + istart -1)*wbrac
                           A(j,3)=-wbrac*((catdef(isvd + istart -1))**2.)
                           B(j)=wbrac-1.
                           w4 = wmin(isvd + istart -1)
                           j = j + 1 
                        endif
                     endif
                  endif
               end do

               j = j -1 
               mp = j
               ap => a (1:j,:)
               bp => b (1:j)
               ap(j,1) = catdef(nbdep)
               ap(j,2) = 0.
               ap(j,3) = 0.
               bp (j) = -1.

               call svdcmp(ap,mp,3,w,v)

               sdmax=0.
               do j=1,3
                  if(w(j).gt.sdmax)sdmax=w(j)
               end do

               sdmin=sdmax*1.0e-6
               do j=1,3
                  if(w(j).lt.sdmin)w(j)=0.
               end do

               call svbksb(ap,w,v,mp,3,bp,ans)

               arw1 = real(ans(1))
               arw2 = real(ans(2))
               arw3 = real(ans(3))
               arw4 = y0
               
!c wmin=arw4+(1.-arw4)*(1.+arw1*catdef(idep))
!c     /(1.+arw2*catdef(idep)+arw3*catdef(idep)*catdef(idep))
!c we want to check the roots of the denominator
               
               adjust=0         
               flag=0
               
               call curve2(arw1,arw2,arw3,arw4,cdcr1,WPWET,flag)
               
               IF (FLAG.EQ.1) THEN
                  !            WRITE(*,*)'Curve2 problem in the catchment:pfc=',pfc
                  
                  arw1 = 9999.
                  arw2 = 9999.
                  arw3 = 9999.
                  arw4 = 9999.
                  
                  flag=0
               end if
               deallocate (A,  B )
               NULLIFY    (AP, BP)
            end if
         endif
         
         if(present(dbg_unit)) then
             write (dbg_unit,*) ars1,ars2,ars3
             write (dbg_unit,*) arw1,arw2,arw3,arw4 
         endif

         if (bug) write(*,*) 'wmin adjustment ok'
         
!c**** SHAPE PARAMETER ADJUSTMENT: with a straight if coeskew > 0.25
!c                                 with 2 segments if not

         if (bug) write(*,*) 'STARTING SHAPE'

         x3=catdef(nbdepl)
         w3=aa(nbdepl)
         x1=0.

         if (coeskew .lt. 0.25) then
            w1=0.1
            loc2=20
            do idep=1,nbdepl
               if (catdef(idep) .gt. ref1) then
                  loc2=idep
                  goto 45
               endif
            enddo
 45         x2=catdef(loc2)
            w2=aabis(loc2)
            ara1 = (w1-w2)/(x1-x2)
            ara2 = w1-ara1*x1
            ara3 = (w2-w3)/(x2-x3)
            ara4 = w2-ara3*x2
         else
            w1=1.
            x2=x1
            w2=w1
            ara3 = (w2-w3)/(x2-x3)
            ara4 = w2-ara3*x2
            ara1 = ara3
            ara2 = ara4       
         endif

         if (bug) write(*,*) 'x1,w1,x2,w2,x3,w3',x1,w1,x2,w2,x3,w3

!**** RMSE checking: on ar1, ar2, swsrf2 and rzeq

         do idep=1,nbdepl
            if(catdef(idep) <= cdcr1) then
            nar1(idep)=AMIN1(1.,AMAX1(0.,(1.+ars1*catdef(idep)) &
                /(1.+ars2*catdef(idep)                          &
                +ars3*catdef(idep)*catdef(idep))))                    
                                                                 
           nwm=AMIN1(1.,AMAX1(0.,arw4+(1.-arw4)*                &
                (1.+arw1*catdef(idep))                          &
                /(1.+arw2*catdef(idep)                          &
                +arw3*catdef(idep)*catdef(idep))))

!c we have to first determine if there is one or two segments
            if (ara1 .ne. ara3) then
               cdi=(ara4-ara2)/(ara1-ara3)
            else
               cdi=0.
            endif

            if (catdef(idep) .ge. cdi) then
               shape=ara3*catdef(idep)+ara4
            else
               shape=ara1*catdef(idep)+ara2
            endif
	    shape =AMIN1(40.,shape)
            area1=exp(-shape*(1.-nwm))*(shape*(1.-nwm)+1.)

!c the threshold for truncation problems is higher than the "usual"
!c E-8 to E-10, because it plays together with the uncertainties coming 
!c from the approximation of the parameters nwm, nar1 and shape.
            if (area1 .ge. 1.-1.E-8) then
		nar1(idep)=1.
		nar2(idep)=0.
		nar3(idep)=0.
		nmean2(idep)=0.  
		nmean3=0.       
                neq(idep)=1.
             else
                
                if (nwm .gt. wpwet) then
                   nar2(idep)=1.-nar1(idep)
                else
                   nar2(idep)=AMAX1(0.,((shape*(wpwet-nwm)+1.)        &
                       *exp(-shape*(wpwet-nwm))                       &
     			- (shape*(1.-nwm)+1.)*exp(-shape*(1.-nwm)))   &
     			* (1.-nar1(idep))/(1.-area1))                  
               endif                                                        
                                                                       
               nar3(idep)=1.-nar1(idep)-nar2(idep)                          
                                                                       
               if (nar3(idep) .lt. 1.E-8) then ! for nwm le wpwet           
                                                                       
                  nmean2(idep)=AMAX1(0.,AMIN1(1.,(nwm + 2./shape +    &
                       shape*exp(-shape*(1.-nwm))*                    &
     			(nwm+nwm/shape-1.-2./shape-2./(shape*shape))) &
     			/(1.-area1)))
                   nmean3=0.

                else

!c WARNING: I think the two values below are false. 
!c But it is never used in this context, because nwm > wpwet !!
                   nmean2(idep)=AMAX1(0.,AMIN1(1.,-shape*(exp(-shape*&
                       (wpwet-nwm))* (nwm*wpwet                      &
                       +nwm/shape-wpwet*wpwet                        &
                       -2.*wpwet/shape-2./(shape*shape))             &
                       - exp(-shape*(1.-nwm))*                       &
     			(nwm+nwm/shape-1.-2./shape-2./(shape*shape)))& 
     			* (1.-nar1(idep))/(1.-area1) / (nar2(idep)+1.e-20)))        
                                                                      
                   nmean3=AMAX1(0.,AMIN1(1.,(nwm+2./shape +          &
                       shape*exp(-shape*(wpwet-nwm))*                &
     			(nwm*wpwet+nwm/shape-wpwet                   &
                       *wpwet-2.*wpwet/shape                         &
                       -2./(shape*shape))) * (1.-nar1(idep))         &
                       /(1.-area1)/(nar3(idep) + 1.e-20)))                       
               endif                                                       
                                                                      
               neq(idep)=nar1(idep)+nar2(idep)*nmean2(idep)          &
                    +nar3(idep)*nmean3
                
                if (area1 .ge. 1.-1.E-5) then
                   nmean2(idep)=1.  
                   nmean3=0.       
                   neq(idep)=1.
                endif

             endif
             endif
          enddo

          if (bug) write(*,*) 'shape adjustment ok'
!c
!c RMSE

!c ERR1
          icount=0
          iref=0
          sum=0.
          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             tabact(i)=0.
             tabfit(i)=0.
             endif
          enddo

          do i=1,nbdepl
            if(catdef(i) <= cdcr1) then 
             if (catdef(i) .gt. lim) then
               icount=icount+1
               sum=sum+ar1(i)
               tabfit(icount)=nar1(i)
               tabact(icount)=ar1(i)
             endif
             endif
          enddo

	  if(icount.gt.1) then
          sum=sum/icount
          call RMSE(tabact,tabfit,icount,err1)
          taberr1=err1
          normerr1=err1/sum
	  endif
!c ERR2
          icount=0
          iref=0
          sum=0.
          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             tabact(i)=0.
             tabfit(i)=0.
             endif
          enddo

          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             if (catdef(i) .gt. lim) then
               icount=icount+1
               sum=sum+ar2(i)
               tabfit(icount)=nar2(i)
               tabact(icount)=ar2(i)
             endif
             endif
          enddo

	  if(icount.gt.1) then
          sum=sum/icount
          call RMSE(tabact,tabfit,icount,err2)
          taberr2=err2
          normerr2=err2/sum
	  endif

!c ERR3
          icount=0
          iref=0
          sum=0.
          do i=1,nbdep
             if(catdef(i) <= cdcr1) then
             tabact(i)=0.
             tabfit(i)=0.
             endif
          enddo

          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             if (catdef(i) .gt. lim) then
               icount=icount+1
               sum=sum+swsrf2(i)
               tabfit(icount)=nmean2(i)
               tabact(icount)=swsrf2(i)
             endif
             endif
          enddo

	  if(icount.gt.1) then
          sum=sum/icount
          call RMSE(tabact,tabfit,icount,err3)
          taberr3=err3
          normerr3=err3/sum
	  endif
!c ERR4
          icount=0
          iref=0
          sum=0.
          do i=1,nbdepl
             tabact(i)=0.
             tabfit(i)=0.
          enddo

          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             if (catdef(i) .gt. lim) then
               icount=icount+1
               sum=sum+rzeq(i)
               tabfit(icount)=neq(i)
               tabact(icount)=rzeq(i)
             endif
             endif
          enddo

	  if(icount.gt.1) then 
          sum=sum/icount
          call RMSE(tabact,tabfit,icount,err4)
          taberr4=err4
          normerr4=err4/sum
       endif

     contains

  ! ******************************************************************

       !c
       SUBROUTINE CURVE1(ars1,ars2,ars3,cdcr2,flag)
         REAL ars1,ars2,ars3,y,x,yp,cdcr2
         INTEGER i,flag
         !c
         yp=1.
         if (abs(ars1+ars2+ars3).le.1.e25) then
            do i=0,CEILING(cdcr2)
               x=float(i)
               if(x > cdcr2) x = cdcr2
               y=(1.+ars1*x)/(1.+ars2*x+ars3*x*x + 1.e-20)
               if((y.gt.0.0).and.(((yp -y) .lt. -1.e-4).or.(y.gt.1.)))then
                  flag=1
                  goto 99
               endif
               yp=y
            end do
99          continue
         else
            flag=1
         endif

    end SUBROUTINE CURVE1

    
    ! ******************************************************************
    
    SUBROUTINE CURVE2(arw1,arw2,arw3,arw4,cdcr1,WPWET,flag)
      REAL arw1,arw2,arw3,arw4,y,x,yp,cdcr1, wpwet
      INTEGER i,flag
!c
      yp=1.
      if (abs(arw1+arw2+arw3+arw4).le.1.e25) then
         do i=0,CEILING(cdcr1)
            x=float(i)
            if(x > cdcr1) x = cdcr1
            y=arw4+(1.-arw4)*(1.+arw1*x)/(1.+arw2*x+arw3*x*x + 1.e-20)
            if ((y .lt. wpwet).or.((yp -y) .lt. -1.e-4).or.(y.gt.1.)) then
               flag=1
               goto 99
            endif
            yp=y
         end do
99       continue
      else
         flag=1
      endif
    end SUBROUTINE CURVE2
! -----------------------------------------------------------------------------------


      SUBROUTINE svbksb(u,w,v,m,n,b,x) 
        implicit none
        INTEGER m,mp,n,np,NMAX 
        REAL*8 b(m),u(m,n),v(n,n),w(n),x(n) 
        PARAMETER (NMAX=500)  !Maximum anticipated value of n
        !------------------------------------------------------------------------------------------- 
        ! Solves A · X = B for a vector X, where A is specified by the arrays u, w, v as returned by 
        ! svdcmp. m and n are the dimensions of a, and will be equal for square matrices. b(1:m) is 
        ! the input right-hand side. x(1:n) is the output solution vector. No input quantities are 
        ! destroyed, so the routine may be called sequentially with different b’s. 
        !-------------------------------------------------------------------------------------------

        INTEGER i,j,jj 
        REAL*8 s,tmp(NMAX) 
        do j=1,n !Calculate UTB. 
           s=0. 
           if(w(j).ne.0.)then !Nonzero result only if wj is nonzero. 
              do i=1,m 
                 s=s+u(i,j)*b(i) 
              end do
              s=s/(w(j) + 1.d-20) !This is the divide by wj . 
           endif
           tmp(j)=s 
        end do
        do j=1,n !Matrix multiply by V to get answer. 
           s=0. 
           do jj=1,n 
              s=s+v(j,jj)*tmp(jj) 
           end do
           x(j)=s 
        end do
        return 
      END SUBROUTINE svbksb

!---------------------------------------------------------------------

      SUBROUTINE svdcmp(a,m,n,w,v) 
        implicit none
        INTEGER m,n,NMAX 
        REAL*8, intent (inout)  :: a(m,n)
        REAL*8, intent (out) :: v(n,n),w(n) 
        PARAMETER (NMAX=500)  !Maximum anticipated value of n. 
        !-------------------------------------------------------------------------------------- 
        ! Given a matrix A(1:m,1:n), this routine computes its singular value decomposition, 
        ! A = U · W · Vt. The matrix U replaces A on output. The diagonal matrix of singular 
        ! values W is output as a vector W(1:n). The matrix V (not the transpose Vt) is output 
        ! as V(1:n,1:n). 
        !--------------------------------------------------------------------------------------

        INTEGER i,its,j,jj,k,l,nm 
        REAL*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX) 
        real*8, parameter :: EPS=epsilon(1.0d0)
        g=0.d0  !Householder reduction to bidiagonal form. 
        scale=0.d0 
        anorm=0.d0 
        c =0.d0 
        f =0.d0 
        g =0.d0 
        h =0.d0 
        s =0.d0 
        x =0.d0 
        y =0.d0 
        z =0.d0 
        rv1=0.d0 
        w = 0.d0 
        v = 0.d0 
        do i=1,n 
           l=i+1 
           rv1(i)=scale*g 
           g=0.d0 
           s=0.d0 
           scale=0.d0 
           if(i.le.m)then 
              do k=i,m 
                 scale=scale+abs(a(k,i)) 
              end do
              if(scale.ne.0.d0)then 
                 do k=i,m 
                    a(k,i)=a(k,i)/scale 
                    s=s+a(k,i)*a(k,i) 
                 end do
                 f=a(i,i) 
                 g=-dsign(dsqrt(s),f) 
                 h=f*g-s 
                 a(i,i)=f-g 
                 do j=l,n 
                    s=0.d0 
                    do k=i,m 
                       s=s+a(k,i)*a(k,j) 
                    end do
                    f=s/h 
                    do k=i,m 
                       a(k,j)=a(k,j)+f*a(k,i) 
                    end do
                 end do
                 do k=i,m 
                    a(k,i)=scale*a(k,i) 
                 end do
              endif
           endif
           w(i)=scale *g 
           g=0.d0 
           s=0.d0 
           scale=0.d0 
           if((i.le.m).and.(i.ne.n))then 
              do k=l,n 
                 scale=scale+abs(a(i,k)) 
              end do
              if(scale.ne.0.d0)then 
                 do k=l,n 
                    a(i,k)=a(i,k)/scale 
                    s=s+a(i,k)*a(i,k) 
                 end do
                 f=a(i,l) 
                 g=-sign(sqrt(s),f) 
                 h=f*g-s 
                 a(i,l)=f-g 
                 do k=l,n 
                    rv1(k)=a(i,k)/h 
                 end do
                 do j=l,m 
                    s=0.d0 
                    do k=l,n 
                       s=s+a(j,k)*a(i,k) 
                    end do
                    do k=l,n 
                       a(j,k)=a(j,k)+s*rv1(k) 
                    end do
                 end do
                 do k=l,n 
                    a(i,k)=scale*a(i,k) 
                 end do
              endif
           endif
           anorm=max(anorm,(abs(w(i))+abs(rv1(i)))) 
        end do !do i=1,n
        
        do i=n,1,-1 !Accumulation of right-hand transformations. 
           if(i.lt.n)then 
              if(g.ne.0.d0)then 
                 do j=l,n       !Double division to avoid possible underflow. 
                    v(j,i)=(a(i,j)/a(i,l))/g 
                 end do
                 do j=l,n 
                    s=0.d0 
                    do k=l,n 
                       s=s+a(i,k)*v(k,j) 
                    end do
                    do k=l,n 
                       v(k,j)=v(k,j)+s*v(k,i) 
                    end do
                 end do
              endif
              do j=l,n 
                 v(i,j)=0.d0 
                 v(j,i)=0.d0 
              end do
           endif
           v(i,i)=1.d0
           g=rv1(i) 
           l=i 
        end do
        
        do i=min(m,n),1,-1 !Accumulation of left-hand transformations. 
           l=i+1 
           g=w(i) 
           do j=l,n 
              a(i,j)=0.d0 
           end do
           if(g.ne.0.d0)then 
              g=1.d0/g 
              do j=l,n 
                 s=0.d0 
                 do k=l,m 
                    s=s+a(k,i)*a(k,j) 
                 end do
                 f=(s/a(i,i))*g 
                 do k=i,m 
                    a(k,j)=a(k,j)+f*a(k,i) 
                 end do
              end do
              do j=i,m 
                 a(j,i)=a(j,i)*g 
              end do
           else
              do j= i,m 
                 a(j,i)=0.d0 
              end do
           endif
           a(i,i)=a(i,i)+1.d0 
        end do
        
        do k=n,1,-1 !Diagonalization of the bidiagonal form: Loop over 
           !singular values, and over allowed iterations. 
           do its=1,30 
              do l=k,1,-1 !Test for splitting. 
                 nm=l-1 !Note that rv1(1) is always zero.
                 if( abs(rv1(l)) <= EPS*anorm ) goto 2 
                 if( abs(w(nm) ) <= EPS*anorm ) goto 1  
              end do
1             c=0.d0 !Cancellation of rv1(l), if l > 1. 
              s=1.d0 
              do i=l,k 
                 f=s*rv1(i) 
                 rv1(i)=c*rv1(i) 
                 if( abs(f) <= EPS*anorm ) goto 2 
                 g=w(i) 
                 h=pythag(f,g) 
                 w(i)=h 
                 h=1.d0/h 
                 c= (g*h) 
                 s=-(f*h) 
                 do j=1,m 
                    y=a(j,nm) 
                    z=a(j,i) 
                    a(j,nm)=(y*c)+(z*s) 
                    a(j,i)=-(y*s)+(z*c) 
                 end do
              end do
2             z=w(k) 
              if(l.eq.k)then   !Convergence. 
                 if(z.lt.0.d0)then !Singular value is made nonnegative. 
                    w(k)=-z 
                    do j=1,n 
                       v(j,k)=-v(j,k) 
                    end do
                 endif
                 goto 3 
              endif
              if(its.eq.30) print *, 'no convergence in svdcmp' 
 !             if(its.ge.4)  print *, 'its = ',its
              x=w(l) !Shift from bottom 2-by-2 minor. 
              nm=k-1 
              y=w(nm) 
              g=rv1(nm) 
              h=rv1(k) 
              f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y) 
              g=pythag(f,1.d0) 
              f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x 
              c=1.d0 !Next QR transformation: 
              s=1.d0 
              do j=l,nm 
                 i=j+1 
                 g=rv1(i) 
                 y=w(i) 
                 h=s*g 
                 g=c*g 
                 z=pythag(f,h) 
                 rv1(j)=z 
                 c=f/z 
                 s=h/z 
                 f= (x*c)+(g*s) 
                 g=-(x*s)+(g*c) 
                 h=y*s 
                 y=y*c 
                 do jj=1,n 
                    x=v(jj,j) 
                    z=v(jj,i) 
                    v(jj,j)= (x*c)+(z*s) 
                    v(jj,i)=-(x*s)+(z*c) 
                 end do
                 z=pythag(f,h) 
                 w(j)=z !Rotation can be arbitrary if z = 0. 
                 if(z.ne.0.d0)then 
                    z=1.d0/z 
                    c=f*z 
                    s=h*z 
                 endif
                 f= (c*g)+(s*y) 
                 x=-(s*g)+(c*y) 
                 do jj=1,m 
                    y=a(jj,j) 
                    z=a(jj,i) 
                    a(jj,j)= (y*c)+(z*s) 
                    a(jj,i)=-(y*s)+(z*c) 
                 end do
              end do !j=l;nm 
              rv1(l)=0.d0 
              rv1(k)=f 
              w(k)=x 
           end do !its=1,30
3          continue 
        end do !k=n,1,-1 
        return 
      END SUBROUTINE svdcmp

!
! ________________________________________________________________________________
!     
      REAL*8 FUNCTION pythag(a,b) 
        REAL*8 a,b 
        !Computes sqrt(a**2 + b**2) without destructive underflow or overflow.
        REAL*8 absa,absb 
        absa=abs(a) 
        absb=abs(b) 
        if(absa.gt.absb)then 
           pythag=absa*sqrt(1.+(absb/absa)**2) 
        else 
           if(absb.eq.0.)then 
              pythag=0. 
           else
              pythag=absb*sqrt(1.+(absa/absb)**2) 
           endif
        endif
        return 
      END FUNCTION pythag
!
! ________________________________________________________________________________
!     

      SUBROUTINE savgol(c,np,nl,nr,ld,m)
      implicit none
      INTEGER ld,m,nl,np,nr,MMAX 
      real c(np) 
      PARAMETER (MMAX=6)
!-------------------------------------------------------------------------------------------- 
!USES lubksb,ludcmp given below. 
!Returns in c(1:np), in wrap-around order (see reference) consistent with the argument respns 
!in routine convlv, a set of Savitzky-Golay filter coefficients. nl is the number of leftward 
!(past) data points used, while nr is the number of rightward (future) data points, making 
!the total number of data points used nl +nr+1. ld is the order of the derivative desired 
!(e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also 
!equal to the highest conserved moment; usual values are m = 2 or m = 4. 
!--------------------------------------------------------------------------------------------
INTEGER d,icode,imj,ipj,j,k,kk,mm,indx(MMAX+1) 
real fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX  & 
	  .or.nl+nr.lt.m) pause ' Bad args in savgol.' 
	do ipj=0,2*m        !Set up the normal equations of the desired leastsquares fit. 
    	sum=0. 
    if(ipj.eq.0) sum=1. 
    do k=1,nr 
      sum=sum+dfloat(k)**ipj 
    end do 
    do k=1,nl 
      sum=sum+dfloat(-k)**ipj 
    end do 
    mm=min(ipj,2*m-ipj) 
    do imj=-mm,mm,2 
      a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum 
    end do 
  end do

  call ludcmp(a,m+1,MMAX+1,indx,d,icode)    !Solve them: LU decomposition. 

  do j=1,m+1 
    b(j)=0. 
  end do 
  b(ld+1)=1.      !Right-hand side vector is unit vector, depending on which derivative we want. 

  call lubksb(a,m+1,MMAX+1,indx,b)   !Backsubstitute, giving one row of the inverse matrix. 

  do kk=1,np                         !Zero the output array (it may be bigger than the number 
    c(kk)=0.                         !of coefficients).  
  end do 
  do k=-nl,nr                        !Each Savitzky-Golay coefficient is the dot product 
    sum=b(1)                         !of powers of an integer with the inverse matrix row. 
    fac=1. 
    do mm=1,m 
      fac=fac*k 
      sum=sum+b(mm+1)*fac 
    end do 
    kk=mod(np-k,np)+1                !Store in wrap-around order. 
    c(kk)=sum 
  end do
  return 
END SUBROUTINE savgol

!***************************************************************
!* Given an N x N matrix A, this routine replaces it by the LU *
!* decomposition of a rowwise permutation of itself. A and N   *
!* are input. INDX is an output vector which records the row   *
!* permutation effected by the partial pivoting; D is output   *
!* as -1 or 1, depending on whether the number of row inter-   *
!* changes was even or odd, respectively. This routine is used *
!* in combination with LUBKSB to solve linear equations or to  *
!* invert a matrix. Return code is 1, if matrix is singular.   *
!***************************************************************
 Subroutine LUDCMP(A,N,NP,INDX,D,CODE)
INTEGER, PARAMETER :: NMAX=100
REAL, PARAMETER :: TINY=1E-12
 real  AMAX,DUM, SUM, A(NP,NP),VV(NMAX)
 INTEGER CODE, D, INDX(N),NP,N,I,J,K,IMAX

 D=1; CODE=0

 DO I=1,N
   AMAX=0.
   DO J=1,N
     IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = 1. / AMAX
 END DO ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*ABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(ABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.N) THEN
     DUM = 1. / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

 RETURN
 END Subroutine LUDCMP


!******************************************************************
!* Solves the set of N linear equations A . X = B.  Here A is     *
!* input, not as the matrix A but rather as its LU decomposition, *
!* determined by the routine LUDCMP. INDX is input as the permuta-*
!* tion vector returned by LUDCMP. B is input as the right-hand   *
!* side vector B, and returns with the solution vector X. A, N and*
!* INDX are not modified by this routine and can be used for suc- *
!* cessive calls with different right-hand sides. This routine is *
!* also efficient for plain matrix inversion.                     *
!******************************************************************
 Subroutine LUBKSB(A,N,NP,INDX,B)
 INTEGER :: II,I,J,LL,N,NP
 real  SUM, A(NP,NP),B(N)
 INTEGER INDX(N)

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN
END Subroutine LUBKSB
         
  END SUBROUTINE SAT_PARAM
!



! ******************************************************************

      subroutine tgen (                         &
          TOPMEAN,TOPVAR,TOPSKEW,               &
          STO,ACO,COESKEW)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                         c
! The difference between tgen4 and tgen3 is that tgen4 deals with arrays  c
! of topmean, topvar and topskew and 2-dim arrays of st and ac.           c
!                                                                         c
! This routine determine the theoretical gamma distribution for the       c
! soil-topographic indexes (Sivapalan et al., 1987), knowing the three    c
! first moments, the min and the max of the observed topographic indexes  c
! in a given catchment.                                                   c
!                                                                         c
! Routine from Dave Wolock.                                               c
! Modified by Agnes (11-06-98): we don't use min and max anymore, and     c
! this strongly improves the behavior for negative skewnesses. It also    c
! improves in general the matching of the moments.                        c
!                                                                         c
! We also add a correction on the skewness to have gamma distributions    c
! that start and end from the x-axis. It is based on the fact that if     c
! TOPETA=1, the gamma is an exponential distribution, and if TOPETA<1,    c
! then the gamma distribution increases towards the infinite when x       c
! decreases towards 0.                                                    c
! To eliminate some numerical pb due to teh discretization of the gamma   c
! distribution, we choose skewness=MAX(MIN(1.9, skewness),-1.6)           c
!                                                                         c
! WE MAY NEED TO COMPUTE IN DOUBLE RESOLUTION !!!! BECAUSE OF THE SMALL   c
! BIN WIDTH
!                                                                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

      IMPLICIT NONE

      real, parameter :: VALMAX=50.
      REAL, intent (in) :: TOPMEAN,TOPVAR,TOPSKEW
      REAL, intent (out) :: COESKEW
      REAL, dimension (NAR), intent (out) :: STO,ACO

      INTEGER I
      REAL ST(NAR),AC(NAR)
      REAL TOPETA,TOPLAM,TOPSCAL,GAMLN,SCALE,ACLN
      real cumac, cum2,cum3

!-------------------------------------------------------------------------

! topmean is the mean of the ln(a/tanB) distribution
! topvar is the variance (2nd moment centerd around the mean) of the ...
! topskew is the skew (3rd moment centerd around the mean) of the ...
! compute the coefficient of skew or skewness (coeskew)

         COESKEW=TOPSKEW/TOPVAR**1.5
         if (coeskew .ge. 0.) then
            COESKEW=AMAX1(0.005, AMIN1(1.9, COESKEW))
         else
            COESKEW=AMAX1(-1.6, AMIN1(-0.005, COESKEW))
         endif

! compute the gamma parameters, eta (topeta) and lambda (toplam), and topscal
! which is the translation parameter

         TOPETA=4./COESKEW**2
         TOPLAM=SQRT(TOPETA)/SQRT(TOPVAR)
         TOPSCAL=TOPMEAN-TOPETA/TOPLAM

! evaluate the gamma function

         CALL GAMMLN(TOPETA,GAMLN)

         CUMAC=0.0

! compute the frequency distribution of ln(a/tanB)
! st(i) are the values of ln(a/tanB)
! ac(i) are the relative frequency values (they should sum to 1)

         DO I=1,NAR
         
            ST(I)=(FLOAT(I)-0.95)*(VALMAX-TOPSCAL)/FLOAT(NAR)+TOPSCAL
            SCALE=ST(I)-TOPSCAL

! below is the logarithmic form of the gamma distribution; this is required 
! because the numerical estimate of the logarithm of the gamma function 
! is more stable than the one of the gamma function.
          
            ACLN=TOPETA*ALOG(TOPLAM)+(TOPETA-1.)*ALOG(SCALE)  &
                -TOPLAM*SCALE-GAMLN
                      
            IF(ACLN.LT.-10.) THEN
               AC(I)=0.
            ELSE
               AC(I)=EXP(ACLN)
            ENDIF

            CUMAC=CUMAC+AC(I)

         ENDDO

! we want the relative frequencies to sum 1.

         IF (CUMAC.eq.0.) THEN
!            write(*,*) 'distrib sum=',CUMAC
            stop
         endif
         CUM2=0.
         DO I=1,NAR
            AC(I) = AC(I) / CUMAC
            CUM2=CUM2+AC(I)
         ENDDO
      
! if the real distribution of the topographic indices is negativeley skewed, 
! we symetrize the gamma distribution (depending on coeskew**2 and always 
! positively skewed), centering on topmean, which preserves topmean and
! topvar, and re-establishes a negative skewness.

         IF (COESKEW.LT.0.) then

            do i=1,nar
               STO(I)=2.*TOPMEAN-ST(I)
               ACO(I)=AC(I)

            enddo
         ELSE
!            if (n .eq. idmax) then
!               write(*,*) 'last catchment'
!            endif
            do i=1,nar
               STO(I)=ST(-I+NAR+1)
               ACO(I)=AC(-I+NAR+1)
            enddo
         ENDIF

!         sum=0.
!         do i=1,nar
!            sum=sum+sto(i)*aco(i)
!         end do

!         sum=0.
!         do i=1,nar
!            sum=sum+aco(i)
!         end do


    END subroutine tgen

  
  ! ********************************************************************

    SUBROUTINE GAMMLN(XX,GAMLN)
      
      implicit none
      DOUBLE PRECISION :: COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      REAL, intent(in) :: XX 
      REAL, intent(out) :: GAMLN
      integer :: j

      DATA COF /76.18009173D0,-86.50532033D0,24.01409822D0,    &
         -1.231739516D0,.120858003D-2,-.536382D-5/
      STP = 2.50662827465D0
      HALF= 0.5D0
      ONE = 1.0D0
      FPF = 5.5D0
      
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE

      DO  J=1,6
         X=X+ONE
         SER=SER+COF(J)/X
      END DO

      GAMLN=TMP+LOG(STP*SER)
      
    END SUBROUTINE GAMMLN
  
  ! ********************************************************************

    SUBROUTINE FUNCIDEP(                                         &
                         NAR0,IDEP,                              &!I
                         BEE,PSIS,POROS,COND,RZDEP,WPWET,        &!I
                         VALX,PX,COESKEW,TIMEAN,SUMA,            &!I
                         CATDEF,AR1,WMIN,AA,AABIS,               &!O
                         AR2,AR3,SWSRF2,SWSRF3,RZEQ)             

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c This program returns the eight parameters for the areal fractioning  c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer, intent (in) :: NAR0,idep 
      REAL, intent (in) :: BEE, PSIS, POROS, COND, RZDEP, WPWET, COESKEW
      REAL, intent (inout) ::  VALX(NAR), PX(NAR),TIMEAN,SUMA
!      logical, intent(in) :: bug
      real, dimension (nbdep), intent (inout) :: CATDEF,AR1,WMIN,AA,   &
           AABIS,AR2,AR3,SWSRF2,SWSRF3,RZEQ 
      INTEGER :: width, nref, nind,nmax,indmin,locmax,shift,ord,locmin,ordref
      integer :: indimax10,indmin0,k,n,n1,n2
      real dx,zbar

      real test,term1,term2,sum
      real zdep(nar),locdef(nar),wrz(nar),frcunsat
      real valtest(nbdep,nar),ptest(nbdep,nar),denstest(nbdep,nar)
      real dtest(nbdep,nar),cump
      real x1,x2,y1,y2,wa,wb
      real densaux(nar),densaux2(nar),densmax,aux10
      real :: dz, sumdef
!c-------------------------------------------------------------------------

!c integral(f(x)dx)=1. for a pdf
!c here px=f(x)dx
      dx=valx(1)-valx(2)

      if (bug) write(*,*) 'IDEP=',IDEP,' dx=',dx

!c the loops over idmax and nbdep are initiated in sta_params4.f

      zbar=float(idep-10)*slice ! zdep in meters
            
!c**** Compute array of water table depths:
      do k=1,nar0
         term1=(1/gnu)*(valx(k)-timean)
         zdep(k)=AMAX1(0.,zbar-term1)
      enddo
            
!c variable change must be reflected in dx
      dz=dx/gnu

      if (bug) write(*,*) 'funcidep: ok1'
      
!c**** Compute array of moisture deficits:
      do k=1,nar0
         term1=(psis-zdep(k))/psis
         term1=term1**(1.-1./bee)
         term2=-psis*(bee/(bee-1.))*(term1-1.)
         locdef(k)=zdep(k)-term2
      enddo

!c**** Add deficits to produce catdef:
      sumdef=0.
      do k=1,nar0
         sumdef=sumdef+locdef(k)*px(k)
      enddo
      catdef(idep)=poros*1000.*sumdef/suma

      if (bug) write(*,*) 'funcidep: ok2'

!c**** Compute array of root zone moisture (degree of wetness in root zone):
      do k=1,nar0              
         term1=((psis-zdep(k))/psis)              &
             **(1.-1./bee)
        if(zdep(k).le.0.) then
           wrz(k)=1.
        elseif(zdep(k)-rzdep.lt.0.) then
           term2=(-psis/zdep(k))*(bee/(bee-1.))  &
                *(term1-1.)
           frcunsat=zdep(k)/rzdep
           wrz(k)=frcunsat*term2+(1.-frcunsat)*1.
        else
           term2=((psis-zdep(k)+rzdep)           &
                /psis)**(1.-1./bee)
           wrz(k)=(-psis/rzdep)*(bee/            &
                (bee-1.))*(term1-term2)
         endif

      enddo

      if (bug) write(*,*) 'funcidep: ok3'

!c**** compute the densities and dx
!c**** we use a usefull property that is due to the construction of the 
!c**** gamma distribution in tgen3.f : this distribution is continuous, 
!c**** with decreasing values on ln(a/tanb) when n goes from 1 to nar0

!c first we gather in the same bin all the bins with values ge 1 
      nref=1
      nind=1
      ptest(idep,1)=0.
      do k=1,nar0
         if (wrz(k) .eq. 1.) then
            nref=nref+1
            ptest(idep,1) = ptest(idep,1) + px(k)
         endif
      enddo
      if (nref .gt. 1) then
         nind=2
         valtest(idep,1)=1.
      endif
      nmax=nar0-nref+nind
      if (bug) write(*,*) 'nmax,nind,nar0,nref=',nmax,nind,nar0,nref
      
!c definition of the probabilities ptest
      if (nmax .eq. 1) then     ! all the bins have values ge 1 
         dtest(idep,1) = 0.0001
         ptest(idep,1) = 1.
      else                      ! distribution in ar2/ar3
         do n=0,nmax-nind
            valtest(idep,nind+n)=wrz(nref+n)
            ptest(idep,nind+n)=px(nref+n)
         enddo
         
!c we have to define dtest, the size of each bin
         if (nmax .eq. 2) then
            dtest(idep,2) = valtest(idep,1)-valtest(idep,2)
            dtest(idep,1) = dtest(idep,2)/2.
         else                   ! nmax .gt. 2
            do n=2,nmax-1                     
               dtest(idep,n)=(valtest(idep,n-1)-valtest(idep,n+1))/2.   
            enddo
               dtest(idep,1) = dtest(idep,2)/2.
            dtest(idep,nmax) = dtest(idep,nmax-1)
         endif
      endif

      if (bug) write(*,*) 'funcidep: ok4'

!c we can now define the probability density: denstest=ptest/dtest
!c where ptest is the probability and dtest the size of the bin
      do n=1,nmax
         if (ptest(idep,n) .eq. 0.) then
            denstest(idep,n)=0.
         else
            denstest(idep,n)=ptest(idep,n)/dtest(idep,n)
         endif
      enddo

      if (bug) write(*,*) 'funcidep: ok5'

!c NOW we can estimate the parameters for the approximated distrib
!c from the actual distrib

!c 1. AR1=saturated area and AR2 and AR3 + averages of the RZ wetness 
!c    in the different fractions

      ar1(idep)=0.
      ar2(idep)=0.
      ar3(idep)=0.
      swsrf3(idep)=0.
      swsrf2(idep)=0.
      rzeq(idep)=0.
    
      if(valtest(idep,1).eq.1.) ar1(idep)=dtest(idep,1)*denstest(idep,1)
      
      if (nmax .gt. 1) then 
         do n=nind,nmax           
            if (valtest(idep,n) .lt. wpwet) then
               ar3(idep)=ar3(idep)+denstest(idep,n)*dtest(idep,n)
               swsrf3(idep)=swsrf3(idep)+valtest(idep,n)*         &
                    denstest(idep,n)*dtest(idep,n)
            else
               ar2(idep)=ar2(idep)+denstest(idep,n)*dtest(idep,n)
               swsrf2(idep)=swsrf2(idep)+valtest(idep,n)*         &
                   denstest(idep,n)*dtest(idep,n)
            endif
         enddo
      endif
         
      test=ar1(idep)+ar2(idep)+ar3(idep)
      if (test .gt. 1.+1.e-5 .or. test .lt. 1.-1.e-5) then
!         write(*,*) 'PROBLEM at depth ',zbar
!         write(*,*) '  ar1+ar2+ar3=',test
!         write(*,*) '  ar1=',ar1(idep),' ar2=',ar2(idep),' ar3=', &
!             ar3(idep)
      endif
         
      ar1(idep)=ar1(idep)/test
      ar2(idep)=ar2(idep)/test
      ar3(idep)=ar3(idep)/test
      if (ar2(idep) .ne. 0.) swsrf2(idep)=swsrf2(idep)/ar2(idep)
      if (ar3(idep) .ne. 0.) swsrf3(idep)=swsrf3(idep)/ar3(idep)

      rzeq(idep)=ar1(idep)+ar2(idep)*swsrf2(idep)+ar3(idep)*swsrf3(idep)
      
      if (bug) write(*,*) 'funcidep: ok6'

!c 2. Maximum density -> shape parameter 
!c                    -> wmin 

      locmax=3
      shift=15
      ordref=1
      do n=1,nmax
         densaux2(n)=denstest(idep,n)
      enddo
         
      if (nmax .ge. shift*2) then
               
!c we start with sliding mean to facilitate the search for the maximum
        
         ord=MIN(ordref,nmax/shift)
	 
         call smtot(densaux2,nmax,ord,densaux) 
!	 print *,nmax,ord,shift,densaux(shift-14),shift-14,size(densaux)
         do n=nmax,shift,-1
            if (densaux(n) .gt. densaux(n-1) .and.     &
                densaux(n) .gt. densaux(n-2) .and.     &
                densaux(n) .gt. densaux(n-3) .and.     &
                densaux(n) .gt. densaux(n-4) .and.     &
                densaux(n) .gt. densaux(n-5) .and.     &
                densaux(n) .gt. densaux(n-6) .and.     &
                densaux(n) .gt. densaux(n-7) .and.     &
                densaux(n) .gt. densaux(n-8) .and.     &
                densaux(n) .gt. densaux(n-9) .and.     &
                densaux(n) .gt. densaux(n-10) .and.    &
                densaux(n) .gt. densaux(n-11) .and.    &
                densaux(n) .gt. densaux(n-12) .and.    &
                densaux(n) .gt. densaux(n-13) .and.    &
                densaux(n) .gt. densaux(n-14))then ! .and.    &
!                densaux(n) .gt. densaux(n-15)) then
               locmax=n
               goto 30
            endif
         enddo
         
      else

         aux10=-9999.
         indimax10=3
         do n=1,nmax
            if (densaux2(n) .gt. aux10) then
               aux10=densaux2(n)
               indimax10=n
            endif
         enddo
         locmax=MAX(3,indimax10)

      endif      ! if (nmax .ge. shift+1) 
                  
 30   densmax=denstest(idep,locmax)
      aa(idep)=exp(1.)*densmax

      if (bug) write(*,*) 'funcidep: ok7'

!c WMIN=lowest value where the density is strictly gt densmax/100.

      indmin=1
      indmin0=0
      do n=1,nmax
         if (denstest(idep,n) .gt. 0.) indmin0=n
         if (denstest(idep,n) .gt. densmax/100. .and.    &
             valtest(idep,n) .lt. valtest(idep,locmax)) indmin=n
      enddo
      if (indmin .eq.0) indmin=indmin0

      if (indmin .le. 2) then
         wmin(idep) = 0.99999
      else
         x1=valtest(idep,indmin)
         wmin(idep)=x1
      endif

      if (bug) write(*,*) 'funcidep: ok8; first wmin=',wmin(idep)

!c for negative or low coeskew the previous wmin doesn't give good results...
!c wmin is higher !!!

      if (coeskew .lt. 1. ) then

         if (locmax .gt. 3 .and. indmin .ge. locmax+4) then
            n2=MAX(locmax+1,(indmin-locmax)/2+locmax)
            x2=valtest(idep,n2)
            y2=denstest(idep,n2)
            n1=locmax
            x1=valtest(idep,n1)
            y1=denstest(idep,n1)
            wa=(y2-y1)/(x2-x1)
            wb=y1-wa*x1
            wmin(idep)=AMAX1(wmin(idep),-wb/wa)
         endif

!c wmin is even higher in some cases !!!
         if (coeskew .lt. 0.2 ) wmin(idep)=wmin(idep)+0.01
         
      endif  

      if (bug) write(*,*) 'funcidep: ok9; 2nd wmin=',wmin(idep)
      
      if (valtest(idep,locmax) .le. wmin(idep)) then ! doesn't make sense
         wmin(idep)=valtest(idep,locmax)-dx
      endif
      aabis(idep)=1./(valtest(idep,locmax)-wmin(idep)+1.e-20)

      if (bug) write(*,*) 'funcidep: ok10'

    END SUBROUTINE FUNCIDEP
  
  ! ********************************************************************

      SUBROUTINE FUNCZBAR(                                     &   
                         NAR0,ZBAR,                            &
                         BEE,PSIS,POROS,COND,RZDEP,WPWET,      &
                         VALX,PX,COESKEW,TIMEAN,SUMA,          &
                         CATDEF,WMIN)                          

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c This program returns the eight parameters for the areal fractioning  c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none                   
      INTEGER , intent (in) :: NAR0
      integer nref,nind,nmax,indmin,locmax,shift,ord,locmin,ordref
      integer indimax10,indmin0
      REAL, intent (in) :: BEE, PSIS, POROS, COND, RZDEP, WPWET, COESKEW
      REAL, intent (inout) ::  VALX(NAR), PX(NAR),TIMEAN,SUMA,zbar
      real, intent (inout) :: catdef,wmin

      REAL  dx,dz,sumdef
      real term1,term2
      real zdep(nar),locdef(nar),wrz(nar),frcunsat
      real valtest(nar),ptest(nar),denstest(nar),dtest(nar)
      real x1,x2,y1,y2,wa,wb
      integer n1,n2,k,n
      real densaux(nar),densaux2(nar),densmax,aux10

!c-------------------------------------------------------------------------
!c integral(f(x)dx)=1. for a pdf
!c here px=f(x)dx
      dx=valx(1)-valx(2)
            
!c**** Compute array of water table depths:
      do k=1,nar0
         term1=(1/gnu)*(valx(k)-timean)
         zdep(k)=AMAX1(0.,zbar-term1)
      enddo

!c variable change must be reflected in dx
      dz=dx/gnu
      
!c**** Compute array of moisture deficits:
      do k=1,nar0
         term1=(psis-zdep(k))/psis
         term1=term1**(1.-1./bee)
         term2=-psis*(bee/(bee-1.))*(term1-1.)
         locdef(k)=zdep(k)-term2
      enddo

!c**** Add deficits to produce catdef:
      sumdef=0.
      do k=1,nar0
         sumdef=sumdef+locdef(k)*px(k)
      enddo
      catdef=poros*1000.*sumdef/suma

!c**** Compute array of root zone moisture (degree of wetness in root zone):
      do k=1,nar0              
         term1=((psis-zdep(k))/psis)  &
             **(1.-1./bee)
         if(zdep(k).le.0.) then
            wrz(k)=1.
         elseif(zdep(k)-rzdep.lt.0.) then
            term2=(-psis/zdep(k))*(bee/(bee-1.))   &
                *(term1-1.)
            frcunsat=zdep(k)/rzdep
            wrz(k)=frcunsat*term2+(1.-frcunsat)*1.
         else
            term2=((psis-zdep(k)+rzdep)     &
                /psis)**(1.-1./bee)
            wrz(k)=(-psis/rzdep)*(bee/      &
                (bee-1.))*(term1-term2)
         endif
      enddo

!c**** compute the densities and dx
!c**** we use a usefull property that is due to the construction of the 
!c**** gamma distribution in tgen3.f : this distribution is continuous, 
!c**** with decreasing values on ln(a/tanb) when n goes from 1 to nar0
!c first we gather in the same bin all the bins with values ge 1 
      nref=1
      nind=1
      ptest(1)=0.
      do k=1,nar0
         if (wrz(k) .eq. 1.) then
            nref=nref+1
            ptest(1) = ptest(1) + px(k)
         endif
      enddo
      if (nref .gt. 1) then
         nind=2
         valtest(1)=1.
      endif
      nmax=nar0-nref+nind
      
!c definition of the probabilities ptest
      if (nmax .eq. 1) then     ! all the bins have values ge 1 
         dtest(1) = 0.0001
         ptest(1) = 1.
      else                      ! distribution in ar2/ar3
         do n=0,nmax-nind
            valtest(nind+n)=wrz(nref+n)
            ptest(nind+n)=px(nref+n)
         enddo
         
!c we have to define dtest, the size of each bin
         if (nmax .eq. 2) then
            dtest(2) = valtest(1)-valtest(2)
            dtest(1) = dtest(2)/2.
         else                   ! nmax .gt. 2
            do n=2,nmax-1
               dtest(n)=(valtest(n-1)-valtest(n+1))/2.            
            enddo             
            dtest(1) = dtest(2)/2.
            dtest(nmax) = dtest(nmax-1)
         endif
      endif

!c we can now define the probability density: denstest=ptest/dtest
!c where ptest is the probability and dtest the size of the bin
      do n=1,nmax
         if (ptest(n) .eq. 0.) then
            denstest(n)=0.
         else
            denstest(n)=ptest(n)/dtest(n)
         endif
      enddo

!c NOW we can estimate the parameters for the approximated distrib
!c from the actual distrib

!c 2. Maximum density -> shape parameter 
!c                    -> wmin 

      locmax=3
      shift=15
      ordref=1
      do n=1,nmax
         densaux2(n)=denstest(n)
      enddo

      if (nmax .ge. shift*2) then
               
!c we start with sliding mean to facilitate the search for the maximum
         
         ord=MIN(ordref,nmax/shift)
         call smtot(densaux2,nmax,ord,densaux)

         do n=nmax,shift,-1
            if (densaux(n) .gt. densaux(n-1) .and.         &
                densaux(n) .gt. densaux(n-2) .and.         &
                densaux(n) .gt. densaux(n-3) .and.         &
                densaux(n) .gt. densaux(n-4) .and.         &
                densaux(n) .gt. densaux(n-5) .and.         &
                densaux(n) .gt. densaux(n-6) .and.         &
                densaux(n) .gt. densaux(n-7) .and.         &
                densaux(n) .gt. densaux(n-8) .and.         &
                densaux(n) .gt. densaux(n-9) .and.         &
                densaux(n) .gt. densaux(n-10) .and.        &
                densaux(n) .gt. densaux(n-11) .and.        &
                densaux(n) .gt. densaux(n-12) .and.        &
                densaux(n) .gt. densaux(n-13) .and.        &
                densaux(n) .gt. densaux(n-14) .and.        &
                densaux(n) .gt. densaux(n-15)) then
               locmax=n
               goto 30
            endif
         enddo

      else

         aux10=-9999.
         indimax10=3
         do n=1,nmax
            if (densaux2(n) .gt. aux10) then
               aux10=densaux2(n)
               indimax10=n
            endif
         enddo
         locmax=MAX(3,indimax10)

      endif      ! if (nmax .ge. shift+1) 
         
 30   densmax=denstest(locmax)

!c WMIN=lowest value where the density is strictly gt densmax/100.

      indmin=1
      indmin0=0
      do n=1,nmax
         if (denstest(n) .gt. 0.) indmin0=n
         if (denstest(n) .gt. densmax/100. .and.         &
             valtest(n) .lt. valtest(locmax)) indmin=n
      enddo
      if (indmin .eq. 0) indmin=indmin0

      if (indmin .le. 2) then
         wmin = 0.99999
      else
         x1=valtest(indmin)
         wmin=x1
      endif

!c for negative or low coeskew the previous wmin doesn't give good results...
!c wmin is higher !!!

      if (coeskew .lt. 1. ) then
         
         if (locmax .gt. 3 .and. indmin .ge. locmax+4) then
            
            n2=MAX(locmax+1,(indmin-locmax)/2+locmax)
            x2=valtest(n2)
            y2=denstest(n2)
            n1=locmax
            x1=valtest(n1)
            y1=denstest(n1)
            wa=(y2-y1)/(x2-x1)
            wb=y1-wa*x1
            wmin=AMAX1(wmin,-wb/wa)
         endif

!c wmin is even higher in some cases !!!
         if (coeskew .lt. 0.2 ) wmin=wmin+0.01
         
      endif  

    END SUBROUTINE FUNCZBAR

! ******************************************************************

       SUBROUTINE RMSE(XX,YY,LEN,ERROR)

!c---------------------------------------------------------------------------
!c Computes the root-mean square error ERROR between two one-dimensional
!c random variables XX and YY of same length LEN
!c---------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, intent (in) :: LEN
      REAL, intent (in) ::  XX(LEN),YY(LEN)
      REAL, intent (out) :: ERROR
      INTEGER :: I

!c---------------------------------------------------------------------------     
      error=0.
      do i=1,len
         if(abs(xx(i)-yy(i)) >=1.e-10) then
         error=error+(xx(i)-yy(i))*(xx(i)-yy(i))
         endif
      enddo
      error=SQRT(error/float(len))

    END SUBROUTINE RMSE

! ******************************************************************
       SUBROUTINE SMTOT(XX,LEN,ORD,YY)

!c---------------------------------------------------------------------------
!c Runs a sliding average of order ORD through the one-dimensional array XX
!c of length LEN and returns the smoothed YY
!!c---------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, intent (in) ::  LEN
      INTEGER :: ORD,WIDTH,i,ini,n,end
      REAL, intent (in) :: XX(NAR)
      REAL, intent (out) :: YY(NAR)

!c---------------------------------------------------------------------------     
      do i=1,nar
         yy(i)=0.
      enddo

      width=ord*2+1
      if (width .gt. len/2) then
         write(*,*) 'the order for the sliding average is too large !!!'
         write(*,*) 'regard with the length of the array to be smoothed'
         stop
      endif

      do i=1,len
         ini=MAX(1,i-ord)
         end=MIN(len,i+ord)
         yy(i)=0.
         do n=ini,end
            yy(i)=yy(i)+xx(n)
         enddo
         yy(i)=yy(i)/(end-ini+1)
      enddo

    END SUBROUTINE SMTOT



end module CLSM_param_routines
