!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module CLSMJ32_parmsMod
!BOP
!
! !MODULE: CLSMJ32_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read Catchment 
!  parameter datasets.
!
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  Catchment LSM  data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  20 Nov 2012: K. Arsenault: Updated for Catchment LSM parameters.
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_gfracMod
  use LDT_albedoMod
  use LDT_paramDataMod
  use LDT_logMod
  use CLSM_util, only : init_geos2lis_mapping, LDT_g5map, G52LIS, LISv2g
  use CLSM_param_routines, only : create_CLSM_parameters
  use mod_HWSD_STATSGO2_texture, ONLY : &
       derive_CLSM_HWSD_soiltypes
  

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: catchmentParms_init_J32    !allocates memory for required structures
  public :: catchmentParms_writeHeader_J32
  public :: catchmentParms_writeData_J32

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: CLSMJ32_struc

  type, public :: clsmJ32_type_dec

     character*100 :: albnirfile
     character*100 :: albvisfile
     real          :: catchparms_gridDesc(20)
     character*50  :: catchparms_proj
     character*50  :: catchparms_gridtransform
!  - Catchment (J3.2):
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

     type(LDT_paramEntry) :: porosity       ! Soil porosity
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
     type(LDT_paramEntry) :: albnirdir   ! Albedo NIR direct scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albnirdif   ! Albedo NIR diffuse scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albvisdir   ! Albedo VIS direct scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albvisdif   ! Albedo VIS diffuse scale factor (CLSM F2.5)

  end type clsmJ32_type_dec

  type(clsmJ32_type_dec), allocatable :: CLSMJ32_struc(:)

contains

!BOP
! 
! !ROUTINE: catchmentParms_init_J32
! \label{catchmentParms_init_J32}
! 
! !INTERFACE:
  subroutine catchmentParms_init_J32
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: &! LDT_catchparmsOptChecks, &
                       LDT_gridOptChecks


! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the Catchment LSM Parameter datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[catchmentParms](\ref{catchmentParms}) \newline
!    calls the registry to invoke the Catchment parms setup methods. 
!  \end{description}
!
!EOP
    implicit none
    integer      :: n, c, r, ntiles, i, glpnr, glpnc
    integer      :: rc
    character*50 :: catchparms_proj
    real         :: catchparms_gridDesc(LDT_rc%nnest, 20)
    real         :: param_grid(20)
    real, pointer, dimension (:) :: ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4, &
         ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3, tsa1, tsa2,tsb1, tsb2, GNU,  &
         ATAU2, BTAU2, ATAU, BTAU, BEE, POROS, WPWET, PSIS, KS, SOILDEPTH
    integer, pointer, dimension (:) :: SOIL_TOP, SOIL_COM

! ________________________________________________________

    write(LDT_logunit,*)" - - - - - - - - - Catchment LSM Parameters - - - - - - - - - - - -"
    
    allocate(CLSMJ32_struc(LDT_rc%nnest))

    ! (1) initialize LDT_g5map
    ! ------------------------

    call init_geos2lis_mapping 
    
    ! (2) Derive soil types, atau and btau
    ! ------------------------------------

    NTILES = LDT_g5map%NT_GEOS
 
    allocate (BEE       (1:NTILES))
    allocate (POROS     (1:NTILES))
    allocate (WPWET     (1:NTILES))
    allocate (PSIS      (1:NTILES))
    allocate (KS        (1:NTILES))
    allocate (SOILDEPTH (1:NTILES))
    allocate (ATAU2     (1:NTILES))
    allocate (BTAU2     (1:NTILES))
    allocate (ATAU      (1:NTILES))
    allocate (BTAU      (1:NTILES))
    allocate (SOIL_TOP  (1:NTILES))
    allocate (SOIL_COM  (1:NTILES))

    call derive_CLSM_HWSD_soiltypes (NTILES, BEE, POROS, WPWET, PSIS, KS, SOILDEPTH, &
         ATAU2, BTAU2, ATAU, BTAU, SOIL_TOP, SOIL_COM)


    ! (3) Derive CLSM model - the parameters that were in ar.new, bf.dat, ts.dat 
    ! ---------------------------------------------------------------------------

    allocate (GNU  (1:NTILES))
    allocate (ARS1 (1:NTILES))
    allocate (ARS2 (1:NTILES))
    allocate (ARS3 (1:NTILES))
    allocate (ARA1 (1:NTILES))
    allocate (ARA2 (1:NTILES))
    allocate (ARA3 (1:NTILES))
    allocate (ARA4 (1:NTILES))
    allocate (ARW1 (1:NTILES))
    allocate (ARW2 (1:NTILES))
    allocate (ARW3 (1:NTILES))
    allocate (ARW4 (1:NTILES))
    allocate (BF1  (1:NTILES))
    allocate (BF2  (1:NTILES))
    allocate (BF3  (1:NTILES))
    allocate (TSA1 (1:NTILES))
    allocate (TSA2 (1:NTILES))
    allocate (TSB1 (1:NTILES))
    allocate (TSB2 (1:NTILES))

    call create_CLSM_parameters (NTILES, BEE, POROS, WPWET, PSIS, KS, SOILDEPTH, & 
         SOIL_TOP, SOIL_COM, GNU, ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4,           &
         ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3, tsa1, tsa2,tsb1, tsb2)
    

    ! (4) write out gridded arrays in LIS input file
    ! ----------------------------------------------

    write(LDT_logunit,'(A60)')'[CLSM catchmentParms_init_J32] writing LIS input file ...'

    do n=1,LDT_rc%nnest  

       param_grid(:) = LDT_rc%mask_gridDesc(n,:)
       glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
       glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9))  + 1

       call set_param_attribs(CLSMJ32_struc(n)%porosity,"POROSITY", &
            full_name="CLSM soil porosity")

       call set_param_attribs(CLSMJ32_struc(n)%bexp,"BEXP", &
            full_name="CLSM Bexp Clapp-Hornberger parameter")

       call set_param_attribs(CLSMJ32_struc(n)%psisat,"PSISAT", &
            full_name="CLSM saturated soil moisture potential")

       call set_param_attribs(CLSMJ32_struc(n)%wpwet,"WPWET", &
            full_name="CLSM wilting point wetness")

       call set_param_attribs(CLSMJ32_struc(n)%ksat,"KSAT",&
            units="ms-1",full_name="CLSM saturated hydraulic conductivity")

       call set_param_attribs(CLSMJ32_struc(n)%gnu,"GNUCLSM",&
            units="m-1", &
            full_name="CLSM vertical transm. decay term")

       call set_param_attribs(CLSMJ32_struc(n)%ars1,"ARS1CLSM",&
            units="m2kg-1", &
            full_name="CLSM (ARS) wetness parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%ara1,"ARA1CLSM", &
            units="m2kg-1", &
            full_name="CLSM (ARA) topographic shape parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%arw1,"ARW1CLSM",&
            units="m2kg-1", &
            full_name="CLSM (ARW) minimum theta parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%bf1,"BF1CLSM",&
            units="kgm-4", &
            full_name="CLSM (BF) baseflow topographic parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%tsa1,"TSA1CLSM",&
            full_name="CLSM (TS) water transfer parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%atau,"ATAUCLSM",&
            full_name="CLSM (TAU) topographic tau parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%bdrckdpth,"BEDROCKDEPTH",&
            units="mm", &
            full_name="CLSM depth to bedrock" )

       call set_param_attribs(CLSMJ32_struc(n)%albnirdir,"ALBNIRDIR",&
            vlevels=12, &
            full_name="CLSM alb near-IR (direct) scale factor" )

       call set_param_attribs(CLSMJ32_struc(n)%albvisdir,"ALBVISDIR",&
            vlevels=12, &
            full_name="CLSM alb visible (direct) scale factor" )

       call set_param_attribs(CLSMJ32_struc(n)%albnirdif,"ALBNIRDIF",&
            vlevels=12, &
            full_name="CLSM alb near-IR (diffuse) scale factor" )

       call set_param_attribs(CLSMJ32_struc(n)%albvisdif,"ALBVISDIF",&
            vlevels=12, &
            full_name="CLSM alb near-IR (diffuse) scale factor" )

       allocate(CLSMJ32_struc(n)%ksat%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%ksat%num_bins))
       
       allocate(CLSMJ32_struc(n)%psisat%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%psisat%num_bins))
       
       allocate(CLSMJ32_struc(n)%bexp%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%bexp%num_bins))
       
       allocate(CLSMJ32_struc(n)%wpwet%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%wpwet%num_bins))
       allocate(CLSMJ32_struc(n)%bdrckdpth%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%bdrckdpth%num_bins))
       allocate(CLSMJ32_struc(n)%gnu%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%gnu%vlevels))

       allocate(CLSMJ32_struc(n)%albnirdir%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%albnirdir%vlevels))
       allocate(CLSMJ32_struc(n)%albvisdir%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%albvisdir%vlevels))
       allocate(CLSMJ32_struc(n)%albnirdif%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%albnirdif%vlevels))
       allocate(CLSMJ32_struc(n)%albvisdif%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%albvisdif%vlevels))
       
        ! Fill in derived parameter entries:
        ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "ARS2CLSM", &
            "Catchment wetness parameters","m2kg-1", &
            CLSMJ32_struc(n)%ars1, &
            CLSMJ32_struc(n)%ars2 )
       
       call populate_param_attribs( "ARS3CLSM", &
            "Catchment wetness parameters","m2kg-1", &
            CLSMJ32_struc(n)%ars1, &
            CLSMJ32_struc(n)%ars3 )
       
       allocate(CLSMJ32_struc(n)%ars1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%ars1%vlevels))
       allocate(CLSMJ32_struc(n)%ars2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%ars2%vlevels))
       allocate(CLSMJ32_struc(n)%ars3%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%ars3%vlevels))
       
       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "ARA2CLSM", &
            "Catchment topographic shape parameters","m2kg-1", &
            CLSMJ32_struc(n)%ara1, &
            CLSMJ32_struc(n)%ara2 )
       
       call populate_param_attribs( "ARA3CLSM", &
            "Catchment topographic shape parameters","m2kg-1", &
            CLSMJ32_struc(n)%ara1, &
            CLSMJ32_struc(n)%ara3 )
       
       call populate_param_attribs( "ARA4CLSM", &
            "Catchment topographic shape parameters","m2kg-1", &
            CLSMJ32_struc(n)%ara1, &
            CLSMJ32_struc(n)%ara4 )
       
       allocate(CLSMJ32_struc(n)%ara1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%ara1%vlevels))
       allocate(CLSMJ32_struc(n)%ara2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%ara2%vlevels))
       allocate(CLSMJ32_struc(n)%ara3%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%ara3%vlevels))
       allocate(CLSMJ32_struc(n)%ara4%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%ara4%vlevels))

       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "ARW2CLSM", &
            "Catchment minimum theta parameters","m2kg-1", &
            CLSMJ32_struc(n)%arw1, &
            CLSMJ32_struc(n)%arw2 )
       
       call populate_param_attribs( "ARW3CLSM", &
            "Catchment minimum theta parameters","m2kg-1", &
            CLSMJ32_struc(n)%arw1, &
            CLSMJ32_struc(n)%arw3 )
       
       call populate_param_attribs( "ARW4CLSM", &
            "Catchment minimum theta parameters","m2kg-1", &
            CLSMJ32_struc(n)%arw1, &
            CLSMJ32_struc(n)%arw4 )          
       
       allocate(CLSMJ32_struc(n)%arw1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%arw1%vlevels))
       allocate(CLSMJ32_struc(n)%arw2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%arw2%vlevels))
       allocate(CLSMJ32_struc(n)%arw3%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%arw3%vlevels))
       allocate(CLSMJ32_struc(n)%arw4%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%arw4%vlevels))
       
       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "BF2CLSM", &
            "Catchment baseflow topographic parameters","m", &
            CLSMJ32_struc(n)%bf1, &
            CLSMJ32_struc(n)%bf2 )
       
       call populate_param_attribs( "BF3CLSM", &
            "Catchment baseflow topographic parameters","log(m)", &
            CLSMJ32_struc(n)%bf1, &
            CLSMJ32_struc(n)%bf3 )
       
       allocate(CLSMJ32_struc(n)%bf1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%bf1%vlevels))
       allocate(CLSMJ32_struc(n)%bf2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%bf2%vlevels))
       allocate(CLSMJ32_struc(n)%bf3%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%bf3%vlevels))
       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "TSA2CLSM", &
            "Catchment water transfer parameters","-", &
            CLSMJ32_struc(n)%tsa1, &
            CLSMJ32_struc(n)%tsa2 )
       
       call populate_param_attribs( "TSB1CLSM", &
            "Catchment water transfer parameters","-", &
            CLSMJ32_struc(n)%tsa1, &
            CLSMJ32_struc(n)%tsb1 )
       
       call populate_param_attribs( "TSB2CLSM", &
            "Catchment water transfer parameters","-", &
            CLSMJ32_struc(n)%tsa1, &
            CLSMJ32_struc(n)%tsb2 )
       
       allocate(CLSMJ32_struc(n)%tsa1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%tsa1%vlevels))
       allocate(CLSMJ32_struc(n)%tsa2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%tsa2%vlevels))
       allocate(CLSMJ32_struc(n)%tsb1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%tsb1%vlevels))
       allocate(CLSMJ32_struc(n)%tsb2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%tsb2%vlevels))
       
       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "BTAUCLSM", &
            "Catchment topographic tau parameters","-", &
            CLSMJ32_struc(n)%atau, &
            CLSMJ32_struc(n)%btau )
       
       allocate(CLSMJ32_struc(n)%atau%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%atau%vlevels))
       allocate(CLSMJ32_struc(n)%btau%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%btau%vlevels))
    enddo
    write(LDT_logunit,'(A60)')'[CLSM catchmentParms_init_J32] writing LIS input file2...'

!    call ESMF_ConfigFindLabel(LDT_config,"CLSMJ32 top soil layer depth:",rc=rc)
!    do n=1,LDT_rc%nnest
!       call ESMF_ConfigGetAttribute(LDT_config,CLSMJ32_struc(n)%dzsfcrd,rc=rc)
!       call LDT_verify(rc,'CLSMJ32 top soil layer depth: not specified')
!    enddo

    call ESMF_ConfigGetAttribute(LDT_config,catchparms_proj,&
         label="CLSMJ32 map projection:",rc=rc)
    call LDT_verify(rc,'CLSMJ32 map projection: option not specified in the config file')
    CLSMJ32_struc(:)%catchparms_proj = catchparms_proj

    call ESMF_ConfigFindLabel(LDT_config,"CLSMJ32 spatial transform:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            CLSMJ32_struc(n)%catchparms_gridtransform,&
            rc=rc)
       call LDT_verify(rc,'CLSMJ32 transform: option not specified in the config file')
    enddo
    
    call LDT_readDomainConfigSpecs("CLSMJ32",catchparms_proj,&
         catchparms_gridDesc)
    
    do n=1,LDT_rc%nnest
       CLSMJ32_struc(n)%catchparms_gridDesc = catchparms_gridDesc(n,:)
    enddo

!-- Read in Catchment Parameter Datasets:

    do n=1,LDT_rc%nnest

       call LDT_gridOptChecks( n, "CLSMJ32", &
            CLSMJ32_struc(n)%catchparms_gridtransform, &
            CLSMJ32_struc(n)%catchparms_proj, &
            CLSMJ32_struc(n)%catchparms_gridDesc(9) )

       CLSMJ32_struc(n)%ars1%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (ars1))
       CLSMJ32_struc(n)%ars2%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (ars2))
       CLSMJ32_struc(n)%ars3%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (ars3))
       CLSMJ32_struc(n)%ara1%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (ara1))
       CLSMJ32_struc(n)%ara2%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (ara2))
       CLSMJ32_struc(n)%ara3%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (ara3))
       CLSMJ32_struc(n)%ara4%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (ara4))
       CLSMJ32_struc(n)%arw1%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (arw1))
       CLSMJ32_struc(n)%arw2%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (arw2))
       CLSMJ32_struc(n)%arw3%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (arw3))
       CLSMJ32_struc(n)%arw4%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (arw4))
       CLSMJ32_struc(n)%bf1%value (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (bf1 ))
       CLSMJ32_struc(n)%bf2%value (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (bf2 ))
       CLSMJ32_struc(n)%bf3%value (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (bf3 ))
       CLSMJ32_struc(n)%tsa1%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (tsa1))
       CLSMJ32_struc(n)%tsa2%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (tsa2))
       CLSMJ32_struc(n)%tsb1%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (tsb1))
       CLSMJ32_struc(n)%tsb2%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (tsb2))
       CLSMJ32_struc(n)%atau%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (atau))
       CLSMJ32_struc(n)%btau%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (btau))
       CLSMJ32_struc(n)%psisat%value   (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (psis))
       CLSMJ32_struc(n)%bexp%value     (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (bee ))
       CLSMJ32_struc(n)%wpwet%value    (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (wpwet))
       CLSMJ32_struc(n)%ksat%value     (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (Ks))
       CLSMJ32_struc(n)%bdrckdpth%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (soildepth))
       CLSMJ32_struc(n)%porosity%value (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (poros))
       CLSMJ32_struc(n)%gnu%value      (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (gnu))
       
       call populate_param_attribs( "ALBNIRDIFF", &
           "Alb near-IR diffuse scale factor", "-",  &
            CLSMJ32_struc(n)%albnirdir, &
            CLSMJ32_struc(n)%albnirdif )
       
       call populate_param_attribs( "ALBVISDIFF", &
            "Alb visible diffuse scale factor", "-",  &
            CLSMJ32_struc(n)%albvisdir, &
            CLSMJ32_struc(n)%albvisdif )

    enddo

    write(LDT_logunit,'(A60)')'[CLSM catchmentParms_init_J32] writing LIS input file3 ...'

    call ESMF_ConfigFindLabel(LDT_config,"Albedo NIR factor file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            CLSMJ32_struc(n)%albnirfile,rc=rc)
       call LDT_verify(rc,"Albedo NIR factor file: not specified")
    enddo

  ! Albedo VIS:
    call ESMF_ConfigFindLabel(LDT_config,"Albedo VIS factor file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            CLSMJ32_struc(n)%albvisfile,rc=rc)
       call LDT_verify(rc,"Albedo VIS factor file: not specified")
    enddo
    
    do n=1,LDT_rc%nnest
       write(LDT_logunit,*) 'Reading '//trim(CLSMJ32_struc(n)%albnirfile)
       call read_CLSMF25_albnir(&
            n,CLSMJ32_struc(n)%albnirdir%value, &   ! Direct
            CLSMJ32_struc(n)%albnirdif%value, &   ! Diffuse
            LDT_LSMparam_struc(n)%landmask%value )
       write(LDT_logunit,*) 'Done reading '//&
            trim(CLSMJ32_struc(n)%albnirfile)
       
       write(LDT_logunit,*) 'Reading '//trim(CLSMJ32_struc(n)%albvisfile)
       call read_CLSMF25_albvis(&
            n,CLSMJ32_struc(n)%albvisdir%value, &   ! Direct
            CLSMJ32_struc(n)%albvisdif%value, &   ! Diffuse
            LDT_LSMparam_struc(n)%landmask%value )
       write(LDT_logunit,*) 'Done reading '//trim(CLSMJ32_struc(n)%albvisfile)

       LDT_rc%monthlyData(n) = .true.
       LDT_albedo_struc(n)%albInterval = "monthly"
       LDT_gfrac_struc(n)%gfracInterval = "monthly"
    enddo
    write(LDT_logunit,'(A60)')'[CLSM catchmentParms_init_J32] leaving ...'

  end subroutine catchmentParms_init_J32

  subroutine catchmentParms_writeHeader_J32(n,ftn,dimID,monthID)

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
       if(LDT_gfrac_struc(n)%gfracInterval.eq."monthly") then !monthly
          t_dimID(3) = monthID
       endif
    end if

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%gnu)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%ars1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%ars2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%ars3)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%ara1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%ara2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%ara3)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%ara4)
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%arw1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%arw2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%arw3)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%arw4)
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%bf1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%bf2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%bf3)
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%tsa1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%tsa2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%tsb1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%tsb2)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%atau)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMJ32_struc(n)%btau)

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMJ32_struc(n)%psisat)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMJ32_struc(n)%ksat)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMJ32_struc(n)%bexp)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMJ32_struc(n)%wpwet)
!    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
!         LDT_LSMparam_struc(n)%quartz)
!    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
!         LDT_LSMparam_struc(n)%soildepth)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMJ32_struc(n)%bdrckdpth)

 !- Albedo NIR scale factors:
    if( LDT_albedo_struc(n)%albInterval.eq."monthly" ) t_dimID(3) = monthID
    call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
         CLSMJ32_struc(n)%albnirdif)
    
 !- Albedo VIS scale factors:
    if( LDT_albedo_struc(n)%albInterval.eq."monthly" ) t_dimID(3) = monthID
    call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
         CLSMJ32_struc(n)%albvisdif)

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ALBEDO_DATA_INTERVAL", &
         LDT_albedo_struc(n)%albInterval))
    
  end subroutine catchmentParms_writeHeader_J32

  subroutine catchmentParms_writeData_J32(n,ftn)

    integer   :: n 
    integer   :: ftn

    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%gnu)

    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%ars1)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%ars2)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%ars3)

    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%ara1)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%ara2)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%ara3)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%ara4)

    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%arw1)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%arw2)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%arw3)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%arw4)

    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%bf1)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%bf2)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%bf3)

    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%tsa1)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%tsa2)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%tsb1)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%tsb2)

    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%atau)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%btau)

    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%psisat)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%ksat)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%bexp)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%wpwet)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%porosity)  ! overwriting LDT_soil_struc.F90 created variable

!    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%quartz)
!    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%soildepth)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%bdrckdpth)

    
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%albnirdir)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%albnirdif)
    
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%albvisdir)
    call LDT_writeNETCDFdata(n,ftn,CLSMJ32_struc(n)%albvisdif)

  end subroutine catchmentParms_writeData_J32


!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name, vlevels, &
                 units, full_name )

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name
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
   paramEntry%source = "CLSMJ3.2"
   paramEntry%units = trim(unit_temp)
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

  ! ---------------------------------------

  SUBROUTINE create_CLSM_PARAMS

    implicit none

    real, dimension (:), allocatable :: TOPMEAN, TOPVAR, TOPSKEW
      
!         ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4,ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,   &
!         tsa1, tsa2,tsb1, tsb2, 

  end SUBROUTINE create_CLSM_PARAMS

end module CLSMJ32_parmsMod
