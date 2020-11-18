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
  use CLSM_param_routines, only : create_CLSM_parameters,    &
       clsm_type_dec, CLSMJ32_struc,                         &
       set_param_attribs => set_CLSM_param_attribs, jpl_canoph
  use mod_HWSD_STATSGO2_texture, ONLY : &
       derive_CLSM_HWSD_soiltypes
  use LDT_ClimateBCsReader, ONLY : ClimateBCsReader
  
  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: catchmentParms_init_J32    !allocates memory for required structures

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------

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
         ATAU2, BTAU2, ATAU, BTAU, BEE, POROS, WPWET, PSIS, KS, SOILDEPTH, Z2CH
    integer, pointer, dimension (:) :: SOIL_TOP, SOIL_COM
    character*20           :: albInterval, source, proj
    type (ClimateBCsReader):: bcr
! ________________________________________________________

    write(LDT_logunit,*)" - - - - - - - - - Catchment LSM Parameters - - - - - - - - - - - -"
    
    allocate(CLSMJ32_struc(LDT_rc%nnest))

    ! (1) initialize LDT_g5map
    ! ------------------------

    if (.not.LDT_g5map%init) call init_geos2lis_mapping 
    
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
    allocate (Z2CH      (1:NTILES))

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
    call jpl_canoph (ntiles, z2ch)

    ! (4) write out gridded arrays in LIS input file
    ! ----------------------------------------------

    write(LDT_logunit,'(A60)')'[CLSM catchmentParms_init_J32] writing LIS input file ...'

    do n=1,LDT_rc%nnest  

       param_grid(:) = LDT_rc%mask_gridDesc(n,:)
       glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
       glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9))  + 1

       call set_param_attribs(CLSMJ32_struc(n)%z2,"Z2CH","CLSMJ3.2", &
            full_name="CLSM vegetation height", units="m")

       call set_param_attribs(CLSMJ32_struc(n)%porosity,"POROSITY","CLSMJ3.2", &
            full_name="CLSM soil porosity")

       call set_param_attribs(CLSMJ32_struc(n)%bexp,"BEXP", "CLSMJ3.2",&
            full_name="CLSM Bexp Clapp-Hornberger parameter")

       call set_param_attribs(CLSMJ32_struc(n)%psisat,"PSISAT", "CLSMJ3.2",&
            full_name="CLSM saturated soil moisture potential")

       call set_param_attribs(CLSMJ32_struc(n)%wpwet,"WPWET", "CLSMJ3.2",&
            full_name="CLSM wilting point wetness")

       call set_param_attribs(CLSMJ32_struc(n)%ksat,"KSAT","CLSMJ3.2",&
            units="ms-1",full_name="CLSM saturated hydraulic conductivity")

       call set_param_attribs(CLSMJ32_struc(n)%gnu,"GNUCLSM","CLSMJ3.2",&
            units="m-1", &
            full_name="CLSM vertical transm. decay term")

       call set_param_attribs(CLSMJ32_struc(n)%ars1,"ARS1CLSM","CLSMJ3.2",&
            units="m2kg-1", &
            full_name="CLSM (ARS) wetness parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%ara1,"ARA1CLSM", "CLSMJ3.2",&
            units="m2kg-1", &
            full_name="CLSM (ARA) topographic shape parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%arw1,"ARW1CLSM","CLSMJ3.2",&
            units="m2kg-1", &
            full_name="CLSM (ARW) minimum theta parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%bf1,"BF1CLSM","CLSMJ3.2",&
            units="kgm-4", &
            full_name="CLSM (BF) baseflow topographic parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%tsa1,"TSA1CLSM","CLSMJ3.2",&
            full_name="CLSM (TS) water transfer parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%atau,"ATAUCLSM","CLSMJ3.2",&
            full_name="CLSM (TAU) topographic tau parameters" )

       call set_param_attribs(CLSMJ32_struc(n)%bdrckdpth,"BEDROCKDEPTH","CLSMJ3.2",&
            units="mm", &
            full_name="CLSM depth to bedrock" )

       call set_param_attribs(CLSMJ32_struc(n)%albnirdir,"ALBNIRDIR","CLSMJ3.2",&
            vlevels=12, &
            full_name="CLSM alb near-IR (direct) scale factor" )

       call set_param_attribs(CLSMJ32_struc(n)%albvisdir,"ALBVISDIR","CLSMJ3.2",&
            vlevels=12, &
            full_name="CLSM alb visible (direct) scale factor" )

       call set_param_attribs(CLSMJ32_struc(n)%albnirdif,"ALBNIRDIF","CLSMJ3.2",&
            vlevels=12, &
            full_name="CLSM alb near-IR (diffuse) scale factor" )

       call set_param_attribs(CLSMJ32_struc(n)%albvisdif,"ALBVISDIF","CLSMJ3.2",&
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

       allocate(CLSMJ32_struc(n)%porosity%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%porosity%num_bins))
       
       allocate(CLSMJ32_struc(n)%wpwet%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%wpwet%num_bins))
       allocate(CLSMJ32_struc(n)%bdrckdpth%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%bdrckdpth%num_bins))
       allocate(CLSMJ32_struc(n)%z2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMJ32_struc(n)%z2%num_bins))
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

    call ESMF_ConfigFindLabel(LDT_config,"CLSMJ32 top soil layer depth:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,CLSMJ32_struc(n)%dzsfcrd, label="CLSMJ32 top soil layer depth:", DEFAULT = 0.05, rc=rc)
       call LDT_verify(rc,'CLSMJ32 top soil layer depth: not specified')
    enddo

    call ESMF_ConfigGetAttribute(LDT_config,catchparms_proj,&
         label="CLSMJ32 map projection:",DEFAULT="latlon", RC=RC)
    call LDT_verify(rc,'CLSMJ32 map projection: option not specified in the config file')
    CLSMJ32_struc(:)%catchparms_proj = catchparms_proj

    call ESMF_ConfigFindLabel(LDT_config,"CLSMJ32 spatial transform:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            CLSMJ32_struc(n)%catchparms_gridtransform,label="CLSMJ32 spatial transform:",DEFAULT="none", &
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
       if(CLSMJ32_struc(n)%dzsfcrd == 0.02) then
          CLSMJ32_struc(n)%atau%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (atau2))
          CLSMJ32_struc(n)%btau%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (btau2))
       else
          CLSMJ32_struc(n)%atau%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (atau))
          CLSMJ32_struc(n)%btau%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (btau))          
       endif
       CLSMJ32_struc(n)%psisat%value   (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (psis))
       CLSMJ32_struc(n)%bexp%value     (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (bee ))
       CLSMJ32_struc(n)%wpwet%value    (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (wpwet))
       CLSMJ32_struc(n)%ksat%value     (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (Ks))
       CLSMJ32_struc(n)%bdrckdpth%value(:,:,1) = LISv2g (glpnc,glpnr,G52LIS (soildepth))
       CLSMJ32_struc(n)%porosity%value (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (poros))
       CLSMJ32_struc(n)%gnu%value      (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (gnu))
       CLSMJ32_struc(n)%z2%value       (:,:,1) = LISv2g (glpnc,glpnr,G52LIS (z2ch))
       
       call populate_param_attribs( "ALBNIRDIFF", &
           "Alb near-IR diffuse scale factor", "-",  &
            CLSMJ32_struc(n)%albnirdir, &
            CLSMJ32_struc(n)%albnirdif )
       
       call populate_param_attribs( "ALBVISDIFF", &
            "Alb visible diffuse scale factor", "-",  &
            CLSMJ32_struc(n)%albvisdir, &
            CLSMJ32_struc(n)%albvisdif )

    enddo

!    call ESMF_ConfigGetAttribute(LDT_config,albInterval, label = "Albedo climatology interval:",rc=rc) ; VERIFY_(RC)
!    call ESMF_ConfigGetAttribute(LDT_config,source, label = "Albedo data source:", rc=rc)              ; VERIFY_(RC)
!    call ESMF_ConfigGetAttribute(LDT_config,proj,label = "Albedo map projection:", rc=rc)              ; VERIFY_(RC)
    
!    do n=1,LDT_rc%nnest
!       write(LDT_logunit,*) 'Reading CLSM ALBEDO : '//trim(source)
!       call bcr%readDataset (n, SOURCE, albinterval, proj, &
!            CLSMJ32_struc(n)%albvisdif%value, CLSMJ32_struc(n)%albnirdif%value, &
!            .true., maskarray=LDT_LSMparam_struc(n)%landmask%value(:,:,n))
!       CLSMJ32_struc(n)%albnirdir%value = CLSMJ32_struc(n)%albnirdif%value
!       CLSMJ32_struc(n)%albvisdir%value = CLSMJ32_struc(n)%albvisdif%value
! 
!       LDT_rc%monthlyData(n) = .true.
!       LDT_albedo_struc(n)%albInterval = "monthly"
!       LDT_gfrac_struc(n)%gfracInterval = "monthly"
!    enddo

  end subroutine catchmentParms_init_J32

end module CLSMJ32_parmsMod
