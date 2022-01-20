!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LIS_irrigationMod
!BOP
!
! !MODULE: LIS_irrigationMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  11 Nov 2012: Sujay Kumar; Initial implementation
!  29 May 2019; Jessica Erlingis; Incorporate Wanshu Nie's max/min GVF update
!  10 Dec 2020: Hiroko Beaudoing; Incorporate crop calendar and concurrent
!                                 irrigation schemes
!                                 Made irrig_type_dec public (was private)
!
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod

  implicit none
  
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_irrigation_init
  public :: LIS_irrigation_run
  public :: LIS_irrigation_output
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_irrig_state      !data structure containing irrigation states
  public :: LIS_irrig_struc      !data structure containing irrigation variables
!EOP  

  type, public :: irrig_type_dec
     real               :: outInterval
     character*100      :: models_used
     logical            :: stats_file_open
     character*50       :: cropcalendar
     real               :: veg_thresh       
     real               :: sprinkler_start  !sprinkler start time
     real               :: sprinkler_duration   !sprinkler duration
     real               :: sprinkler_thresh  !sprinkler threshold
     real               :: sprinkler_efcor   !sprinkler efficiency
     real               :: drip_start  !drip start time
     real               :: drip_duration   !drip duration
     real               :: drip_thresh  !drip threshold
     real               :: drip_efcor   !drip efficiency
     real               :: flood_start  !flood start time
     real               :: flood_duration   !flood duration
     real               :: flood_thresh  !flood threshold
     real               :: flood_efcor   !flood efficency
     integer            :: cropseasons
     real,allocatable   :: plantDay(:,:)
     real,allocatable   :: harvestDay(:,:)
     real               :: irrigation_thresh !BZ
     real               :: irrigation_mxsoildpth
     real               :: irrigation_GVFparam1   
     real               :: irrigation_GVFparam2   
     integer            :: irrigation_SourcePartition  
     integer            :: irrigation_GWabstraction
     integer            :: irrigation_dveg
     integer            :: irrigation_schedule
  end type irrig_type_dec

  type(irrig_type_dec),allocatable :: LIS_irrig_struc(:)

  type(ESMF_State),    allocatable :: LIS_irrig_state(:)

contains

!BOP
! 
! !ROUTINE: LIS_irrigation_init
! \label{LIS_irrigation_init}
! 
! !DESCRIPTION:
!
! Allocates memory for data structures used for reading 
! irrigation datasets. The irrigationdepth field is updated by the external
! files. The irrigation water equivalent fields are expected to be set
! by the model. 
! 
! !INTERFACE:
  subroutine LIS_irrigation_init

! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif
   use LIS_timeMgrMod,   only : LIS_registerAlarm, LIS_parseTimeString
!EOP
    integer       :: n
    integer       :: status
    integer       :: rc
    integer       :: ios, nid
    character*100 :: temp
    character*10  :: time
    character*1   :: nestid(2)
    logical       :: file_exists
    
! ___________________________________________________

 !- Read in Config file irrigation inputs:

    call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%irrigation_type,label='Irrigation scheme: ',DEFAULT="none", rc=rc)
    write(LIS_logunit,*) "[INFO] Irrigation scheme : ",LIS_rc%irrigation_type

    if (LIS_rc%irrigation_type .ne. "none") then
       
       allocate(LIS_irrig_state(LIS_rc%nnest))
       allocate(LIS_irrig_struc(LIS_rc%nnest))
              
       do n=1,LIS_rc%nnest

          ! Read in type of irrigation scheme selected (CONCURRENT, SPRINKLER, DRIP,FLOOD):
          call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_type,   &
               label="Irrigation scheme:",default="none",rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation scheme selected:  ",       &
               trim(LIS_rc%irrigation_type)

          ! HKB--added irrigation type specific configurations 
          ! Set trigger check start time [local hour] and duration in lis.config
          ! sprinkler
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%sprinkler_start,     &
               label="Irrigation Sprinkler start time:",default=6.,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation Sprinkler start time: ",                &
               LIS_irrig_struc(n)%sprinkler_start
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%sprinkler_duration,  &
               label="Irrigation Sprinkler duration:",default=4.,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation Sprinkler duration:   ",                &
               LIS_irrig_struc(n)%sprinkler_duration
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%sprinkler_thresh,    &
               label="Irrigation threshold for Sprinkler:",default=0.5,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation threshold for Sprinkler:  ",            &
               LIS_irrig_struc(n)%sprinkler_thresh
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%sprinkler_efcor,     &
               label="Irrigation Sprinkler efficiency:",default=1.,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation Sprinkler efficiency:  ",               &
               LIS_irrig_struc(n)%sprinkler_efcor
          ! drip
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%drip_start,          &
               label="Irrigation Drip start time:",default=6.,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation Drip start time:  ",                    &
               LIS_irrig_struc(n)%drip_start
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%drip_duration,       &
               label="Irrigation Drip duration:",default=1.,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation Drip duration:  ",                      &
               LIS_irrig_struc(n)%drip_duration
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%drip_thresh,         &
               label="Irrigation threshold for Drip:",default=0.5,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation threshold for Drip:  ",                 &
               LIS_irrig_struc(n)%drip_thresh
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%drip_efcor,          &
               label="Irrigation Drip efficiency:",default=1.,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation Drip efficiency:  ",                    &
               LIS_irrig_struc(n)%drip_efcor
          ! flood
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%flood_start,         &
               label="Irrigation Flood start time:",default=9.,rc=rc) 
          write(LIS_logunit,*) "[INFO] Irrigation Flood start time:  ",                   &
               LIS_irrig_struc(n)%flood_start
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%flood_duration,      &
               label="Irrigation Flood duration:",default=18.,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation Flood duration:  ",                     &
               LIS_irrig_struc(n)%flood_duration
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%flood_thresh,        &
               label="Irrigation threshold for Flood:",default=0.5,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation threshold for Flood:  ",                &
               LIS_irrig_struc(n)%flood_thresh
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%flood_efcor,         &
               label="Irrigation Flood efficiency:",default=1.,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation Flood efficiency:  ",                   &
               LIS_irrig_struc(n)%flood_efcor
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%irrigation_schedule, &
               label="IRRIGATION ON SCHEDULE:",default=0,rc=rc)
          write(LIS_logunit,*) "[INFO] IRRIGATION ON SCHEDULE: ",                         &
               LIS_irrig_struc(n)%irrigation_schedule

          ! Dynamic vegetation trigger parameters
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%veg_thresh,          &
               label="THRESHOLD OF DYNAMIC VEGETATION RANGE:", default=1., rc=rc)
          write(LIS_logunit,*) "THRESHOLD OF DYNAMIC VEGETATION RANGE:  ",                &
               LIS_irrig_struc(n)%veg_thresh 
          ! Parameters to control the GVF threshold based on the range of GVF
          ! (shdmax-shdmin) for which sprinkler irrigation is triggered for NoahMP:(WN)
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%irrigation_GVFparam1,&
               label="Irrigation GVF parameter 1:",default = 0., rc=rc)
          write(LIS_logunit,*) "irrigation GVF parameter 1:  ",                           & 
               LIS_irrig_struc(n)%irrigation_GVFparam1
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%irrigation_GVFparam2,&
               label="Irrigation GVF parameter 2:",default = 0.4, rc=rc)
          write(LIS_logunit,*) "irrigation GVF parameter 2:  ",                           &
               LIS_irrig_struc(n)%irrigation_GVFparam2

          ! Max. soil layer depth for irrigation to reach to (available for flood only):
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%irrigation_mxsoildpth,&
               label="Irrigation max soil layer depth:", default=1., rc=rc)
          write(LIS_logunit,*) "Irrigation max soil layer depth:  ",                       &
               LIS_irrig_struc(n)%irrigation_mxsoildpth

          ! Crop calendar options
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%cropcalendar,         &
               label="CROP CALENDAR USE:", default="none",rc=rc)
          write(LIS_logunit,*) "CROP CALENDAR USE:  ",                                     &
               LIS_irrig_struc(n)%cropcalendar
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%cropseasons,          &
               label="Crop seasons:", default=2,rc=rc)
          write(LIS_logunit,*) "Crop seasons:  ",                                          &
               LIS_irrig_struc(n)%cropseasons

          if ( LIS_irrig_struc(n)%cropcalendar .ne. "none" ) then
             allocate(LIS_irrig_struc(n)%plantDay( &
                  LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_irrig_struc(n)%cropseasons))
             allocate(LIS_irrig_struc(n)%harvestDay( &
                  LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_irrig_struc(n)%cropseasons))
             LIS_irrig_struc(n)%plantDay   = 0.0
             LIS_irrig_struc(n)%harvestDay = 0.0
          endif

          ! JE Remove irrigated water from groundwater
          ! Need to add model sanity check here to make sure model contains GW (?)
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%irrigation_GWabstraction,&
               label="Groundwater abstraction for irrigation:",default=0,rc=rc)
          write(LIS_logunit,*) "[INFO] Groundwater abstraction for irrigation:  ",            &
               LIS_irrig_struc(n)%irrigation_GWabstraction
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%irrigation_SourcePartition, &
               label="Irrigation source water partition:",default=0,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation source water partition:  ",                 &
               LIS_irrig_struc(n)%irrigation_SourcePartition
          call ESMF_ConfigGetAttribute(LIS_config,LIS_irrig_struc(n)%irrigation_dveg,         &
               label="Irrigation scheduling based on dynamic vegetation:",default=0,rc=rc)
          write(LIS_logunit,*) "[INFO] Irrigation scheduling based on dynamic vegetation:  ", &
               LIS_irrig_struc(n)%irrigation_dveg
          
       enddo
       
       ! Frequency with which irrigation field is written out:
       call ESMF_ConfigGetAttribute(LIS_config,time,&
            label="Irrigation output interval:",rc=rc)
       call LIS_verify(rc,"Irrigation output interval: not defined")
       write(LIS_logunit,*) "[INFO] Irrigation output interval:  ",time
       
       ! Register irrigation output interval:
       do n=1,LIS_rc%nnest
          call LIS_parseTimeString(time,LIS_irrig_struc(n)%outInterval)
          call LIS_registerAlarm("LIS irrigation output interval",&
               real(LIS_irrig_struc(n)%outInterval), &
               LIS_irrig_struc(n)%outInterval)
          LIS_irrig_struc(n)%models_used = trim(LIS_rc%irrigation_type)
          LIS_irrig_struc(n)%stats_file_open = .true.
       enddo

       do n=1,LIS_rc%nnest
          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid
          
          LIS_irrig_state(n) = ESMF_StateCreate(name="LSM Irrigation State"//&
               nestid(1)//nestid(2), rc=status)
          call LIS_verify(status, &
               "ESMF_StateCreate failed in LIS_irrigation_init")
       enddo

    !- Initiate the irrigation scheme selected in lis.config file:
       call irrigationschemeinit(trim(LIS_rc%irrigation_type)//char(0),&
            LIS_irrig_state)

    endif

  end subroutine LIS_irrigation_init

!BOP
! 
! !ROUTINE: LIS_irrigation_run
! \label{LIS_irrigation_run}
! 
! !INTERFACE:
  subroutine LIS_irrigation_run(n)
! !USES: 
    implicit none

! !ARGUMENTS: 
    integer  :: n 

! !DESCRIPTION:
! This routine runs the specified irrigation model.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP

    if(LIS_rc%irrigation_type.ne."none") then     

       call getirrigationlsmstates(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%irrigation_type)//char(0), n,LIS_irrig_state(n))
       call applyirrigationupdates(trim(LIS_rc%irrigation_type)//char(0),&
            n,LIS_irrig_state(n))
       
    endif

  end subroutine LIS_irrigation_run

!BOP
! 
! !ROUTINE: LIS_irrigation_output
! \label{LIS_irrigation_output}
! 
! !INTERFACE:
  subroutine LIS_irrigation_output(n)
! !USES: 

    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_historyMod, only : LIS_writeModelOutput
    use LIS_fileIOMod,  only : LIS_create_output_directory, &
         LIS_create_output_filename,  &
         LIS_create_stats_filename

! !ARGUMENTS: 
    integer, intent(in)   :: n 

! !DESCRIPTION:
! This routine writes the irrigation model output.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP
    
    logical           :: alarmCheck,open_stats
    character*100     :: outfile, statsfile

    if(LIS_rc%irrigation_type.ne."none") then 
       alarmCheck = LIS_isAlarmRinging(LIS_rc,&
            "LIS irrigation output interval")
       if(alarmCheck) then 
          open_stats = .false. 
          if(LIS_rc%wopt.ne."none") then 
             if(LIS_masterproc) then 
                call LIS_create_output_directory('IRRIGATION')
                if (LIS_irrig_struc(n)%stats_file_open) then
                   call LIS_create_stats_filename(n,statsfile,"IRRIGATION")
                   LIS_irrig_struc(n)%stats_file_open = .false.
                   open_stats = .true.
                endif
             endif

             call LIS_create_output_filename(n,outfile,&
                  model_name ="IRRIGATION")

             call LIS_writeModelOutput(n,outfile,statsfile,              &
                  open_stats,outInterval=LIS_irrig_struc(n)%outInterval, &
                  nsoillayers=1, lyrthk = (/1.0/),                       &
                  nsoillayers2=1,                                        &
                  model_name=LIS_irrig_struc(n)%models_used,group=4)
          endif
       endif
    endif
    
  end subroutine LIS_irrigation_output

end module LIS_irrigationMod
