!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file xDE_sampler.F90
!! This application takes care of doing 1D and 2D parameter space sampling enforcing
!! EFTCAMB viability conditions, as specified by a parameter file.
!! Part of this file is taken from the inidriver file of CAMB.

! hardcoded options:
module hardcoded_options

    use precision

    real(dl), parameter :: log_sampling_value = 1.d-4 ! the value at which the sampler switches from linear to log sampling
    real(dl), parameter :: log_sampling_min   = 1.d-8 ! the assumed minimum log sampling value

end module hardcoded_options

! the program itself:
program IC_sampler

    use IniFile
    use CAMB
    use LambdaGeneral
    use Lensing
    use AMLUtils
    use Transfer
    use constants
    use Bispectrum
    use CAMBmain
    use omp_lib
    use EFTCAMB_stability
    use NonLinear

    !> use the sampling xDE class
    !use EFT_sampler
    use EFTCAMB_cache
    ! use EFTCAMB_cache
#ifdef NAGF95
    use F90_UNIX
#endif
    use hardcoded_options




    implicit none

    Type(CAMBparams) P

    character(LEN=Ini_max_string_len) numstr, VectorFileName, &
        InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
        LensedTotFileName, LensPotentialFileName,ScalarCovFileName, benchmark_buffer
    integer i,j,sample_points, benchmark_count
    character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
        MatterPowerFileNames(max_transfer_redshifts), outroot, version_check
    character(LEN=8) name_1, name_2
    real(dl) output_factor, nmassive

    !> sources parameters
    character(LEN=Ini_max_string_len) S,TransferClFileNames(max_transfer_redshifts)
    Type (TRedWin), pointer :: RedWin

    !> this is part of the sampler parameters.
    real(dl)  t2, param_1_max, param_1_min, param_2_max, param_2_min
    real(dl),allocatable :: t1(:)
    integer param_number, iter, nu_i
    logical do_log_sampling_param_1, do_log_sampling_param_2, success, error
    real(dl) k_max, astart, aend

#ifdef WRITE_FITS
    character(LEN=Ini_max_string_len) FITSfilename
#endif


logical :: DoCounts = .false.

    !> IC MOD START: adding the parameters for the sampling of xDE
    integer, parameter :: n_samples = 10000
    !integer :: n_ic = 100   !< extra boost of IC
    integer :: num_success=0
    integer :: i_sample
    integer :: model_params_num
    real(dl), dimension(:), allocatable :: random_params
    real(dl) :: param1_max, param1_min
    real(dl) :: param2_max, param2_min
    real(dl) :: param3_max, param3_min
    real(dl) :: param4_max, param4_min
    real(dl) :: param5_max, param5_min
    real(dl) :: param6_max, param6_min
    real(dl) :: param7_max, param7_min
    !type(EFTSamplingParameters) :: EFT_sampling_params
    !> output files
    real(dl), dimension(100)    :: output_scale_factor, output_xDE_sample, output_EFTOmegaV, output_phi
    real(dl)                    :: a_now
    real(dl)                    :: delta_omega_i
    !integer                     :: i_ic, n_init_cond

    !> error handling stat
    integer         :: xDE_iostat
    character(256)  :: xDE_iomsg
    integer :: i_num_par

    Type( EFTCAMB_timestep_cache ) :: eft_cache_output
    !> xDE MOD END

    logical bad

    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') stop 'No parameter input file'

    call Ini_Open(InputFile, 1, bad, .false.)
    if (bad) stop 'Error opening parameter file'

    Ini_fail_on_not_found = .false.

    outroot = Ini_Read_String('output_root')
    if (outroot /= '') outroot = trim(outroot) // '_'

    highL_unlensed_cl_template = Ini_Read_String_Default('highL_unlensed_cl_template',highL_unlensed_cl_template)

    call CAMB_SetDefParams(P)

    P%WantScalars = Ini_Read_Logical('get_scalar_cls')
    P%WantVectors = Ini_Read_Logical('get_vector_cls',.false.)
    P%WantTensors = Ini_Read_Logical('get_tensor_cls',.false.)

    !> sources add-on
    P%Want_CMB =  Ini_Read_Logical('want_CMB',.true.)
    P%Want_CMB_lensing =  P%Want_CMB .or. Ini_Read_Logical('want_CMB_lensing',.true.)

    if (P%WantScalars) then
        num_redshiftwindows = Ini_Read_Int('num_redshiftwindows',0)
    else
        num_redshiftwindows = 0
    end if
    limber_windows = Ini_Read_Logical('limber_windows',limber_windows)
    if (limber_windows) limber_phiphi = Ini_Read_Int('limber_phiphi',limber_phiphi)
    if (num_redshiftwindows>0) then
        DoRedshiftLensing = Ini_Read_Logical('DoRedshiftLensing',.false.)
        Kmax_Boost = Ini_Read_Double('Kmax_Boost',Kmax_Boost)
    end if
    Do21cm = Ini_Read_Logical('Do21cm', .false.)
    !> AZ TEST: START
    !write(*,*) "Do21cm", Do21cm
    !> AZ TEST: END
    num_extra_redshiftwindows = 0
    do i=1, num_redshiftwindows
        RedWin => Redshift_w(i)
        call InitRedshiftWindow(RedWin)
        write (numstr,*) i
        numstr=adjustl(numstr)
        RedWin%Redshift = Ini_Read_Double('redshift('//trim(numstr)//')')
        S = Ini_Read_String('redshift_kind('//trim(numstr)//')')
        if (S=='21cm') then
            RedWin%kind = window_21cm
        elseif (S=='counts') then
            RedWin%kind = window_counts
        elseif (S=='lensing') then
            RedWin%kind = window_lensing
        else
            write (*,*) i, 'Error: unknown type of window '//trim(S)
            stop
        end if
        RedWin%a = 1/(1+RedWin%Redshift)
        if (RedWin%kind /= window_21cm) then
            RedWin%sigma = Ini_Read_Double('redshift_sigma('//trim(numstr)//')')
            RedWin%sigma_z = RedWin%sigma
        else
            Do21cm = .true.
            RedWin%sigma = Ini_Read_Double('redshift_sigma_Mhz('//trim(numstr)//')')
            if (RedWin%sigma < 0.003) then
                write(*,*) 'WARNING:Window very narrow.'
                write(*,*) ' --> use transfer functions and transfer_21cm_cl =T ?'
            end if
            !with 21cm widths are in Mhz, make dimensionless scale factor
            RedWin%sigma = RedWin%sigma/(f_21cm/1e6)
            RedWin%sigma_z = RedWin%sigma*(1+RedWin%RedShift)**2
            write(*,*) i,'delta_z = ', RedWin%sigma_z
        end if
        if (RedWin%kind == window_counts) then
            DoCounts = .true.
            RedWin%bias = Ini_Read_Double('redshift_bias('//trim(numstr)//')')
            RedWin%dlog10Ndm = Ini_Read_Double('redshift_dlog10Ndm('//trim(numstr)//')',0.d0)
            if (DoRedshiftLensing) then
                num_extra_redshiftwindows=num_extra_redshiftwindows+1
                RedWin%mag_index = num_extra_redshiftwindows
            end if
        end if
    end do


    if (Do21cm) then
        line_basic = Ini_Read_Logical('line_basic')
        line_distortions = Ini_read_Logical('line_distortions')
        line_extra = Ini_Read_Logical('line_extra')

        line_phot_dipole = Ini_read_Logical('line_phot_dipole')
        line_phot_quadrupole = Ini_Read_Logical('line_phot_quadrupole')
        line_reionization = Ini_Read_Logical('line_reionization')

        use_mK = Ini_read_Logical('use_mK')
        if (DebugMsgs) then
            write (*,*) 'Doing 21cm'
            write (*,*) 'dipole = ',line_phot_dipole, ' quadrupole =', line_phot_quadrupole
        end if
    else
        line_extra = .false.
    end if

    if (DoCounts) then
        counts_density = Ini_read_Logical('counts_density')
        counts_redshift = Ini_read_Logical('counts_redshift')
        counts_radial = Ini_read_Logical('counts_radial')
        counts_evolve = Ini_read_Logical('counts_evolve')
        counts_timedelay = Ini_read_Logical('counts_timedelay')
        counts_ISW = Ini_read_Logical('counts_ISW')
        counts_potential = Ini_read_Logical('counts_potential')
        counts_velocity = Ini_read_Logical('counts_velocity')
        !> AZ MOD: START
        write(*,*) 'counts_density  :', counts_density
        write(*,*) 'counts_redshift :', counts_redshift
        write(*,*) 'counts_radial   :', counts_radial
        write(*,*) 'counts_evolve   :', counts_evolve
        write(*,*) 'counts_timedelay:', counts_timedelay
        write(*,*) 'counts_ISW      :', counts_ISW
        write(*,*) 'counts_potential:', counts_potential
        write(*,*) 'counts_velocity :', counts_velocity
        !> AZ MOD: END
    end if



    P%OutputNormalization=outNone
    output_factor = Ini_Read_Double('CMB_outputscale',1.d0)

    P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

    P%PK_WantTransfer=Ini_Read_Logical('get_transfer')

    AccuracyBoost  = Ini_Read_Double('accuracy_boost',AccuracyBoost)
    lAccuracyBoost = Ini_Read_Real('l_accuracy_boost',lAccuracyBoost)
    HighAccuracyDefault = Ini_Read_Logical('high_accuracy_default',HighAccuracyDefault)

    P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)

    P%DoLensing = .false.
    if (P%WantCls) then
        if (P%WantScalars  .or. P%WantVectors) then
            P%Max_l = Ini_Read_Int('l_max_scalar')
            P%Max_eta_k = Ini_Read_Double('k_eta_max_scalar',P%Max_l*2._dl)
            if (P%WantScalars) then
                P%DoLensing = Ini_Read_Logical('do_lensing',.false.)
                if (P%DoLensing) lensing_method = Ini_Read_Int('lensing_method',1)
            end if
            if (P%WantVectors) then
                if (P%WantScalars .or. P%WantTensors) stop 'Must generate vector modes on their own'
                i = Ini_Read_Int('vector_mode')
                if (i==0) then
                    vec_sig0 = 1
                    Magnetic = 0
                else if (i==1) then
                    Magnetic = -1
                    vec_sig0 = 0
                else
                    stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
                end if
            end if
        end if

        if (P%WantTensors) then
            P%Max_l_tensor = Ini_Read_Int('l_max_tensor')
            P%Max_eta_k_tensor =  Ini_Read_Double('k_eta_max_tensor',Max(500._dl,P%Max_l_tensor*2._dl))
        end if
    endif

    !  Read initial parameters.

    call DarkEnergy_ReadParams(DefIni)

    P%h0     = Ini_Read_Double('hubble')
    P%eft_par_cache%h0  =P%h0
    if (Ini_Read_Logical('use_physical',.false.)) then
        P%omegab = Ini_Read_Double('ombh2')/(P%H0/100)**2
        P%eft_par_cache%omegab =P%omegab
        P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
        P%eft_par_cache%omegac=P%omegac
        P%omegan = Ini_Read_Double('omnuh2')/(P%H0/100)**2
        P%eft_par_cache%omegan= P%omegan
        P%omegav = 1- Ini_Read_Double('omk') - P%omegab-P%omegac - P%omegan
        P%eft_par_cache%omegav=P%omegav
    else
        P%omegab = Ini_Read_Double('omega_baryon')
        P%omegac = Ini_Read_Double('omega_cdm')
        P%omegav = Ini_Read_Double('omega_lambda')
        P%omegan = Ini_Read_Double('omega_neutrino')
    end if

    P%tcmb   = Ini_Read_Double('temp_cmb',COBE_CMBTemp)
    P%yhe    = Ini_Read_Double('helium_fraction',0.24_dl)
    P%Num_Nu_massless  = Ini_Read_Double('massless_neutrinos')

    ! EFTCAMB MOD START: read and allocate EFTCAMB
    ! read the EFTCAMB model selection flags:
    call P%EFTCAMB%EFTCAMB_init_from_file( DefIni )
    if ( P%EFTCAMB%EFTFlag /= 0 ) then
        ! print the EFTCAMB header:
        call P%EFTCAMB%EFTCAMB_print_header()
        ! initialize the output root name:
        P%EFTCAMB%outroot = TRIM( outroot )
        ! initialize the model from file:
        call P%EFTCAMB%EFTCAMB_init_model_from_file( DefIni )
        ! print feedback:
        call P%EFTCAMB%EFTCAMB_print_model_feedback( print_params=.False. )
    end if
    ! EFTCAMB MOD END.

    P%Nu_mass_eigenstates = Ini_Read_Int('nu_mass_eigenstates',1)
    if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'

    numstr = Ini_Read_String('massive_neutrinos')
    read(numstr, *) nmassive
    if (abs(nmassive-nint(nmassive))>1e-6) stop 'massive_neutrinos should now be integer (or integer array)'
    read(numstr,*, end=100, err=100) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
    P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

    if (P%Num_Nu_massive>0) then
        P%share_delta_neff = Ini_Read_Logical('share_delta_neff', .true.)
        numstr = Ini_Read_String('nu_mass_degeneracies')
        if (P%share_delta_neff) then
            if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
        else
            if (numstr=='') stop 'must give degeneracies for each eigenstate if share_delta_neff=F'
            read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
        end if
        numstr = Ini_Read_String('nu_mass_fractions')
        if (numstr=='') then
            if (P%Nu_mass_eigenstates >1) stop 'must give nu_mass_fractions for the eigenstates'
            P%Nu_mass_fractions(1)=1
        else
            read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
        end if
    end if

    !JD 08/13 begin changes for nonlinear lensing of CMB + LSS compatibility
    !P%Transfer%redshifts -> P%Transfer%PK_redshifts and P%Transfer%num_redshifts -> P%Transfer%PK_num_redshifts
    !in the P%WantTransfer loop.
    if (((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) &
        .or. P%PK_WantTransfer) then
        P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
    else
        P%transfer%high_precision = .false.
    endif

    if (P%PK_WantTransfer)  then
        P%WantTransfer  = .true.
        P%transfer%kmax          =  Ini_Read_Double('transfer_kmax')
        P%transfer%k_per_logint  =  Ini_Read_Int('transfer_k_per_logint')
        P%transfer%PK_num_redshifts =  Ini_Read_Int('transfer_num_redshifts')

        transfer_interp_matterpower = Ini_Read_Logical('transfer_interp_matterpower ', transfer_interp_matterpower)
        transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)
        if (P%transfer%PK_num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
        do i=1, P%transfer%PK_num_redshifts
            P%transfer%PK_redshifts(i)  = Ini_Read_Double_Array('transfer_redshift',i,0._dl)
            transferFileNames(i)     = Ini_Read_String_Array('transfer_filename',i)
            MatterPowerFilenames(i)  = Ini_Read_String_Array('transfer_matterpower',i)
            if (TransferFileNames(i) == '') then
                TransferFileNames(i) =  trim(numcat('transfer_',i))//'.dat'
            end if
            if (MatterPowerFilenames(i) == '') then
                MatterPowerFilenames(i) =  trim(numcat('matterpower_',i))//'.dat'
            end if
            if (TransferFileNames(i)/= '') &
                TransferFileNames(i) = trim(outroot)//TransferFileNames(i)
            if (MatterPowerFilenames(i) /= '') &
                MatterPowerFilenames(i)=trim(outroot)//MatterPowerFilenames(i)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts = 0
    end if

    if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) then
        P%WantTransfer  = .true.
        call Transfer_SetForNonlinearLensing(P%Transfer)
    end if

    call Transfer_SortAndIndexRedshifts(P%Transfer)
    !JD 08/13 end changes

    P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)

    Ini_fail_on_not_found = .false.

    DebugParam = Ini_Read_Double('DebugParam',DebugParam)
    ALens = Ini_Read_Double('Alens',Alens)

    call Reionization_ReadParams(P%Reion, DefIni)
    call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors)
    call Recombination_ReadParams(P%Recomb, DefIni)
    if (Ini_HasKey('recombination')) then
        i = Ini_Read_Int('recombination',1)
        if (i/=1) stop 'recombination option deprecated'
    end if

    call Bispectrum_ReadParams(BispectrumParams, DefIni, outroot)

    if (P%WantScalars .or. P%WantTransfer) then
        P%Scalar_initial_condition = Ini_Read_Int('initial_condition',initial_adiabatic)
        if (P%Scalar_initial_condition == initial_vector) then
            P%InitialConditionVector=0
            numstr = Ini_Read_String('initial_vector',.true.)
            read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
        end if
        if (P%Scalar_initial_condition/= initial_adiabatic) use_spline_template = .false.
    end if

    if (P%WantScalars) then
        ScalarFileName = trim(outroot)//Ini_Read_String('scalar_output_file')
        LensedFileName =  trim(outroot) //Ini_Read_String('lensed_output_file')
        LensPotentialFileName =  Ini_Read_String('lens_potential_output_file')
        if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
        ScalarCovFileName =  Ini_Read_String_Default('scalar_covariance_output_file','scalCovCls.dat',.false.)
        if (ScalarCovFileName/='') then
            has_cl_2D_array = .true.
            ScalarCovFileName = concat(outroot,ScalarCovFileName)
        end if
    end if
    if (P%WantTensors) then
        TensorFileName =  trim(outroot) //Ini_Read_String('tensor_output_file')
        if (P%WantScalars)  then
            TotalFileName =  trim(outroot) //Ini_Read_String('total_output_file')
            LensedTotFileName = Ini_Read_String('lensed_total_output_file')
            if (LensedTotFileName/='') LensedTotFileName= trim(outroot) //trim(LensedTotFileName)
        end if
    end if
    if (P%WantVectors) then
        VectorFileName =  trim(outroot) //Ini_Read_String('vector_output_file')
    end if

#ifdef WRITE_FITS
    if (P%WantCls) then
        FITSfilename =  trim(outroot) //Ini_Read_String('FITS_filename',.true.)
        if (FITSfilename /='') then
            inquire(file=FITSfilename, exist=bad)
            if (bad) then
                open(unit=18,file=FITSfilename,status='old')
                close(18,status='delete')
            end if
        end if
    end if
#endif


    Ini_fail_on_not_found = .false.

    !optional parameters controlling the computation

    P%AccuratePolarization = Ini_Read_Logical('accurate_polarization',.true.)
    P%AccurateReionization = Ini_Read_Logical('accurate_reionization',.false.)
    P%AccurateBB = Ini_Read_Logical('accurate_BB',.false.)
    P%DerivedParameters = Ini_Read_Logical('derived_parameters',.true.)

    version_check = Ini_Read_String('version_check')
    if (version_check == '') then
        !tag the output used parameters .ini file with the version of CAMB being used now
        call TNameValueList_Add(DefIni%ReadValues, 'version_check', version)
    else if (version_check /= version) then
        write(*,*) 'WARNING: version_check does not match this CAMB version'
    end if
    !Mess here to fix typo with backwards compatibility
    if (Ini_HasKey('do_late_rad_trunction')) then
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_trunction',.true.)
        if (Ini_HasKey('do_late_rad_truncation')) stop 'check do_late_rad_xxxx'
    else
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_truncation',.true.)
    end if

    if (HighAccuracyDefault) then
        DoTensorNeutrinos = .true.
    else
        DoTensorNeutrinos = Ini_Read_Logical('do_tensor_neutrinos',DoTensorNeutrinos )
    end if
    FeedbackLevel = Ini_Read_Int('feedback_level',FeedbackLevel)

    P%MassiveNuMethod  = Ini_Read_Int('massive_nu_approx',Nu_best)

    ThreadNum      = Ini_Read_Int('number_of_threads',ThreadNum)
    use_spline_template = Ini_Read_Logical('use_spline_template',use_spline_template)

    if (do_bispectrum) then
        lSampleBoost   = 50
    else
        lSampleBoost   = Ini_Read_Double('l_sample_boost',lSampleBoost)
    end if

    ! get the number of parameters:
    param_number = P%EFTCAMB%model%parameter_number

    !write(*,*) 'number of parameters:', param_number
    !pause

    !> xDE MOD START
    !> get the number of parameters in the model
    model_params_num = P%EFTCAMB%model%parameter_number
    allocate(random_params(model_params_num))

    write(*,*) 'number of parameters:', model_params_num
    !pause

    !> here I should also get the bounds for the params
    param1_max = Ini_Read_Double('param1_max', 1.d0)
    if (model_params_num .ge. 2) param2_max = Ini_Read_Double('param2_max', 1.d0)
    if (model_params_num .ge. 3) param3_max = Ini_Read_Double('param3_max', 1.d0)
    if (model_params_num .ge. 4) param4_max = Ini_Read_Double('param4_max', 1.d0)
    if (model_params_num .ge. 5) param5_max = Ini_Read_Double('param5_max', 1.d0)
    if (model_params_num .ge. 6) param6_max = Ini_Read_Double('param6_max', 1.d0)
    if (model_params_num .ge. 7) param7_max = Ini_Read_Double('param7_max', 1.d0)

    param1_min = Ini_Read_Double('param1_min', -1.d0)
    if (model_params_num .ge. 2) param2_min = Ini_Read_Double('param2_min', -1.d0)
    if (model_params_num .ge. 3) param3_min = Ini_Read_Double('param3_min', -1.d0)
    if (model_params_num .ge. 4) param4_min = Ini_Read_Double('param4_min', -1.d0)
    if (model_params_num .ge. 5) param5_min = Ini_Read_Double('param5_min', -1.d0)
    if (model_params_num .ge. 6) param6_min = Ini_Read_Double('param6_min', -1.d0)
    if (model_params_num .ge. 7) param7_min = Ini_Read_Double('param7_min', -1.d0)
    !> xDE MOD END

    call Ini_Close

    if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'

    FeedbackLevel = 2

    !Fill the cache
    grhom = 3*P%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)

    !grhom=3.3379d-11*h0*h0
    grhog = kappa/c**2*4*sigma_boltz/c**3*P%tcmb**4*Mpc**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
    ! grhog=1.4952d-13*tcmb**4
    grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*grhog !7/8*(4/11)^(4/3)*grhog (per neutrino species)
    !grhor=3.3957d-14*tcmb**4

    !correction for fractional number of neutrinos, e.g. 3.04 to give slightly higher T_nu hence rhor
    !for massive Nu_mass_degeneracies parameters account for heating from grhor

    grhornomass=grhor*P%Num_Nu_massless
    grhormass=0
    do nu_i = 1, P%Nu_mass_eigenstates
        grhormass(nu_i)=grhor*P%Nu_mass_degeneracies(nu_i)
    end do
    grhoc=grhom*P%omegac
    grhob=grhom*P%omegab
    grhov=grhom*P%omegav
    grhok=grhom*P%omegak
    !  adotrad gives the relation a(tau) in the radiation era:
    adotrad = sqrt((grhog+grhornomass+sum(grhormass(1:P%Nu_mass_eigenstates)))/3)

    Nnow = P%omegab*(1-P%yhe)*grhom*c**2/kappa/m_H/Mpc**2

    akthom = sigma_thomson*Nnow*Mpc
    !sigma_T * (number density of protons now)

    fHe = P%YHe/(mass_ratio_He_H*(1.d0-P%YHe))  !n_He_tot / n_H_tot

        ! EFTCAMB MOD START: clean up the EFTCAMB parameter cache
        if ( P%EFTCAMB%EFTFlag /= 0 ) then
            call P%eft_par_cache%initialize()
        end if
        ! EFTCAMB MOD END.

        call init_massive_nu(P%omegan /=0)

        ! EFTCAMB MOD START: initialize the EFTCAMB parameter choice
        if ( P%EFTCAMB%EFTFlag /= 0 ) then

            ! 1) parameter cache:
            !    - relative densities:
            P%eft_par_cache%omegac      = P%omegac
            P%eft_par_cache%omegab      = P%omegab
            P%eft_par_cache%omegav      = P%omegav
            P%eft_par_cache%omegak      = P%omegak
            P%eft_par_cache%omegan      = P%omegan
            P%eft_par_cache%omegag      = grhog/grhom
            P%eft_par_cache%omegar      = grhornomass/grhom
            !    - Hubble constant:
            P%eft_par_cache%h0          = P%h0
            P%eft_par_cache%h0_Mpc      = P%h0/c*1000._dl
            !    - densities:
            P%eft_par_cache%grhog       = grhog
            P%eft_par_cache%grhornomass = grhornomass
            P%eft_par_cache%grhoc       = grhoc
            P%eft_par_cache%grhob       = grhob
            P%eft_par_cache%grhov       = grhov
            P%eft_par_cache%grhok       = grhok
            !    - massive neutrinos:
            P%eft_par_cache%Num_Nu_Massive       = P%Num_Nu_Massive
            P%eft_par_cache%Nu_mass_eigenstates  = P%Nu_mass_eigenstates
            allocate( P%eft_par_cache%grhormass(max_nu), P%eft_par_cache%nu_masses(max_nu) )
            P%eft_par_cache%grhormass            = grhormass
            P%eft_par_cache%nu_masses            = nu_masses
            ! 2) now run background initialization:
            !call P%EFTCAMB%model%initialize_background( P%eft_par_cache, P%EFTCAMB%EFTCAMB_feedback_level, success )
            !if ( .not. success ) then
                ! global_error_flag         = 1
                ! global_error_message      = 'EFTCAMB: background solver failed'
                ! error = global_error_flag
                ! 5) final feedback:
                !if ( P%EFTCAMB%EFTCAMB_feedback_level > 1 ) then
                    !write(*,'(a)') '***************************************************************'
                !end if
                !call MpiStop('EFTCAMB: background solver failed')
            !end if
          end if

    !open(unit=1, name=trim(outroot) // 'Stability_Space.dat', action='write')
    ! do the sampling and save to file:
    !allocate(t1(param_number))
    !astart = 0.1_dl
    !aend = 1._dl
    !k_max = 10._dl

    !> xDE MOD START


    write(*,*) " Opening the files."
    !> opening the files for the plots

    open(unit=750, file = trim(outroot) // 'IC_cls_TT.dat',     status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining xDE_cls_TT.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=751, file = trim(outroot) // 'IC_cls_TxW1.dat',   status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining IC_cls_TxW1.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=217, file = trim(outroot) // 'IC_cls_TxW2.dat',  status='unknown', iostat=xDE_iostat, iomsg=xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining xDE_cls_TxW2.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=753, file = trim(outroot) // 'IC_cls_TxW3.dat',   status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg,recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining xDE_cls_TxW3.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=754, file = trim(outroot) // 'IC_cls_W1xW1.dat',  status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining xDE_cls_W1xW1.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=755, file = trim(outroot) // 'IC_cls_W2xW2.dat',  status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining xDE_cls_W2xW2.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=756, file = trim(outroot) // 'IC_cls_W3xW3.dat',  status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining xDE_cls_W3xW3.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=752, file = trim(outroot) // 'IC_delta_Om.dat',   status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining xDE_delta_Om.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=88, file = trim(outroot) // 'IC_EFTOmega.dat',   status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining xDE_EFTOmega.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if


    open(unit=98, file = trim(outroot) // 'scale_factor.dat',   status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining scale_factor.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    open(unit=22, file = trim(outroot) // 'init_cond.dat',      status = 'unknown', iostat = xDE_iostat, iomsg = xDE_iomsg, recl=99999 )
    if (xDE_iostat /= 0 ) then
        write(*,*) 'Opeining init_cond.dat failed with', xDE_iostat, xDE_iomsg
        !pause
    end if

    write(*,*) "files opened"


    !> building the scale factor array
    do i = 1, 100
        output_scale_factor(i) = 1.d-3 + (i-1) * (1.d0 - 1.d-3)/99
        write(98,*) output_scale_factor(i)
    end do

    !> fix a few things here
    astart = 0.01_dl
    aend = 1._dl
    k_max = 10._dl
    i_sample = 1

    write(*,*) "    Sampling the Initial Conditions"

    do while ( num_success < n_samples )! .and. i_sample < 10000000)

        write(*,*) 'Monte Carlo Step:',i_sample,'. Successfull reconstruction n.:', num_success, ' /', n_samples

        !> randomly drawing the parameters for the models (uniform distribution)
        call random_number(random_params)

        !> now adjust the parameters range
        !random_params(1) = param1_min+random_params(1)*(param1_max-param1_min)
        !if ( model_params_num .ge. 2 )  random_params(2) = param2_min+random_params(2)*(param2_max-param2_min)
        !if ( model_params_num .ge. 3 )  random_params(3) = param3_min+random_params(3)*(param3_max-param3_min)
        !if ( model_params_num .ge. 4 )  random_params(4) = param4_min+random_params(4)*(param4_max-param4_min)
        !if ( model_params_num .ge. 5 )  random_params(5) = param5_min+random_params(5)*(param5_max-param5_min)
        !if ( model_params_num .ge. 6 )  random_params(6) = param6_min+random_params(6)*(param6_max-param6_min)
        !if ( model_params_num .ge. 7 )  random_params(7) = param7_min+random_params(7)*(param7_max-param7_min)

        do i_num_par=1, model_params_num
            select case ( i_num_par )
                case (1)
                    random_params(1)=param1_min+random_params(1)*(param1_max-param1_min)
                case(2)
                    random_params(2) = param2_min+random_params(2)*(param2_max-param2_min)
                case(3)
                    random_params(3) = param3_min+random_params(3)*(param3_max-param3_min)
                case(4)
                    random_params(4) = param4_min+random_params(4)*(param4_max-param4_min)
                case(5)
                    random_params(5) = param5_min+random_params(5)*(param5_max-param5_min)
                case(6)
                    random_params(6) = param6_min+random_params(6)*(param6_max-param6_min)
                case(7)
                    random_params(7) = param7_min+random_params(7)*(param7_max-param7_min)
            end select

        end do

        !> then I need to initialize the model with these parameters
        !> intializing model parameters
        call P%EFTCAMB%model%init_model_parameters( random_params )

        !> and then solve the background
        success = .true.
        call P%EFTCAMB%model%initialize_background( P%eft_par_cache, P%EFTCAMB%EFTCAMB_feedback_level, success )

        !> if the theory is stable, call CAMB
        if (success) then

            write(*,*) "Model Parameters:", random_params

            success = .true.
            call EFTCAMB_Stability_Check( success, P%EFTCAMB, P%eft_par_cache, astart, aend, k_max )

            if(success) then

                call eft_cache_output%initialize
                call  P%EFTCAMB%model%compute_background_EFT_functions(1._dl, P%eft_par_cache, eft_cache_output )

                delta_omega_i = eft_cache_output%EFTOmegaV
                write(*,*) 'Delta Omega:',delta_omega_i

                if (abs(delta_omega_i) .le. 5.d-1) then

                write(*,*) 'theory stable!!! Woooo!'
                    num_success = num_success+1

                    !> call CAMB
                    write(*,*) 'Calling CAMB'
                    call CAMB_GetResults(P)

                    !> Re-computing Omega(a)
                    do i = 1, 100
                        a_now = output_scale_factor(i)
                        !> fill the EFT cache
                        call  P%EFTCAMB%model%compute_background_EFT_functions(a_now, P%eft_par_cache, eft_cache_output )
                        output_EFTOmegaV(i) = eft_cache_output%EFTOmegaV
                        if ( i == 100 .and. (abs(output_EFTOmegaV(i) - 1.d0) <1.d-10 )) write(*,*) '1+Omega(1) = 1!!!'
                    end do

                    write(*,*) 'Dumping data in files'
                    !> dump in files
                    write(*,*) '    dumping delta_omega'
                    write(752,*) delta_omega_i                   !< dumping Delta Omega(a)
                    write(*,*) '    dumping Omega(a)'
                    write(88,*) output_EFTOmegaV                !< dumping Omega(a)
                    write(*,*) '    dumping initial conditions'
                    write(22,*) random_params                   !< dumping initial conditions


                    !> dump cls in file
                    !> first write the CMB TT spectra
                    write(*,*) '    dumping Cls TT'
                    write(750,*) output_factor*Cl_scalar(:,1,1)

                    !> Then write the CMB-GNC cross correlations
                    write(*,*) '    dumping Cls TxW1'
                    write(751,*) sqrt(output_factor)*Cl_Scalar_Array(:,1,1,4)
                    write(*,*) '    dumping Cls TxW2'
                    write(217,*) sqrt(output_factor)*Cl_Scalar_Array(:,1,1,5)
                    write(*,*) '    dumping Cls TxW3'
                    write(753,*) sqrt(output_factor)*Cl_Scalar_Array(:,1,1,6)

                    !> Write the GNC self-correlations
                    write(*,*) '    dumping Cls W1xW1'
                    write(754,*) Cl_Scalar_Array(:,1,4,4)
                    write(*,*) '    dumping Cls W2xW2'
                    write(755,*) Cl_Scalar_Array(:,1,5,5)
                    write(*,*) '    dumping Cls W3xW3'
                    write(756,*) Cl_Scalar_Array(:,1,6,6)


                end if ! if delta omega < 0.5


            end if !< if theory is stable

        end if !< if the background is solved successfully


        i_sample = i_sample + 1

    end do !< end do on the IC sampling

    close(750)
    close(751)
    close(752)
    close(88)
    close(22)
    close(217)
    close(753)
    close(754)
    close(755)
    close(756)

!> xDE MOD END

!> xDE MOD END


    call CAMB_cleanup
    ! stop

100 stop 'Must give num_massive number of integer physical neutrinos for each eigenstate'

end program IC_sampler
