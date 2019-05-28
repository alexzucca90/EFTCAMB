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

!> @file 09p3_Designer_GBD.f90
!! This file contains the relevant code for designer GBD models.


!----------------------------------------------------------------------------------------
!> This module contains the relevant code for designer GBD models.

!> @author Bin Hu, Marco Raveri, Simone Peirone

!> @author Alex Zucca

module EFTCAMB_designer_GBD

    use precision
    use IniFile
    use AMLutils
    use equispaced_linear_interpolation_1D
    use EFT_def
    use EFTCAMB_rootfind
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_abstract_model_designer

    !> adding the interpolated function and the reconstructed dark energy
    use EFTCAMB_interpolated_function_1D
    use EFTCAMB_power_law_DE_parametrizations_1D
    use EFTCAMB_hyperbolic_tangent_parametrizations_1D


    !use EFT_sampler

    implicit none

    private

    public EFTCAMB_GBD_designer

    !----------------------------------------------------------------------------------------
    !> This is the designer GBD model. Inherits from the abstract designer model and has the
    !! freedom of defining the expansion history.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_GBD_designer

        ! theory parameters:
        real(dl) :: phi_ini, dphi_ini                                     !< The initial values of the scalar dof and its derivative
        real(dl) :: xi                                                    !< The value of the coupling constant xi


        ! the pure EFT functions model selection flags:
        integer  :: EFTxDE                                                !< Model selection flag for designer GBD w DE. (this will choose only the part 8)

        !> GBD coupling type
        integer :: coupling_type                                          !< which coupling is being used for the reconstruction

        ! the pure EFT functions:
        class( parametrized_function_1D ), allocatable ::  DesGBDxDE     !< The pure EFT function X := rhoDE(a)/rho_crit(a).


        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: EFTOmega       !< The interpolated function Omega (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTc           !< The interpolated funtcion c (and derivatives).

        ! some designer parameters:
        integer  :: designer_num_points = 100000                           !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(10._dl**(-8._dl))           !< log(a start)
        real(dl) :: x_final             = 0.0_dl                          !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBDesignerGBDReadModelSelectionFromFile   !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBDesignerGBDAllocateModelSelection       !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBDesignerGBDInitModelParameters          !< subroutine taht initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBDesignerGBDInitModelParametersFromFile  !< subroutine that reads the parameters of the model from file.
        !procedure :: init_model_parameters_for_sampling => EFTCAMBDesignerGBDInitModelParametersSampling  !< subroutine that initializes the model params for sampling xDE


        ! background solver:
        procedure :: initialize_background           => EFTCAMBDesignerGBDInitBackground               !< subroutine that initializes the background of designer GBD.
        procedure :: solve_designer_equations        => EFTCAMBDesignerGBDSolveDesignerEquations       !< subroutine that solves the designer GBD background equations.

        ! coupling function
        procedure :: omega_phi                       => EFTCAMBDesignerGBDCoupling                      !< function that computes omega(phi) and its derivatives

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBDesignerGBDComputeParametersNumber     !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBDesignerGBDFeedback                    !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBDesignerGBDParameterNames              !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBDesignerGBDParameterNamesLatex         !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBDesignerGBDParameterValues             !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBDesignerGBDBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBDesignerGBDSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBDesignerGBDComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBDesignerGBDComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBDesignerGBDComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBDesignerGBDAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_GBD_designer

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBDesignerGBDReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_GBD_designer)   :: self   !< the base class
        type(TIniFile)                  :: Ini     !< Input ini file

        ! read model selection flags:
        self%EFTxDE             = Ini_Read_Int_File( Ini, 'EFTxDE', 0 )

        !> read coupling type
        self%coupling_type      = Ini_Read_Int_File( Ini, 'GBD_coupling_type', 3)


    end subroutine EFTCAMBDesignerGBDReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBDesignerGBDAllocateModelSelection( self, Ini )

        implicit none

        class(EFTCAMB_GBD_designer)             :: self              !< the base class
        type(TIniFile)                          :: Ini               !< Input ini file
        character, allocatable, dimension(:)    :: param_names       !< an array of strings containing the names of the function parameters
        character, allocatable, dimension(:)    :: param_names_latex !< an array of strings containing the latex names of the function parameters

        ! allocate xDE:
        if ( allocated(self%DesGBDxDE) ) deallocate(self%DesGBDxDE)
        select case ( self%EFTxDE )
            case(0)
                allocate( constant_parametrization_1D::self%DesGBDxDE )
            case(1)
                allocate( power_law_DE_parametrization_1D::self%DesGBDxDE )
                call self%DesGBDxDE%set_param_names(['EFTxDE_wDE'], ['w_{\rm DE}'])
            case(2)
                allocate( hyperbolic_tangent_parametrization_1D::self%DesGBDxDE )
                call self%DesGBDxDE%set_param_names(['EFTxDE_A','EFTxDE_B', 'EFTxDE_C'], ['A','B','C'])
            case(3)
                !allocate( interpolated_function_1D::self%DesGBDxDE )
                !call self%DesGBDxDE%set_param_names(['xDE_filename  '])
                !> here need to modify this
                allocate( interpolated_function_1D::self%DesGBDxDE )
                call self%DesGBDxDE%set_param_names(['xDE_filename  ', 'xDE_num_bins  '])
            case default
                write(*,'(a,I3)') 'No model corresponding to EFTxDE =', self%EFTxDE
                write(*,'(a)')    'Choose EFTxDE < 4.'
        end select

        ! initialize the names:
        call self%DesGBDxDE%set_name( 'EFTx', 'X' )

        ! additional initialization of the function:
        call self%DesGBDxDE%init_func_from_file( Ini )

    end subroutine EFTCAMBDesignerGBDAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBDesignerGBDInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_GBD_designer)                          :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.
        real(dl), dimension(self%parameter_number -4)          :: temp
        integer                                                :: i

        self%phi_ini    = array(1)
        self%dphi_ini   = array(2)
        self%xi         = array(3)
        self%x_initial  = array(4)


        do i = 1, self%parameter_number -1
            temp(i) = array(i+4)
        end do
        call self%DesGBDxDE%init_parameters(temp)

    end subroutine EFTCAMBDesignerGBDInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBDesignerGBDInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_GBD_designer)   :: self   !< the base class
        type(TIniFile)                  :: Ini    !< Input ini file

        !> read phi_ini
        self%phi_ini = Ini_Read_Double_File( Ini, 'GBD_phi_ini', 0._dl)

        !> read dphi_ini
        self%dphi_ini = Ini_Read_Double_File( Ini, 'GBD_dphi_ini', 0._dl)

        !> read xi
        self%xi      = Ini_Read_Double_File( Ini, 'GBD_xi',      1._dl)

        !> read initial scale factor
        self%x_initial          = Ini_Read_Double_File( Ini, 'GBD_a_initial', 1.d-8)
        self%x_initial          = log(self%x_initial)


        !> read x_DE parameters:
        call self%DesGBDxDE%init_from_file( Ini )
    

    end subroutine EFTCAMBDesignerGBDInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes the coupling function and its derivatives given the value of the field
    function EFTCAMBDesignerGBDCoupling( self, phi, deriv ) result(Omega)

        implicit none

        class(EFTCAMB_GBD_designer) :: self     !< the base class

        real(dl) :: phi         !< the value of the field
        integer  :: deriv       !< the derivative of the coupling function that needs to be computed

        real(dl) :: Omega       !< the result, Omega(phi) or its derivatives w.r.t. phi

        if(self%coupling_type == 1) then

            !> Linear coupling

            if (deriv == 0) then
                Omega = 1._dl + self%xi*phi

            else if (deriv == 1) then
                Omega = self%xi

            else if (deriv == 2) then
                Omega = 0._dl

            else if (deriv == 3) then
                Omega = 0._dl

            else
                write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong derivative"
                stop
            end if


        else if (self%coupling_type == 2) then

            !> Quadratic coupling

            if (deriv == 0) then
                Omega = 1._dl + self%xi*phi**2

            else if (deriv == 1) then
                Omega = 2._dl*self%xi*phi

            else if (deriv == 2) then
                Omega = 2.+dl*self%xi

            else if (deriv == 3) then
                Omega = 0._dl

            else
                write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong derivative"
                stop
            end if


        else if (self%coupling_type == 3) then

            !> Exponential coupling

            if (deriv == 0) then
                Omega = exp(self%xi*phi)

            else if (deriv == 1) then
                Omega = self%xi*exp(self%xi*phi)

            else if (deriv == 2) then
                Omega = self%xi**2 * exp(self%xi*phi)

            else if (deriv == 3) then
                Omega = self%xi**3 * exp(self%xi*phi)

            else
                write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong derivative"
                stop
            end if

        else if (self%coupling_type == 4) then

            !> Negative exponential coupling
            if (deriv == 0) then
                Omega = exp(-self%xi*phi)

            else if (deriv == 1) then
                Omega = - self%xi * exp(-self%xi*phi)

            else if (deriv == 2) then
                Omega = self%xi**2 * exp(-self%xi*phi)

            else if (deriv == 3) then
                Omega = -self%xi**3 * exp(-self%xi*phi)

            else
                write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong derivative"
                stop
            end if

        else
            write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong coupling type"
            stop
        end if

    end function EFTCAMBDesignerGBDCoupling

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of designer GBD
    subroutine EFTCAMBDesignerGBDInitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_GBD_designer)                     :: self             !< the base class
        type(EFTCAMB_parameter_cache),  intent(in)      :: params_cache     !< a EFTCAMB parameter cache containing cosmological parameters
        integer,                        intent(in)      :: feedback_level   !< a level of feedback from the background code. 0=none; 1=some; 2=chatty
        logical,                        intent(out)     :: success          !< whether the background initialization succeded or not
        character(LEN=*), optional   , intent(in)    :: outroot        !< the output root for the debug files

        real(dl) :: phi_i, dphi_i
        real(dl) :: init_lambda


        !> some feedback
        if (feedback_level>0) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a)') "EFTCAMB designer GBD background solver"
            write(*,'(a)')
        end if

        !> some temporary output
        write(*,*) "GBD designer: a_ini, dphi_ini, xi:", exp(self%x_initial), self%dphi_ini, self%xi

        !> intial value of Lambda
        init_lambda = - 3._dl * params_cache%h0_Mpc**2 * params_cache%omegav *self%DesGBDxDE%value(exp(self%x_initial))

        !> initialize interpolating functions
        call self%EFTOmega%initialize   ( self%designer_num_points, self%x_initial, self%x_final, null_value=0.d0 )
        call self%EFTLambda%initialize  ( self%designer_num_points, self%x_initial, self%x_final, null_value=init_lambda )
        !call self%EFTLambda%initialize  ( self%designer_num_points, self%x_initial, self%x_final, null_value=0.d0 )
        call self%EFTc%initialize       ( self%designer_num_points, self%x_initial, self%x_final, null_value=0.d0 )


        !> Re-introducing this at some point
        !> debug code:
        if ( DebugEFTCAMB ) then
            print*, 'EFTCAMB DEBUG ( GBD designer ): Printing GBD results'
            call CreateTxtFile( TRIM(outroot)//'debug_designer_GBD_solution.dat', 33 )
            write(33,'(a)') "#  1:a       2:phi         3:dphi       4:ddphi        5:V       6:Omega       7:c       8:Lambda      9:X     10:Xp"
            call self%solve_designer_equations( params_cache, success=success )
            close(33)
            print*, 'EFTCAMB DEBUG ( GBD designer ): GBD results printed'
        end if

        !> solve the background equations and store the solution:
        write(*,*) "solving the equations"
        call self%solve_designer_equations( params_cache, success=success )
        write(*,*) "equations solved: success,", success


        !write(27,*) self%EFTOmega%value(1.d0)

        if ( DebugEFTCAMB ) then
            write(*,'(a)') "EFTCAMB designer GBD background solver. Equation solved."
            write(*,"(L1)") "success:",success
        end if

    end subroutine EFTCAMBDesignerGBDInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the designer GBD background equations.
    subroutine EFTCAMBDesignerGBDSolveDesignerEquations( self, params_cache, success )

        implicit none

        class(EFTCAMB_GBD_designer)                 :: self          !< the base class.
        type(EFTCAMB_parameter_cache), intent(in)   :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                       :: success       !< whether the calculation ended correctly or not

        integer, parameter :: num_eq = 2   !<  Number of equations

        real(dl) :: Omegam_EFT, Omegavac_EFT, OmegaMassiveNu_EFT, OmegaGamma_EFT, OmegaNu_EFT
        real(dl) :: Omegarad_EFT

        !real(dl) :: y(num_eq+1), ydot(num_eq+1)
        real(dl) :: y(num_eq), ydot(num_eq)
        real(dl) :: N, dN, N_fin

        integer :: i, i_y

        !> adding the parameters for dverk
        integer :: ind = 1
        real(dl), dimension(24) :: c = 0.d0
        integer, parameter :: nw = num_eq + 2
        real(dl), dimension(nw,9) :: w
        real(dl) :: tol = 1.d-6

        !> fixing the index
        ind = 1

        success = .True.

        ! 1) Cosmological parameters:
        Omegam_EFT         = params_cache%omegab + params_cache%omegac
        Omegavac_EFT       = params_cache%omegav
        OmegaMassiveNu_EFT = params_cache%omegan
        OmegaGamma_EFT     = params_cache%omegag
        OmegaNu_EFT        = params_cache%omegar

        Omegarad_EFT       = OmegaGamma_EFT + OmegaNu_EFT

        ! 2) Set initial conditions:
        !> need to modify this
        !y(1) = self%x_initial
        !y(2) = self%phi_ini
        !y(3) = self%dphi_ini
        !ydot = 0._dl

        N = self%x_initial
        y(1) = self%phi_ini
        y(2) = self%dphi_ini
        ydot = 0._dl

        ! 3) Solve the equation of motion
        !> defining the time-step
        dN = (self%x_final - self%x_initial)/self%designer_num_points

        if( DebugEFTCAMB ) then
            !> test
            write(*,*) "Solving the equation of motion for phi"
            write(*,*) "dN:",dN
        end if

        !> Loop to fill the interpolation arrays
        do  i = 1, self%EFTOmega%num_points


            !> tests (to be removed in the final version)
            N_fin = N + dN

            !> calling the solver, in this case gl10
            !> substituting with dverk_std
            !call gl10(num_eq+1,derivs, y, dN)
            call dverk_std(num_eq, derivs, N, y, N_fin, tol, ind, c, nw, w)

            !> filling the interpolation arrays
            !write(*,*) 'output at:', exp(N)
            call output(num_eq, N, y, i)

            !> DEBUG
            if ( i == self%EFTOmega%num_points ) then
                write(*,*) 'Omega0:', self%omega_phi(y(1), 0)
            end if


            if (ind<0) then
                write(*,*) "ERROR dverk: ind:", ind
                success = .False.
                return
            end if

        end do

        !> end of the reconstruction
        return

    contains

        subroutine derivs(num,N, y,yprime)

            implicit none

            integer :: num                          !< number of equations
            real(dl), dimension(num) :: y, yprime   !< input status of the system and output derivative
            real(dl) :: N,a,a2                      !< efolds number, scale factor and squared scale factor

            real(dl) :: om, omp, ompp               !< coupling functoin and its derivatives
            real(dl) :: H, adotdot                  !< expansion history parameters

            real(dl) :: rhonu, presnu               !< massive neutrinos background variables
            real(dl) :: rhonu_tot, presnu_tot
            real(dl) :: grhormass_t
            integer  :: nu_i

            real(dl) :: Em,   Er,   Enu,   X        !< normalized energy densities
            real(dl) :: Em_p, Er_p, Enu_p, X_p      !< derivatives of normalized energy densities
            real(dl) :: X_pp                        !< second derivative of the normalized DE density

            associate (phi => y(1), pi=>y(2),  dphi => yprime(1) , ddphi => yprime(2))

                !> convert N in a
                a = exp(N)

                !> compute coupling function and its derivatives
                om      = self%omega_phi(phi, 0)
                omp     = self%omega_phi(phi, 1)
                ompp    = self%omega_phi(phi, 2)

                !> normalized energy densities
                Em = Omegam_EFT   * exp(-3._dl*N)
                Er = OmegaRad_EFT * exp(-4._dl*N)

                !> compute the normalized dark energy density and its derivatives
                X       = self%DesGBDxDE%value(a)
                X_p     = self%DesGBDxDE%first_derivative(a)*a
                X_pp    = self%DesGBDxDE%second_derivative(a)*a**2 +self%DesGBDxDE%first_derivative(a)*a

                !> compute massive neutrinos contribution:
                rhonu_tot   = 0._dl
                presnu_tot  = 0._dl
                Enu         = 0._dl
                Enu_p       = 0._dl

                if ( params_cache%Num_Nu_Massive /= 0) then
                    do nu_i = 1, params_cache%Nu_mass_eigenstates

                        rhonu  = 0._dl
                        presnu = 0._dl
                        grhormass_t= params_cache%grhormass(nu_i)/a**2

                        call params_cache%Nu_background(a*params_cache%nu_masses(nu_i),rhonu,presnu)

                        rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                        presnu_tot = presnu_tot + grhormass_t*presnu

                        Enu    = Enu   + params_cache%grhormass(nu_i)/3._dl/a**4/params_cache%h0_Mpc**2*rhonu
                        Enu_p  = Enu_p - params_cache%grhormass(nu_i)/params_cache%h0_Mpc**2/a**4*(rhonu +presnu)

                    end do
                end if

                !> now the equations of motion

                !> derivative of N=ln(a)
                !dN = 1.d0

                !> phi prime
                dphi = pi

                !> phi prime prime
                ddphi = ((om - 1._dl)*(3._dl * Em + 4._dl * Er - Enu_p) - om * Omegavac_EFT * X_p)/omp/(Em+Er+Enu+Omegavac_EFT*X) &
                        - (1._dl + ompp)*dphi**2/omp + (1._dl+0.5_dl*(3._dl*Em+4._dl*Er-Enu_p-Omegavac_EFT*X_p)/(Em+Er+Enu+Omegavac_EFT*X))*dphi

            end associate

        end subroutine derivs

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that takes the solution of the background GBD equations and stores the values of
        !> the EFT functions.
        subroutine output( num, N,  y, ind )

            implicit none

            integer , intent(in)                    :: num  !< number of equations in the ODE system.
            integer , intent(in)                    :: ind  !< index of the EFT functions interpolation tables to fill.
            real(dl), intent(in) , dimension(num)   :: y    !< input status of the system.
            real(dl), intent(in)                    :: N    !< e-fold number

            real(dl) :: ydot(num)                           !< array of derivatives of the system
            !real(dl) :: EFT_E_gfun, EFT_E_gfunp, EFT_E_gfunpp, EFT_E_gfunppp !< effective dark energy variables (can be removed)

            real(dl) :: om, omp, ompp, omppp         !< coupling function and its derivatives

            !> some massive neutrinos variables
            real(dl) :: rhonu_tot, presnu_tot, presnudot_tot, presnudotdot_tot
            real(dl) :: rhonu, presnu, grhormass_t, presnudot, presnudotdot

            real(dl) :: Em,    Er,    Enu,    X     !< normalized energy densities
            real(dl) :: Em_p,  Er_p,  Enu_p,  X_p   !< derivatives of normalized energy densities
            real(dl) :: Em_pp, Er_pp, Enu_pp, X_pp  !< second derivatives of normalized energy densities
            real(dl) :: Enu_ppp                     !< third derivative of normalized energy density

            integer  :: nu_i
            logical  :: is_open

            !> other parameters
            real(dl) :: a, phi, dphi, ddphi, dddphi
            real(dl) :: adotoa, Hdot, adotdotoa
            real(dl) :: calF, V, Vprime, Vdot
            real(dl) :: Etot

            !> convert N in a
            !N = y(1)
            a = exp(N)

            !> AZ TEST:
            !write(*,*) "GBD designer output. a, X, y :", a, self%DesGBDxDE%value(a), y

            !> extract values of the field
            phi     = y(1)
            dphi    = y(2)

            !> compute coupling function and its derivatives
            om      = self%omega_phi(phi, 0)
            omp     = self%omega_phi(phi, 1)
            ompp    = self%omega_phi(phi, 2)
            omppp   = self%omega_phi(phi, 3)

            !> normalized energy densities
            Em = Omegam_EFT   * exp(-3._dl*N)
            Er = OmegaRad_EFT * exp(-4._dl*N)

            !> compute the normalized dark energy density
            X       = self%DesGBDxDE%value(a)
            X_p     = self%DesGBDxDE%first_derivative(a)*a
            X_pp    = self%DesGBDxDE%second_derivative(a)*a**2 +self%DesGBDxDE%first_derivative(a)*a


            !> compute massive neutrinos contribution:
            rhonu_tot   = 0._dl
            presnu_tot  = 0._dl
            Enu         = 0._dl
            Enu_p       = 0._dl

            if ( params_cache%Num_Nu_Massive /= 0) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates

                    rhonu  = 0._dl
                    presnu = 0._dl
                    grhormass_t= params_cache%grhormass(nu_i)/a**2
                    call params_cache%Nu_background(a*params_cache%nu_masses(nu_i),rhonu,presnu)
                    rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                    presnu_tot = presnu_tot + grhormass_t*presnu

                    Enu     = Enu   + params_cache%grhormass(nu_i)/3._dl/a**4/params_cache%h0_Mpc**2*rhonu
                    Enu_p   = Enu_p - params_cache%grhormass(nu_i)/params_cache%h0_Mpc**2/a**4*(rhonu +presnu)

                end do
            end if

            !> calculate H, Hdot and adotdotoa
            Etot        = Em+Er+Enu+Omegavac_EFT*X
            adotoa      = a * params_cache%h0_Mpc* sqrt(Etot)
            Hdot        = a**2 * params_cache%h0_Mpc**2 * ( (Em+Er+Enu+Omegavac_EFT*X) + &
                          0.5_dl * (-3._dl * Em - 4._dl*Er + Enu_p + Omegavac_EFT*X_p) )
            adotdotoa   = Hdot + adotoa**2

            !> start filling the EFT interpolated functions
            call derivs( num ,N, y, ydot )

            !> extract ddphi
            ddphi = ydot(2)

            !> calculating the potential - needed for the function Lambda
            !> begin with {\cal F} in the notes
            calF = om - dphi**2 / 6._dl + omp * dphi

            !> inverted Friedmann equation
            V = 3._dl*(Em+Er+Enu+Omegavac_EFT*X)*calF - 3.d0*(Em+Er+Enu)
            V = V * params_cache%h0_Mpc**2 ! this is in Mpc-2

            !> get Vprime from the equation of motion for the scalar field
            Vprime = (3._dl * omp * adotdotoa - adotoa**2 *ddphi - (adotoa**2 + adotdotoa)*dphi)/a**2 ! in Mpc-2

            !> compute Vdot - needed for Lambda dot -
            Vdot = Vprime * dphi * adotoa

            !------> calculating dddphi: this requires many steps
            ! 1) Compute everything of massive nu again to get the time derivatives:
            rhonu_tot        = 0._dl
            presnu_tot       = 0._dl
            presnudot_tot    = 0._dl
            presnudotdot_tot = 0._dl

            Enu     = 0._dl
            Enu_p   = 0._dl
            Enu_pp  = 0._dl
            Enu_ppp = 0._dl

            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates

                    rhonu        = 0._dl
                    presnu       = 0._dl
                    presnudot    = 0._dl
                    presnudotdot = 0._dl

                    grhormass_t  = params_cache%grhormass(nu_i)/a**2

                    call params_cache%Nu_background(a*params_cache%nu_masses(nu_i),rhonu,presnu)
                    presnudot = params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),adotoa,presnu)
                    presnudotdot = params_cache%Nu_pidotdot(a*params_cache%nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)

                    rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                    presnu_tot = presnu_tot + grhormass_t*presnu
                    presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
                    presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot &
                                    & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))

                    Enu     = Enu   + params_cache%grhormass(nu_i)/3._dl/a**4/params_cache%h0_Mpc**2*rhonu
                    Enu_p   = Enu_p - params_cache%grhormass(nu_i)/params_cache%h0_Mpc**2/a**4*(rhonu +presnu)
                    Enu_pp  = Enu_pp + 3._dl/params_cache%h0_Mpc**2*params_cache%grhormass(nu_i)/a**4*(rhonu +presnu)&
                                & -grhormass_t*(presnudot -4._dl*adotoa*presnu)/params_cache%h0_Mpc**3/sqrt(Etot)/a**3
                    Enu_ppp = Enu_ppp -9._dl/params_cache%h0_Mpc**2*params_cache%grhormass(nu_i)/a**4*(rhonu +presnu)&
                                & +(3._dl/adotoa/params_cache%h0_Mpc**2/a**2+Hdot/adotoa**3/params_cache%h0_Mpc**2/a**2)&
                                &*grhormass_t*(presnudot -4._dl*adotoa*presnu)&
                                & -grhormass_t*(presnudotdot &
                                & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))/adotoa**2/params_cache%h0_Mpc**2/a**2
                end do
            end if

            !> E tot
            Etot = Em + Er + Enu + Omegavac_EFT*X

            !> check dddphi with the MAPLE output

            !> now compute dddphi
            dddphi = (3._dl*Em + 4._dl * Er - Enu_p)/Etot*dphi                                                                      &
                    -(om - 1._dl)*(3._dl*Em + 4._dl*Er - Enu_p)/Etot*ompp/omp**2 *dphi                                              &
                    +(om-1._dl)/Etot/omp * (- 9._dl*Em -16._dl*Er - Enu_pp)                                                         &
                    -(om-1._dl)/omp/Etot**2*(3._dl*Em+4._dl*Er - Enu_p)*(-3._dl*Em-4._dl*Er+Enu_p + Omegavac_EFT*X_p)               &
                    + Omegavac_EFT*X_p * ompp * dphi / omp**2 / Etot                                                                &
                    - Omegavac_EFT*X_pp / omp / Etot                                                                                &
                    + Omegavac_EFT*X_p * (-3._dl*Em-4._dl*Er + Enu_p + Omegavac_EFT*X_p)/omp/Etot**2                                &
                    +(1._dl+ompp)/omp**2 * dphi**3 * ompp                                                                           &
                    - omppp*dphi**3/omp                                                                                             &
                    -2._dl * (1._dl+ompp)*dphi*ddphi/omp                                                                            &
                    +(0.5_dl * (-9._dl*Em - 16._dl*Er - Enu_pp - Omegavac_EFT*X_pp)/Etot                                            &
                    -0.5_dl*(3._dl*Em+4._dl*Er-Enu_p-Omegavac_EFT*X_p)*(-3._dl*Em-4._dl*Er+Enu_p+Omegavac_EFT*X_p)/Etot**2)*dphi    &
                    +0.5_dl* (5._dl*Em+6._dl*Er +2._dl*Enu+2._dl*Omegavac_EFT*X - Enu_p-Omegavac_EFT*X_p)/Etot*ddphi


            !> Filling the EFT functions:
            !> Omega
            self%EFTOmega%y(ind)    = om-1._dl
            self%EFTOmega%yp(ind)   = omp * dphi / a
            self%EFTOmega%ypp(ind)  = (ompp * dphi**2 + omp*ddphi - omp * dphi)/a**2
            self%EFTOmega%yppp(ind) = (omppp * dphi**3 + 3._dl * ompp * dphi * ddphi - 3._dl* ompp*dphi**2 &
                                        - 3._dl*omp * ddphi + omp * dddphi + 2._dl*omp *dphi)/a**3

            !> c
            self%EFTc%y(ind)    = 0.5_dl*adotoa**2 * dphi**2
            self%EFTc%yp(ind)   = (- dphi**2 + dphi * ddphi)*adotoa**3 + dphi**2 * adotoa * Hdot

            !> Lambda
            self%EFTLambda%y(ind)   = self%EFTc%y(ind)  - V    * a**2
            self%EFTLambda%yp(ind)  = self%EFTc%yp(ind) - Vdot * a**2


            !> Debug info
            if ( DebugEFTCAMB ) then
            !if ( .true. ) then
                inquire( unit=33, opened=is_open )
                if ( is_open ) then
                    !write (33,'(20E15.5)') a, phi, dphi, ddphi,V, self%EFTOmega%y(ind), self%EFTc%y(ind), self%EFTLambda%y(ind), X, X_p/a, &
                    !        & self%EFTOmega%yp(ind), self%EFTOmega%ypp(ind), self%EFTOmega%yppp(ind), adotoa, adotdotoa, Hdot
                    !> the following was used to check the solutions and the pieces of the differential equation
                    !                       0  1   2      3    4   5      6          7
                    write (33,'(20E15.5)') a, phi, dphi, ddphi,V,adotoa,adotdotoa,(1._dl + ompp)/omp,                   &
                                            (1._dl+0.5_dl*(3._dl*Em+4._dl*Er-Enu_p-X_p)/(Em+Er+Enu+X)),                 & !8
                                            ((om - 1._dl)*(3._dl * Em + 4._dl * Er - Enu_p) - om*X_p)/omp/(Em+Er+Enu+X),& !9
                                            X_p, Em, Er, self%EFTOmega%y(ind), self%EFTc%y(ind), self%EFTLambda%y(ind)
                                            !10  11  12    13                   14                  15
                end if
            end if


        end subroutine output

    end subroutine EFTCAMBDesignerGBDSolveDesignerEquations




    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBDesignerGBDComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_GBD_designer)  :: self   !< the base class

        self%parameter_number = 4
        self%parameter_number = self%parameter_number +self%DesGBDxDE%parameter_number

    end subroutine EFTCAMBDesignerGBDComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBDesignerGBDFeedback( self, print_params )

        implicit none

        class(EFTCAMB_GBD_designer)     :: self         !< the base class
        logical, optional               :: print_params !< optional flag that decides whether to print numerical values
                                                     !! of the parameters.

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        !> print model functions informations:
        write(*,*)
        write(*,'(a,I3)')  '   EFTxDE              =', self%EFTxDE
        write(*,*)
        write(*,'(a24,F12.6)') '   phi_ini             =', self%phi_ini
        write(*,'(a24,F12.6)') '   dphi_ini            =', self%dphi_ini
        write(*,'(a24,F12.6)') '   xi                  =', self%xi
        write(*,*)
        write(*,'(a,I3)')      '   coupling type       =', self%coupling_type
        write(*,'(a24,F12.6)') '   a_ini               =', exp(self%x_initial)
        write(*,*)


        call self%DesGBDxDE%feedback( print_params )

    end subroutine EFTCAMBDesignerGBDFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerGBDParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_GBD_designer)     :: self   !< the base class
        integer     , intent(in)        :: i      !< the index of the parameter
        character(*), intent(out)       :: name   !< the output name of the i-th parameter

        !> check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        !> the first parameter is phi_ini:
        else if ( i==1 ) then
            name = TRIM('phi_ini')
            return
        !> the second parameter is dphi_ini
        else if ( i==2 ) then
            name = TRIM('dphi_ini')
            return
        !> the third parameter is xi
        else if ( i==3 ) then
            name = TRIM('xi')
            return
        !> the fourth parameter is x_initial
        else if ( i==4 ) then
            name = TRIM('x_initial')
            return
        !> the other parameters are the w_DE parameters:
        else
            call self%DesGBDxDE%parameter_names( i-3, name )
            return
        end if

    end subroutine EFTCAMBDesignerGBDParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerGBDParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_GBD_designer)     :: self       !< the base class
        integer     , intent(in)        :: i         !< The index of the parameter
        character(*), intent(out)       :: latexname !< the output latex name of the i-th parameter

        !> check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        !> the first parameter is phi_ini:
        else if ( i==1 ) then
            latexname = TRIM('\phi_{\rm ini}')
            return
        !> the second parameter is dphi_ini
        else if ( i==2 ) then
            latexname = TRIM('\phi_{\rm ini}^{\prime}')
            return
        !> the third parameter is xi
        else if ( i==3 ) then
            latexname = TRIM('\xi')
            return
        else if ( i==3 ) then
            latexname = TRIM('x_{\rm ini}')
            return
        !> the other parameters are the w_DE parameters:
        else
            call self%DesGBDxDE%parameter_names_latex( i-4, latexname )
            return
        end if

    end subroutine EFTCAMBDesignerGBDParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerGBDParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_GBD_designer)     :: self   !< the base class
        integer , intent(in)            :: i      !< The index of the parameter
        real(dl), intent(out)           :: value  !< the output value of the i-th parameter

        !> check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        !> the first parameter is phi_ini:
        else if ( i==1 ) then
            value = self%phi_ini
            return
        !> the second parameter is dphi_ini
        else if ( i==2 ) then
            value = self%dphi_ini
            return
        !> the third parameter is xi
        else if ( i==3 ) then
            value = self%xi
            return
        !> the other parameters are the w_DE parameters:
        else if ( i==4 ) then
            value = self%x_initial
            return
        else
            call self%DesGBDxDE%parameter_value( i-4, value )
            return
        end if

    end subroutine EFTCAMBDesignerGBDParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBDesignerGBDBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_GBD_designer)                   :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout)  :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout)  :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        if ( a .ge. 1 ) then
            x = 0.d0
        else
            x   = log(a)
        end if

        call self%EFTOmega%precompute( x, ind, mu )

        eft_cache%EFTOmegaV    = self%EFTOmega%value( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaP    = self%EFTOmega%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPP   = self%EFTOmega%second_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPPP  = self%EFTOmega%third_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTc         = self%EFTc%value( x, index=ind, coeff=mu )
        eft_cache%EFTcdot      = self%EFTc%first_derivative( x, index=ind, coeff=mu )

        !> if a < a_ini then the Lambda funciton is related to the potential and it's not constant.
        !  so we need the following tweak.
        ! NOTE: this tweak assumes that phi = 0 before x_initial!!!!
        if ( x .ge. self%x_initial) then
            eft_cache%EFTLambda    = self%EFTLambda%value( x, index=ind, coeff=mu )
            eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( x, index=ind, coeff=mu )
        else
            eft_cache%EFTLambda = -3._dl * eft_par_cache%h0_Mpc**2 * a**2 * eft_par_cache%omegav *self%DesGBDxDE%value(a)
            eft_cache%EFTLambdadot = -3._dl * eft_par_cache%h0_Mpc**2 * a**3 * eft_par_cache%omegav *self%DesGBDxDE%first_derivative(a) * eft_cache%adotoa
        end if

    end subroutine EFTCAMBDesignerGBDBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBDesignerGBDSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_GBD_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTGamma1V  = 0._dl
        eft_cache%EFTGamma1P  = 0._dl
        eft_cache%EFTGamma2V  = 0._dl
        eft_cache%EFTGamma2P  = 0._dl
        eft_cache%EFTGamma3V  = 0._dl
        eft_cache%EFTGamma3P  = 0._dl
        eft_cache%EFTGamma4V  = 0._dl
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl
        eft_cache%EFTGamma5V  = 0._dl
        eft_cache%EFTGamma5P  = 0._dl
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBDesignerGBDSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBDesignerGBDComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_GBD_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBDesignerGBDComputeDtauda               !< the output dtauda

        real(dl) :: temp



        !> this is just a trick
        if (a .le. 1.d-10) then
            !a = 1.d-10
            temp = eft_cache%grhoa2 !+ 3._dl*eft_par_cache%h0_Mpc**2 * self%DesGBDxDE%value(1.d-10) * a**4
        else
            temp = eft_cache%grhoa2 + 3._dl*eft_par_cache%h0_Mpc**2 * eft_par_cache%omegav * self%DesGBDxDE%value(a) * a**4
        end if
        EFTCAMBDesignerGBDComputeDtauda = sqrt(3._dl/temp)

!write(*,*) "dtauda: temp=",temp

    end function EFTCAMBDesignerGBDComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBDesignerGBDComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_GBD_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        !eft_cache%grhov_t = eft_par_cache%grhov*self%DesGBDwDE%integral(a)
        !> Adapting the Friedmann equation
        eft_cache%grhov_t = 3._dl * eft_par_cache%h0_Mpc**2 * eft_par_cache%omegav *self%DesGBDxDE%value(a) * a**2

        !> for some reason the Friedmann equation is missing the radiation contribution
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )

    end subroutine EFTCAMBDesignerGBDComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBDesignerGBDComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_GBD_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        !> Dark energy pressure from the reconstructed dark energy density
        eft_cache%gpiv_t    = eft_par_cache%omegav*(-a**3 * self%DesGBDxDE%first_derivative(a) * eft_par_cache%h0_Mpc**2 &
                            & - 3._dl*eft_par_cache%h0_Mpc**2 * a**2 * self%DesGBDxDE%value(a))
        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )

        !> replacing the following
        !eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
        !    & +eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%DesGBDwDE%value(a) +1.5_dl*self%DesGBDwDE%value(a)**2 -0.5_dl*a*self%DesGBDwDE%first_derivative(a) ) &
        !    & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

        eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
            & +3._dl*eft_par_cache%h0_Mpc**2 * a**2* eft_cache%adotoa * (-7._dl/6._dl*self%DesGBDxDE%first_derivative(a)*a + a**2 * self%DesGBDxDE%second_derivative(a)/6._dl +2._dl*self%DesGBDxDE%value(a)/3._dl) * eft_par_cache%omegav &
            & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

    end subroutine EFTCAMBDesignerGBDComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBDesignerGBDAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_GBD_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBDesignerGBDAdditionalModelStability          !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        !> nothing to add right now.
        EFTCAMBDesignerGBDAdditionalModelStability = .True.


    end function EFTCAMBDesignerGBDAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model for the sampling of the DE density
    !subroutine EFTCAMBDesignerGBDInitModelParametersSampling( self, array, sampling_params )
    !
    !    implicit none
    !
    !    class(EFTCAMB_GBD_designer)                             :: self             !< the base class
    !    real(dl), dimension(self%parameter_number), intent(in)  :: array            !< array of parameters
    !    type(EFTSamplingParameters)                             :: sampling_params  !< the class sampling parameters that contains the DE density
    !
    !    self%phi_ini    = array(1)
    !    self%dphi_ini   = array(2)
    !    self%xi         = array(3)
    !    self%x_initial  = log(10._dl**array(4))
    !
    !    !> Now I have to initialize the DE density
    !    call self%DesGBDxDE%init_from_code( sampling_params%xDE_scale_factor_means, sampling_params%xDE_sample )
    !
    !
    !end subroutine EFTCAMBDesignerGBDInitModelParametersSampling

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_designer_GBD

!----------------------------------------------------------------------------------------
