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

!> @file 04p9_interpolated_parametrizations_1D.f90
!! This file contains the interpolation of a function of the scale factor.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the interpolation parametrization,
!! inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone
!> @author Alex Zucca

module EFTCAMB_interpolated_function_1D

    use precision
    use AMLutils
    use GBD_Utils
    use IniFile
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none


    private

    public interpolated_function_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: interpolated_function_1D

        ! Interpolation function variables
        character(len=EFT_names_max_length)             :: file, file_name, num_bins_string
        logical                                         :: is_function_loaded = .false.

        integer :: num_points
        real(dl) :: grid_width


        real(dl), allocatable, dimension(:) :: x              !< array containing the values of x.
        real(dl), allocatable, dimension(:) :: y              !< array containing the values of the function \f$ y_i=f(x_i) \f$.
        real(dl), allocatable, dimension(:) :: yp             !< array containing the values of the function derivative \f$ yp_i= \frac{d f(x_i)}{dx} \f$.
        real(dl), allocatable, dimension(:) :: ypp            !< array containing the values of the function second derivative \f$ ypp_i= \frac{d^2 f(x_i)}{dx^2} \f$.
        real(dl), allocatable, dimension(:) :: yppp           !< array containing the values of the function third derivative \f$ yppp_i= \frac{d^3 f(x_i)}{dx^3} \f$.
        real(dl), allocatable, dimension(:) :: yint           !< array containing the values of the function w DE integral \f$ yint_i= \exp\left(-3\int_1^{x_i} \frac{1+f(x)}{x} \, dx \right) \f$.

        !> other variables for the interpolation
        real(dl), allocatable, dimension(:) :: y2             !< array containing the values of the function \f$ y_i=f(x_i) \f$.
        real(dl), allocatable, dimension(:) :: yp2            !< array containing the values of the function derivative \f$ yp_i= \frac{d f(x_i)}{dx} \f$.
        real(dl), allocatable, dimension(:) :: ypp2           !< array containing the values of the function second derivative \f$ ypp_i= \frac{d^2 f(x_i)}{dx^2} \f$.
        real(dl), allocatable, dimension(:) :: yppp2          !< array containing the values of the function third derivative \f$ yppp_i= \frac{d^3 f(x_i)}{dx^3} \f$.
        real(dl), allocatable, dimension(:) :: yint2

    contains

        ! utility functions
        procedure :: set_param_number    => InterpolatedFunction1DSetParamNumber   !< subroutine that sets the number of parameters of the interpolated function.
        procedure :: init_parameters     => InterpolatedFunction1DInitParams       !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value     => InterpolatedFunction1DParameterValues  !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback            => InterpolatedFunction1DFeedback         !< subroutine that prints to screen the informations about the function.
        procedure :: init_func_from_file => InterpolatedFunction1DInitFromFile     !< subroutine that initializes a few parameters from files
        procedure :: set_param_names     => InterpolatedFunction1DSetParamNames    !< subroutine that sets the file_name

        ! Initialization here
        procedure :: initialize_function => InterpolatedFunction1DInitialization   !< subroutine that initializes the interpolation
        procedure :: init_from_code      => InterpolatedFunction1DInitFromCode

        ! evaluation procedures
        procedure :: value               => InterpolatedFunction1DValue            !< function that returns the value of the interpolation.
        procedure :: first_derivative    => InterpolatedFunction1DFirstDerivative  !< function that returns the first derivative of the interpolation.
        procedure :: second_derivative   => InterpolatedFunction1DSecondDerivative !< function that returns the second derivative of the interpolation.
        procedure :: third_derivative    => InterpolatedFunction1DThirdDerivative  !< function that returns the third derivative of the interpolation.
        procedure :: integral            => InterpolatedFunction1DIntegral         !< function that returns the strange integral that we need for w_DE.



    end type interpolated_function_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the interpolated function
    ! ---------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
!> Subroutine the sets the number of parameter of the interpolation of the function
subroutine InterpolatedFunction1DSetParamNumber(self)
    implicit none
    class(interpolated_function_1D) :: self

    self%parameter_number = 0 !no parameters here, just an interpolation

end subroutine InterpolatedFunction1DSetParamNumber


! ---------------------------------------------------------------------------------------------
!> Subroutine that initializes the function parameters based on the values found in an input array.
subroutine InterpolatedFunction1DInitParams( self, array )

    implicit none

    class(interpolated_function_1D) :: self
    real(dl), dimension(self%parameter_number), intent(in) :: array

    ! That's it, nothing to do..
end subroutine InterpolatedFunction1DInitParams



! ---------------------------------------------------------------------------------------------
!> Subroutine that returns the value of the function i-th parameter.
subroutine InterpolatedFunction1DParameterValues( self, i, value )

    implicit none

    class(interpolated_function_1D) :: self
    integer,intent(in)              :: i
    real(dl), intent(out)           :: value

    value = 0.d0

end subroutine InterpolatedFunction1DParameterValues


! ---------------------------------------------------------------------------------------------
!> Subroutine that prints to screen the informations about the function.
subroutine InterpolatedFunction1DFeedback( self, print_params )

    implicit none

    class(interpolated_function_1D) :: self
    logical, optional               :: print_params
    logical                         :: print_params_temp

    if ( present(print_params) ) then
        print_params_temp = print_params
    else
        print_params_temp = .True.
    end if

    !> some useless stuff
    !write(*,*) self%file_name
    !pause

    write(*,*)     'Interpolated Function: ', self%name
    if ( print_params_temp ) then
        write(*,*) TRIM(self%file_name), ' = ', TRIM(self%file)
    end if

end subroutine InterpolatedFunction1DFeedback


! ---------------------------------------------------------------------------------------------
!> Subroutine that reads a Ini file looking for name of the file containing the function.
subroutine InterpolatedFunction1DInitFromFile( self, Ini )

    implicit none

    class(interpolated_function_1D)     :: self   !< the base class
    type(TIniFile)                      :: Ini    !< Input ini file
    character(len=EFT_names_max_length) :: file_name, file, n_bins_string
    integer :: i

    !> intializing the number of points
    n_bins_string = self%num_bins_string
    self%num_points = Ini_Read_Int_File( Ini, TRIM(n_bins_string), 40 )

    !> file containing the table of the function
    file_name = self%file_name
    file = Ini_Read_String_File(Ini, TRIM(file_name))

    self%file = TRIM(file)


    ! allocate the other vectors:
    if ( allocated(self%x)    ) deallocate( self%x    )
    if ( allocated(self%y)    ) deallocate( self%y    )
    if ( allocated(self%yp)   ) deallocate( self%yp   )
    if ( allocated(self%ypp)  ) deallocate( self%ypp  )
    if ( allocated(self%yppp) ) deallocate( self%yppp )
    if ( allocated(self%yint) ) deallocate( self%yint )

    !> spline arrays
    if ( allocated(self%y2)    ) deallocate( self%y2    )
    if ( allocated(self%yp2)   ) deallocate( self%yp2   )
    if ( allocated(self%ypp2)  ) deallocate( self%ypp2  )
    if ( allocated(self%yppp2) ) deallocate( self%yppp2 )
    if ( allocated(self%yint2) ) deallocate( self%yint2 )

    allocate( self%x( self%num_points )    )
    allocate( self%y( self%num_points )    )
    allocate( self%yp( self%num_points )   )
    allocate( self%ypp( self%num_points )  )
    allocate( self%yppp( self%num_points ) )
    allocate( self%yint( self%num_points ) )

    !> spline arrays
    allocate( self%y2( self%num_points )    )
    allocate( self%yp2( self%num_points )   )
    allocate( self%ypp2( self%num_points )  )
    allocate( self%yppp2( self%num_points ) )
    allocate( self%yint2( self%num_points ) )


    !> Open the file
    open(unit=7, file = TRIM(self%file), status = "unknown")

    do i  = 1, self%num_points
        write(*,*) i
        ! read a, f(a) from file and fill the array
        read(7,*) self%x(i), self%y(i)
    end do

    self%grid_width = abs((self%x(self%num_points) - self%x(1))/self%num_points)

    self%is_function_loaded = .false.

    !initialize the function parameters from the vector:
    !call self%init_parameters( parameters )

end subroutine InterpolatedFunction1DInitFromFile





! ---------------------------------------------------------------------------------------------
!> subroutine that sets the name of the file that we want to load
subroutine InterpolatedFunction1DSetParamNames(self, param_names, param_names_latex)

    implicit none

    class(interpolated_function_1D)                   :: self              !< the base class
    character(*), intent(in), dimension(:)            :: param_names       !< an array of strings containing the names of the function parameters
    character(*), intent(in), dimension(:), optional  :: param_names_latex

    self%file_name = param_names(1)
    self%num_bins_string = param_names(2)

end subroutine InterpolatedFunction1DSetParamNames




! ---------------------------------------------------------------------------------------------
!> subroutine that initializes the spline arrays of the function
subroutine InterpolatedFunction1DInitialization(self)

    implicit none

    class(interpolated_function_1D) :: self
    integer                         :: i

    real(dl)                        :: a_start, a_end, a !< output parameters
    real(dl)                        :: f, fprime, fprimeprime

    real(dl), dimension(self%num_points) :: spline_workspace

!> TEST:
!write(*,*) " initializing xDE interp."
!pause
    ! initialize the calculation:
    call splini( spline_workspace, self%num_points )

    ! compute the numerical first derivative:
    call splder( self%y, self%yp, self%num_points, spline_workspace )
    self%yp = self%yp/self%grid_width

    ! compute the numerical second derivative:
    call splder( self%yp, self%ypp, self%num_points, spline_workspace )
    self%ypp = self%ypp/self%grid_width

    ! compute the numerical third derivative:
    call splder( self%ypp, self%yppp, self%num_points, spline_workspace )
    self%yppp = self%yppp/self%grid_width

    !> intialize the spline interpolations
    call GBD_spline_double( self%x, self%y,     self%num_points, self%y2    )
    call GBD_spline_double( self%x, self%yp,    self%num_points, self%yp2   )
    call GBD_spline_double( self%x, self%ypp,   self%num_points, self%ypp2  )
    call GBD_spline_double( self%x, self%yppp,  self%num_points, self%yppp2 )
    call GBD_spline_double( self%x, self%yint,  self%num_points, self%yint2 )

    self%is_function_loaded = .true.

end subroutine InterpolatedFunction1DInitialization


! ---------------------------------------------------------------------------------------------
!> Function that returns the value of the function in the scale factor.
function InterpolatedFunction1DValue( self, x, eft_cache )

    implicit none

    class(interpolated_function_1D)                    :: self        !< the base class
    real(dl), intent(in)                               :: x           !< the input scale factor
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache   !< the optional input EFTCAMB cache
    real(dl) :: InterpolatedFunction1DValue                           !< the output value

    if (.not. self%is_function_loaded) call self%initialize_function

    InterpolatedFunction1DValue = GBD_spline_val(x, self%x, self%y,self%y2, self%num_points)

end function InterpolatedFunction1DValue


! ---------------------------------------------------------------------------------------------
!> Function that returns the value of the first derivative, wrt scale factor, of the function.
function InterpolatedFunction1DFirstDerivative(self, x, eft_cache )

    implicit none

    class(interpolated_function_1D)                    :: self
    real(dl), intent(in)                               :: x
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache
    real(dl) :: InterpolatedFunction1DFirstDerivative

    if (.not. self%is_function_loaded) call self%initialize_function

    InterpolatedFunction1DFirstDerivative = GBD_spline_val(x, self%x, self%yp,self%yp2, self%num_points)

end function InterpolatedFunction1DFirstDerivative


! ---------------------------------------------------------------------------------------------
!> Function that returns the second derivative of the function.
function InterpolatedFunction1DSecondDerivative( self, x, eft_cache )

    implicit none

    class(interpolated_function_1D)                    :: self
    real(dl), intent(in)                               :: x
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache
    real(dl) :: InterpolatedFunction1DSecondDerivative

    if (.not. self%is_function_loaded) call self%initialize_function

    InterpolatedFunction1DSecondDerivative = GBD_spline_val( x, self%x,  self%ypp, self%ypp2, self%num_points )

end function InterpolatedFunction1DSecondDerivative



! ---------------------------------------------------------------------------------------------
!> Function that returns the second derivative of the function.
function InterpolatedFunction1DThirdDerivative( self, x, eft_cache )

    implicit none

    class(interpolated_function_1D)                    :: self
    real(dl), intent(in)                               :: x
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache
    real(dl) :: InterpolatedFunction1DThirdDerivative

    if (.not. self%is_function_loaded) call self%initialize_function

    InterpolatedFunction1DThirdDerivative = GBD_spline_val( x, self%x, self%yppp, self%yppp2, self%num_points )

end function InterpolatedFunction1DThirdDerivative


! ---------------------------------------------------------------------------------------------
!> Function that computes the integral for rho_DE
function InterpolatedFunction1DIntegral( self, x, eft_cache )

    implicit none

    class(interpolated_function_1D) :: self
    real(dl), intent(in) :: x
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache
    real(dl) :: InterpolatedFunction1DIntegral

    !other variables
    real(dl) :: int
    real(dl) :: tol = 1.d-6

    if (.not. self%is_function_loaded) call self%initialize_function

    int = GBD_spline_val( x, self%x, self%yint, self%yint2, self%num_points )

    InterpolatedFunction1DIntegral = x**2.d0 * Exp(3.d0*int)

end function InterpolatedFunction1DIntegral

! ---------------------------------------------------------------------------------------------
!> Subroutine that reads a Ini file looking for name of the file containing the function.
subroutine InterpolatedFunction1DInitFromCode( self, x, y )

    implicit none

    class(interpolated_function_1D)     :: self   !< the base class
    real(dl), dimension(:), intent(in) :: x,y

    !> intializing the number of points
    self%num_points = size(x)

    !> TEST
    !write(*,*) "initializing xDE from CODE"
    !write(*,*) 'x:',x
    !write(*,*) 'y:',y
    !pause

    ! allocate the other vectors:
    if ( allocated(self%x)    ) deallocate( self%x    )
    if ( allocated(self%y)    ) deallocate( self%y    )
    if ( allocated(self%yp)   ) deallocate( self%yp   )
    if ( allocated(self%ypp)  ) deallocate( self%ypp  )
    if ( allocated(self%yppp) ) deallocate( self%yppp )
    if ( allocated(self%yint) ) deallocate( self%yint )

    !> spline arrays
    if ( allocated(self%y2)    ) deallocate( self%y2    )
    if ( allocated(self%yp2)   ) deallocate( self%yp2   )
    if ( allocated(self%ypp2)  ) deallocate( self%ypp2  )
    if ( allocated(self%yppp2) ) deallocate( self%yppp2 )
    if ( allocated(self%yint2) ) deallocate( self%yint2 )

    allocate( self%x( self%num_points )    )
    allocate( self%y( self%num_points )    )
    allocate( self%yp( self%num_points )   )
    allocate( self%ypp( self%num_points )  )
    allocate( self%yppp( self%num_points ) )
    allocate( self%yint( self%num_points ) )

    !> spline arrays
    allocate( self%y2( self%num_points )    )
    allocate( self%yp2( self%num_points )   )
    allocate( self%ypp2( self%num_points )  )
    allocate( self%yppp2( self%num_points ) )
    allocate( self%yint2( self%num_points ) )

    !> Now pass the parameters
    self%x = x
    self%y = y

    self%grid_width = abs((self%x(self%num_points) - self%x(1))/self%num_points)
    self%is_function_loaded = .false.

end subroutine InterpolatedFunction1DInitFromCode



! ---------------------------------------------------------------------------------------------

end module EFTCAMB_interpolated_function_1D
!----------------------------------------------------------------------------------------
