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

!> @file 04p10_reconstructed_fit_parametrizations_1D.f90
!! This file contains the definition of the reconstructed fit parametrization for the
!! dark energy density reconstruction, inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of power law parametrization for the DE density.

!> @author Bin Hu, Marco Raveri, Simone Peirone

!> @author Alex Zucca: azucca@sfu.ca


module EFTCAMB_power_law_DE_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public power_law_DE_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: power_law_DE_parametrization_1D

        real(dl) :: wDE

    contains

        ! utility functions:
        procedure :: set_param_number      => PowerLawDEParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
        procedure :: init_parameters       => PowerLawDEParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => PowerLawDEParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => PowerLawDEParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => PowerLawDEParametrized1DValue               !< function that returns the value of the Taylor expansion.
        procedure :: first_derivative      => PowerLawDEParametrized1DFirstDerivative     !< function that returns the first derivative of the Taylor expansion.
        procedure :: second_derivative     => PowerLawDEParametrized1DSecondDerivative    !< function that returns the second derivative of the Taylor expansion.
        procedure :: third_derivative      => PowerLawDEParametrized1DThirdDerivative     !< function that returns the third derivative of the Taylor expansion.
        procedure :: integral              => PowerLawDEParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type power_law_DE_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor expansion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine PowerLawDEParametrized1DSetParamNumber( self )

        implicit none

        class(power_law_DE_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 1

    end subroutine PowerLawDEParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine PowerLawDEParametrized1DInitParams( self, array )

        implicit none

        class(power_law_DE_parametrization_1D)                      :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)      :: array  !< input array with the values of the parameters.

        self%wDE    = array(1)


    end subroutine PowerLawDEParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine PowerLawDEParametrized1DParameterValues( self, i, value )

        implicit none

        class(power_law_DE_parametrization_1D) :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%wDE
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine PowerLawDEParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine PowerLawDEParametrized1DFeedback( self, print_params )

        implicit none

        class(power_law_DE_parametrization_1D)  :: self         !< the base class
        logical, optional                       :: print_params !< optional flag that decised whether to print numerical values
                                                         !! of the parameters.

        integer                                     :: i
        real(dl)                                    :: param_value
        character(len=EFT_names_max_length)         :: param_name
        logical                                     :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,*)     'Power Law DE parametrization: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine PowerLawDEParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function PowerLawDEParametrized1DValue( self, x, eft_cache )

        implicit none

        class(power_law_DE_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional  :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawDEParametrized1DValue                        !< the output value

        !> value of the function
        PowerLawDEParametrized1DValue = x**(-3._dl * (1._dl + self%wDE))

        if (x == 0.d0) then
            PowerLawDEParametrized1DValue = 1.d3
        end if

    end function PowerLawDEParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function PowerLawDEParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(power_law_DE_parametrization_1D)                  :: self        !< the base class
        real(dl), intent(in)                                    :: x           !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional      :: eft_cache   !< the optional input EFTCAMB cache
        real(dl) :: PowerLawDEParametrized1DFirstDerivative                    !< the output value


        !> firt derivative
        PowerLawDEParametrized1DFirstDerivative = -(3._dl + 3._dl * self%wDE) * x**(-3._dl *  self%wDE-4._dl)


    end function PowerLawDEParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function PowerLawDEParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(power_law_DE_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional  :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawDEParametrized1DSecondDerivative      !< the output value

        !> second derivative
        PowerLawDEParametrized1DSecondDerivative =  (3._dl + 3._dl *self%wDE) * (3._dl * self%wDE + 4._dl ) *  x**(-3._dl *self%wDE - 5._dl)

    end function PowerLawDEParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function PowerLawDEParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(power_law_DE_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                                        :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional          :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawDEParametrized1DThirdDerivative                      !< the output value


        !> third derivative
        PowerLawDEParametrized1DThirdDerivative = - (3._dl  + 3._dl * self%wDE) * (3._dl * self%wDE + 4._dl ) * &
                                                    (3._dl * self%wDE + 5._dl) * x**(- 3._dl *self%wDE - 6._dl)


    end function PowerLawDEParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function PowerLawDEParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(power_law_DE_parametrization_1D)              :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional  :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PowerLawDEParametrized1DIntegral                     !< the output value

        write(*,*) "WARNING: the function Power_Law_DE%Integral is not supposed to be called. Returning 0"
        PowerLawDEParametrized1DIntegral = 0._dl

    end function PowerLawDEParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_power_law_DE_parametrizations_1D

!----------------------------------------------------------------------------------------
