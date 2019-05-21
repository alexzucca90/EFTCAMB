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

!> @file 04p10_hyperbolic_tangent_parametrizations_1D.f90
!! This file contains the definition of the reconstructed fit parametrization for the
!! dark energy density reconstruction, inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Taylor expansion parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

!> @author Alex Zucca: azucca@sfu.ca


module EFTCAMB_hyperbolic_tangent_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public hyperbolic_tangent_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: hyperbolic_tangent_parametrization_1D

        real(dl) :: A
        real(dl) :: B
        real(dl) :: C



    contains

        ! utility functions:
        procedure :: set_param_number      => HyperbolicTangentParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
        procedure :: init_parameters       => HyperbolicTangentParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => HyperbolicTangentParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => HyperbolicTangentParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => HyperbolicTangentParametrized1DValue               !< function that returns the value of the Taylor expansion.
        procedure :: first_derivative      => HyperbolicTangentParametrized1DFirstDerivative     !< function that returns the first derivative of the Taylor expansion.
        procedure :: second_derivative     => HyperbolicTangentParametrized1DSecondDerivative    !< function that returns the second derivative of the Taylor expansion.
        procedure :: third_derivative      => HyperbolicTangentParametrized1DThirdDerivative     !< function that returns the third derivative of the Taylor expansion.
        procedure :: integral              => HyperbolicTangentParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type hyperbolic_tangent_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor expansion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine HyperbolicTangentParametrized1DSetParamNumber( self )

        implicit none

        class(hyperbolic_tangent_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 3

    end subroutine HyperbolicTangentParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine HyperbolicTangentParametrized1DInitParams( self, array )

        implicit none

        class(hyperbolic_tangent_parametrization_1D)             :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%A     = array(1)
        self%B     = array(2)
        self%C     = array(3)

    end subroutine HyperbolicTangentParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine HyperbolicTangentParametrized1DParameterValues( self, i, value )

        implicit none

        class(hyperbolic_tangent_parametrization_1D)       :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%A
            case(2)
                value = self%B
            case(3)
                value = self%C
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop( 'EFTCAMB error' )
        end select

    end subroutine HyperbolicTangentParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine HyperbolicTangentParametrized1DFeedback( self, print_params )

        implicit none

        class(hyperbolic_tangent_parametrization_1D) :: self         !< the base class
        logical, optional                           :: print_params !< optional flag that decised whether to print numerical values
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

        write(*,*)     'Reconstruction Fit parametrization: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine HyperbolicTangentParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function HyperbolicTangentParametrized1DValue( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional  :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentParametrized1DValue                 !< the output value

        !> value of the function
        HyperbolicTangentParametrized1DValue = self%A*tanh(self%B*(x-self%C)) + 1._dl - self%A*tanh(self%B*(1._dl-self%C))

    end function HyperbolicTangentParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function HyperbolicTangentParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_parametrization_1D)                    :: self      !< the base class
        real(dl), intent(in)                                            :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional              :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentParametrized1DFirstDerivative                    !< the output value

        !> first derivative
        HyperbolicTangentParametrized1DFirstDerivative = self%A * self%B * (1._dl - tanh(self%B*(x-self%C))**2)




    end function HyperbolicTangentParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function HyperbolicTangentParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentParametrized1DSecondDerivative      !< the output value



        !> second derivative
        HyperbolicTangentParametrized1DSecondDerivative = -2._dl*self%A*self%B**2*tanh(self%B*(x-self%C))*(1._dl-tanh(self%B*(x-self%C))**2)


    end function HyperbolicTangentParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function HyperbolicTangentParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional  :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentParametrized1DThirdDerivative       !< the output value

        !> third derivative
        HyperbolicTangentParametrized1DThirdDerivative = -2._dl*self%A*self%B**3*(1._dl-tanh(self%B*(x-self%C))**2)**2 &
                                & +4._dl*self%A*self%B**3*tanh(self%B*(x-self%C))**2*(1._dl-tanh(self%B*(x-self%C))**2)


    end function HyperbolicTangentParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function HyperbolicTangentParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentParametrized1DIntegral              !< the output value

        !> Integral
        HyperbolicTangentParametrized1DIntegral = 0._dl

    end function HyperbolicTangentParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_hyperbolic_tangent_parametrizations_1D

!----------------------------------------------------------------------------------------
