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

!> @file 02_GBD_utils.f90
!! This file contains various generic algorithms that are useful to EFTCAMB.


!----------------------------------------------------------------------------------------
!> This module contains various generic algorithms that are useful to EFTCAMB.

!> @author Bin Hu, Marco Raveri
!> @author Alex Zucca

module GBD_Utils
    use precision
    implicit none

contains

subroutine GBD_spline_double(x,y,n,d2)
    integer, intent(in) :: n
    integer, parameter :: dp=KIND(1.d0)
    real(dp), intent(in) :: x(n), y(n)
    real(dp), intent(out) :: d2(n)
    real(dp), dimension(:), allocatable :: u
    integer i
    real(dp) xp,sig,xxdiv,d1l,d1r

    allocate(u(1:n-1))

    d2(1)=0._dp
    u(1)=0._dp

    d1r= (y(2)-y(1))/(x(2)-x(1))
    do i=2,n-1
        d1l=d1r
        d1r=(y(i+1)-y(i))/(x(i+1)-x(i))
        xxdiv=1._dp/(x(i+1)-x(i-1))
        sig=(x(i)-x(i-1))*xxdiv
        xp=1._dp/(sig*d2(i-1)+2._dp)
        d2(i)=(sig-1._dp)*xp
        u(i)=(6._dp*(d1r-d1l)*xxdiv-sig*u(i-1))*xp
    end do

    d2(n)=0._dp

    do i=n-1,1,-1
        d2(i)=d2(i)*d2(i+1)+u(i)
    end do

    deallocate(u)
end subroutine GBD_spline_double

! General routine for cubic spline interpolation (see NR)
real(dl) function GBD_spline_val(x, xv, yv, y2, n)

    integer, intent(in) :: n
    real(dl), intent(in) :: x
    real(dl), intent(in) :: xv(n), yv(n), y2(n)

    integer :: kh,kl,kn
    real(dl) :: h,a,b,c,d

    ! Extrapolate if value is above of below interval
    if(x < xv(1)) then
        h = xv(2) - xv(1)
        a = (yv(2) - yv(1)) / h
        GBD_spline_val = (a - h * y2(2) / 6) * (x - xv(1)) + yv(1)
    else if(x > xv(n)) then
        h = xv(n) - xv(n-1)
        a = (yv(n) - yv(n-1)) / h
        GBD_spline_val = (a + h * y2(n-1) / 6) * (x - xv(n)) + yv(n)
    else
    ! Bisection to find correct interval
        kh = n
        kl = 1
        do while(kh - kl > 1)
            kn = (kh + kl) / 2
            if(xv(kn) > x) then
                kh = kn
            else
                kl = kn
            end if
        end do

        ! Set up constants (a la NR)
        h = xv(kh) - xv(kl)

        a = (xv(kh) - x) / h
        b = (x - xv(kl)) / h
        c = (a**3 - a)* h**2 / 6
        d = (b**3 - b)* h**2 / 6

        GBD_spline_val = (a*yv(kl) + b*yv(kh) + c*y2(kl) + d*y2(kh))

    end if
end function GBD_spline_val


! Adapted from NR to compute derivatives

! General routine for cubic spline interpolation (see NR)
real(dl) function GBD_spline_1der(x, xv, yv, y2, n)

    integer, intent(in) :: n
    real(dl), intent(in) :: x
    real(dl), intent(in) :: xv(n), yv(n), y2(n)

    integer :: kh,kl,kn
    real(dl) :: h,a,b,c,d

    ! Extrapolate if value is above of below interval
    ! AZ: I NEED TO MODIFY THIS
    if(x < xv(1)) then
        h = xv(2) - xv(1)
        a = (yv(2) - yv(1)) / h
        !GBD_spline_1der = (a - h * y2(2) / 6) * (x - xv(1)) + yv(1)
        GBD_spline_1der =  - h * y2(2) / 6
        !WRITE(*,*) "spline_der: extrapolation not implemented yet"
        !write(*,*) x, xv(1)
        !stop
    else if(x > xv(n)) then
        h = xv(n) - xv(n-1)
        !a = (yv(n) - yv(n-1)) / h
        !GBD_spline_1der = (a + h * y2(n-1) / 6) * (x - xv(n)) + yv(n)
        GBD_spline_1der =  h * y2(n-1) / 6
        !WRITE(*,*) "spline_1der: extrapolation not implemented yet"
        !write(*,*) x, xv(n)
        !stop
    else
    ! Bisection to find correct interval
        kh = n
        kl = 1
        do while(kh - kl > 1)
            kn = (kh + kl) / 2
            if(xv(kn) > x) then
                kh = kn
            else
                kl = kn
            end if
        end do

        ! Set up constants (a la NR)
        h = xv(kh) - xv(kl)

        a = (xv(kh) - x) / h
        b = (x - xv(kl)) / h
        c = (a**3 - a)* h**2 / 6
        d = (b**3 - b)* h**2 / 6



        !this is the derivative
        GBD_spline_1der = (yv(kh) - yv(kl))/h - (3.d0*a**2 -1.d0)/6.d0 * h * y2(kl) +(3.d0*b**2 -1.d0)/6.d0*&
                h * y2(kh)

    end if

end function GBD_spline_1der

real(dl) function GBD_spline_2der(x, xv, yv, y2, n)

    integer, intent(in) :: n
    real(dl), intent(in) :: x
    real(dl), intent(in) :: xv(n), yv(n), y2(n)

    integer :: kh,kl,kn
    real(dl) :: h,a,b,c,d

    ! Extrapolate if value is above of below interval
    ! AZ: I NEED TO MODIFY THIS
    if(x < xv(1)) then
        h = xv(2) - xv(1)
        a = (yv(2) - yv(1)) / h
        !GBD_spline_2der = (a - h * y2(2) / 6) * (x - xv(1)) + yv(1)
        GBD_spline_2der = 0.d0
        !WRITE(*,*) "spline_2der: extrapolation not implemented yet"
        !write(*,*) x, xv(1)
        !stop
    else if(x > xv(n)) then
        h = xv(n) - xv(n-1)
        a = (yv(n) - yv(n-1)) / h
        !GBD_spline_2der = (a + h * y2(n-1) / 6) * (x - xv(n)) + yv(n)
        GBD_spline_2der = 0.d0
        !WRITE(*,*) "spline_2der: extrapolation not implemented yet"
        !write(*,*) x, xv(n)
        !stop
    else
        ! Bisection to find correct interval
        kh = n
        kl = 1
        do while(kh - kl > 1)
            kn = (kh + kl) / 2
            if(xv(kn) > x) then
                kh = kn
            else
                kl = kn
            end if
        end do

        ! Set up constants (a la NR)
        h = xv(kh) - xv(kl)

        a = (xv(kh) - x) / h
        b = (x - xv(kl)) / h
        c = (a**3 - a)* h**2 / 6
        d = (b**3 - b)* h**2 / 6

        !this is the derivative
        !GBD_spline_1der = (yv(kh) - yv(kl))/h - (3.d0*a**2 -1.d0)/6.d0 * h * y2(kl) +(3.d0*b**2 -1.d0)/6.d0*&
        !h * y2(kh)

        ! this is the second derivative
        GBD_spline_2der = a*y2(kl)+b*y2(kh)

    end if

end function GBD_spline_2der



function GBD_rombint(f,a,b,tol)
!use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
    implicit none
    integer, parameter :: MAXITER=20
    integer, parameter :: MAXJ=5
    dimension g(MAXJ+1)
    real(dl) f
    external f
    real(dl) :: GBD_rombint
    real(dl), intent(in) :: a,b,tol
    integer :: nint, i, k, jmax, j
    real(dl) :: h, gmax, error, g, g0, g1, fourj
    !

    h=0.5d0*(b-a)
    gmax=h*(f(a)+f(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
        go to 40
    !  Calculate next trapezoidal rule approximation to integral.
    g0=0._dl
    do 20 k=1,nint
        g0=g0+f(a+(k+k-1)*h)
20  continue
    g0=0.5d0*g(1)+h*g0
    h=0.5d0*h
    nint=nint+nint
    jmax=min(i,MAXJ)
    fourj=1._dl
    do 30 j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4._dl*fourj
        g1=g0+(g0-g(j))/(fourj-1._dl)
        g(j)=g0
        g0=g1
30  continue
    if (abs(g0).gt.tol) then
        error=1._dl-gmax/g0
    else
        error=gmax
    end if
    gmax=g0
    g(jmax+1)=g0
    go to 10
40  GBD_rombint=g0
    if (i.gt.MAXITER.and.abs(error).gt.tol)  then
        write(*,*) 'Warning: GBD_Rombint failed to converge; '
        write (*,*)'integral, error, tol:', GBD_rombint,error, tol
    end if

end function GBD_rombint

end module GBD_Utils
