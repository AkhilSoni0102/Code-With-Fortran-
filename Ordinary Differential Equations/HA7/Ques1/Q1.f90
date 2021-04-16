!The radial Schr¨odinger equation for the central potential V (r) is given by
!                       [d^2/dt^2 + (2/r)d/dr]R(r)] + 2*µ/h^2[E + V(r) - (l(l+1)h^2)/2*µ*r^2]R(r) = 0
! Here, µ is the reduced mass of the system, l is the orbital-angular momentum quantum
!number, and R(r) is the radial wave function. The above equation, in atomic units,
!for the ground state (l=0) of the hydrogen atom can be written as
!                       [d^2/dt^2 + (2/r)d/dr]R(r)] + 2[E_0 + 1/r]R(r) = 0
!Write a Fortran program to solve the above equation using RK4 to find E0 with
!following starting values: 
!                       R(r = 0.0005) = 0.000001, R0 (r = 0.0005) = −1000.0. 
!The r grid will be from 0.0005 unit to 5 unit with 10000 points. The code will be for a
!range of E values, −0.6 ≤ E ≤ −0.4, with ∆E=0.01. For finding the correct value of
!E, plot R(r) and the radial distribution function, |rR(r)|^2 , against r and check their convergence with respect to E. At the E where both R(r) and |rR(r)|^2
!behave properly against r, that will be the answer. You will submit one zip file, q1.zip, containing the
! following files: the code, and 5 different output files containing the data for R(r) and |rR(r)|^2 vs r. 
!You will show the results for E= -0.60, -0.55, -0.50, -0.45 and -0.40.

module precision_use
    implicit none 
    INTEGER, PARAMETER:: dp = SELECTED_REAL_KIND(12)
end module precision_use

program Ques_1
    use precision_use
    implicit none
    real(kind = dp) :: R(2)
    real(kind = dp) :: t, h, E, end
    integer :: N, i
    N = 10000
    t = 0.0005
    end = 5.0
    h = (end-t)/N
    R(1) = 0.000001 
    R(2) = -1000.0 
    E = -0.60
    OPEN(unit = 10, file = 'E1.txt')
    do i = 1,N
        R = RK4(t, R)
        WRITE(10,*) t, R(1), abs(t*t*R(1)*R(1))
        t = t + h
    end do
    CLOSE(10)
    
    R(1) = 0.000001 
    R(2) = -1000.0 
    E = -0.55
    t = 0.0005
    
    OPEN(unit = 10, file = 'E2.txt')
    do i = 1,N
        R = RK4(t, R)
        WRITE(10,*) t, R(1), abs(t*t*R(1)*R(1))
        t = t + h
    end do
    CLOSE(10)
    
    R(1) = 0.000001 
    R(2) = -1000.0 
    E = -0.50
    t = 0.0005
    OPEN(unit = 10, file = 'E3.txt')
    do i = 1,N
        R = RK4(t, R)
        WRITE(10,*) t, R(1), abs(t*t*R(1)*R(1))
        t = t + h
    end do
    CLOSE(10)
    
    R(1) = 0.000001 
    R(2) = -1000.0 
    E = -0.45
    t = 0.0005
    OPEN(unit = 10, file = 'E4.txt')
    do i = 1,N
        R = RK4(t, R)
        WRITE(10,*) t, R(1), abs(t*t*R(1)*R(1))
        t = t + h
    end do
    CLOSE(10)
    
    R(1) = 0.000001 
    R(2) = -1000.0 
    E = -0.40
    t = 0.0005
    OPEN(unit = 10, file = 'E5.txt')
    do i = 1,N
        R = RK4(t, R)
        WRITE(10,*) t, R(1), abs(t*t*R(1)*R(1))
        t = t + h
    end do
    CLOSE(10)
    
contains
    function f(x, R_0) result (R_1)
        real(kind = dp) :: x
        real(kind = dp) :: R_0(2), R_1(2)
        R_1(1) = R_0(2) 
        R_1(2) = (-2./x)*R_0(2) - ((2./x) + 2.*E)*R_0(1)
    end function

    function RK4(x, R_0) result (R_1)
        real(kind = dp) :: x
        real(kind = dp) :: R_0(2), R_1(2)
        real(kind = dp) :: s1(2), s2(2), s3(2), s4(2)
        s1 = h*f(x, R_0)
        s2 = h*f(x + h/2, R_0 + s1/2)
        s3 = h*f(x + h/2, R_0 + s2/2)
        s4 = h*f(x + h, R_0 + s3)
        R_1 = R_0 + ((s1 + 2*s2 + 2*s3 + s4)/(6.))
    end function
end program Ques_1