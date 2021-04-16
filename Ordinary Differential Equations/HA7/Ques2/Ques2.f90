!The energy of the simple harmonic oscillator is:
!                                                E = p^2/2m + 0.5*k*x^2
!Here m is the mass, k is the spring force constant, and p = mx˙ is the momentum. 
!Let us set m = k=1, so the angular frequency, ω, and period, T, are given by ω = sqrt(k/m) = 1, and T = 2*pi/ω = 2*pi Hence, 2E = p^2 + x^2 is a constant, and a plot of p vs x will be a circle of radius √2E.
!We will solve the Newton’s equation of motion ¨x = −x which will be written as two first order differential equations
!                                               x˙ = p and p = −x. 
!(1) Numerically integrate these equations using the Euler and RK4 methods.
!Use initial conditions, x =1 and p =0. Hence 2E =1. Use a time step h = 0.02T
!Use 200 time steps in total. Plot 2E vs t/period, x vs p and x vs t/period.
!In the figure below, results obtained using Euler’s, RK2 and RK4 methods are shown.
!But you need to show the results only for Euler’s and RK4 methods. You
!will submit a zip file, q2.zip, containing the code(s), and six output files, 3 for each method.


module precision_use
    implicit none 
    INTEGER, PARAMETER:: dp = SELECTED_REAL_KIND(12)
end module precision_use

program Ques_2
    use precision_use 
    implicit none
    real(kind = dp) :: xp(2)
    real(kind = dp) :: t, h, pi
    integer :: i 

    pi = acos(-1.0)
    t = 0.0
    h = 0.02*2.*pi
    xp(1) = 1 
    xp(2) = 0 
    OPEN(unit=20, file = 'Euler.txt')
    do i = 1,600
        xp = Euler(xp)
        WRITE(20,*) t/(2.*pi), xp(1), xp(2), xp(1)*xp(1) + xp(2)*xp(2)
        t = t + h
    end do
    CLOSE(20)
    t = 0.0
    xp(1) = 1 
    xp(2) = 0
    OPEN(unit=10, file = 'RK4.txt')
    do i = 1,600
        xp = RK4(xp)
        WRITE(10,*) t/(2.*pi), xp(1), xp(2), xp(1)*xp(1) + xp(2)*xp(2)
        t = t + h
    end do
    CLOSE(10)

contains

    function f_dash(xp_0) result (xp_1)
        real(kind = dp) :: xp_0(2), xp_1(2)

        xp_1(1) = xp_0(2) 
        xp_1(2) = -xp_0(1)

    end function

    function RK4(xp_0) result (xp_1)
        real(kind = dp) :: xp_0(2), xp_1(2)
        real(kind = dp) :: s1(2), s2(2), s3(2), s4(2)
        
        s1 = h*f_dash(xp_0)
        s2 = h*f_dash(xp_0 + s1/2)
        s3 = h*f_dash(xp_0 + s2/2)
        s4 = h*f_dash(xp_0 + s3)

        xp_1 = xp_0 + ((s1 + 2*s2 + 2*s3 + s4)/(6.))

    end function

    function Euler(xp_0) result (xp_1)
        real(kind = dp) :: xp_0(2), xp_1(2)
        xp_1 = xp_0 + h*f_dash(xp_0)
    end function
end program Ques_2