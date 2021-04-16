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