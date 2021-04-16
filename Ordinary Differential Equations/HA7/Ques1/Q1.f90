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