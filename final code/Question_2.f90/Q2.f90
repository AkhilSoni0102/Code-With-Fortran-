program Question_2
    implicit none
    real :: Y(2)
    real :: a, b, time, h, T, pi, E
    integer :: N, i 

    pi = acos(-1.0)
    T = 2*pi
    
    N = 200

    a = 0
    time = a
    b = 8*pi
    h = (b-a)/N

    Y(1) = 1 ! x(0)
    Y(2) = 0 ! p(0)

    OPEN(unit=10, file = 'rk4')
    do i = 1,N
        Y = RK4_method(time, Y)
        E = Y(1)*2 + Y(2)*2
        WRITE(10,*) Y(1), Y(2), E, time/T
        
        time = time + h
    end do
    CLOSE(10)
    
    time = a
    OPEN(unit=20, file = 'euler')
    do i = 1,N
        Y = Euler_method(time, Y)
        E = Y(1)*2 + Y(2)*2
        WRITE(20,*) Y(1), Y(2), E, time/T
        
        time = time + h
    end do
    CLOSE(20)

    contains
    function f(x, Yvec) result (fvec)
        real :: x
        real :: Yvec(2), fvec(2)

        fvec(1) = Yvec(2) 
        fvec(2) = -Yvec(1)
    end function

    function RK4_method(x, Y_n) result (Y_nplus1)
        real :: x
        real :: Y_n(2), Y_nplus1(2)
        real :: k1(2), k2(2), k3(2), k4(2)
        
        k1 = h*f(x, Y_n)
        k2 = h*f(x + h/2, Y_n + k1/2)
        k3 = h*f(x + h/2, Y_n + k2/2)
        k4 = h*f(x + h, Y_n + k3)

        Y_nplus1 = Y_n + (k1 + 2*k2 + 2*k3 + k4)/6

    end function

    function Euler_method(x, Y_n) result (Y_nplus1)
        real :: x
        real :: Y_n(2), Y_nplus1(2)
        
        Y_nplus1 = Y_n + h*f(x, Y_n)

    end function

end program Question_2