module Constant
    implicit none
    integer,parameter:: nt = 5000
    integer,parameter:: istep=100
    integer,parameter:: n=256
    real,parameter:: pi = acos(-1.0)
    real,parameter:: alpha = 20.0
    real,parameter:: x_0 = -0.5
    real,parameter:: p_0 = 50.0
    real,parameter:: Mass_of_e = 14500.0
    real,parameter:: dx = 0.02
    real,parameter:: dt = 0.1
    complex,parameter:: ai = (0.0,1.0)
    character*64 File
end module Constant


program TDSE
    use Constant
    implicit none
    integer:: i, it
    real:: x(n), vx(n), psi_psi(n)
    real:: wnorm, x_min
    complex:: psi_0(n), psi_1(n), ftmp(n), d2y(n), psi_2(n)
    x_min = -2.0
    do i = 1, n
        x(i) = x_min + (i-1) * dx
    enddo

    open(unit = 10, file = 'Potential.out')
    do i = 1, n
        if(x(i) < 0.0)then
            vx(i) = 0.0
        else
            vx(i) = 0.1
        endif
        write(10,*) x(i), vx(i)
    enddo
    close(10)

    wnorm = (2.0 * alpha / pi)**0.25

    open(unit=11,file='wave1.out')
    do i = 1, n
        psi_0(i) = wnorm * (exp(ai * p_0 * (x(i) - x_0))) * exp(-alpha * (x(i) - x_0)**2)
        psi_psi = psi_0(i) * conjg(psi_0(i))
        write(11,*) x(i), psi_psi
    enddo
    close(11)

    do i=1,n
        ftmp(i) = psi_0(i)
    end do
    call sd2y(x, n, ftmp, d2y, dx)

    open(unit=12,file='wave2.out')
    do i = 1, n
        psi_1(i) = psi_0(i) - ai * dt * (-0.5 * d2y(i) / Mass_of_e + vx(i) * psi_0(i))
        write(12,*) x(i), abs(psi_1(i))**2
    end do
    close(12)

    do it = 2, nt
        do i= 1, n
            ftmp(i) = psi_1(i)
        end do

        call sd2y(x, n, ftmp, d2y, dx)
        do i= 1, n
            psi_2(i) = psi_0(i) - 2.0 * ai * dt * (-0.5 * d2y(i) / Mass_of_e + vx(i) * psi_1(i))
            psi_psi = psi_2(i) * conjg(psi_2(i))
            if (mod(it, istep) == 0) then
                write(File, 13) it, p_0
                13 format('Wave.', i4.4, '.', f5.2)
                open(unit=14, file = File, form = 'formatted', status = 'unknown')
                write(14,*) x(i), psi_psi
            endif
        end do

        do  i= 1, n
            psi_0(i) = psi_1(i)
            psi_1(i) = psi_2(i)
        end do
    end do
end

subroutine sd2y(x, n, fn, d2y, dx)
    implicit none
    integer:: n, i
    real:: dx, x(n)
    complex:: fn(n) , d2y(n), d2x(n), d1x(n)
    x = x

    d2y = 0
    do i = 2, n-1
        d2y(i) = (fn(i+1) - 2 * fn(i) + fn(i-1)) / dx**2 !Central Difference 
    enddo
    d2y(1) = (fn(3) - 2.0 * fn(2) + fn(1)) / dx**2  !Forward Difference
    d2y(n) = (fn(n) - 2.0 * fn(n-1) + fn(n-2)) / dx**2  !Backward Difference
    
    d1x=0.0
    do i=2,n-1
        d1x(i) = (fn(i+1) - fn(i-1)) / (2.0 * dx) !Central Difference 
    enddo
    d1x(1) = (-3.0 * fn(1) + 4.0 * fn(2) - fn(3)) / (2.0 * dx) !Forward Difference
    d1x(n) = (-3.0 * fn(n) + 4.0 * fn(n - 1) - fn(n - 2)) / (-2.0 * dx) !Backward Difference
    
    d2x=0.0
    do i=2,n-1
        d2x(i)=(d1x(i+1)-d1x(i-1))/(2.0*dx) !Central Difference 
    enddo
    d2x(1)=(-3.0*d1x(1)+4.0*d1x(2)-d1x(3))/(2.0*dx)  !Forward Difference
    d2x(n)=(-3.0*d1x(n)+4.0*d1x(n-1)-d1x(n-2))/(-2.0*dx)  !Backward Difference

    return
end
