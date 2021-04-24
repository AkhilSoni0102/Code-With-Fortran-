module Constant
    implicit none
    integer,parameter:: nt = 5000
    integer,parameter:: steps = 100
    integer,parameter:: n=256
    real,parameter:: pi = acos(-1.0)
    real,parameter:: alpha = 20.0
    real,parameter:: x_0 = -0.5
    real,parameter:: p_0 = 50.0
    real,parameter:: Mass_of_e = 14500.0
    real,parameter:: dx = 0.02
    real,parameter:: dt = 0.1
    complex,parameter:: ai = (0.0,1.0)
    character*64 outran
end module Constant

program TDSE
    use Constant
    implicit none
    integer:: i, t
    real:: x(n), vx(n), psi_psi(n), akx(n), akx_2(n)
    real:: wnorm, x_min, dk, p_max , p_min 
    complex:: psi_0(n), psi_1(n), ftmp(n), psi_2(n)
    x_min=-2.0

    do i = 1, n
        x(i) = x_min + (i - 1) * dx
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

    open(unit = 11, file = 'Wave1.out')
    wnorm = (2.0 * alpha / pi)**0.25

    do i = 1, n
        psi_0(i) = wnorm * (exp(ai * p_0 * (x(i) - x_0))) * exp(-alpha * (x(i) - x_0)**2)
        psi_psi = psi_0(i) * conjg(psi_0(i))
        write(11,*) x(i), psi_psi
    enddo
    close(11)

    dk = 2.0 * pi / (float(n) * dx)
    p_max = pi / dx
    p_min = (-1)*p_max 
    do i = 1, (n / 2) - 1
        akx(i) = 2.0 * pi * float(i - 1) / (float(n) * dx)
    enddo

    do i = n / 2, n
        akx(i) = 2.0 * pi * float(i - 1 - n) / (float(n) * dx)
    enddo

    do i = 1, n
        akx_2(i) = akx(i) * akx(i)
    enddo

    
    do i = 1, n
        ftmp(i) = psi_0(i)
    end do

    call Descrete_FT(ftmp, x, akx, akx_2)
    open(unit=12,file='Wave2.out')
    do i = 1, n
        psi_1(i) = psi_0(i) - ai * dt * (-0.5 * ftmp(i) / Mass_of_e + vx(i) * psi_0(i))
        write(12,*) x(i), abs(psi_1(i))**2
    end do
    close(12)

    do t = 2, nt
        do i = 1, n
            ftmp(i) = psi_1(i)
        end do
        call Descrete_FT(ftmp, x, akx, akx_2)
        do i = 1, n
            psi_2(i) = psi_0(i) - 2.0 * ai * dt * (-0.5 * ftmp(i) / Mass_of_e + vx(i) * psi_1(i))
            psi_psi = psi_2(i) * conjg(psi_2(i))
            if (mod(t, steps) == 0) then
                write(outran,1161) t, p_0
                1161 format('wave.', i4.4,'.', f5.2)
                open(unit = 16, file = outran, form = 'formatted', status = 'unknown')
                write(16,*) x(i) , psi_psi
            endif
        end do
        do  i = 1, n
            psi_0(i) = psi_1(i)
            psi_1(i) = psi_2(i)
        end do
    end do
end
subroutine Descrete_FT(psi, x_Grid, akx, akx_2)
    use Constant
    implicit none
    integer:: i, j
    real:: akx(n), akx_2(n), x_Grid(n)
    complex phi(n),psi(n)
    complex e
    e = (0.0, 1.0)
    do i = 1, n
        phi(i) = 0.0
        do j = 1, n
            phi(i) = phi(i) + psi(j) * exp(-e * akx(i) * x_Grid(j))
        enddo
    enddo

    phi = phi / sqrt(float(n))
    do i = 1, n 
        phi(i) = (-1) * akx_2(i) * phi(i)
    end do
    do j = 1, n 
        psi(j) = 0.0
        do i=1,n
            psi(j) = psi(j) + phi(i) * exp(e * akx(i) * x_Grid(j))
        enddo
    enddo
    psi = psi / sqrt(float(n))
    return
end