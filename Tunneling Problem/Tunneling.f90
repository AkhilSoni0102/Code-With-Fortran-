program AkhilSoni_180122004
    implicit none
    real, dimension(:), allocatable:: Potential, x, prob_density
    complex, dimension(:), allocatable:: psi
    real:: x_max, dx, E, k, max_Val, min_Val, P, trans_probability
    integer:: nx, j, l
    complex:: i
    i = complex(0., 1.)
    nx = 1000
    x_max = 10.0
    dx = 0.01
    allocate(Potential(nx), x(nx), psi(nx), prob_density(nx))

    do j = 1, nx
        x(j) = (j-1)*dx
    enddo
    do j = 1, nx
        if(x(j) < 4) then
            Potential(j) = 0.0
        else if(x(j) <= 5) then
            Potential(j) = 9.0
        else
            Potential(j) = 0.0
        endif
    enddo
    open(unit = 10, file = 'trans.txt')
        E = 1.00
        do l = 1, 251
            k = sqrt(2*E)
            psi(1) = complex(1.,0.)
            psi(2) = EXP(-i*k*dx)
            do j = 2, nx-1
                psi(j+1)=(2-2 * (E-Potential(j)) * dx * dx) * psi(j)-psi(j-1)
            enddo
            prob_density = abs(psi)**2
            max_Val = 0.0
            min_Val = 100000.0
            do j = 600, nx
                max_Val = max(max_Val, prob_density(j))
                min_Val = min(min_Val, prob_density(j))
            enddo
            P = (max_Val + min_Val) / 2.0
            trans_probability = 2.0 / (1 + P)
            write(10, *) E, trans_probability
            if (l == 81) then
            open(unit = 11, file = 'prob_dens.txt')
                do j = 1, nx
                    write(11, *) x(j), prob_density(j)
                enddo
            close(11)
            endif
            E = E + 0.10
        enddo
    open(10)
end