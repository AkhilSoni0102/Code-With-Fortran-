PROGRAM Akhil_Soni_180122004
    IMPLICIT NONE
    real, PARAMETER:: Pi = acos(-1.0)
    real, ALLOCATABLE:: potential(:), x(:)
    double PRECISION, ALLOCATABLE:: h_mat(:,:), eigen_value(:), W(:)
    double precision:: energy_difference
    real:: x_max, x_min, dx, h_bar, t_constant, b, k, x_0, c, mu
    integer:: nx, l_W, i, i_fail
    write(*, fmt = '(/A)', advance = 'NO') "Enter the no. of grid points"
    read*, nx
    x_max = 10.0
    x_min = 0.0
    dx = (x_max - x_min)/(nx - 1)
    
    ALLOCATE(potential(nx), h_mat(nx,nx), eigen_value(nx), x(nx), W(64*nx))
    do i = 1, nx
        x(i) = x_min + (i-1)*dx
    end do
    h_bar = 1.00
    k = 0.08
    x_0 = 5.00
    b = 0.06
    c = 1.37
    mu = 4668 
    
    t_constant = (h_bar*h_bar) / (2.0 * mu * dx*dx)
    do i = 1, nx
        potential(i) = 0.5 * k * (x(i)-x_0)*(x(i)-x_0) + b * EXP(-c * (x(i)-x_0)*(x(i)-x_0) ) 
    end do
    
    h_mat = 0.0
    do i = 1, nx-1
        h_mat(i, i) = potential(i) + 2.0 * t_constant
        h_mat(i, i+1) = -t_constant
        h_mat(i+1, i) = -t_constant
    end do 
    h_mat(nx, nx) = potential(nx) + 2.0 * t_constant
    
    l_W = 64*nx
    CALL dsyev('V', 'U', nx, h_mat, nx, eigen_value, W, l_W, i_fail)

    open(unit = 10, file = 'potential.txt')
        do i = 1, nx 
            write(10, 20) x(i), potential(i)
            20 format(f0.6, 3x, f0.6)
        end do
    close(10)

    open(unit = 20, file = 'eigenstate.txt')
        do i = 1, nx 
            WRITE(20, 30) x(i), h_mat(i, 1)/sqrt(dx)
            30 FORMAT(f0.6, 3x, f0.6)
        ENDDO
    CLOSE(20)

    energy_difference = 6579689.75 * (eigen_value(2) - eigen_value(1))
    write(*, *) "The difference : ", energy_difference, "GHz"

END PROGRAM Akhil_Soni_180122004