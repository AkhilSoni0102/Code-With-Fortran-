!                                       Chemical Kinetics
!In this problem, you will solve a chemical kinetics problem.
!Consider the following set of reactions
!                                           A ↔ B
!                                           B ↔ 2C
!                                           C ↔ D.
!For these elementary reactions, the relevant equations are
!d[A]/dt = −kab[A] + kba[B]
!d[B]/dt = +kab[A] − kba[B] − kbc[B] + kcb[C]^2
!d[C]/dt = +2kbc[B] − 2kcb[C]^2 − kcd[C]
!d[D]/dt = +kcd[C]

!Write a Fortran program to propagate this set of equations, using Euler’s and 4thorder Runge-Kutta methods. Plot the concentrations of all the species as a function
!of time. The code will be modular, i.e. separate routines (modules), e.g. constants are
!declared in a module, the Euler’s in a subroutine etc. Use arrays for concentrations
!• Use the following : kab = 1, kba = 3, kbc = 4.2, kcb = 7.3, kcd = 0.4 and dt = 0.0001,
!tmax = 60, [A]0 = 5, [B]0=[C]0=[D]0=0.0


MODULE precision_use
    ! dp = double precision_use
    INTEGER, PARAMETER:: dp = SELECTED_REAL_KIND(12)
END MODULE precision_use

MODULE  Reaction_Constants
    use precision_use
    implicit none 
    integer, parameter:: n = 4
    real(kind = dp), parameter:: kab = 1, kba = 3, kbc = 4.2, kcb = 7.3, kcd = 0.4
end module Reaction_Constants

PROGRAM kinetics
    use precision_use
    use Reaction_Constants
    implicit none 
    REAL(KIND = dp):: concentration(n)
    REAL(KIND = dp):: t, dt, tmax
    INTEGER:: nt, i

    t = 0
    dt = 0.0001
    tmax = 60;
    ! Number of iteration
    nt = int(tmax/dt);
    concentration = 0
    concentration(1) = 5.0      

    OPEN(unit = 10, file = "output.txt")
    DO i = 1,nt 
        t = t + dt;
        ! CALL Euler(concentration, dt)
        ! CALL RK2(concentration, dt)
        CALL RK4(concentration, dt)
        WRITE(10, *) t, concentration
    ENDDO
    CLOSE(10)

    CONTAINS
    SUBROUTINE Euler(concentration, dt)
        USE precision_use
        USE Reaction_Constants
        IMPLICIT NONE

        REAL(KIND = dp):: concentration(n), Func(n)
        REAL(KIND = dp):: dt
        CALL DbyDT(concentration, Func)
        concentration = concentration + dt * Func
    END SUBROUTINE Euler

    SUBROUTINE RK2(concentration, dt)
        USE precision_use
        IMPLICIT NONE

        REAL(KIND = dp):: concentration(n), sl(n), sr(n)
        REAL(KIND = dp):: dt

        ! RK2 Propagation
        CALL DbyDT(concentration, sl)
        CALL DbyDT(concentration + dt * sl, sr)
        concentration = concentration + dt * (sl + sr)/2.0_dp
    END SUBROUTINE RK2

    SUBROUTINE RK4(concentration, dt)
        USE precision_use
        IMPLICIT NONE

        REAL(KIND = dp):: concentration(n), s1(n), s2(n), s3(n), s4(n)
        REAL(KIND = dp):: dt

        ! RK2 Propagation
        CALL DbyDT(concentration, s1)
        CALL DbyDT(concentration + dt * s1/2.0_dp, s2)
        CALL DbyDT(concentration + dt * s2/2.0_dp, s3)
        CALL DbyDT(concentration + dt * s3, s4)
        concentration = concentration + dt * (s1 + 2*s2 + 2*s3 + s4)/6.0_dp
    END SUBROUTINE RK4

    SUBROUTINE DbyDT(concentration, Func)
        USE precision_use
        USE Reaction_Constants
        IMPLICIT NONE

        REAL(KIND = dp):: concentration(n), Func(n)

        Func(1) = -kab * concentration(1) + kba * concentration(2)
        Func(2) = kab * concentration(1) - kba * concentration(2) - kbc * concentration(2) + kcb * concentration(3)**2
        Func(3) = 2 * kbc * concentration(2) - 2 * kcb * concentration(3)**2 - kcd * concentration(3)
        Func(4) = kcd * concentration(3)
    END SUBROUTINE DbyDT

END PROGRAM kinetics

