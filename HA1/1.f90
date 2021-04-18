!Write a Fortran code which reads in your name and your date of birth. The code will
!print the following: ”dd.mm.yyyy is the date of birth of xxxx”

PROGRAM Akhil_1
    IMPLICIT NONE
    INTEGER :: D, M, Y
    CHARACTER Name*20
    READ*, Name
    Read*, D, M, Y
    PRINT*, (D),".",(M),".",Y, "is the date of birth of", Name
end
    