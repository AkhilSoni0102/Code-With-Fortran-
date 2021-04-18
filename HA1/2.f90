!Write a Fortran code which reads in two complex numbers, and then multiplies those
!two. Print the answer

PROGRAM Akhil_2
IMPLICIT NONE
INTEGER :: a,b,c,d
READ*, a, b, c, d
print*, a*c - b*d, " ", a*d + b*c
END