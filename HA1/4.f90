!Evaluate the following when A = 3, B = 5 and C = 4. A, B and C are all integer variables.
!  A/B×C, A×B/C, 3×AB, B+A/C, and A/B+C
PROGRAM Akhil_4
    IMPLICIT NONE
    INTEGER :: A = 3, B = 5, C = 4
    print*, A/B*C
    print *, A*B/C
    print*, 3*(A**B)
    print *, B + A/C
    print *, A/B+C
END
    