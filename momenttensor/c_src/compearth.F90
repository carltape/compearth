      MODULE COMPEARTH_MODULE 
         INTERFACE !COMPEARTH_INTERFACE
            INTEGER(C_INT) FUNCTION compearth_angleMT(n, M1, M2, theta) &
                           BIND(C, NAME='compearth_angleMT')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: n
            REAL(C_DOUBLE), INTENT(IN) :: M1(6*n), M2(6*n)
            REAL(C_DOUBLE), INTENT(OUT) :: theta(n)
            END FUNCTION

            SUBROUTINE compearth_beta2u(n, beta, u) &
                       BIND(C, NAME='compearth_beta2u')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: n
            REAL(C_DOUBLE), INTENT(IN) :: beta(n)
            REAL(C_DOUBLE), INTENT(OUT) :: u(n)
            END SUBROUTINE

            SUBROUTINE compearth_gamma2v(n, gama, v) &
                       BIND(C, NAME='compearth_gamma2v')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: n
            REAL(C_DOUBLE), INTENT(IN) :: gama(n)
            REAL(C_DOUBLE), INTENT(OUT) :: v(n)
            END SUBROUTINE

            REAL(C_DOUBLE) FUNCTION compearth_matlab_fangle(n, va, vb) &
                           BIND(C, NAME='compearth_matlab_fangle')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN), VALUE :: n
            REAL(C_DOUBLE), INTENT(IN) :: va(n), vb(n)
            END FUNCTION

         END INTERFACE !COMPEARTH_INTERFACE 
      END MODULE
