! Crank-Nicolson method

MODULE MethodModule
  CONTAINS
    SUBROUTINE gauss_seidel(A, b, N)
      REAL *8, ALLOCATABLE, INTENT(in):: A(:,:)
      REAL *8, ALLOCATABLE, INTENT(inout):: b(:)
      INTEGER, INTENT(in):: N
      
      INTEGER:: i, j, iter
      REAL *8:: aux
      REAL *8, ALLOCATABLE:: x(:)
      
      ALLOCATE(x(N))
      
      DO i=1, N
        x(i) = 1.0
      ENDDO
      
      DO iter=5*N, 1, -1
        DO i=1, N
          aux = 0.0
          DO j=1, N
            IF (j .NE. i) THEN
              aux = aux + (A(i,j) * x(j))
            ENDIF
            x(i) = (b(i) - aux) / A(i,i)
          ENDDO
        ENDDO
      ENDDO
      
      DO i=1, N
        b(i) = x(i)
      ENDDO
      DEALLOCATE(x)
      RETURN
    END SUBROUTINE gauss_seidel

    REAL *8 FUNCTION infinity_norm(nu, an, N)
      REAL *8, ALLOCATABLE, INTENT(in):: nu(:)
      REAL *8, ALLOCATABLE, INTENT(in):: an(:)
      INTEGER, INTENT(in):: N

      INTEGER:: i
      REAL *8, ALLOCATABLE:: l(:)

      ALLOCATE(l(N))

      DO i=1, N
        l(i) = SQRT((nu(i) - an(i)) * (nu(i) - an(i)))
      ENDDO

      infinity_norm = MAXVAL(l)
    END FUNCTION infinity_norm

    SUBROUTINE initial_value(arr, dx, N)
      REAL *8, ALLOCATABLE, INTENT(inout):: arr(:)
      REAL *8, INTENT(in):: dx
      INTEGER, INTENT(in):: N

      INTEGER:: i, j

      DO i=1, N
        j = i - 1
        arr(i) = analytic_solution((j * dx), 0.0)
      ENDDO

      RETURN
    END SUBROUTINE initial_value
END MODULE MethodModule


MODULE ArrayModule
  CONTAINS
    SUBROUTINE right_hand_updt(nu, N, s)
      REAL *8, ALLOCATABLE, INTENT(inout):: nu(:)
      INTEGER, INTENT(in):: N
      REAL *8, INTENT(in):: s

      REAL *8, ALLOCATABLE:: x(:)
      INTEGER:: i

      ALLOCATE(x(N))

      x(1) = 0.0
      nu(1) = 0.0
      
      DO i=2, (N-1)
        x(i) = (s * nu(i-1)) + ((1.0 - 2.0*s) * nu(i)) + (s * nu(i+1))
      ENDDO
      
      x(N) = 0.0
      nu(N) = 0.0

      DO i=1, N
        nu(i) = x(i)
      ENDDO

      DEALLOCATE(x)
      RETURN
    END SUBROUTINE right_hand_updt

    SUBROUTINE fill_values(an, N, dx, ta)
      REAL *8, ALLOCATABLE, INTENT(inout):: an(:)
      INTEGER, INTENT(in):: N
      REAL *8, INTENT(in):: dx, ta

      INTEGER:: i, j

      DO i=1, N
        j = i - 1
        an(i) = analytic_solution((j * dx), ta)
      ENDDO

      RETURN
    END SUBROUTINE fill_values

    SUBROUTINE create_diag(A, N, s)
      REAL *8, ALLOCATABLE, INTENT(inout):: A(:,:)
      REAL *8, INTENT(in):: s
      INTEGER, INTENT(in):: N

      INTEGER:: i, j

      DO i=1, N
        DO j=1, N
          IF (i .EQ. j) THEN
            A(i, j) = (1.0 + 2.0*s)
          ELSE IF ((j .EQ. i-1) .OR. (j .EQ. i+1)) THEN
            A(i, j) = -s
          ELSE
            A(i, j) = 0.0
          ENDIF
        ENDDO
      ENDDO

      DO i=1, N
        A(1,i) = 0.0
        A(N, i) = 0.0
      ENDDO
      
      A(1,1) = 1.0
      A(N, N) = 1.0

      RETURN
    END SUBROUTINE create_diag
END MODULE ArrayModule



PROGRAM crank_nicolson
  USE MethodModule
  USE ArrayModule

  REAL *8:: L, dx, dt, s
  REAL *8, ALLOCATABLE:: A(:,:)
  REAL *8, ALLOCATABLE:: nu(:)
  REAL *8, ALLOCATABLE:: an(:)
  REAL *8, ALLOCATABLE:: x(:)
  REAL *8, ALLOCATABLE:: erro(:)
  INTEGER:: N, o, ts, i
  

  DO o=2, 7
    N = 2 ** o

    ts = 100
    dt = 0.001
    L = 2.0 * ACOS(DBLE(-1.0))
    dx = L / (N - 1)
    s = dt / (2.0 * (dx * dx))

    ALLOCATE(nu(N), an(N), x(N), A(N, N), erro(ts))

    CALL create_diag(A, N, s)
    CALL initial_value(nu, dx, N)
    CALL initial_value(an, dx, N)

    DO i=1, ts-1
      erro(i) = infinity_norm(nu, an, N)
      CALL fill_values(an, N, dx, i * dt)
      CALL right_hand_updt(nu, N, s)
      CALL gauss_seidel(A, nu, N)
    ENDDO

    PRINT *, MAXVAL(erro)

    DEALLOCATE(nu, an, x, A, erro)

  ENDDO
END PROGRAM crank_nicolson



REAL *8 FUNCTION analytic_solution(x, t)
  REAL *8, INTENT(in):: x, t
  
  analytic_solution = SIN(x) * EXP(-t)
END FUNCTION analytic_solution

