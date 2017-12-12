! Crank-Nicolson method

MODULE MethodModule
  CONTAINS
    SUBROUTINE thomas_algorithm(a, b, c, nu, N)
      REAL *8, ALLOCATABLE, INTENT(in):: a(:), b(:), c(:)
      REAL *8, ALLOCATABLE, INTENT(inout):: nu(:)
      INTEGER, INTENT(in):: N
      
      
      REAL *8, ALLOCATABLE:: x(:), bb(:)
      REAL *8:: t
      INTEGER:: i
      
      ALLOCATE(x(N), bb(N))
            
      DO i=1, N
        bb(i) = b(i)
      ENDDO
      
      DO i=2, N
        t = a(i) / bb(i-1)
        bb(i) = bb(i) - (C(i-1) * t)
        nu(i) = nu(i) - (nu(i-1) * t)
      ENDDO
      
      x(N) = nu(N) / bb(N)
      DO i=N-1, 1, -1
        x(i) = (nu(i) - c(i) * x(i+1)) / bb(i)
      ENDDO
      
      DO i=1, N
        nu(i) = x(i)
      ENDDO
      
      DEALLOCATE(x, bb)
      RETURN
    END SUBROUTINE thomas_algorithm

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

    SUBROUTINE create_diag(a, b, c, N, s)
      REAL *8, ALLOCATABLE, INTENT(inout):: a(:), b(:), c(:)
      REAL *8, INTENT(in):: s
      INTEGER, INTENT(in):: N

      INTEGER:: i

      DO i=1, N
        a(i) = -s
        b(i) = (1.0 + 2.0*s)
        c(i) = -s
      ENDDO
      
      a(1) = 0.0
      b(1) = 1.0
      c(1) = 0.0
      
      a(N) = 0.0
      b(N) = 1.0
      c(N) = 0.0
      
      RETURN
    END SUBROUTINE create_diag
END MODULE ArrayModule



PROGRAM crank_nicolson
  USE MethodModule
  USE ArrayModule

  REAL *8:: L, dx, dt, s
  REAL *8, ALLOCATABLE:: a(:)
  REAL *8, ALLOCATABLE:: b(:)
  REAL *8, ALLOCATABLE:: c(:)
  REAL *8, ALLOCATABLE:: nu(:)
  REAL *8, ALLOCATABLE:: an(:)
  REAL *8, ALLOCATABLE:: x(:)
  REAL *8, ALLOCATABLE:: erro(:)
  INTEGER:: N, o, ts, i
  

  DO o=2, 14
    N = 2 ** o

    ts = 100
    dt = 0.000001
    L = 2.0 * ACOS(DBLE(-1.0))
    dx = L / (N - 1)
    s = dt / (2.0 * (dx * dx))

    ALLOCATE(nu(N), an(N), x(N), a(N), b(N), c(N), erro(ts))

    CALL create_diag(a, b, c, N, s)
    CALL initial_value(nu, dx, N)
    CALL initial_value(an, dx, N)

    DO i=1, ts-1
      erro(i) = infinity_norm(nu, an, N)
      CALL fill_values(an, N, dx, i * dt)
      CALL right_hand_updt(nu, N, s)
      CALL thomas_algorithm(a, b, c, nu, N)
    ENDDO

    PRINT *, MAXVAL(erro)

    DEALLOCATE(nu, an, x, a, b, c, erro)

  ENDDO
END PROGRAM crank_nicolson



REAL *8 FUNCTION analytic_solution(x, t)
  REAL *8, INTENT(in):: x, t
  
  analytic_solution = SIN(x) * EXP(-t)
END FUNCTION analytic_solution

