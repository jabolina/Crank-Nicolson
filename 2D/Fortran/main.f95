! 2-D diffusion equation

MODULE method_module
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


    SUBROUTINE initial_value(f, N, h)
      REAL *8, ALLOCATABLE, INTENT(inout):: f(:,:)
      REAL *8, INTENT(in):: h
      INTEGER, INTENT(in):: N

      INTEGER:: i, j, k, l

      DO i=1, N
        k = i - 1

        f(i, 1) = 0.0
        f(i, N) = 0.0

        f(1, i) = 0.0
        f(N, i) = 0.0

        IF ((i .NE. 1) .AND. (i .NE. N)) THEN
          DO j=2, N-1
            l = j - 1
            f(i, j) = analytic_solution((k*h), (l*h), 0.0)
          ENDDO
        ENDIF
      ENDDO

      RETURN
    END SUBROUTINE initial_value


    REAL *8 FUNCTION infinity_norm(nu, an, N)
      REAL *8, ALLOCATABLE, INTENT(in):: nu(:,:), an(:,:)
      INTEGER, INTENT(in):: N

      REAL *8, ALLOCATABLE:: eps(:,:)
      INTEGER:: i, j

      ALLOCATE(eps(N,N))

      DO i=1, N
        DO j=1, N
          eps(i,j) = SQRT((nu(i,j)-an(i,j)) * (nu(i,j)-an(i,j)))
        ENDDO
      ENDDO

      infinity_norm = MAXVAL(eps)
    END FUNCTION infinity_norm


    SUBROUTINE updt_analytic(an, N, h, ta)
      REAL *8, ALLOCATABLE, INTENT(inout):: an(:,:)
      REAL *8, INTENT(in):: h, ta
      INTEGER, INTENT(in):: N

      INTEGER:: y, x, i, j

      DO i=1, N
        y = i - 1
        DO j=1, N
          x = j - 1
          an(i,j) = analytic_solution((x*h), (y*h), ta)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE updt_analytic
END MODULE method_module


PROGRAM bidimensional
  USE method_module

  REAL *8, ALLOCATABLE:: nu(:,:), an(:,:)
  REAL *8, ALLOCATABLE:: a(:), b(:), c(:), r(:), eps(:)
  REAL *8:: L, h, dt, s, at
  INTEGER:: N, o, st, i, j, k

  DO o=2, 8
    N = 2 ** o
    L = 1.0
    st = 100
    dt = 0.00001
    h = L / (N - 1)
    s = dt / (2.0 * (h * h))

    ALLOCATE(nu(N,N), an(N,N))
    ALLOCATE(a(N), b(N), c(N), r(N), eps(st))

    CALL create_diag(a, b, c, N, s)

    CALL initial_value(nu, N, h)

    DO k=1, st
      at = k - 1
      CALL updt_analytic(an, N, h, (at*dt))
      eps(k) = infinity_norm(nu, an, N)

      DO i=1, N
        r(1) = 0.0
        DO j=2, N-1
          r(j) = (s*nu(i,j-1)) + ((1.0 - 2.0*s)*nu(i,j)) + (s*nu(i,j+1))
        ENDDO
        r(N) = 0.0

        CALL thomas_algorithm(a, b, c, r, N)

        DO j=1, N
          nu(i, j) = r(j)
        ENDDO
      ENDDO

      DO j=1, N
        r(1) = 0.0
        DO i=2, N-1
          r(i) = (s*nu(i-1,j)) + ((1.0 - 2.0*s)*nu(i,j)) + (s*nu(i+1,j))
        ENDDO
        r(N) = 0.0

        CALL thomas_algorithm(a, b, c, r, N)

        DO i=1, N
          nu(i, j) = r(i)
        ENDDO
      ENDDO

    ENDDO

    PRINT *, MAXVAL(eps)

    DEALLOCATE(a, b, c, r, eps)
    DEALLOCATE(nu, an)
  ENDDO

END PROGRAM bidimensional

!#########################################################################################
REAL *8 FUNCTION analytic_solution(x, y, t)
  REAL *8, INTENT(in):: x, y, t
  REAL *8:: pi = ACOS(DBLE(-1.0))

  analytic_solution = SIN(pi*x) * SIN(pi*y) * EXP(-2.0*pi*t)
END FUNCTION analytic_solution

