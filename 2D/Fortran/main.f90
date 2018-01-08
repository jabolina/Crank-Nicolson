! 2D diffusion equation

module functions_module
  contains
    subroutine thomas_algorithm(a, b, c, r, x, s, e)
      real *8, allocatable, intent(in):: a(:), c(:)
      real *8, allocatable, intent(inout):: x(:), b(:), r(:)
      integer, intent(in):: s, e

      integer:: i

      do i=s+1, e
        b(i) = b(i) - a(i)/b(i-1)*c(i-1)
        r(i) = r(i) - a(i)/b(i-1)*r(i-1)
      end do

      x(e) = r(e)/b(e)
      do i=e-1,s,-1
        x(i) = (r(i)-c(i)*x(i+1))/b(i)
      end do

      return
    end subroutine thomas_algorithm


    subroutine initial_value(f, n, h)
      real *8, allocatable, intent(inout):: f(:,:)
      real *8, intent(in):: h
      integer, intent(in):: n

      integer:: i, j

      do i=0, n
        do j=0, n
          f(i, j) = analytic_solution((j*h), (i*h), 0.0d0)
        end do
      end do

      return
    end subroutine initial_value


    real *8 function infinity_norm(nu, an, n)
      real *8, allocatable, intent(in):: nu(:,:), an(:,:)
      integer, intent(in):: n

      real *8, allocatable:: eps(:,:)
      real *8:: maximum
      integer:: i, j

      allocate(eps(0:n,0:n))

      do i=0, n
        do j=0, n
          eps(i,j) = SQRT((nu(i,j) - an(i,j)) * (nu(i,j) - an(i,j)))
        end do
      end do

      maximum = eps(0, 0)

      do i=0, n
        do j=0, n
          if (eps(i, j) > maximum) then
            maximum = eps(i, j)
          end if
        end do
      end do

      deallocate(eps)

      infinity_norm = maximum
    end function infinity_norm


    subroutine updt_analytic(an, n, h, ta)
      real *8, allocatable, intent(inout):: an(:,:)
      real *8, intent(in):: h, ta
      integer, intent(in):: n

      integer:: y, x

      do y=0, n
        do x=0, n
          an(y,x) = analytic_solution((x*h), (y*h), ta)
        end do
      end do

      return
    end subroutine updt_analytic


    subroutine save_analytic(an, n)
      real *8, allocatable, intent(in):: an(:,:)
      integer, intent(in):: n

      integer:: i, j

      open(unit=2, file="analytic.dat")

      do i=0, n
        do j=0, n
          write (2,*) i, j, an(i, j)
        end do
      end do

      return
    end subroutine save_analytic


    subroutine save_numeric(nu, n)
      real *8, allocatable, intent(in):: nu(:,:)
      integer, intent(in):: n

      integer:: i, j

      open(unit=2, file="numeric.dat")

      do i=0, n
        do j=0, n
          write (2,*) i, j, nu(i, j)
        end do
      end do

      return
    end subroutine save_numeric


    subroutine crank_nicolson(nu, n, s)
      real *8, allocatable, intent(inout):: nu(:,:)
      real *8, intent(in):: s
      integer, intent(in):: n

      real *8, allocatable:: a(:), b(:), c(:), r(:), x(:)
      real *8, allocatable:: z(:,:)
      integer:: i, j

      allocate(z(0:n, 0:n))

      do j=0, n
        do i=0, n
          z(i, j) = 0.0d0
        end do
      end do

      do j=1, n-1
        allocate(a(1:n-1), b(1:n-1), c(1:n-1), r(1:n-1), x(1:n-1))

        do i=1, n-1
          a(i) = -s
          b(i) = (1.0d0 + 2.0d0*s)
          c(i) = -s
          r(i) = (s*nu(i, j-1)) + ((1.0d0 - 2.0d0*s)*nu(i, j)) + (s*nu(i, j+1))
        end do

        call thomas_algorithm(a, b, c, r, x, 1, n-1)

        z(0, j) = 0.0d0
        do i=1, n-1
          z(i, j) = x(i)
        end do
        z(n, j) = 0.0d0

        deallocate(a, b, c, r, x)
      end do

      do i=1, n-1
        allocate(a(1:n-1), b(1:n-1), c(1:n-1), r(1:n-1), x(1:n-1))

        do j=1, n-1
          a(j) = -s
          b(j) = (1.0d0 + 2.0d0*s)
          c(j) = -s
          r(j) = (s*z(i-1, j)) + ((1.0d0 - 2.0d0*s)*z(i, j)) + (s*z(i+1, j))
        end do

        call thomas_algorithm(a, b, c, r, x, 1, n-1)

        nu(i, 0) = 0.0d0
        do j=1, n-1
          nu(i, j) = x(j)
        end do
        nu(i, n) = 0.0d0

        deallocate(a, b, c, r, x)
      end do

      deallocate(z)

      return
    end subroutine crank_nicolson

end module functions_module

program diffusion
  use functions_module

  integer:: n, o, k, st, i, j
  real *8:: d, dt, l, s
  real *8, allocatable:: eps(:)
  real *8, allocatable:: an(:,:), nu(:,:)

  l = 1.0d0
  st = 100

  open(unit=3, file="error.dat")
  do o=2, 6
    n = 2 ** o
    d = l / (dfloat(n))
    dt = 0.0001d0
    s = dt / (2.0d0 * (d * d))

    allocate(nu(0:n, 0:n), an(0:n, 0:n))
    allocate(eps(0:st))

    call initial_value(nu, n, d)
    call updt_analytic(an, n, d, 0.0d0)
    eps(0) = infinity_norm(nu, an, n)

    do k=1, st-1
      call save_analytic(an, n)
      call save_numeric(nu, n)
      call updt_analytic(an, n, d, dfloat(k)*dt)

      call crank_nicolson(nu, n, s)

      eps(k) = infinity_norm(nu, an, n)
    end do

    write (3,*) MAXVAL(eps)
    deallocate(an, nu)
    deallocate(eps)
  end do
  close(3)

end program diffusion


!###############################################################################
real *8 function analytic_solution(x, y, t)
  real *8, intent(in):: x, y, t
  real *8:: pi = ACOS(DBLE(-1.0))

  analytic_solution = SIN(pi*x) * SIN(pi*y) * EXP(-2.0d0*(pi*pi)*t)
end function analytic_solution
