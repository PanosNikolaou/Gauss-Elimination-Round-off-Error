module tools

  contains

  subroutine printmtx(a)
  use types

    real(kind=rk), dimension (:,:), intent(in) :: a
    integer :: i

        do i = 1, ubound (a,1)
            print *, "row numbers", i, "/", a(i,:)
        end do
    end subroutine printmtx

  subroutine gauss(a,b,x)
  use types

    implicit none

    real (kind=rk), dimension (:,:), intent (in),allocatable :: a
    real (kind=rk), dimension (:), intent (out),allocatable :: x
    real (kind=rk), dimension (:), intent (in),allocatable :: b

    real (kind=rk), dimension (:,:), allocatable :: acopy

    integer :: i, j, n

    if (.not. allocated(a) .or. .not.(allocated(b))) then
        print *, "inputs to gem are not allocated..."
        return
    end if

    if ((ubound(a,1) /= ubound(a,2))) then
        print*, "matrix a is not square..."
        return
    end if

    if ((ubound(a,1) /= ubound(b,1))) then
        print*, "vector b and matrix a have incompatible lengths..."
        return
    end if

    n = ubound(a,1)

    if (.not. allocated (x)) then
        allocate(x(n))
    end if

    allocate(acopy(ubound(a,1), ubound(a,1)+1))
    acopy(:, 1:ubound(a,1)) = a
    acopy(:, ubound(a,1)+1) = b

    do j=2, ubound(a,1)
       do i = j, ubound(a,1)
         acopy(i, :) = acopy(i,:) - acopy(j-1,:)/acopy(j-1,j-1)*acopy(i,j-1)
       end do
    end do

    forall (i=n:1:-1) x(i) = ( acopy(i,n+1) - sum(acopy(i,i+1:n)*x(i+1:n)) ) / acopy(i,i)

    !call printmtx(acopy)

  end subroutine gauss

  subroutine swap_reals(a, b)
  use types
    real (kind=rk), intent(inout) :: a, b
    real (kind=rk) :: tmp

    tmp = a
    a = b
    b = tmp
  end subroutine swap_reals

  subroutine gen_hilbert(n, A)
  use types
    integer, intent(in) :: n
    integer :: i,j

    real (kind=rk), dimension (:,:), intent (out), allocatable :: A

    allocate(A(n,n))
    do i=1,n
        do j=1,n
            A(i,j) = 1.0_rk/(i+j-1)
        end do
    end do
  end subroutine gen_hilbert

  subroutine inverse(a,c,n)
    implicit none
    integer n
    real (kind=rk) ::  a(n,n), c(n,n)
    real (kind=rk) ::  L(n,n), U(n,n), b(n), d(n), x(n)
    real (kind=rk) ::  coeff
    integer i, j, k

    L=0.0
    U=0.0
    b=0.0

    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
                do j=k+1,n
                    a(i,j) = a(i,j)-coeff*a(k,j)
                end do
        end do
    end do

    do i=1,n
        L(i,i) = 1.0
    end do

    do j=1,n
        do i=1,j
            U(i,j) = a(i,j)
        end do
    end do

    do k=1,n
        b(k)=1.0
        d(1) = b(1)

        do i=2,n
            d(i)=b(i)
                do j=1,i-1
                    d(i) = d(i) - L(i,j)*d(j)
                end do
        end do

        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
            x(i) = d(i)
                do j=n,i+1,-1
                    x(i)=x(i)-U(i,j)*x(j)
                end do
            x(i) = x(i)/u(i,i)
        end do

        do i=1,n
            c(i,k) = x(i)
        end do
        b(k)=0.0
    end do
  end subroutine inverse

end module tools
