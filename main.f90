program gem
use tools
use types

    implicit none
    real (kind=rk), dimension (:,:), allocatable :: a,c
    real (kind=rk), dimension (:), allocatable :: x
    real (kind=rk), dimension (:), allocatable :: r
    real (kind=rk), dimension (:), allocatable :: b

    real (kind=rk) :: error

    integer, dimension(:), allocatable :: seed

    integer :: n

    integer :: i,IUNIT

    integer :: values(1:100), k=12

    print*, "Enter the size of the array:"
    read*, n

    allocate ( a(n,n) )
    allocate ( x(n) )
    allocate ( b(n) )
    allocate ( r(n) )
    allocate ( c(n,n) )

    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:) = values(100)
    call random_seed(put=seed)

    call random_number(a)
    call random_number(b)

!-------------------- testing matrix ------------------------------------

!    a = transpose(reshape((/1.0,  6.0, -5.0, 2.0, &
!                            2.0, -3.0,  12.0, 4.0, &
!                            2.0, -6.0,  1.0,  7.0, &
!                            4.0,  4.0, -4.0,  4.0  /), &
!                            (/4,4/)))
!
!    b = (/1.0,  5.0,  1.0,  5.0/)

!    r:   1.4730971128608925       0.19553805774278207       0.26968503937007876      -0.14895013123359568
!    x:   1.4730971128608923       0.19553805774278216       0.26968503937007876      -0.14895013123359577
!------------------------------------------------------------------------

    call gen_hilbert(n,a)

    call gauss(a,b,x)

    call inverse(a,c,n)

    r = matmul(c,b)

    open(NEWUNIT=IUNIT,FILE='data.dat',FORM="FORMATTED")

    do i=1,n
        error = abs(r(i) - x(i))
        write(IUNIT,*) error,i
        print*,x(i),r(i)
    end do

    close(IUNIT)

    call SYSTEM('gnuplot -p data_plot.plt')

end program gem
