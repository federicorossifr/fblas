module fblas
    implicit none
    
    type Matrix 
        integer :: r,c
        real, allocatable :: data(:,:)
    end type Matrix

contains


function cMatrix(R,C) result(mat)
    implicit none
    integer,intent(in) :: R,C
    type(Matrix) :: mat
    mat%r = R
    mat%c = C
    allocate(mat%data(R,C))
    call fill_random(mat)
end function cMatrix


subroutine matrixMultiply(m1,m2,res)
    implicit none
    type(Matrix),intent(in) :: m1,m2
    type(Matrix) :: res

    integer i,j,l,R,C,K
    R = m1%r
    C = m1%c
    K = m2%c

    do i = 0, R-1
        do j = 0, K-1
            res%data(i,j) = 0

            do l = 0, C-1
                res%data(i,j) =  res%data(i,j) + m1%data(i,l)*m2%data(l,j)
            end do
        end do
    end do
end subroutine matrixMultiply

subroutine matrixConvolution(m,filter,stride,res)
    implicit none
    type(Matrix),intent(in) :: m,filter
    type(integer),intent(in) :: stride
    type(Matrix) :: res

    integer i,j,k,l,R,C,fR,fC
    real tmp
    
    R = res%r
    C = res%c
    fR = filter%r
    fC = filter%c
    do i = 0, R-1
        do j = 0, C-1 
            tmp = 0

            do k = 0, fR
                do l = 0, fC
                    tmp = tmp + m%data(i*stride + k,j*stride + l) * filter%data(k,l)
                end do
            end do

            res(i,j) = tmp;

        end do
    end do
end subroutine

function dotProduct(m1,m2) result(dot)
    use, intrinsic:: iso_fortran_env, only: stderr=>error_unit
    implicit none
    type(Matrix),intent(in) :: m1,m2
    real :: dot

    integer i,j,R,C

    R = m1%r
    C = m1%c
    if( (m1%R .ne. m2%R) .OR. (m1%C .ne. m2%C)) then
        write(stderr,*) 'Shape mismatch in dotProduct' 
        stop -1
    endif

    dot = 0
    do i = 0, R-1
        do j = 0, C-1
            dot = dot + m1%data(i,j) * m2%data(i,j)
        end do
    end do
    
end function dotProduct

subroutine print_mat(mat)
    implicit none
    type(Matrix), intent(in) :: mat

    integer i,j,R,C
    R = mat%r
    C = mat%c
    print *,"Shape:(",R,",",C,")"
    do i = 0, R-1
        do j = 0, C-1
            write(*, fmt="(f8.2,A) ", advance="no") (mat%data(i,j)),' '
        end do
        print *,''
    end do

end subroutine print_mat

subroutine write_mat(mat,unit)
    implicit none
        integer, intent(in) :: unit
        type(Matrix), intent(in) :: mat            
        integer i,j,R,C

        R = mat%R
        C = mat%C
        do i = 0, R-1
            do j = 0, C-1
                write(unit, *) mat%data(i,j)
            end do
        end do
end subroutine write_mat

subroutine fill_random(mat)
    implicit none
    type(Matrix) :: mat
    real tmp
    integer i,j

    do i=0, mat%R - 1
        do j = 0, mat%C - 1
            call random_number(tmp)
            mat%data(i,j) = tmp
        end do
    end do
end subroutine fill_random

end module fblas