program dot_product
    use fblas

    implicit none
    
    integer R,C,K,F
    parameter (R=16,C=16,K=16,F=3)
    real dot
    type(Matrix) m1,m2,m3

    m1 = cMatrix(R,C)
    m2 = cMatrix(C,K)
    m3 = cMatrix(R,K)
    f  = cMatrix(F,F)
    call fill_random(m1)
    call fill_random(m2)
    call fill_random(f)
    call matrixMultiply(m1,m2,m3)

    open(2,file="mat.txt")
    call write_mat(m3,2)
    close(2)

    dot = dotProduct(m1,m2)
    print *,dot
end program dot_product




