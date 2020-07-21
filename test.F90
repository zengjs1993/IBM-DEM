program main
implicit none
    integer,dimension(1:3,1:3,1:3)::u
    integer,dimension(1:2)::idx
    integer::i,j,k
    do i=1,3
    do j=1,3
    do k=1,3
        u(i,j,k) = 9*i+3*j+k
    enddo
    enddo
    enddo
    u = u**2
    write(*,*)u
    !idx = 1:2
    
end
