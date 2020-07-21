!-----------------------------------------------
subroutine CloudMaxVelocityMagnitude(vmax)
use vars
use ZjsVar
implicit none
    integer::i
    real*8::vmax,vi

    vmax = 0d0
    do i=1,np
        vi = part(i)%v(1)**2+part(i)%v(2)**2+part(i)%v(3)**2
        if (vi>vmax)then
            vmax = vi
        endif
    enddo
    vmax = sqrt(vmax)
return
endsubroutine
!-----------------------------------------------
subroutine CrossMultiplication(a,b,d)
use vars
use ZjsVar
implicit none
    real*8,dimension(1:3)::a,b,c,d
    c(1) =  a(2)*b(3)-a(3)*b(2)
    c(2) = -a(1)*b(3)+a(3)*b(1)
    c(3) =  a(1)*b(2)-a(2)*b(1)
    d = c
return
endsubroutine
!-----------------------------------------------
subroutine TensorMultiplication(A,B,D)
use vars
use ZjsVar
implicit none
    real*8,dimension(1:9)::A,B,C,D
    C(1) = A(1)*B(1)+A(2)*B(4)+A(3)*B(7)
    C(2) = A(1)*B(2)+A(2)*B(5)+A(3)*B(8)
    C(3) = A(1)*B(3)+A(2)*B(6)+A(3)*B(9)
    C(4) = A(4)*B(1)+A(5)*B(4)+A(6)*B(7)
    C(5) = A(4)*B(2)+A(5)*B(5)+A(6)*B(8)
    C(6) = A(4)*B(3)+A(5)*B(6)+A(6)*B(9)
    C(7) = A(7)*B(1)+A(8)*B(4)+A(9)*B(7)
    C(8) = A(7)*B(2)+A(8)*B(5)+A(9)*B(8)
    C(9) = A(7)*B(3)+A(8)*B(6)+A(9)*B(9)
    D = C
return
endsubroutine
!-----------------------------------------------
subroutine TensorVectorMultiplication(A,b,d)
use vars
use ZjsVar
implicit none
    real*8,dimension(1:3)::b,c,d
    real*8,dimension(1:9)::A
    c(1) = A(1)*b(1)+A(2)*b(2)+A(3)*b(3)
    c(2) = A(4)*b(1)+A(5)*b(2)+A(6)*b(3)
    c(3) = A(7)*b(1)+A(8)*b(2)+A(9)*b(3)
    d = c
return
endsubroutine
!-----------------------------------------------
subroutine InverseOrthogonalMatrix(A,B)
use vars
use ZjsVar
implicit none
    real*8,dimension(1:9)::A,B
    B(1) = A(1); B(2) = A(4); B(3) = A(7)
    B(4) = A(2); B(5) = A(5); B(6) = A(8)
    B(7) = A(3); B(8) = A(6); B(9) = A(9)
return
endsubroutine
!-----------------------------------------------
