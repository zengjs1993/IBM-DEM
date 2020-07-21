!----------------------------------------------------------
subroutine DeltaFunction3D(xcp1,xcp2,xcp3,weight)
use vars
implicit none
    real*8::weight,r1,r2,r3,w1,w2,w3
    real*8::xcp1,xcp2,xcp3
    r1 = abs(xcp1/dx)
    r2 = abs(xcp2/dy)
    r3 = abs(xcp3/dz)
    if (r1<=0.5d0)then
        w1 = (1+sqrt(1-3*r1*r1))/3
    else if (r1<=1.5d0)then
        w1 = (5-3*r1-sqrt(1-3*(1-r1)**2))/6
    else
        w1 =0
    endif
    if (r2<=0.5d0)then
        w2 = (1+sqrt(1-3*r2*r2))/3
    else if (r2<=1.5d0)then
        w2 = (5-3*r2-sqrt(1-3*(1-r2)**2))/6
    else
        w2 =0
    endif
    if (r3<=0.5d0)then
        w3 = (1+sqrt(1-3*r3*r3))/3
    else if (r3<=1.5d0)then
        w3 = (5-3*r3-sqrt(1-3*(1-r3)**2))/6
    else
        w3 =0
    endif
    weight = w1*w2*w3/(dx*dy*dz)
return
end subroutine
!----------------------------------------------------------
