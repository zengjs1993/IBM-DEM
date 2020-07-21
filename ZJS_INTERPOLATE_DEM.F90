!-----------------------------------------------
subroutine GetWeight(xp1,xp2,xp3,xc1,xc2,xc3,weight)
use vars
use ZjsVar
implicit none 
    real*8::xp1,xp2,xp3,xc1,xc2,xc3,deltax1,deltax2,deltax3,weight(1:3)
    integer::flag
    deltax1 = (xp1-xc1)
    deltax2 = (xp2-xc2)
    deltax3 = (xp3-xc3)
        
    call DeltaFunction3D(deltax1+0.5d0*di,deltax2,deltax3,weight(1))
    call DeltaFunction3D(deltax1,deltax2+0.5d0*dj,deltax3,weight(2))
    call DeltaFunction3D(deltax1,deltax2,deltax3+0.5d0*dk,weight(3))
    !weight = (1-deltax1/di)*(1-deltax2/dj)*(1-deltax3/dk)/(di*dj*dk)
return
endsubroutine
!-----------------------------------------------
subroutine GetWeightSumForParcels()
use vars
use ZjsVar
implicit none
return
endsubroutine
!-----------------------------------------------
subroutine LocalWeightNormalization()
use vars
use ZjsVar
implicit none
return
endsubroutine
!-----------------------------------------------
