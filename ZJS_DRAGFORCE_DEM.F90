!-----------------------------------------------
subroutine CalculateDragForce(i,denf,alphaf_local,velf,beta)
use vars
use ZjsVar
implicit none
    integer::i,j
    real*8,dimension(1:3)::velf,beta    
    real*8::denf,alphaf_local

    real*8::f,Cd,Rep,diameter,Ur,spg

    real*8::zeta,fbase

    spg = part(i)%rho/denf
    diameter = 2d0*part(i)%rad
    do j=1,3
        Ur = abs(velf(j)-part(i)%v(j))
        Rep = Ur*diameter/viscf
        Rep = max(1d0,Rep)
        zeta = 3.7-0.65*exp(-0.5*(1.5-log10(Rep))**2)
        if (Rep<1000) then
            fbase = 1+1.0d0/6d0*Rep**0.687
        else
            fbase = 0.0183*Rep
        endif
        if (alphaf_local<0.8)then
            f = 8.33*(1-alphaf_local)/alphaf_local+0.0972*Rep
        else
            f = fbase*alphaf_local**(-zeta)
        endif
        Cd = 24d0/Rep*f
        beta(j) = 0.75d0*(Cd*Ur)/(diameter*spg)
    enddo
return
endsubroutine
!-----------------------------------------------  
