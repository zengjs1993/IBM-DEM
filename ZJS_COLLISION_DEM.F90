!-----------------------------------------------
subroutine Collision(i,j,flag)
use vars
use ZjsVar
implicit none
    integer::i,j,flag
    real*8::distSqr,xij1,xij2,xij3
    xij1 = part(i)%x(1)-part(j)%x(1)
    xij2 = part(i)%x(2)-part(j)%x(2)
    xij3 = part(i)%x(3)-part(j)%x(3)
    distSqr = xij1**2+xij2**2+xij3**2
    
    if (distSqr < (collisionDist*(1d0+collisionBuffer))**2)then
        flag = 1
    else
        flag = 0
    endif
return
endsubroutine
!-----------------------------------------------  
subroutine CalculateCollisionForce()
use vars
use ZjsVar
implicit none
    integer::i,j,ii,jj
    real*8::dist,rij,deltan,kn,gamman,eta,etan,meff
    real*8::vnij
    real*8::r1,r2,r3,n1,n2,n3,f1,f2,f3

    eta = -2*dlog(ResistCoeff)/sqrt(3.1415926536**2+dlog(ResistCoeff))
    do ii=1,npLocal
    part(ii)%fcol(:) = 0
    do jj=1,neighPartNo(ii)
        i = ii
        j = neighPartIdx(jj,i)
        r1 = part(i)%x(1)-part(j)%x(1)
        r2 = part(i)%x(2)-part(j)%x(2)
        r3 = part(i)%x(3)-part(j)%x(3)
        dist = sqrt(r1**2+r2**2+r3**2)
        n1 = r1/dist
        n2 = r2/dist
        n3 = r3/dist
        meff = 2*part(i)%m*part(j)%m/(part(i)%m+part(j)%m)
        rij = part(i)%radeff+part(j)%radeff
        deltan = rij-dist
        kn = knCoeff/(0.1d0*di)*meff
        etan = eta*sqrt(meff*kn)

        f1 = 0; f2 = 0; f3 = 0;
        vnij = n1*(part(i)%v(1)-part(j)%v(1)) + &
               n2*(part(i)%v(2)-part(j)%v(2)) + &
               n3*(part(i)%v(3)-part(j)%v(3))
        if (deltan>0)then
            f1 = n1*(kn*deltan-etan*vnij)
            f2 = n2*(kn*deltan-etan*vnij)
            f3 = n3*(kn*deltan-etan*vnij)

            !open(cpuLocal,file=logFile,access="append")
            !write(cpuLocal,*) "collision: ",part(i)%x,part(j)%x
            !write(cpuLocal,*) part(i)%MSIdx,rij,deltan,f1,f2,f3
            !close(cpuLocal)
        endif
            !open(cpuLocal,file=logFile,access="append")
            !write(cpuLocal,*) part(i)%x,part(j)%x
            !write(cpuLocal,*) rij,deltan,f1,f2,f3
            !close(cpuLocal)
        part(i)%fcol(1) = part(i)%fcol(1) + f1
        part(i)%fcol(2) = part(i)%fcol(2) + f2
        part(i)%fcol(3) = part(i)%fcol(3) + f3

    enddo
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine CalculateSolidBoundaryCollisionForce()
use vars
use ZjsVar
implicit none
    integer::i,j,k,cellIdx1,cellIdx2,cellIdx3,flag
    real*8::dist,rij,deltan,kn,gamman,eta,etan,meff
    real*8::vnij
    real*8::xij1,xij2,xij3,r1,r2,r3,f1,f2,f3,fn1,fn2,fn3
    real*8::p1,p2,p3,n1,n2,n3,cp1,cp2,cp3
    
    do i=1,6
    do j=1,numWallProbe(i)
        wallProbe(j,i)%force(1:3) = 0d0
    enddo
    enddo

    eta = -2*dlog(ResistCoeff)/sqrt(3.1415926536**2+dlog(ResistCoeff))
    do i=1,npLocal
        cellIdx1 = map1Dto3D(1,part(i)%cellIdx)
        cellIdx2 = map1Dto3D(2,part(i)%cellIdx)
        cellIdx3 = map1Dto3D(3,part(i)%cellIdx)
        xij1 = part(i)%x(1); xij2 = part(i)%x(2); xij3 = part(i)%x(3)
        if (cellIdx1==SG1.or.cellIdx1==EG1.or.&
            cellIdx2==SG2.or.cellIdx2==EG2.or.&
            cellIdx3==SG3.or.cellIdx3==EG3)then
        do j=1,6
            if (bType(j)==1)then
                p1 = pointBound(1,j); p2 = pointBound(2,j); p3 = pointBound(3,j)
                n1 = normBound(1,j); n2 = normBound(2,j); n3 = normBound(3,j)
                r1 = xij1-p1; r2 = xij2-p2; r3 = xij3-p3
                dist = r1*n1+r2*n2+r3*n3
                cp1 = xij1-dist*n1
                cp2 = xij2-dist*n2
                cp3 = xij3-dist*n3
                meff = part(i)%m
                rij = part(i)%radeff
                deltan = rij-dist
                kn = knCoeff/(0.1d0*di)*meff
                etan = eta*sqrt(meff*kn)

                f1 = 0; f2 = 0; f3 = 0
                vnij = n1*part(i)%v(1)+n2*part(i)%v(2)+n3*part(i)%v(3)
                if (deltan>0)then
                    f1 = n1*(kn*deltan-etan*vnij)
                    f2 = n2*(kn*deltan-etan*vnij)
                    f3 = n3*(kn*deltan-etan*vnij)
                    fn1 = n1*kn*deltan
                    fn2 = n2*kn*deltan
                    fn3 = n3*kn*deltan
                    part(i)%fcol(1) = part(i)%fcol(1) + f1
                    part(i)%fcol(2) = part(i)%fcol(2) + f2 
                    part(i)%fcol(3) = part(i)%fcol(3) + f3

                    ! search wall probes
                    do k=1,numWallProbe(j)
                        call WallProbeCollision(cp1,cp2,cp3,wallProbe(k,j),flag)
                        if (flag==1)then
                            wallProbe(k,j)%force(1) = wallProbe(k,j)%force(1) + fn1 
                            wallProbe(k,j)%force(2) = wallProbe(k,j)%force(2) + fn2 
                            wallProbe(k,j)%force(3) = wallProbe(k,j)%force(3) + fn3
                            
                            open(cpuLocal,file=logFile,access="append")
                            write(cpuLocal,*)i,j,f1,f2,f3
                            close(cpuLocal)
                        endif
                    enddo
                endif
            endif
        enddo
        endif
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine WallProbeCollision(cp1,cp2,cp3,prob,flag)
use ZjsVar
implicit none
real*8::cp1,cp2,cp3,ni1,ni2,ni3,nj1,nj2,nj3,thetai,theta,pi
integer::flag,i,idxi,idxj
type(WallColProbe)::prob
theta = 0d0
do i=1,4
    idxi = i; idxj = i+1
    if (idxj>4) idxj = idxj-4
    ni1 = prob%points(1,idxi)-cp1
    ni2 = prob%points(2,idxi)-cp2
    ni3 = prob%points(3,idxi)-cp3
    nj1 = prob%points(1,idxj)-cp1
    nj2 = prob%points(2,idxj)-cp2
    nj3 = prob%points(3,idxj)-cp3
    
    thetai = (ni1*nj1+ni2*nj2+ni3*nj3)/ &
        (sqrt(ni1**2+ni2**2+ni3**2)*sqrt(nj1**2+nj2**2+nj3**2))
    thetai = acos(thetai)
    theta = theta+thetai
enddo
pi = acos(-1.0d0)
if (abs(theta-2*pi)>0.001*pi)then
    flag = 0
else
    flag = 1
endif
return
endsubroutine
!-----------------------------------------------
