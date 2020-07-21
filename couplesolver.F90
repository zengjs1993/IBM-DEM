!----------------------------------------------------------
subroutine EvolutionMS()
use vars
use ZjsVar
implicit none
    integer::i
    real*8::ratio,vTmp(1:3),omegaTmp(1:3)
    real*8::g
    g = 9.81d0/1
    !call CalculateCollisionForceMS()    
    !call CalculateSolidBoundaryCollisionForceMS()    
    do i=1,nMSLocal
        ratio = (MS(i)%rho-rhof)/MS(i)%rho
        !MS(i)%v(1) = MS(i)%v(1) + MS(i)%fdrag(1)/MS(i)%m/ratio*rhof*dtDEM  + MS(i)%fcol(1)/MS(i)%m/ratio*dtDEM
        !MS(i)%v(2) = MS(i)%v(2) + MS(i)%fdrag(2)/MS(i)%m/ratio*rhof*dtDEM  + MS(i)%fcol(2)/MS(i)%m/ratio*dtDEM
        !MS(i)%v(3) = MS(i)%v(3) + MS(i)%fdrag(3)/MS(i)%m/ratio*rhof*dtDEM  + MS(i)%fcol(3)/MS(i)%m/ratio*dtDEM - g*dtDEM
        !MS(i)%omega(1) = MS(i)%omega(1) + MS(i)%fmom(1)/MS(i)%I*dtDEM/ratio*rhof
        !MS(i)%omega(2) = MS(i)%omega(2) + MS(i)%fmom(2)/MS(i)%I*dtDEM/ratio*rhof
        !MS(i)%omega(3) = MS(i)%omega(3) + MS(i)%fmom(3)/MS(i)%I*dtDEM/ratio*rhof
        vTmp = MS(i)%v
        omegaTmp = MS(i)%omega
        MS(i)%v(1) = MS(i)%v(1) + (1-ratio)*(MS(i)%v(1) - MS(i)%vold(1)) + MS(i)%fdrag(1)/MS(i)%m*rhof*dtDEM  + MS(i)%fcol(1)/MS(i)%m*dt
        MS(i)%v(2) = MS(i)%v(2) + (1-ratio)*(MS(i)%v(2) - MS(i)%vold(2)) + MS(i)%fdrag(2)/MS(i)%m*rhof*dtDEM  + MS(i)%fcol(2)/MS(i)%m*dt
        MS(i)%v(3) = MS(i)%v(3) + (1-ratio)*(MS(i)%v(3) - MS(i)%vold(3)) + MS(i)%fdrag(3)/MS(i)%m*rhof*dtDEM  + MS(i)%fcol(3)/MS(i)%m*dt - ratio*g*dt
        MS(i)%omega(1) = MS(i)%omega(1) + (1-ratio)*(MS(i)%omega(1) - MS(i)%omegaold(1)) + MS(i)%fmom(1)/MS(i)%I*dtDEM*rhof
        MS(i)%omega(2) = MS(i)%omega(2) + (1-ratio)*(MS(i)%omega(2) - MS(i)%omegaold(2)) + MS(i)%fmom(2)/MS(i)%I*dtDEM*rhof
        MS(i)%omega(3) = MS(i)%omega(3) + (1-ratio)*(MS(i)%omega(3) - MS(i)%omegaold(3)) + MS(i)%fmom(3)/MS(i)%I*dtDEM*rhof
        !MS(i)%v(1) = MS(i)%v(1) + MS(i)%fdrag(1)/MS(i)%m*rhof*dtDEM  + MS(i)%fcol(1)/MS(i)%m*dt
        !MS(i)%v(2) = MS(i)%v(2) + MS(i)%fdrag(2)/MS(i)%m*rhof*dtDEM  + MS(i)%fcol(2)/MS(i)%m*dt
        !MS(i)%v(3) = MS(i)%v(3) + MS(i)%fdrag(3)/MS(i)%m*rhof*dtDEM  + MS(i)%fcol(3)/MS(i)%m*dt - ratio*g*dt
        !MS(i)%omega(1) = MS(i)%omega(1) + MS(i)%fmom(1)/MS(i)%I*dtDEM*rhof
        !MS(i)%omega(2) = MS(i)%omega(2) + MS(i)%fmom(2)/MS(i)%I*dtDEM*rhof
        !MS(i)%omega(3) = MS(i)%omega(3) + MS(i)%fmom(3)/MS(i)%I*dtDEM*rhof
        !MS(i)%omega = 0d0
        
        MS(i)%vold = vTmp
        MS(i)%omegaold = omegaTmp

        MS(i)%x(1) = MS(i)%x(1) + MS(i)%v(1)*dtDEM
        MS(i)%x(2) = MS(i)%x(2) + MS(i)%v(2)*dtDEM
        MS(i)%x(3) = MS(i)%x(3) + MS(i)%v(3)*dtDEM
        call CalculateRotateMatrix(MS(i)%omega*dtDEM,MS(i)%rotMat)
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine ResetVelocityAndPositionLagrangianPoints()
use vars
use ZjsVar
implicit none
    integer::i,msIdx
    do i=1,npLocal
        msIdx = part(i)%MSIdx
        part(i)%x(1) = part(i)%x(1) + MSGlo(msIdx)%v(1)*dtDEM
        part(i)%x(2) = part(i)%x(2) + MSGlo(msIdx)%v(2)*dtDEM
        part(i)%x(3) = part(i)%x(3) + MSGlo(msIdx)%v(3)*dtDEM
        call RotateLagrangianPoints(part(i)%x,MSGlo(msIdx)%x,MSGlo(msIdx)%rotMat)
        call CrossMultiplication(MSGlo(msIdx)%omega,part(i)%x-MSGlo(msIdx)%x,part(i)%v)
        part(i)%v = part(i)%v + MSGlo(msIdx)%v
        !part(i)%v = MSGlo(msIdx)%v
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine RotateLagrangianPoints(xp,xcell,mat)
use vars
use ZjsVar
implicit none
    real*8,dimension(1:3)::xp,xcell,xtmp
    real*8,dimension(1:9)::mat
    xtmp = xp-xcell
    xp(1) = mat(1)*xtmp(1)+mat(2)*xtmp(2)+mat(3)*xtmp(3)
    xp(2) = mat(4)*xtmp(1)+mat(5)*xtmp(2)+mat(6)*xtmp(3)
    xp(3) = mat(7)*xtmp(1)+mat(8)*xtmp(2)+mat(9)*xtmp(3)
    xp = xp+xcell
return
endsubroutine
!----------------------------------------------------------
subroutine CalculateRotateMatrix(theta,mat)
use vars
use ZjsVar
implicit none
    real*8,dimension(1:3)::theta,norm
    real*8,dimension(1:9)::mat
    real*8::mag,ux,uy,uz      
    mat = 0d0                 
    mag = sqrt(theta(1)**2+theta(2)**2+theta(3)**2)
    mag = max(mag,1d-10)
    ux = theta(1)/mag; uy = theta(2)/mag; uz = theta(3)/mag
    mat(1) = (1-cos(mag))*ux*ux + cos(mag)
    mat(2) = (1-cos(mag))*ux*uy - sin(mag)*uz
    mat(3) = (1-cos(mag))*ux*uz + sin(mag)*uy
    mat(4) = (1-cos(mag))*uy*ux + sin(mag)*uz
    mat(5) = (1-cos(mag))*uy*uy + cos(mag)
    mat(6) = (1-cos(mag))*uy*uz - sin(mag)*ux
    mat(7) = (1-cos(mag))*uz*ux - sin(mag)*uy
    mat(8) = (1-cos(mag))*uz*uy + sin(mag)*ux
    mat(9) = (1-cos(mag))*uz*uz + cos(mag)
return
endsubroutine
!----------------------------------------------------------
subroutine CalculateCollisionForceMS()
use vars
use ZjsVar
    integer::i,j,ii,jj
    real*8::dist,rij,deltan,kn,gamman,eta,etan,meff
    real*8::vnij
    real*8::r1,r2,r3,n1,n2,n3,f1,f2,f3

    eta = -2*dlog(ResistCoeff)/sqrt(3.1415926536**2+dlog(ResistCoeff))
    do ii=1,nMSLocal
    MS(ii)%fcol(:) = 0
    do jj=ii+1,nMS
        i = ii
        j = jj
        r1 = MS(i)%x(1)-MS(j)%x(1)
        r2 = MS(i)%x(2)-MS(j)%x(2)
        r3 = MS(i)%x(3)-MS(j)%x(3)
        dist = sqrt(r1**2+r2**2+r3**2)
        n1 = r1/dist
        n2 = r2/dist
        n3 = r3/dist
        meff = 2*MS(i)%m*MS(j)%m/(MS(i)%m+MS(j)%m)
        rij = MS(i)%rad+MS(j)%rad+dx
        deltan = rij-dist
        kn = knCoeff/(0.1d0*di)*meff
        etan = eta*sqrt(meff*kn)

        f1 = 0; f2 = 0; f3 = 0;
        vnij = n1*(MS(i)%v(1)-MS(j)%v(1)) + &
               n2*(MS(i)%v(2)-MS(j)%v(2)) + &
               n3*(MS(i)%v(3)-MS(j)%v(3))
        etan = 0d0
        if (deltan>0)then
            f1 = n1*(kn*deltan-etan*vnij)
            f2 = n2*(kn*deltan-etan*vnij)
            f3 = n3*(kn*deltan-etan*vnij)
        endif
        MS(i)%fcol(1) = MS(i)%fcol(1) + f1
        MS(i)%fcol(2) = MS(i)%fcol(2) + f2
        MS(i)%fcol(3) = MS(i)%fcol(3) + f3
    enddo
    enddo
    
return
endsubroutine
!----------------------------------------------------------
subroutine CalculateSolidBoundaryCollisionForceMS()
use vars
use ZjsVar
    integer::ii,i,j,k,cellIdx1,cellIdx2,cellIdx3,flag
    real*8::dist,rij,deltan,kn,gamman,eta,etan,meff
    real*8::vnij
    real*8::xij1,xij2,xij3,r1,r2,r3,f1,f2,f3,fn1,fn2,fn3
    real*8::p1,p2,p3,n1,n2,n3,cp1,cp2,cp3
    
    knCoeff = 1d0*4
    eta = -2*dlog(ResistCoeff)/sqrt(3.1415926536**2+dlog(ResistCoeff))
    do ii=1,nMSLocal
        i = ii
        xij1 = MS(i)%x(1); xij2 = MS(i)%x(2); xij3 = MS(i)%x(3)
        do j=1,6
            if (bType(j)==1)then
                p1 = pointBound(1,j); p2 = pointBound(2,j); p3 = pointBound(3,j)
                n1 = normBound(1,j); n2 = normBound(2,j); n3 = normBound(3,j)
                r1 = xij1-p1; r2 = xij2-p2; r3 = xij3-p3
                dist = r1*n1+r2*n2+r3*n3
                cp1 = xij1-dist*n1
                cp2 = xij2-dist*n2
                cp3 = xij3-dist*n3
                meff = MS(i)%m
                rij = MS(i)%rad+dx
                deltan = rij-dist
                kn = knCoeff/(0.1d0*di)*meff
                etan = eta*sqrt(meff*kn)

                f1 = 0; f2 = 0; f3 = 0
                vnij = n1*MS(i)%v(1)+n2*MS(i)%v(2)+n3*MS(i)%v(3)
                !etan = 0d0
                if (deltan>0)then
                    f1 = n1*(kn*deltan-etan*vnij)
                    f2 = n2*(kn*deltan-etan*vnij)
                    f3 = n3*(kn*deltan-etan*vnij)
                    fn1 = n1*kn*deltan
                    fn2 = n2*kn*deltan
                    fn3 = n3*kn*deltan
                    MS(i)%fcol(1) = MS(i)%fcol(1) + f1
                    MS(i)%fcol(2) = MS(i)%fcol(2) + f2 
                    MS(i)%fcol(3) = MS(i)%fcol(3) + f3

                endif
            endif
        enddo
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine CloudEvolutionMS()
use vars
use ZjsVar
    integer::ierr,i,ncpuTot
    ncpuTot = cpunx*cpuny*cpunz

    !call RebuildAllListOrNot()
    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call CalculateCollisionForce()
    call CalculateSolidBoundaryCollisionForce()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)


    if (ncputot>1)call AssembleForceMS()
    if (ncputot==1)then
        call AssembleLocalForceMS()
        nMSlocal = 1
        MS(1)%fdrag = MSGlo(1)%fdrag
        MS(1)%fmom = MSGlo(1)%fmom
        MS(1)%fcol = MSGlo(1)%fcol
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call EvolutionMS()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (ncputot==1)MSGlo(1) = MS(1)
    
    if (ncputot>1)call UpdateMSInfoInGhostRegion()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (ncputot>1)call RebuildListMS()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "step",step,"After update ms info in ghost region"
    !do i=1,nMSLocal
    !    write(cpuLocal,*)"idx, x, v, omega", MS(i)%GlobalIdx,MS(i)%x,MS(i)%v,MS(i)%omega
    !    write(cpuLocal,*)"drag, mom", MS(i)%fdrag,MS(i)%fmom
    !    write(cpuLocal,*)"rot matrix", MS(i)%rotMat
    !enddo
    !write(cpuLocal,*)" " 
    !close(cpuLocal)

    call ResetVelocityAndPositionLagrangianPoints()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call AdjustParcelPosition()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    
    call UpdateParcelInfoInGhostRegion()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    
    call RebuildAllListOrNot()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call UpdateLocalCellWeightListDEM()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)  

    call UpdateParcelInfoInGhostRegion()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   

    call UpdateGhostCellWeightListDEM()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
    
    !open(cpuLocal,file=logFile,access="append")
    !do i=1,npLocal
    !    write(cpuLocal,*)"part Idx, x, v", part(i)%GlobalIdx,part(i)%x,part(i)%v,part(i)%fdrag
    !enddo
    !close(cpuLocal)
return
endsubroutine
!----------------------------------------------------------
