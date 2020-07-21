!----------------------------------------------------------
subroutine EvolutionMS()
use vars
use ZjsVar
implicit none
    integer::i
    real*8::ratio,vTmp(1:3),omegaTmp(1:3),angMomTmp(1:3)
    real*8::ITInv(1:9) ! Inverse matrix of inertia tensor
    real*8::RMInv(1:9) ! Inverse matrix of rotation matrix
    real*8::TMInvTmp(1:9),TMInv(1:9) ! Inertia Tensor matrix inverse
    real*8::g(1:3)
    g = 0d0
    g(3) = -9.81d0
    do i=1,nMSLocal
        ratio = (MS(i)%rho-rhof)/MS(i)%rho
        MS(i)%rotMatGloOld = MS(i)%rotMatGlo
        vTmp = MS(i)%v
        omegaTmp = MS(i)%omega
        angMomTmp = MS(i)%AngMom

        MS(i)%v = MS(i)%v + MS(i)%fhydro*dtDEM + MS(i)%fcol/MS(i)%m*dtDEM + ratio*g*dtDEM

        call CalculateRotateMatrix(MS(i)%omega*dtDEM,MS(i)%rotMat)
        ITInv(1:9) = MS(i)%InertiaTensor
        ITInv(1) = 1d0/ITInv(1); ITInv(5) = 1d0/ITInv(5); ITInv(9) = 1d0/ITInv(9)
        call TensorMultiplication(MS(i)%rotMatGlo,MS(i)%rotMat,MS(i)%rotMatGlo)
        call InverseOrthogonalMatrix(MS(i)%rotMatGlo,RMInv)
        call TensorMultiplication(MS(i)%rotMatGlo,ITInv,TMInvTmp)
        call TensorMultiplication(TMInvTmp,RMInv,TMInv)
        MS(i)%AngMom = MS(i)%AngMom + (MS(i)%fmomHydro + MS(i)%fmomCol)*dtDEM
        call TensorVectorMultiplication(TMInv,MS(i)%AngMom,MS(i)%omega)

        !MS(i)%vold = vTmp
        !MS(i)%omegaold = omegaTmp
        !MS(i)%AngMomOld = angMomTmp

        MS(i)%x(1) = MS(i)%x(1) + MS(i)%v(1)*dtDEM
        MS(i)%x(2) = MS(i)%x(2) + MS(i)%v(2)*dtDEM
        MS(i)%x(3) = MS(i)%x(3) + MS(i)%v(3)*dtDEM
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
subroutine CloudEvolutionStepMS()
use vars
use ZjsVar
    integer::ierr,i,ncpuTot
    ncpuTot = cpunx*cpuny*cpunz

    !call RebuildAllListOrNot()
    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call CalculateCollisionForce()
    call CalculateSolidBoundaryCollisionForce()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 112"
    !close(cpuLocal)
    if (ncputot>1)call AssembleForceMS()
    if (ncputot==1)then
        call AssembleLocalForceMS()
        nMSlocal = 1
        MS(1)%fdrag = MSGlo(1)%fdrag
        MS(1)%fmom = MSGlo(1)%fmom
        MS(1)%fcol = MSGlo(1)%fcol
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 126"
    !close(cpuLocal)
    call EvolutionMS()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (ncputot==1)MSGlo(1) = MS(1)
    
    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 133"
    !close(cpuLocal)
    if (ncputot>1)call UpdateMSInfoInGhostRegion()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 139"
    !close(cpuLocal)
    if (ncputot>1)call RebuildListMS()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 145"
    !close(cpuLocal)
    call ResetVelocityAndPositionLagrangianPoints()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 151"
    !close(cpuLocal)
    call AdjustParcelPosition()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    
    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 157"
    !close(cpuLocal)
    call UpdateParcelInfoInGhostRegion()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    
    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 163"
    !close(cpuLocal)
    call RebuildAllListOrNot()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 169"
    !close(cpuLocal)
    call UpdateLocalCellWeightListDEM()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)  

    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 175"
    !close(cpuLocal)
    call UpdateParcelInfoInGhostRegion()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   

    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "line 181"
    !close(cpuLocal)
    call UpdateGhostCellWeightListDEM()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
    
return
endsubroutine
!----------------------------------------------------------
subroutine CloudEvolutionMS()
use vars
use ZjsVar
    integer::i,ierr
    real*8::ratio
    do i=1,nMSLocal
        MS(i)%fhydro = 0d0
        MS(i)%fmomHydro = 0d0
    enddo
    do i=1,nMSLocal
        ratio = (MS(i)%rho-rhof)/MS(i)%rho
        MS(i)%fhydro = (MS(i)%v-MS(i)%vold)/dt*(1d0-ratio)
        MS(i)%fmomHydro = (MS(i)%AngMom-MS(i)%AngMomOld)/dt*(1d0-ratio)
        MS(i)%vold = MS(i)%v
        MS(i)%AngMomOld = MS(i)%AngMom
    enddo
    call AssembleForceMS()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    do i=1,nMSLocal
        MS(i)%fhydro = MS(i)%fhydro + MS(i)%fdrag/MS(i)%m
        MS(i)%fmomHydro = MS(i)%fmomhydro + MS(i)%fmom
    enddo
        
    do i=1,coupleInterval
        call CloudEvolutionStepMS()
    enddo
return
endsubroutine
!----------------------------------------------------------
