!-----------------------------------------------
subroutine CompressMS()
use vars
use ZjsVar
implicit none
    integer::i
    real*8,dimension(1:3)::xi
    
    nMSlocal = 0
    do i=1,nMS
        xi = MS(i)%x 
        if (xi(1)>=xStart.and.xi(1)<xEnd.and.xi(2)>=yStart.and.xi(2)<yEnd.and.xi(3)>=zStart.and.xi(3)<zEnd)then
            nMSlocal = nMSlocal+1
            MS(nMSlocal) = MS(i)
        endif
    enddo
    nMS = nMSlocal
return
endsubroutine
!-----------------------------------------------
subroutine RebuildListMS()
use vars
use ZjsVar
implicit none
    integer::ierr
    call CompressMS()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    nMSSend(1:6) = 0
    nMSRecv(1:6) = 0
    call PutMSIntoDomainMS(1,nMs,0)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call RebuildListForGhostRegionMS()
return
endsubroutine
!-----------------------------------------------
subroutine PutMSIntoDomainMS(nMSStart,nMSEnd,mode) 
use vars
use ZjsVar
implicit none
    integer::mode,nMSStart,nMSEnd,i,j,dim
    real*8,dimension(1:6)::boundCoord
    real*8,dimension(1:3)::centerCoord
    boundCoord(1) = xStart
    boundCoord(2) = xEnd
    boundCoord(3) = yStart
    boundCoord(4) = yEnd
    boundCoord(5) = zStart
    boundCoord(6) = zEnd
    centerCoord(1) = 0.5d0*(xStart+xEnd)
    centerCoord(2) = 0.5d0*(yStart+yEnd)
    centerCoord(3) = 0.5d0*(zStart+zEnd)

    do i=nMSStart,nMSEnd
    do j=2*mode+1,6
        dim = (j-1)/2+1
        if (MS(i)%x(dim)<=centerCoord(dim).and.mod(j,2)==1)then
            nMSSend(j) = nMSSend(j) + 1
            MSSend(nMSSend(j),j) = MS(i)
            MSSendIdx(nMSSend(j),j) = i
        endif
        if (MS(i)%x(dim)>centerCoord(dim).and.mod(j,2)==0)then
            nMSSend(j) = nMSSend(j) + 1
            MSSend(nMSSend(j),j) = MS(i)
            MSSendIdx(nMSSend(j),j) = i
        endif

        !open(cpuLocal,file=logFile,access="append")
        !write(cpuLocal,*) "Put MS into domain",i,j,dim,MS(i)%x(1),MS(i)%x(2),MS(i)%x(3),centerCoord(1),centerCoord(2),centerCoord(3)
        !write(cpuLocal,*) nMSSend
        !write(cpuLocal,*) nMSRecv
        !close(cpuLocal)

    enddo
    enddo
return 
endsubroutine
!-----------------------------------------------
subroutine RebuildListForGhostRegionMS()
use vars
use ZjsVar
implicit none
    integer::idx,i,nSend,nRecv,dim
    integer::idxSend,idxRecv,cpuSend,cpuRecv,ierr

    do idx=1,6
        if (idx==1.or.idx==3.or.idx==5)then
            idxSend = idx
            idxRecv = idx+1
        else
            idxSend = idx
            idxRecv = idx-1
        endif    
        
        if (idx==1.or.idx==2)then
            dim = 1
        else if (idx==3.or.idx==4)then
            dim = 2
        else
            dim = 3
        endif
        
        ! send parcels
        cpuSend = cpuLocal
        nSend = nMSSend(idx)
        if (bType(idx)==0.and.cpunx*cpuny*cpunz>1) then ! period
            do i=1,nSend
                MSSend(i,idxSend)%x(dim) = MSSend(i,idxSend)%x(dim) + coordShift(idxSend)
            enddo
        endif
        cpuSend = neighCpu(idxRecv) ! date: 20181211
        cpuRecv = neighCpu(idxSend) ! date: 20181211
        call MyMPISendAndRecvMS(cpuSend,cpuRecv,idxSend,idxRecv)

        ! put new parcels into cells
        nRecv = nMSRecv(idxRecv)
        call PutMSIntoDomainMS(nMS+1,nMS+nRecv,dim)
        startIdxMSRecv(idxRecv) = nMS+1
        nMS = nMS+nRecv
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    enddo
    do i=1,nMS
        MSGlo(MS(i)%globalIdx) = MS(i)
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine UpdateMSInfoInGhostRegion()
use vars
use ZjsVar
implicit none
    integer::idx,i,j,nSend,dim,ierr
    integer::idxSend,idxRecv,cpuSend,cpuRecv

    do idx=1,6
        if (idx==1.or.idx==3.or.idx==5)then
            idxSend = idx
            idxRecv = idx+1
        else
            idxSend = idx
            idxRecv = idx-1
        endif    
        ! send parcels
        nSend = nMSSend(idx)
        if (idx==1.or.idx==2)then
            dim = 1
        else if (idx==3.or.idx==4)then
            dim = 2
        else
            dim = 3
        endif
        do i=1,nSend
            !MSSend(i,idx)%x(:) = MS(MSSendIdx(i,idx))%x(:)
            !MSSend(i,idx)%v(:) = MS(MSSendIdx(i,idx))%v(:)
            !MSSend(i,idx)%omega = MS(MSSendIdx(i,idx))%omega
            !MSSend(i,idx)%fcol  = MS(MSSendIdx(i,idx))%fcol
            !MSSend(i,idx)%fmom  = MS(MSSendIdx(i,idx))%fmom
            !MSSend(i,idx)%fdrag  = MS(MSSendIdx(i,idx))%fdrag
            !MSSend(i,idx)%xold  = MS(MSSendIdx(i,idx))%xold
            !MSSend(i,idx)%vold  = MS(MSSendIdx(i,idx))%vold
            !MSSend(i,idx)%omegaold  = MS(MSSendIdx(i,idx))%omegaold
            !MSSend(i,idx)%rotMat  = MS(MSSendIdx(i,idx))%rotMat
            MSSend(i,idx) = MS(MSSendIdx(i,idx))
        enddo
        if (bType(idx)==0.and.cpunx*cpuny*cpunz>1) then ! period
            do i=1,nSend
                MSSend(i,idxSend)%x(dim) = MSSend(i,idxSend)%x(dim) + coordShift(idxSend)
            enddo
        endif
        cpuSend = neighCpu(idxRecv)
        cpuRecv = neighCpu(idxSend)
        call MyMPIUpdateVelocityAndPositionMS(cpuSend,cpuRecv,idxSend,idxRecv)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    enddo
    do i=1,nMS
        MSGlo(MS(i)%globalIdx) = MS(i)
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine MyMPISendAndRecvMS(cpuSend,cpuRecv,idxSend,idxRecv)
use vars
use ZjsVar
implicit none
    integer::cpuSend,cpuRecv,idxSend,idxRecv
    integer::request,info,request1,request2,info1,info2,ierr
    integer::req(1:2),stat(MPI_STATUS_SIZE,1:2)
    integer::nSend,nRecv,nSendInt,nSendReal,nRecvInt,nRecvReal
    integer::intOffset,realOffset,pIdx
    integer::i
    integer::numMSSendInt,numMSSendReal
    integer::STATUS(MPI_STATUS_SIZE),STATUS2(MPI_STATUS_SIZE),STATUS1(MPI_STATUS_SIZE)

    numMSSendInt = 2; numMSSendReal = 41+33+9
    nSend=1; nSendInt=1; nSendReal=1; nRecv=1; nRecvInt=1; nRecvReal=1
    nMSRecv(idxRecv) = 0
    if (cpuRecv/=-1)then
        nSend = nMSSend(idxSend)
        nSendInt = numMSSendInt*nSend
        nSendReal = numMSSendReal*nSend
        do i=1,nSend
            intOffset = numMSSendInt*(i-1)
            realOffset = numMSSendReal*(i-1)
            sendIntBuff(intOffset+1) = MSSend(i,idxSend)%globalIdx
            sendIntBuff(intOffset+2) = MSSend(i,idxSend)%cpuIdx
            sendRealBuff(realOffset+1) = MSSend(i,idxSend)%I
            sendRealBuff(realOffset+2) = MSSend(i,idxSend)%rad
            sendRealBuff(realOffset+3) = MSSend(i,idxSend)%rho
            sendRealBuff(realOffset+4) = MSSend(i,idxSend)%vol
            sendRealBuff(realOffset+5) = MSSend(i,idxSend)%m
            sendRealBuff(realOffset+6:realOffset+8) = MSSend(i,idxSend)%x
            sendRealBuff(realOffset+9:realOffset+11) = MSSend(i,idxSend)%v
            sendRealBuff(realOffset+12:realOffset+14) = MSSend(i,idxSend)%omega
            sendRealBuff(realOffset+15:realOffset+17) = MSSend(i,idxSend)%fcol
            sendRealBuff(realOffset+18:realOffset+20) = MSSend(i,idxSend)%fmom
            sendRealBuff(realOffset+21:realOffset+23) = MSSend(i,idxSend)%fdrag
            sendRealBuff(realOffset+24:realOffset+26) = MSSend(i,idxSend)%xold
            sendRealBuff(realOffset+27:realOffset+35) = MSSend(i,idxSend)%rotMat
            sendRealBuff(realOffset+36:realOffset+38) = MSSend(i,idxSend)%omegaold
            sendRealBuff(realOffset+39:realOffset+41) = MSSend(i,idxSend)%vold
            sendRealBuff(realOffset+42:realOffset+44) = MSSend(i,idxSend)%AngMom
            sendRealBuff(realOffset+45:realOffset+47) = MSSend(i,idxSend)%AngMomOld
            sendRealBuff(realOffset+48:realOffset+56) = MSSend(i,idxSend)%rotMatGlo
            sendRealBuff(realOffset+57:realOffset+65) = MSSend(i,idxSend)%rotMatGloOld
            sendRealBuff(realOffset+66:realOffset+74) = MSSend(i,idxSend)%InertiaTensor
            sendRealBuff(realOffset+75:realOffset+77) = MSSend(i,idxSend)%fhydro
            sendRealBuff(realOffset+78:realOffset+80) = MSSend(i,idxSend)%fmomHydro
            sendRealBuff(realOffset+81:realOffset+83) = MSSend(i,idxSend)%fmomCol
        enddo
    endif
    call MPI_SENDRECV(nMSSend(idxSend),1,MPI_INTEGER,CPUTransferIdx(idxSend),111,&
                      nMSRecv(idxRecv),1,MPI_INTEGER,CPUTransferIdx(idxRecv),111,MPI_COMM_WORLD,STATUS,info)
    if (cpuSend==-1) nMSRecv(idxRecv) = 0
    !if (cpuRecv/=-1)then
    !    call MPI_ISEND(nMSSend(idxSend),1,MPI_INTEGER,cpuRecv,111,MPI_COMM_WORLD,request,info)
    !    !call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuSend/=-1)then
    !    call MPI_IRECV(nMSRecv(idxRecv),1,MPI_INTEGER,cpuSend,111,MPI_COMM_WORLD,request,info)
    !    call MPI_WAIT(request,STATUS,info)
    !endif
    
    if (cpuSend/=-1)then
        nRecv = nMSRecv(idxRecv)
        nRecvInt = nRecv*numMSSendInt
        nRecvReal = nRecv*numMSSendReal
    endif
    call MPI_SENDRECV(sendIntBuff,nSendInt,MPI_INTEGER,CPUTransferIdx(idxSend),222,&
                      recvIntBuff,nRecvInt,MPI_INTEGER,CPUTransferIdx(idxRecv),222,MPI_COMM_WORLD,STATUS,info)
    call MPI_SENDRECV(sendRealBuff,nSendReal,MPI_DOUBLE_PRECISION,CPUTransferIdx(idxSend),333,&
                      recvRealBuff,nRecvReal,MPI_DOUBLE_PRECISION,CPUTransferIdx(idxRecv),333,MPI_COMM_WORLD,STATUS,info)
    !if (cpuRecv/=-1)then
    !    call MPI_ISEND(sendIntBuff,nSendInt,MPI_INTEGER,cpuRecv,222,MPI_COMM_WORLD,request,info)
    !    !call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuSend/=-1)then
    !    call MPI_IRECV(recvIntBuff,nRecvInt,MPI_INTEGER,cpuSend,222,MPI_COMM_WORLD,request,info)
    !    call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuRecv/=-1)then
    !    call MPI_ISEND(sendRealBuff,nSendReal,MPI_DOUBLE_PRECISION,cpuRecv,333,MPI_COMM_WORLD,request,info)
    !    !call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuSend/=-1)then
    !    call MPI_IRECV(recvRealBuff,nRecvReal,MPI_DOUBLE_PRECISION,cpuSend,333,MPI_COMM_WORLD,request,info)
    !    call MPI_WAIT(request,STATUS,info)
    !endif
    if (cpuSend/=-1)then
        do i=1,nRecv
            pIdx = nMS+i
            intOffset = numMSSendInt*(i-1)
            realOffset = numMSSendReal*(i-1)
            MS(pIdx)%globalIdx = recvIntBuff(intOffset+1)
            MS(pIdx)%cpuIdx    = recvIntBuff(intOffset+2)
            MS(pIdx)%I         = recvRealBuff(realOffset+1)
            MS(pIdx)%rad       = recvRealBuff(realOffset+2)
            MS(pIdx)%rho       = recvRealBuff(realOffset+3)
            MS(pIdx)%vol       = recvRealBuff(realOffset+4)
            MS(pIdx)%m         = recvRealBuff(realOffset+5)
            MS(pIdx)%x         = recvRealBuff(realOffset+6:realOffset+8)
            MS(pIdx)%v         = recvRealBuff(realOffset+9:realOffset+11)
            MS(pIdx)%omega     = recvRealBuff(realOffset+12:realOffset+14)
            MS(pIdx)%fcol      = recvRealBuff(realOffset+15:realOffset+17)        
            MS(pIdx)%fmom      = recvRealBuff(realOffset+18:realOffset+20)        
            MS(pIdx)%fdrag     = recvRealBuff(realOffset+21:realOffset+23)        
            MS(pIdx)%xold      = recvRealBuff(realOffset+24:realOffset+26)        
            MS(pIdx)%rotMat    = recvRealBuff(realOffset+27:realOffset+35)        
            MS(pIdx)%omegaold  = recvRealBuff(realOffset+36:realOffset+38)        
            MS(pIdx)%vold      = recvRealBuff(realOffset+39:realOffset+41)        
            MS(pIdx)%AngMom    = recvRealBuff(realOffset+42:realOffset+44)        
            MS(pIdx)%AngMomOld = recvRealBuff(realOffset+45:realOffset+47)        
            MS(pIdx)%rotMatGlo = recvRealBuff(realOffset+48:realOffset+56)        
            MS(pIdx)%rotMatGloOld = recvRealBuff(realOffset+57:realOffset+65)        
            MS(pIdx)%InertiaTensor= recvRealBuff(realOffset+66:realOffset+74)        
            MS(pIdx)%fhydro       = recvRealBuff(realOffset+75:realOffset+77)        
            MS(pIdx)%fmomHydro    = recvRealBuff(realOffset+78:realOffset+80)        
            MS(pIdx)%fmomCol      = recvRealBuff(realOffset+81:realOffset+83)        
        enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
return
endsubroutine
!-----------------------------------------------
subroutine MyMPIUpdateVelocityAndPositionMS(cpuSend,cpuRecv,idxSend,idxRecv)
use vars
use ZjsVar
implicit none
    integer::cpuSend,cpuRecv,idxSend,idxRecv
    integer::request,info,request1,request2,info1,info2,ierr
    integer::nSend,nRecv,nSendInt,nSendReal,nRecvInt,nRecvReal
    integer::intOffset,realOffset,pIdx
    integer::i,numMSSendReal
    integer::STATUS(MPI_STATUS_SIZE)
    nRecv = 0; numMSSendReal = 36+33+9
    nSend=1; nSendInt=1; nSendReal=1; nRecv=1; nRecvInt=1; nRecvReal=1
    if (cpuRecv/=-1)then
        nSend = nMSSend(idxSend)
        nSendReal = numMSSendReal*nSend
        do i=1,nSend
            pIdx = startIdxMSRecv(idxSend)+i-1
            realOffSet = numMSSendReal*(i-1)
            sendRealBuff(realOffset+1:realOffset+3) = MSSend(i,idxSend)%x
            sendRealBuff(realOffset+4:realOffset+6) = MSSend(i,idxSend)%v
            sendRealBuff(realOffset+7:realOffset+9) = MSSend(i,idxSend)%omega
            sendRealBuff(realOffset+10:realOffset+12) = MSSend(i,idxSend)%fcol
            sendRealBuff(realOffset+13:realOffset+15) = MSSend(i,idxSend)%fmom
            sendRealBuff(realOffset+16:realOffset+18) = MSSend(i,idxSend)%fdrag
            sendRealBuff(realOffset+19:realOffset+21) = MSSend(i,idxSend)%xold
            sendRealBuff(realOffset+22:realOffset+30) = MSSend(i,idxSend)%rotMat
            sendRealBuff(realOffset+31:realOffset+33) = MSSend(i,idxSend)%omegaold
            sendRealBuff(realOffset+34:realOffset+36) = MSSend(i,idxSend)%vold
            sendRealBuff(realOffset+37:realOffset+39) = MSSend(i,idxSend)%AngMom
            sendRealBuff(realOffset+40:realOffset+42) = MSSend(i,idxSend)%AngMomOld
            sendRealBuff(realOffset+43:realOffset+51) = MSSend(i,idxSend)%rotMatGlo
            sendRealBuff(realOffset+52:realOffset+60) = MSSend(i,idxSend)%rotMatGloOld
            sendRealBuff(realOffset+61:realOffset+69) = MSSend(i,idxSend)%InertiaTensor
            sendRealBuff(realOffset+70:realOffset+72) = MSSend(i,idxSend)%fhydro
            sendRealBuff(realOffset+73:realOffset+75) = MSSend(i,idxSend)%fmomHydro
            sendRealBuff(realOffset+76:realOffset+78) = MSSend(i,idxSend)%fmomCol
        enddo
    endif
    if (cpuSend/=-1)then
        nRecv = nMSRecv(idxRecv)
        nRecvReal = nRecv*numMSSendReal
    endif
    
    call MPI_SENDRECV(sendRealBuff,nSendReal,MPI_DOUBLE_PRECISION,CPUTransferIdx(idxSend),111,&
                      recvRealBuff,nRecvReal,MPI_DOUBLE_PRECISION,CPUTransferIdx(idxRecv),111,MPI_COMM_WORLD,STATUS,info)
    !if (cpuRecv/=-1)then
    !    call MPI_ISEND(sendRealBuff,nSendReal,MPI_DOUBLE_PRECISION,cpuRecv,111,MPI_COMM_WORLD,request,info)
    !    !call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuSend/=-1)then
    !    call MPI_IRECV(recvRealBuff,nRecvReal,MPI_DOUBLE_PRECISION,cpuSend,111,MPI_COMM_WORLD,request,info)
    !    call MPI_WAIT(request,STATUS,info)
    !endif
    if (cpuSend/=-1)then
        do i=1,nRecv
            pIdx = startIdxMSRecv(idxRecv)+i-1
            realOffset = numMSSendReal*(i-1)
            MS(pIdx)%x         = recvRealBuff(realOffset+1:realOffset+3)
            MS(pIdx)%v         = recvRealBuff(realOffset+4:realOffset+6)
            MS(pIdx)%omega     = recvRealBuff(realOffset+7:realOffset+9)
            MS(pIdx)%fcol      = recvRealBuff(realOffset+10:realOffset+12)
            MS(pIdx)%fmom      = recvRealBuff(realOffset+13:realOffset+15)
            MS(pIdx)%fdrag     = recvRealBuff(realOffset+16:realOffset+18)
            MS(pIdx)%xold      = recvRealBuff(realOffset+19:realOffset+21)
            MS(pIdx)%rotMat    = recvRealBuff(realOffset+22:realOffset+30)
            MS(pIdx)%omegaold  = recvRealBuff(realOffset+31:realOffset+33)
            MS(pIdx)%vold      = recvRealBuff(realOffset+34:realOffset+36)
            MS(pIdx)%AngMom    = recvRealBuff(realOffset+37:realOffset+39)
            MS(pIdx)%AngMomOld = recvRealBuff(realOffset+40:realOffset+42)
            MS(pIdx)%rotMatGlo = recvRealBuff(realOffset+43:realOffset+51)
            MS(pIdx)%rotMatGloOld = recvRealBuff(realOffset+52:realOffset+60)
            MS(pIdx)%InertiaTensor= recvRealBuff(realOffset+61:realOffset+69)
            MS(pIdx)%fhydro       = recvRealBuff(realOffset+70:realOffset+72)
            MS(pIdx)%fmomHydro    = recvRealBuff(realOffset+73:realOffset+75)
            MS(pIdx)%fmomCol      = recvRealBuff(realOffset+76:realOffset+78)
        enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
return
endsubroutine
!-----------------------------------------------
subroutine AssembleForceMS()
use vars
use ZjsVar
    integer::idx,i,nSend,nRecv,dim
    integer::idxSend,idxRecv,cpuSend,cpuRecv,ierr
    call AssembleLocalForceMS()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !open(cpuLocal,file=logFile,access="append")
    !write(cpuLocal,*) "step ",step,nMS
    !write(cpuLocal,*) nMSSend(1),nMSSend(2),nMSRecv(1),nMSRecv(2)
    !close(cpuLocal)
    do idx=6,1,-1
        if (idx==1.or.idx==3.or.idx==5)then
            idxSend = idx
            idxRecv = idx+1
        else
            idxSend = idx
            idxRecv = idx-1
        endif    
        ! send parcels
        cpuSend = cpuLocal
        nSend = nMSSend(idx)
        if (idx==1.or.idx==2)then
            dim = 1
        else if (idx==3.or.idx==4)then
            dim = 2
        else
            dim = 3
        endif
        cpuSend = neighCpu(idxRecv) 
        cpuRecv = neighCpu(idxSend) 
        call MyMPIAssembleForceMS(cpuSend,cpuRecv,idxSend,idxRecv)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    enddo

    do i=1,nMSlocal
        MS(i)%fcol = MSGlo(MS(i)%GlobalIdx)%fcol
        MS(i)%fmom = MSGlo(MS(i)%GlobalIdx)%fmom
        MS(i)%fmomCol = MSGlo(MS(i)%GlobalIdx)%fmomCol
        MS(i)%fdrag = MSGlo(MS(i)%GlobalIdx)%fdrag
    enddo
    !do i=1,nMSGlo
    !    open(cpuLocal,file=logFile,access="append")
    !    write(cpuLocal,*) "Tot mom and force:", MSGlo(i)%fmom,MSGlo(i)%fdrag
    !    close(cpuLocal)
    !enddo
return
endsubroutine
!-----------------------------------------------
subroutine MyMPIAssembleForceMS(cpuSend,cpuRecv,idxSend,idxRecv)
use vars
use ZjsVar
implicit none
    integer::cpuSend,cpuRecv,idxSend,idxRecv
    integer::request,info,request1,request2,info1,info2,ierr
    integer::nSend,nRecv,nSendInt,nSendReal,nRecvInt,nRecvReal
    integer::intOffset,realOffset,pIdx,MSIdx
    integer::i,numMSSendReal
    integer::STATUS(MPI_STATUS_SIZE)
    nRecv = 0; numMSSendReal = 18+3
    nSend=1; nSendInt=1; nSendReal=1; nRecv=1; nRecvInt=1; nRecvReal=1
    if (cpuRecv/=-1)then
        nSend = nMSRecv(idxSend)
        nSendReal = numMSSendReal*nSend
        do i=1,nSend
            pIdx = startIdxMSRecv(idxSend)+i-1
            MSIdx = MS(pIdx)%globalIdx
            realOffSet = numMSSendReal*(i-1)
            sendRealBuff(realOffset+1:realOffset+3) = MSGlo(MSIdx)%x
            sendRealBuff(realOffset+4:realOffset+6) = MSGlo(MSIdx)%v
            sendRealBuff(realOffset+7:realOffset+9) = MSGlo(MSIdx)%omega
            sendRealBuff(realOffset+10:realOffset+12) = MSGlo(MSIdx)%fcol
            sendRealBuff(realOffset+13:realOffset+15) = MSGlo(MSIdx)%fmom
            sendRealBuff(realOffset+16:realOffset+18) = MSGlo(MSIdx)%fdrag
            sendRealBuff(realOffset+19:realOffset+21) = MSGlo(MSIdx)%fmomCol
        enddo
    endif
    if (cpuSend/=-1)then
        nRecv = nMSSend(idxRecv)
        nRecvReal = nRecv*numMSSendReal
    endif
    
    call MPI_SENDRECV(sendRealBuff,nSendReal,MPI_DOUBLE_PRECISION,CPUTransferIdx(idxSend),111,&
                      recvRealBuff,nRecvReal,MPI_DOUBLE_PRECISION,CPUTransferIdx(idxRecv),111,MPI_COMM_WORLD,STATUS,info)

    !if (cpuRecv/=-1)then
    !    call MPI_ISEND(sendRealBuff,nSendReal,MPI_DOUBLE_PRECISION,cpuRecv,111,MPI_COMM_WORLD,request,info)
    !    !call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuSend/=-1)then
    !    call MPI_IRECV(recvRealBuff,nRecvReal,MPI_DOUBLE_PRECISION,cpuSend,111,MPI_COMM_WORLD,request,info)
    !    call MPI_WAIT(request,STATUS,info)
    !endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (cpuSend/=-1)then
        do i=1,nRecv
            MSIdx = MSSend(i,idxRecv)%globalIdx
            realOffset = numMSSendReal*(i-1)
            MSGlo(MSIdx)%fcol  = MSGlo(MSIdx)%fcol  + recvRealBuff(realOffset+10:realOffset+12)
            MSGlo(MSIdx)%fmom  = MSGlo(MSIdx)%fmom  + recvRealBuff(realOffset+13:realOffset+15)
            MSGlo(MSIdx)%fdrag = MSGlo(MSIdx)%fdrag + recvRealBuff(realOffset+16:realOffset+18)
            MSGlo(MSIdx)%fmomCol = MSGlo(MSIdx)%fmomCol + recvRealBuff(realOffset+19:realOffset+21)
        enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
return
endsubroutine
!-----------------------------------------------
subroutine AssembleLocalForceMS()
use vars
use ZjsVar
implicit none
   integer::i
   real*8::momTmp(1:3),momCol(1:3),rpc(1:3)
   do i=1,nMSglo
      MSGlo(i)%fdrag = 0d0
      MSGlo(i)%fcol = 0d0
      MSGlo(i)%fmom = 0d0
      MSGlo(i)%fmomCol = 0d0
   enddo
   do i=1,npLocal
      MSGlo(part(i)%MSIdx)%fdrag = MSGlo(part(i)%MSIdx)%fdrag + part(i)%fdrag*part(i)%vol*rhof
      MSGlo(part(i)%MSIdx)%fcol = MSGlo(part(i)%MSIdx)%fcol + part(i)%fcol
      rpc = part(i)%x-MSGlo(part(i)%MSIdx)%x
      call CrossMultiplication(rpc,part(i)%fdrag*part(i)%vol,momTmp)
      call CrossMultiplication(rpc,part(i)%fcol,momCol)
      MSGlo(part(i)%MSIdx)%fmom = MSGlo(part(i)%MSIdx)%fmom + momTmp
      MSGlo(part(i)%MSIdx)%fmomCol = MSGlo(part(i)%MSIdx)%fmomCol + momCol
   enddo
   !open(cpuIdx,file=logFile,access="append")
   !write(cpuIdx,*) "MS coords and mom:",MSGlo(1)%x,MSGlo(1)%fdrag,MSGlo(1)%fmom
   !close(cpuIdx)
return
endsubroutine
!-----------------------------------------------


