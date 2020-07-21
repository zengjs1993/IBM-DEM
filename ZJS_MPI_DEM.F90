!-----------------------------------------------
subroutine MyMPISendAndRecv(cpuSend,cpuRecv,idxSend,idxRecv)
use vars
use ZjsVar
implicit none
    integer::cpuSend,cpuRecv,idxSend,idxRecv
    integer::request,info,request1,request2,info1,info2,ierr
    integer::req(1:2),stat(MPI_STATUS_SIZE,1:2)
    integer::nSend,nRecv,nSendInt,nSendReal,nRecvInt,nRecvReal
    integer::intOffset,realOffset,pIdx
    integer::i
    integer::STATUS(MPI_STATUS_SIZE),STATUS2(MPI_STATUS_SIZE),STATUS1(MPI_STATUS_SIZE)

    npRecv(idxRecv) = 0
    nSend=1; nSendInt=1; nSendReal=1; nRecv=1; nRecvInt=1; nRecvReal=1
    if (cpuRecv/=-1)then
        nSend = npSend(idxSend)
        nSendInt = numSendInt*nSend
        nSendReal = numSendReal*nSend
        do i=1,nSend
            intOffset = numSendInt*(i-1)
            realOffset = numSendReal*(i-1)
            sendIntBuff(intOffset+1) = partSend(i,idxSend)%globalIdx
            sendIntBuff(intOffset+2) = partSend(i,idxSend)%localIdx
            sendIntBuff(intOffset+3) = partSend(i,idxSend)%cellIdx
            sendIntBuff(intOffset+4) = partSend(i,idxSend)%pType
            sendIntBuff(intOffset+5) = partSend(i,idxSend)%MSIdx
            sendRealBuff(realOffset+1) = partSend(i,idxSend)%radeff
            sendRealBuff(realOffset+2) = partSend(i,idxSend)%rad
            sendRealBuff(realOffset+3) = partSend(i,idxSend)%rho
            sendRealBuff(realOffset+4) = partSend(i,idxSend)%vol
            sendRealBuff(realOffset+5) = partSend(i,idxSend)%m
            sendRealBuff(realOffset+6:realOffset+8) = partSend(i,idxSend)%x
            sendRealBuff(realOffset+9:realOffset+11) = partSend(i,idxSend)%v
            sendRealBuff(realOffset+12:realOffset+14) = partSend(i,idxSend)%fcol
            sendRealBuff(realOffset+15:realOffset+17) = partSend(i,idxSend)%fdrag
        enddo
    endif
    
    call MPI_SENDRECV(npSend(idxSend),1,MPI_INTEGER,CPUTransferIdx(idxSend),111,&
                      npRecv(idxRecv),1,MPI_INTEGER,CPUTransferIdx(idxRecv),111,MPI_COMM_WORLD,STATUS,info)
    if (cpuSend==-1) npRecv(idxRecv) = 0
    !if (cpuRecv/=-1)then
    !    call MPI_ISEND(npSend(idxSend),1,MPI_INTEGER,cpuRecv,111,MPI_COMM_WORLD,request,info)
    !    !call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuSend/=-1)then
    !    call MPI_IRECV(npRecv(idxRecv),1,MPI_INTEGER,cpuSend,111,MPI_COMM_WORLD,request,info)
    !    call MPI_WAIT(request,STATUS,info)
    !endif
    
    if (cpuSend/=-1)then
        nRecv = npRecv(idxRecv)
        nRecvInt = nRecv*numSendInt
        nRecvReal = nRecv*numSendReal
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
            pIdx = np+i
            intOffset = numSendInt*(i-1)
            realOffset = numSendReal*(i-1)
            part(pIdx)%globalIdx = recvIntBuff(intOffset+1)
            part(pIdx)%localIdx  = recvIntBuff(intOffset+2)
            part(pIdx)%cellIdx   = recvIntBuff(intOffset+3)
            part(pIdx)%pType     = recvIntBuff(intOffset+4)
            part(pIdx)%MSIdx     = recvIntBuff(intOffset+5)
            part(pIdx)%radeff    = recvRealBuff(realOffset+1)
            part(pIdx)%rad       = recvRealBuff(realOffset+2)
            part(pIdx)%rho       = recvRealBuff(realOffset+3)
            part(pIdx)%vol       = recvRealBuff(realOffset+4)
            part(pIdx)%m         = recvRealBuff(realOffset+5)
            part(pIdx)%x         = recvRealBuff(realOffset+6:realOffset+8)
            part(pIdx)%v         = recvRealBuff(realOffset+9:realOffset+11)
            part(pIdx)%fcol      = recvRealBuff(realOffset+12:realOffset+14)
            part(pIdx)%fdrag     = recvRealBuff(realOffset+15:realOffset+17)        
        enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
return
endsubroutine
!-----------------------------------------------
subroutine MyMPIUpdateVelocityAndPosition(cpuSend,cpuRecv,idxSend,idxRecv)
use vars
use ZjsVar
implicit none
    integer::cpuSend,cpuRecv,idxSend,idxRecv
    integer::request,info,request1,request2,info1,info2,ierr
    integer::nSend,nRecv,nSendInt,nSendReal,nRecvInt,nRecvReal
    integer::intOffset,realOffset,pIdx
    integer::i
    integer::STATUS(MPI_STATUS_SIZE)
    nSend=1; nSendInt=1; nSendReal=1; nRecv=1; nRecvInt=1; nRecvReal=1
    if (cpuRecv/=-1)then
        nSend = npSend(idxSend)
        nSendInt = 1*nSend
        nSendReal = 10*nSend
        do i=1,nSend
            pIdx = startIdxRecv(idxSend)+i-1
            realOffSet = 10*(i-1)
            intOffSet  = 1*(i-1)
            sendIntBuff(intOffset+1) = partSend(i,idxSend)%MSIdx
            sendRealBuff(realOffset+1:realOffset+3) = partSend(i,idxSend)%x(1:3)
            sendRealBuff(realOffset+4:realOffset+6) = partSend(i,idxSend)%v(1:3)
            sendRealBuff(realOffset+7)              = partSend(i,idxSend)%wsum
            sendRealBuff(realOffset+8:realOffset+10)= partSend(i,idxSend)%fdrag(1:3)
        enddo
    endif
    if (cpuSend/=-1)then
        nRecv = npRecv(idxRecv)
        nRecvInt = nRecv*1
        nRecvReal = nRecv*10
    endif
    !do i=1,nSend
    !    realOffset = 6*(i-1)
    !    write(cpuLocal,*) sendRealBuff(realOffset+1),sendRealBuff(realOffset+2),sendRealBuff(realOffset+3)
    !enddo
    call MPI_SENDRECV(sendIntBuff,nSendInt,MPI_INTEGER,CPUTransferIdx(idxSend),222,&
                      recvIntBuff,nRecvInt,MPI_INTEGER,CPUTransferIdx(idxRecv),222,MPI_COMM_WORLD,STATUS,info)
    call MPI_SENDRECV(sendRealBuff,nSendReal,MPI_DOUBLE_PRECISION,CPUTransferIdx(idxSend),111,&
                      recvRealBuff,nRecvReal,MPI_DOUBLE_PRECISION,CPUTransferIdx(idxRecv),111,MPI_COMM_WORLD,STATUS,info)
    !if (cpuRecv/=-1)then
    !    call MPI_ISEND(sendIntBuff,nSendInt,MPI_INTEGER,cpuRecv,222,MPI_COMM_WORLD,request,info)
    !    !call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuSend/=-1)then
    !    call MPI_IRECV(recvIntBuff,nRecvInt,MPI_INTEGER,cpuSend,222,MPI_COMM_WORLD,request,info)
    !    call MPI_WAIT(request,STATUS,info)
    !endif
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
            pIdx = startIdxRecv(idxRecv)+i-1
            realOffset = 10*(i-1)
            intOffset = 1*(i-1)
            part(pIdx)%MSIdx     = recvIntBuff(intOffset+1)
            part(pIdx)%x         = recvRealBuff(realOffset+1:realOffset+3)
            part(pIdx)%v         = recvRealBuff(realOffset+4:realOffset+6)
            part(pIdx)%wsum      = recvRealBuff(realOffset+7)
            part(pIdx)%fdrag     = recvRealBuff(realOffset+8:realOffset+10)
        enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
return
endsubroutine
!-----------------------------------------------
subroutine UpdateFullInfoForSolidOutletBoundary(dim,idx)
use vars
use ZjsVar
implicit none
    integer::idx,dim,nSend,nRecv,i,lapsize
    real*8::stepsize,signlabel
    if (dim==1)then
        stepsize = di; lapsize = haloX
    elseif (dim==2)then
        stepsize = dj; lapsize = haloY
    else
        stepsize = dk; lapsize = haloZ
    endif
    if (mod(idx,2)==0)then
        signlabel = 1d0
    else
        signlabel = -1d0
    endif
    nRecv = npSend(idx)
    do i=1,nRecv
        part(np+i) = part(partSendIdx(i,idx))
    enddo
    npRecv(idx) = nRecv
    if (bType(idx)==1)then ! solid
        do i=1,nRecv
            part(np+i)%x(dim) = 2*coordMirror(idx)-part(np+i)%x(dim)
            part(np+i)%v(dim) = -part(np+i)%v(dim)
            part(np+i)%pType = 1
        enddo
    endif
    if (bType(idx)==2)then ! outlet
        do i=1,nRecv
            !part(np+i)%x(dim) = part(np+i)%x(dim) + signlabel*stepsize*lapsize
            part(np+i)%pType = 2
            part(np+i)%x(dim) = 2*coordMirror(idx)-part(np+i)%x(dim)
        enddo
    endif
    !if (step>=4853)then
    !    open(cpuLocal,file=logFile,access="append")
    !    write(cpuLocal,*) "recving no:",nRecv,"side:",idx,"dim:",dim
    !    do i=1,nRecv
    !        write(cpuLocal,*)     i,np+i,partSendIdx(i,idx), &
    !            part(partSendIdx(i,idx))%x(1),part(partSendIdx(i,idx))%x(2),part(partSendIdx(i,idx))%x(3) &
    !            ,part(np+i)%x(1),part(np+i)%x(2),part(np+i)%x(3)
    !    enddo
    !endif
return 
endsubroutine
!-----------------------------------------------
subroutine UpdatePartInfoForSolidOutletBoundary(dim,idx)
use vars
use ZjsVar
implicit none
    integer::idx,dim,nSend,nRecv,i,lapsize,pIdx,idxRecv
    real*8::stepsize,signlabel
!    open(cpuLocal,file=logFile,access="append")
!    write(cpuLocal,*) "side",idx,"recv no",npSend(1),npSend(2),npSend(3),npSend(4),npSend(5),npSend(6)
    if (dim==1)then
        stepsize = di; lapsize = haloX
    elseif (dim==2)then
        stepsize = dj; lapsize = haloY
    else
        stepsize = dk; lapsize = haloZ
    endif
    if (mod(idx,2)==0)then
        signlabel = 1d0
    else
        signlabel = -1d0
    endif
    nRecv = npSend(idx)
    idxRecv = idx
    do i=1,nRecv
        pIdx = startIdxRecv(idxRecv)+i-1
        part(pIdx)%x(:) = part(partSendIdx(i,idx))%x(:)
        part(pIdx)%v(:) = part(partSendIdx(i,idx))%v(:)
        part(pIdx)%fdrag(:) = part(partSendIdx(i,idx))%fdrag(:)
        part(pIdx)%wsum = part(partSendIdx(i,idx))%wsum
!        write(cpuLocal,*) i, part(partSendIdx(i,idx))%wsum, part(pIdx)%wsum
    enddo
    if (bType(idx)==1)then ! solid
        do i=1,nRecv
            pIdx = startIdxRecv(idxRecv)+i-1
            part(pIdx)%x(dim) = 2*coordMirror(idx)-part(pIdx)%x(dim)
            part(pIdx)%v(dim) = -part(pIdx)%v(dim)
            part(pIdx)%pType = 1
        enddo
    endif
    if (bType(idx)==2)then ! outlet
        do i=1,nRecv
            pIdx = startIdxRecv(idxRecv)+i-1
            !part(pIdx)%x(dim) = part(pIdx)%x(dim) + signlabel*stepsize*lapsize
            part(pIdx)%pType = 2
            part(pIdx)%x(dim) = 2*coordMirror(idx)-part(pIdx)%x(dim)
        enddo
    endif
return 
endsubroutine
!-----------------------------------------------
!-----------------------------------------------
