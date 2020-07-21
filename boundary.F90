!----------------------------------------------------------
subroutine BuildMPIRecvAndSendIndexForFluid()
use vars
implicit none
    integer::sendS1,sendS2,sendS3,sendE1,sendE2,sendE3
    integer::recvS1,recvS2,recvS3,recvE1,recvE2,recvE3
    integer::i,j,k,ii,inc1,inc2,inc3,cnt,idx
    
    sendS1 = -1; sendS2 = -1; sendS3 = -1
    sendE1 = -1; sendE2 = -1; sendE3 = -1
    recvS1 = -1; recvS2 = -1; recvS3 = -1
    recvE1 = -1; recvE2 = -1; recvE3 = -1
    do i=1,6
        if (i==1)then
            sendS1 = sg1; sendE1 = sg1+lap-1; sendS2 = sg2; sendE2 = eg2; sendS3 = sg3; sendE3 = eg3
        endif
        if (i==2)then
            sendS1 = eg1-lap+1; sendE1 = eg1; sendS2 = sg2; sendE2 = eg2; sendS3 = sg3; sendE3 = eg3
        endif
        if (i==3)then
            sendS1 = lg1; sendE1 = ug1; sendS2 = sg2; sendE2 = sg2+lap-1; sendS3 = sg3; sendE3 = eg3
        endif
        if (i==4)then
            sendS1 = lg1; sendE1 = ug1; sendS2 = eg2-lap+1; sendE2 = eg2; sendS3 = sg3; sendE3 = eg3
        endif
        if (i==5)then
            sendS1 = lg1; sendE1 = ug1; sendS2 = lg2; sendE2 = ug2; sendS3 = sg3; sendE3 = sg3+lap-1
        endif
        if (i==6)then
            sendS1 = lg1; sendE1 = ug1; sendS2 = lg2; sendE2 = ug2; sendS3 = eg3-lap+1; sendE3 = eg3
        endif

        inc1 = 1; inc2 = 1; inc3 = 1
        if (i==1)then
            if (bw==1)then
                recvS1 = sg1-1; recvE1 = lg1; recvS2 = sg2; recvE2 = eg2; recvS3 = sg3; recvE3 = eg3; inc1 = -1
            else
                recvS1 = lg1; recvE1 = sg1-1; recvS2 = sg2; recvE2 = eg2; recvS3 = sg3; recvE3 = eg3
            endif
        endif
        if (i==2)then
            if (be==1)then
                recvS1 = ug1; recvE1 = eg1+1; recvS2 = sg2; recvE2 = eg2; recvS3 = sg3; recvE3 = eg3; inc1 = -1
            else
                recvS1 = eg1+1; recvE1 = ug1; recvS2 = sg2; recvE2 = eg2; recvS3 = sg3; recvE3 = eg3
            endif
        endif
        if (i==3)then
            if (bb==1)then
                recvS1 = lg1; recvE1 = ug1; recvS2 = sg2-1; recvE2 = lg2; recvS3 = sg3; recvE3 = eg3; inc2 = -1
            else
                recvS1 = lg1; recvE1 = ug1; recvS2 = lg2; recvE2 = sg2-1; recvS3 = sg3; recvE3 = eg3
            endif
        endif
        if (i==4)then
            if (bf==1)then
                recvS1 = lg1; recvE1 = ug1; recvS2 = ug2; recvE2 = eg2+1; recvS3 = sg3; recvE3 = eg3; inc2 = -1
            else
                recvS1 = lg1; recvE1 = ug1; recvS2 = eg2+1; recvE2 = ug2; recvS3 = sg3; recvE3 = eg3
            endif
        endif
        if (i==5)then
            if (bs==1)then
                recvS1 = lg1; recvE1 = ug1; recvS2 = lg2; recvE2 = ug2; recvS3 = sg3-1; recvE3 = lg3; inc3 = -1
            else
                recvS1 = lg1; recvE1 = ug1; recvS2 = lg2; recvE2 = ug2; recvS3 = lg3; recvE3 = sg3-1
            endif
        endif
        if (i==6)then
            if (bn==1)then
                recvS1 = lg1; recvE1 = ug1; recvS2 = lg2; recvE2 = ug2; recvS3 = ug3; recvE3 = eg3+1; inc3 = -1
            else
                recvS1 = lg1; recvE1 = ug1; recvS2 = lg2; recvE2 = ug2; recvS3 = eg3+1; recvE3 = ug3
            endif
        endif
    
        cnt = 0
        do k=sendS3,sendE3
        do j=sendS2,sendE2
        do ii=sendS1,sendE1
            idx = cellIdx3Dto1D(ii,j,k)
            cnt = cnt+1
            sendIdxFluid(cnt,i) = idx
        enddo
        enddo
        enddo
        sendCntFluid(i) = cnt

        cnt = 0
        do k=recvS3,recvE3,inc3
        do j=recvS2,recvE2,inc2
        do ii=recvS1,recvE1,inc1
            idx = cellIdx3Dto1D(ii,j,k)
            cnt = cnt + 1
            recvIdxFluid(cnt,i) = idx
        enddo
        enddo
        enddo
        recvCntFluid(i) = cnt
        ! debug mode
        !open(cpuIdx,file=logFile,access="append")
        !write(cpuIdx,*) "side ",i,sendCntFluid(i),recvCntFluid(i)
        !write(cpuIdx,*) sendS1,sendE1,sendS2,sendE2,sendS3,sendE3
        !write(cpuIdx,*) recvS1,recvE1,recvS2,recvE2,recvS3,recvE3
        !write(cpuIdx,*) inc1,inc2,inc3
        !write(cpuIdx,*) "send cell indice:"
        !do ii=1,sendCntFluid(i)
        !    write(cpuIdx,"(4I5)",advance="no") sendIdxFluid(ii,i), &
        !        cellIdx1Dto3D(1,sendIdxFluid(ii,i)),cellIdx1Dto3D(2,sendIdxFluid(ii,i)),cellIdx1Dto3D(3,sendIdxFluid(ii,i))
        !    write(cpuIdx,*) " "
        !enddo
        !write(cpuIdx,*) "Recv cell indice:"
        !do ii=1,recvCntFluid(i)
        !    write(cpuIdx,"(5I5)",advance="no") i,recvIdxFluid(ii,i), &
        !        cellIdx1Dto3D(1,recvIdxFluid(ii,i)),cellIdx1Dto3D(2,recvIdxFluid(ii,i)),cellIdx1Dto3D(3,recvIdxFluid(ii,i))
        !    write(cpuIdx,*) " "
        !enddo
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine FillGhostRegionForFluid()
use vars
implicit none
    integer::i,dim,sideSend,sideRecv,cpuSend,cpuRecv,ierr
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i=1,6
        dim = (i-1)/2+1
        sideSend = i
        if (mod(i,2)==0)then
            sideRecv = i-1
        else
            sideRecv = i+1
        endif
        cpuSend = CPUNeighIdx(sideRecv)
        cpuRecv = CPUNeighIdx(sideSend)
        
        call MPISendAndRecvFluidInfo(sideSend,sideRecv,cpuSend,cpuRecv)
        if (cpuSend==-1)then
            call DealWithBoundaryFluid(dim,sideRecv)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !write(ii,"(I10)") i 
        !call OutputFluidVariables(trim(adjustl(ii)))
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine MPISendAndRecvFluidInfo(sideSend,sideRecv,cpuSend,cpuRecv)
use vars
implicit none
    integer::sideSend,sideRecv,cpuSend,cpuRecv
    integer::request
    integer::info
    integer::nSend,nRecv
    integer::STATUS(MPI_STATUS_SIZE)
    integer::i,offset,idx1,idx2,idx3,cellIdxTmp,nVar
    nSend = sendCntFluid(sideSend)
    nRecv = recvCntFluid(sideRecv)
    nVar = 9
    do i=1,nSend
        offset = nVar*(i-1)
        cellIdxTmp = sendIdxFluid(i,sideSend)
        idx1 = cellIdx1Dto3D(1,cellIdxTmp)
        idx2 = cellIdx1Dto3D(2,cellIdxTmp)
        idx3 = cellIdx1Dto3D(3,cellIdxTmp)
        sendBuffFluid(offset+1) = u(idx1,idx2,idx3)
        sendBuffFluid(offset+2) = v(idx1,idx2,idx3)
        sendBuffFluid(offset+3) = w(idx1,idx2,idx3)
        sendBuffFluid(offset+4) = P(idx1,idx2,idx3)
        sendBuffFluid(offset+5) = omegaZ(idx1,idx2,idx3)
        sendBuffFluid(offset+6) = Pcorr(idx1,idx2,idx3)
        sendBuffFluid(offset+7) = errorU(idx1,idx2,idx3)
        sendBuffFluid(offset+8) = errorV(idx1,idx2,idx3)
        sendBuffFluid(offset+9) = errorP(idx1,idx2,idx3)
    enddo

    call MPI_SENDRECV(sendBuffFluid,nVar*nSend,MPI_DOUBLE,CPUTransferIdx(sideSend),111,&
                      recvBuffFluid,nVar*nRecv,MPI_DOUBLE,CPUTransferIdx(sideRecv),111,MPI_COMM_WORLD,STATUS,info)
    !if (cpuRecv/=-1)then
    !    call MPI_ISEND(sendBuffFluid,nVar*nSend,MPI_DOUBLE,cpuRecv,111,MPI_COMM_WORLD,request,info)
    !    !call MPI_WAIT(request,STATUS,info)
    !endif
    !if (cpuSend/=-1)then
    !    call MPI_IRECV(recvBuffFluid,nVar*nRecv,MPI_DOUBLE,cpuSend,111,MPI_COMM_WORLD,request,info)
    !    call MPI_WAIT(request,STATUS,info)
    !endif
    !call MPI_WAITALL(2,request,STATUS,info)
    if (cpuSend/=-1)then
        do i=1,nRecv
            offset = nVar*(i-1)
            cellIdxTmp = recvIdxFluid(i,sideRecv)
            idx1 = cellIdx1Dto3D(1,cellIdxTmp)
            idx2 = cellIdx1Dto3D(2,cellIdxTmp)
            idx3 = cellIdx1Dto3D(3,cellIdxTmp)
            u(idx1,idx2,idx3) = recvBuffFluid(offset+1)
            v(idx1,idx2,idx3) = recvBuffFluid(offset+2)
            w(idx1,idx2,idx3) = recvBuffFluid(offset+3)
            P(idx1,idx2,idx3) = recvBuffFluid(offset+4)
            omegaZ(idx1,idx2,idx3) = recvBuffFluid(offset+5)
            Pcorr(idx1,idx2,idx3) = recvBuffFluid(offset+6)
            errorU(idx1,idx2,idx3) = recvBuffFluid(offset+7)
            errorV(idx1,idx2,idx3) = recvBuffFluid(offset+8)
            errorP(idx1,idx2,idx3) = recvBuffFluid(offset+9)
        enddo
    endif

return
endsubroutine
!----------------------------------------------------------
subroutine DealWithBoundaryFluid(dim,idx)
use vars
implicit none
    integer dim,idx,i,cellIdxTmp,nRecv
    integer idxS1,idxS2,idxS3,idxR1,idxR2,idxR3
    nRecv = sendCntFluid(idx)
    do i=1,nRecv
        cellIdxTmp = sendIdxFluid(i,idx)
        idxS1 = cellIdx1Dto3D(1,cellIdxTmp)
        idxS2 = cellIdx1Dto3D(2,cellIdxTmp)
        idxS3 = cellIdx1Dto3D(3,cellIdxTmp)
        cellIdxTmp = RecvIdxFluid(i,idx)
        idxR1 = cellIdx1Dto3D(1,cellIdxTmp)
        idxR2 = cellIdx1Dto3D(2,cellIdxTmp)
        idxR3 = cellIdx1Dto3D(3,cellIdxTmp)
         
        u(idxR1,idxR2,idxR3) = u(idxS1,idxS2,idxS3)
        v(idxR1,idxR2,idxR3) = v(idxS1,idxS2,idxS3)
        w(idxR1,idxR2,idxR3) = w(idxS1,idxS2,idxS3)
        P(idxR1,idxR2,idxR3) = P(idxS1,idxS2,idxS3)
        omegaZ(idxR1,idxR2,idxR3) = omegaZ(idxS1,idxS2,idxS3)
        Pcorr(idxR1,idxR2,idxR3) = Pcorr(idxS1,idxS2,idxS3)
        errorU(idxR1,idxR2,idxR3) = errorU(idxS1,idxS2,idxS3)
        errorV(idxR1,idxR2,idxR3) = errorV(idxS1,idxS2,idxS3)
        errorP(idxR1,idxR2,idxR3) = errorP(idxS1,idxS2,idxS3)
        if (bType(idx)==1)then
            if (idx==6)then
                u(idxR1,idxR2,idxR3) = -u(idxR1,idxR2,idxR3)
            else
                u(idxR1,idxR2,idxR3) = -u(idxR1,idxR2,idxR3)
            endif
            v(idxR1,idxR2,idxR3) = -v(idxR1,idxR2,idxR3) 
            w(idxR1,idxR2,idxR3) = -w(idxR1,idxR2,idxR3) 
        endif
    enddo
return 
endsubroutine
!----------------------------------------------------------
