!-----------------------------------------------
subroutine CompressDEM()
use vars
use ZjsVar
implicit none
    integer::i,np_updated
    integer::idxx,idxy,idxz
    real*8,dimension(1:3)::xi
    np_updated = 0
    do i=1,np
        if (part(i)%pType==-1)then
        xi = part(i)%x
        idxx = floor((xi(1)-xStart)/di)+SG1
        idxy = floor((xi(2)-yStart)/dj)+SG2
        idxz = floor((xi(3)-zStart)/dk)+SG3
        if(idxx>=SG1.and.idxx<=EG1.and.idxy>=SG2.and.idxy<=EG2.and.idxz>=SG3.and.idxz<=EG3)then
        !if (xi(1)>=xStart.and.xi(1)<xEnd.and.xi(2)>=yStart.and.xi(2)<yEnd.and.xi(3)>=zStart.and.xi(3)<zEnd)then
            np_updated = np_updated + 1
            part(np_updated) = part(i)
            part(np_updated)%pType = -1
        endif
        endif
    enddo
    np = np_updated
    npLocal = np
return
endsubroutine
!-----------------------------------------------
subroutine PutParcelsIntoCellsDEM(npStart,npEnd,mode)
use vars
use ZjsVar
implicit none
    integer::mode
    integer::i,npStart,npEnd
    integer::idxx,idxy,idxz,cellIdx1D
    integer::nSendTmp
    real*8,dimension(1:3)::xi

    do i=npStart,npEnd
        xi = part(i)%x
        idxx = floor((xi(1)-xStart)/di)+SG1
        idxy = floor((xi(2)-yStart)/dj)+SG2
        idxz = floor((xi(3)-zStart)/dk)+SG3
        !part(i)%globalIdx = cpuLocal
        cellIdx1D = -1
        if(idxx<LG1.and.idxx>UG1.and.idxy<LG2.and.idxy>UG2.and.idxz<LG3.and.idxz>UG3)then
            write(*,*) "Parcel exceeds the region!!!"
        else
            cellIdx1D = map3Dto1D(idxx,idxy,idxz)
            npPartInCell(cellIdx1D) = npPartInCell(cellIdx1D) + 1
            partIdxInCell(npPartInCell(cellIdx1D),cellIdx1D) = i 
            part(i)%localIdx = i
            part(i)%cellIdx = cellIdx1D
            if (mode==0)then
                part(i)%pType = -1
            endif
        endif
        ! LG~SG-1,SG~SG+NGG-1,EG-NGG+1~EG,EG+1~EG+UG
        if (mode<=0)then
        if (idxx<=SG1+haloX-1.and.idxx>=SG1)then
            npSend(1) = npSend(1)+1
            partSend(npSend(1),1) = part(i)
            partSendIdx(npSend(1),1) = i
        endif
        if (idxx>=EG1-haloX+1.and.idxx<=EG1)then
            npSend(2) = npSend(2)+1
            partSend(npSend(2),2) = part(i)
            partSendIdx(npSend(2),2) = i
        endif
        endif
        if (mode<=1)then
        if (idxy<=SG2+haloY-1.and.idxy>=SG2)then
            npSend(3) = npSend(3)+1
            partSend(npSend(3),3) = part(i)
            partSendIdx(npSend(3),3) = i
        endif
        if (idxy>=EG2-haloY+1.and.idxy<=EG2)then
            npSend(4) = npSend(4)+1
            partSend(npSend(4),4) = part(i)
            partSendIdx(npSend(4),4) = i
        endif
        endif
        if (mode<=2)then
        if (idxz<=SG3+haloZ-1.and.idxz>=SG3)then
            npSend(5) = npSend(5)+1
            partSend(npSend(5),5) = part(i) 
            partSendIdx(npSend(5),5) = i
        endif
        if (idxz>=EG3-haloZ+1.and.idxz<=EG3)then
            npSend(6) = npSend(6)+1
            partSend(npSend(6),6) = part(i)
            partSendIdx(npSend(6),6) = i
        endif
        endif
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine RebuildParcelInCellListDEM()
use vars
use ZjsVar
implicit none
    real*8,dimension(1:3)::xi
    integer::i,cellIdx,idxx,idxy,idxz,ierr
    ! LG~SG-1,SG~SG+NGG-1,EG-NGG+1~EG,EG+1~EG+UG
    call CompressDEM()
    !if (step>=4853)then
    !    call OutputCloudDEM("Compressed")
    !endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    npSend(1:6) = 0
    npRecv(1:6) = 0
    npPartInCell(1:ni*nj*nk) = 0
    ! Inner region
    call PutParcelsIntoCellsDEM(1,np,0)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    ! Ghost region
    call RebuildParcelInCellListForGhostRegionDEM()
return
endsubroutine
!-----------------------------------------------
subroutine RebuildParcelInCellListForGhostRegionDEM()
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
        ! send parcels
        cpuSend = cpuLocal
        nSend = npSend(idx)
        if (idx==1.or.idx==2)then
            dim = 1
        else if (idx==3.or.idx==4)then
            dim = 2
        else
            dim = 3
        endif
        if (bType(idx)==0) then ! period
            do i=1,nSend
                partSend(i,idxSend)%x(dim) = partSend(i,idxSend)%x(dim) + coordShift(idxSend)
            enddo
        endif
        !    cpuRecv = neighCpu(idx)
        !elseif (bType(idx)==1) then ! solid
        !    do i=1,nSend
        !        partSend(i,idxSend)%x(dim) = 2*coordMirror(idxSend)-partSend(i,idxSend)%x(dim)
        !        partSend(i,idxSend)%v(dim) = -partSend(i,idxSend)%v(dim)
        !    enddo
        !    cpuRecv = cpuLocal
        !else ! inner
        !    cpuRecv = neighCpu(idx)
        !endif
        cpuSend = neighCpu(idxRecv) ! date: 20181211
        cpuRecv = neighCpu(idxSend) ! date: 20181211
        call MyMPISendAndRecv(cpuSend,cpuRecv,idxSend,idxRecv)
        if (cpuSend==-1)then
            call UpdateFullInfoForSolidOutletBoundary(dim,idxRecv)
        endif

        ! put new parcels into cells
        nRecv = npRecv(idxRecv)
        call PutParcelsIntoCellsDEM(np+1,np+nRecv,dim)
        startIdxRecv(idxRecv) = np+1
        np = np+nRecv
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine UpdateParcelInfoInGhostRegion()
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
        cpuSend = cpuLocal
        nSend = npSend(idx)
        if (idx==1.or.idx==2)then
            dim = 1
        else if (idx==3.or.idx==4)then
            dim = 2
        else
            dim = 3
        endif
        do i=1,nSend
            partSend(i,idx)%MSIdx = part(partSendIdx(i,idx))%MSIdx
            partSend(i,idx)%x(:)  = part(partSendIdx(i,idx))%x(:)
            partSend(i,idx)%v(:)  = part(partSendIdx(i,idx))%v(:)
            partSend(i,idx)%wsum  = part(partSendIdx(i,idx))%wsum
            partSend(i,idx)%fdrag  = part(partSendIdx(i,idx))%fdrag
        enddo
        if (bType(idx)==0) then ! period
            do i=1,nSend
                partSend(i,idxSend)%x(dim) = partSend(i,idxSend)%x(dim) + coordShift(idxSend)
            enddo
        endif
        !        partSend(i,idxSend)%pType = 0
        !    enddo
        !    cpuRecv = neighCpu(idx)
        !elseif (bType(idx)==1) then ! solid
        !    do i=1,nSend
        !        partSend(i,idxSend)%x(dim) = 2*coordMirror(idxSend)-partSend(i,idxSend)%x(dim)
        !        partSend(i,idxSend)%v(dim) = -partSend(i,idxSend)%v(dim)
        !        partSend(i,idxSend)%pType = 1
        !    enddo
        !    cpuRecv = cpuLocal
        !else ! inner
        !    cpuRecv = neighCpu(idx)
        !    do i=1,nSend
        !        partSend(1,idxSend)%pType = -1
        !    enddo
        !endif
        cpuSend = neighCpu(idxRecv)
        cpuRecv = neighCpu(idxSend)
        call MyMPIUpdateVelocityAndPosition(cpuSend,cpuRecv,idxSend,idxRecv)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (cpuSend==-1)then
            call UpdatePartInfoForSolidOutletBoundary(dim,idxRecv)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ! put new parcels into cells
        !call PutParcelsIntoCellsDEM(np+1,np+nSend)
        !np = np+nSend
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine RebuildParcelCollisionListDEM()
use vars
use ZjsVar
implicit none
    integer::i,j,k,cellIdxTmp,npCellTmp,ii,jj,flag

    open(cpuLocal,file=logFile,access="APPEND")
    !write(cpuLocal,*) "step part ms idx:",step
    !do i=1,np
    !    write(cpuLocal,*) part(i)%MSIdx
    !enddo
    !write(cpuLocal,*) "Time_cnt: ",timeCntDEM
    collisionBufferRecord = collisionBuffer
    do i=1,npLocal
        neighPartNo(i) = 0
        ii = i
        do j=1,maxNeighCellNo
            cellIdxTmp = neighCellIdxForCell(j,part(i)%cellIdx)
            if (cellIdxTmp/=-1)then
            npCellTmp = npPartInCell(cellIdxTmp)
            do k=1,npCelltmp
                jj = partIdxInCell(k,cellIdxTmp)
                if (part(ii)%MSIdx/=part(jj)%MSIdx)then
                    call Collision(ii,jj,flag)
                    if (flag>0.and.part(jj)%pType==-1.and.ii/=jj)then
                        !write(cpuLocal,*) jj,part(ii)%x,part(jj)%x,part(ii)%MSIdx,part(jj)%MSIdx

                        neighPartNo(ii) = neighPartNo(ii) + 1
                        neighPartIdx(neighPartNo(ii),ii) = jj
                    endif
                endif
            enddo
            endif
        enddo
    enddo
    
    close(cpuLocal)
return
endsubroutine
!-----------------------------------------------
subroutine RebuildNeighborCellListForParcelsDEM()
use vars
use ZjsVar
implicit none
    integer::i,j,k,ii,jj,cellIdxTmp,idx1,idx2,idx3,flag
    real*8::xp1,xp2,xp3,xc1,xc2,xc3,weight
!    open(cpuLocal,file=logFile,access="APPEND")
    do i=1,np
        xp1 = part(i)%x(1)
        xp2 = part(i)%x(2)
        xp3 = part(i)%x(3)
        neighCellNo(i) = 0
        do j=1,maxNeighCellNo
            cellIdxTmp = neighCellIdxForCell(j,part(i)%cellIdx)
            if (cellIdxTmp/=-1)then
                neighCellNo(i) = neighCellNo(i) + 1
                neighCellIdx(neighCellNo(i),i) = cellIdxTmp
            endif
        enddo
    enddo
    !write(cpuLocal,*) "NEIGHBOR UPDATED..."
    !do i=1,np
    !write(cpuLocal,*) part(i)%cellIdx
    !do j=1,27
    !    write(cpuLocal,"(I5)",advance="no") neighCellIdx(j,i)
    !enddo
    !write(cpuLocal,*) " "
    !do j=1,27
    !    write(cpuLocal,"(I5)",advance="no") neighCellIdxForCell(j,part(i)%cellIdx)
    !enddo
    !write(cpuLocal,*) " "
    !enddo
return
endsubroutine
!-----------------------------------------------
subroutine UpdateLocalCellWeightListDEM()
use vars
use ZjsVar
implicit none
    integer::i,j,cellIdxTmp,cellNoTmp,idx1,idx2,idx3,flag
    real*8::xp1,xp2,xp3,xc1,xc2,xc3,weight(1:3)
!    open(cpuLocal,file=logFile,access="append") 
    do i=1,npLocal
        xp1 = part(i)%x(1)
        xp2 = part(i)%x(2)
        xp3 = part(i)%x(3)
        cellNoTmp = neighCellNo(i)
        do j=1,cellNoTmp
            cellIdxTmp = neighCellIdx(j,i)
            idx1 = map1Dto3D(1,cellIdxTmp)
            idx2 = map1Dto3D(2,cellIdxTmp)
            idx3 = map1Dto3D(3,cellIdxTmp)
            xc1 = xc(idx1,idx2,idx3)
            xc2 = yc(idx1,idx2,idx3)
            xc3 = zc(idx1,idx2,idx3)
            call GetWeight(xp1,xp2,xp3,xc1,xc2,xc3,weight)
            neighCellWeight(j,i,1) = weight(1)
            neighCellWeight(j,i,2) = weight(2)
            neighCellWeight(j,i,3) = weight(3)
        enddo
    enddo

    call GetWeightSumForParcels()
    !call LocalWeightNormalization()

    !if (step==0.and.cpuLocal==0)then
    !    open(cpuLocal,file="./DEM/log/1.log",access="append")
    !    do i=1,npLocal
    !    cellNoTmp = neighCellNo(i)
    !    write(cpuLocal,*) "particle ",part(i)%globalIdx,"neigh cell num",cellNoTmp,"Coords:",part(i)%x(1),part(i)%x(2),part(i)%x(3),part(i)%wsum
    !    do j=1,cellNoTmp
    !        cellIdxTmp = neighCellIdx(j,i)
    !        weight = neighCellWeight(j,i,:)
    !        idx1 = map1Dto3D(1,cellIdxTmp)
    !        idx2 = map1Dto3D(2,cellIdxTmp)
    !        idx3 = map1Dto3D(3,cellIdxTmp)
    !        write(cpuLocal,*) idx1,idx2,idx3,weight(1),weight(2),weight(3)
    !    enddo
    !    enddo
    !    close(cpuLocal)
    !endif
return
endsubroutine
!-----------------------------------------------
subroutine UpdateGhostCellWeightListDEM()
use vars
use ZjsVar
implicit none
    integer::i,j,cellIdxTmp,cellNoTmp,idx1,idx2,idx3,flag
    real*8::xp1,xp2,xp3,xc1,xc2,xc3,weight(1:3),wsumI
    do i=npLocal+1,np
        xp1 = part(i)%x(1)
        xp2 = part(i)%x(2)
        xp3 = part(i)%x(3)
        cellNoTmp = neighCellNo(i)
        do j=1,cellNoTmp
            cellIdxTmp = neighCellIdx(j,i)
            idx1 = map1Dto3D(1,cellIdxTmp)
            idx2 = map1Dto3D(2,cellIdxTmp)
            idx3 = map1Dto3D(3,cellIdxTmp)
            xc1 = xc(idx1,idx2,idx3)
            xc2 = yc(idx1,idx2,idx3)
            xc3 = zc(idx1,idx2,idx3)
            call GetWeight(xp1,xp2,xp3,xc1,xc2,xc3,weight)
            neighCellWeight(j,i,1) = weight(1)
            neighCellWeight(j,i,2) = weight(2)
            neighCellWeight(j,i,3) = weight(3)

        enddo
    enddo
    !if (step==0.and.cpuLocal==0)then
    !    open(cpuLocal,file="./DEM/log/1.log",access="append")
    !    write(cpuLocal,*) "Here comes the ghost particles"
    !    do i=npLocal+1,np
    !    cellNoTmp = neighCellNo(i)
    !    write(cpuLocal,*) "particle ",part(i)%globalIdx,"neigh cell num",cellNoTmp,"Coords:",part(i)%x(1),part(i)%x(2),part(i)%x(3),part(i)%wsum
    !    do j=1,cellNoTmp
    !        cellIdxTmp = neighCellIdx(j,i)
    !        weight = neighCellWeight(j,i,:)
    !        idx1 = map1Dto3D(1,cellIdxTmp)
    !        idx2 = map1Dto3D(2,cellIdxTmp)
    !        idx3 = map1Dto3D(3,cellIdxTmp)
    !        write(cpuLocal,*) idx1,idx2,idx3,weight(1),weight(2),weight(3)
    !    enddo
    !    enddo
    !    close(cpuLocal)
    !endif
return
endsubroutine
!-----------------------------------------------
