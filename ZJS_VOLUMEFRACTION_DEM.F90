!-----------------------------------------------
subroutine UpdateFluidVelocityOnLagrangianPoints()
use vars
use ZjsVar
implicit none
    integer::i,j,cellIdx,idx1,idx2,idx3,ierr
    real*8::ulp,vlp,wlp,weight(1:3),relaxfactor
    relaxfactor = 1.0d0
    do i=1,npLocal
        !part(i)%fdrag(:) = 0d0
        ulp = 0d0; vlp = 0d0; wlp = 0d0
        do j=1,neighCellNo(i)
            cellIdx = neighCellIdx(j,i)
            weight(1) = neighCellWeight(j,i,1)
            weight(2) = neighCellWeight(j,i,2)
            weight(3) = neighCellWeight(j,i,3)
            idx1 = map1Dto3D(1,cellIdx)
            idx2 = map1Dto3D(2,cellIdx)
            idx3 = map1Dto3D(3,cellIdx)
            ulp = ulp + weight(1)*u(idx1,idx2,idx3)*dx*dy*dz
            vlp = vlp + weight(2)*v(idx1,idx2,idx3)*dx*dy*dz
            wlp = wlp + weight(3)*w(idx1,idx2,idx3)*dx*dy*dz
        enddo
        part(i)%vgap(1) = ulp-part(i)%v(1)
        part(i)%vgap(2) = vlp-part(i)%v(2)
        part(i)%vgap(3) = wlp-part(i)%v(3)
        part(i)%fdrag(1) = part(i)%fdrag(1) + relaxfactor*(ulp-part(i)%v(1))/dt
        part(i)%fdrag(2) = part(i)%fdrag(2) + relaxfactor*(vlp-part(i)%v(2))/dt
        part(i)%fdrag(3) = part(i)%fdrag(3) + relaxfactor*(wlp-part(i)%v(3))/dt
        !open(cpuIdx,file=logFile,access="append")
        !write(cpuIdx,*)"part ",i, part(i)%fdrag,ulp,vlp,wlp
        !close(cpuIdx)
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
return
endsubroutine
!-----------------------------------------------  
subroutine UpdateMomentumSourceTerm()
use vars
use ZjsVar
implicit none
    integer::i,j,cellIdx,idx1,idx2,idx3,ierr
    real*8::weight(1:3)
    fxOld(:,:,:) = fx(:,:,:); fyOld(:,:,:) = fy(:,:,:); fzOld(:,:,:) = fz(:,:,:);
    fx(:,:,:) = 0d0
    fy(:,:,:) = 0d0
    fz(:,:,:) = 0d0
    do i=1,np
    do j=1,neighCellNo(i)
        cellIdx = neighCellIdx(j,i)
        if (cellIdx/=-1)then
            weight(1) = neighCellWeight(j,i,1)
            weight(2) = neighCellWeight(j,i,2)
            weight(3) = neighCellWeight(j,i,3)
            idx1 = map1Dto3D(1,cellIdx)
            idx2 = map1Dto3D(2,cellIdx)
            idx3 = map1Dto3D(3,cellIdx)
            
            fx(idx1,idx2,idx3) = fx(idx1,idx2,idx3) - part(i)%fdrag(1)*part(i)%vol*weight(1)
            fy(idx1,idx2,idx3) = fy(idx1,idx2,idx3) - part(i)%fdrag(2)*part(i)%vol*weight(2)
            fz(idx1,idx2,idx3) = fz(idx1,idx2,idx3) - part(i)%fdrag(3)*part(i)%vol*weight(3)
        endif
    enddo
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
return
endsubroutine
!-----------------------------------------------
