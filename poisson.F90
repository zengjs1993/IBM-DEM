!----------------------------------------------------------
subroutine InitPressureCorrectionMatrix()
use vars
#include<petsc/finclude/petscksp.h>
use petscksp
implicit none
    integer::nGrid,ierr
    nGrid = nTotX*nTotY*nTotZ
    pConstIdx = 1 !nGrid/2
    call MatCreate(PETSC_COMM_WORLD,pMat,ierr)
    call MatSetSizes(pMat,PETSC_DECIDE,PETSC_DECIDE,nGrid,nGrid,ierr)
    call MatSetType(pMat,MATMPIAIJ,ierr)
    call MatMPIAIJSetPreallocation(pMat,7,PETSC_NULL_INTEGER,7,PETSC_NULL_INTEGER,ierr)
    call MatSetFromOptions(pMat,ierr)
    call MatSetup(pMat,ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,nGrid,pCorrVec,ierr)
    call VecSetFromOptions(pCorrVec,ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,nGrid,rhsP,ierr)
    call VecSetFromOptions(rhsP,ierr)

    call SetGlobalToLocalIdx()
    call VecCreateMPI(PETSC_COMM_WORLD,nx*ny*nz,nGrid,pVecLocal,ierr)
    call ISCreateGeneral(PETSC_COMM_WORLD,nx*ny*nz,fromGlobalIdx,PETSC_COPY_VALUES,fromIS,ierr)
    call ISCreateGeneral(PETSC_COMM_WORLD,nx*ny*nz,toLocalIdx,PETSC_COPY_VALUES,toIS,ierr)
    call VecScatterCreate(pCorrVec,fromIS,pVecLocal,toIS,scatter,ierr)
return
end subroutine
!----------------------------------------------------------
subroutine SetGlobalToLocalIdx()
use vars
implicit none
    integer::i,j,k,cnt
    cnt = 0
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        fromGlobalIdx(cnt) = (k-sg3)*nx*ny+(j-sg2)*nx+(i-sg1)+GlobalStartCellIndex1D(cpuIdx+1)
        toLocalIdx(cnt) = (k-sg3)*nx*ny+(j-sg2)*nx+(i-sg1)+GlobalStartCellIndex1D(cpuIdx+1)
        cnt = cnt+1
    enddo
    enddo
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine SetPressureCorrectionJacobian()
use vars
#include<petsc/finclude/petscksp.h>
use petsc
use petscksp
implicit none
    PetscInt::i,j,k,l,idxC,dim,idxNeigh
    PetscErrorCode::ierr
    PetscInt::ione
    PetscScalar::value,vc,dl(1:3)


    !call MatGetOwnershipRange(pMat,IStart,IEnd,ierr)
    !open(cpuIdx,file=logFile,access="append")
    !write(cpuIdx,*) IStart,IEnd
    !close(cpuIdx)

    !ione = 1
    !do II=IStart,Iend-1
    !    i =  
    !    do l=1,6
    !        
    !    enddo
    !enddo


    ione = 1
    dl(1) = dx; dl(2) = dy; dl(3) = dz
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        idxC = (k-sg3)*nx*ny+(j-sg2)*nx+(i-sg1)+GlobalStartCellIndex1D(cpuIdx+1)
        if (idxC/=pConstIdx-1)then
            vc = 0d0
            do l=1,6
                dim = (l+1)/2
                idxNeigh = GlobalNeighCellIndex(l,i,j,k)
                if (idxNeigh/=-1)then
                    value = -1d0/(dl(dim)**2)
                    call MatSetValues(pMat,ione,idxC,ione,idxNeigh-1,value,ADD_VALUES,ierr)
                    vc = vc + value
                endif
            enddo
            vc = -vc
            call MatSetValues(pMat,ione,idxC,ione,idxC,vc,ADD_VALUES,ierr)
        endif
    enddo
    enddo
    enddo

    ! pressure const grid
    if(cpuIdx==0)then
        value = 1d0
        call MatSetValues(pMat,ione,pConstIdx-1,ione,pConstIdx-1,value,ADD_VALUES,ierr)
    endif
    !if (cpuIdx==0) write(*,*)"Assembling Matrix"
    call MatAssemblyBegin(pMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(pMat,MAT_FINAL_ASSEMBLY,ierr)
return
endsubroutine
!------------------------------------------------------------
subroutine SetPressureVector()
#include<petsc/finclude/petscvec.h>
use petscvec
use vars
implicit none
    PetscErrorCode::ierr
    PetscScalar::value
    PetscInt::i,j,k,idx
    !one = 1d0
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        idx = (k-sg3)*nx*ny+(j-sg2)*nx+(i-sg1)+GlobalStartCellIndex1D(cpuIdx+1)
        value = velDiv(i,j,k)
        call VecSetValues(rhsP,1,idx,value,INSERT_VALUES,ierr)
    enddo
    enddo
    enddo
    if (cpuIdx==0)then
        call VecSetValues(rhsP,1,pConstIdx-1,0d0,INSERT_VALUES,ierr)
    endif
    call VecAssemblyBegin(rhsP,ierr)
    call VecAssemblyEnd(rhsP,ierr)
return
endsubroutine
!-----------------------------------------------------------
subroutine SetKSP()
#include<petsc/finclude/petscksp.h>
use petscksp
use vars
implicit none
    integer::ierr
    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
    call KSPSetOperators(ksp,pMat,pMat,ierr)
    call KSPSetFromOptions(ksp,ierr)
return 
endsubroutine
!-----------------------------------------------------------
subroutine PoissonSolver()
#include<petsc/finclude/petscksp.h>
use petscksp
use vars
implicit none
    integer::ierr
    call KSPSolve(ksp,rhsP,pCorrVec,ierr)
    !if (cpuIdx==0)write(*,*)"Linear system solved"
    call GetLocalArray()
    !call VecView(pCorrVec,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call MatMult(pMat,pCorrVec,rhsP,ierr)
    !call VecView(rhsP,PETSC_VIEWER_STDOUT_WORLD,ierr)
return
endsubroutine
!------------------------------------------------------------
subroutine GetLocalArray()
use vars
implicit none
    integer::ierr,i,j,k
    call VecScatterBegin(scatter,pCorrVec,pVecLocal,INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(scatter,pCorrVec,pVecLocal,INSERT_VALUES,SCATTER_FORWARD,ierr)
    !call VecView(pCorrVec,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !call VecView(pVecLocal,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call VecGetArrayF90(pVecLocal,pCorr1D,ierr)
    !open(cpuIdx,file=logFile,access="append")
    !write(cpuIdx,*) "Start output..."
    !close(cpuIdx)
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        Pcorr(i,j,k) = pCorr1D((k-sg3)*nx*ny+(j-sg2)*nx+(i-sg1)+1)
        !open(cpuIdx,file=logFile,access="append")
        !write(cpuIdx,*) Pcorr(i,j,k)
        !close(cpuIdx)
    enddo
    enddo
    enddo
    !open(cpuIdx,file=logFile,access="append")
    !write(cpuIdx,*) "End output..."
    !close(cpuIdx)
return
endsubroutine
!------------------------------------------------------------
