program main
#include<petsc/finclude/petscksp.h>
use petscksp
use vars
use ZjsVar
implicit none
    integer ierr
    integer i,nn
    integer flag
    character(len=20) stepStr
    real*8::pi,lenX,lenY,lenZ

    pi = acos(-1.0d0)
    flag = 0
    ! Initialization
    call MPI_INIT(ierr)
    nn = 50*1
    lenX = 5d-3*1 !1.67d-3*2
    lenY = lenX; lenZ = lenX
    call PartitionInitialization(4,1,4,240,32,240,9d-3,1.2d-3,9d-3)
    call GridInitialization()
    if (cpuIdx==0.and.flag==1) write(*,*) "Allocation completed"

    call SetCellIndexMap()
    call BuildMPIRecvAndSendIndexForFluid()
    if (cpuIdx==0.and.flag==1) write(*,*) "Recv and Send indice build-up completed"

    call OutputGridAndPartitionInformation()
    if (cpuIdx==0.and.flag==1) write(*,*) "Output grid and partition infos completed"
    
    call InitializationDEM()

    call FluidInitialization()
    if (cpuIdx==0.and.flag==1) write(*,*) "Fluid initialization completed"

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    if (cpuIdx==0.and.flag==1) write(*,*) "PETSC initialization completed"

    call InitPressureCorrectionMatrix()
    if (cpuIdx==0.and.flag==1) write(*,*) "Matrix initialization completed"

    call SetPressureCorrectionJacobian()
    call SetKSP()
    if (cpuIdx==0.and.flag==1) write(*,*) "Poisson Jacobian build-up completed"
    
    call OutputFluidTecplot(0)
    do i=1,20000
    time = time+dt
    step = step+1
    !call ThirdOrderLowStorageRungeKutta()
    call FirstOrderForwardEuler()
    call PostProcessing()
    
    !call CloudEvolutionMS()
    !call PostProcessing()
    if (mod(i,100)==0.or.i==1) then
        write(stepStr,*) i
        call OutputFluidTecplot(i)
        call OutputCloudTecplot(adjustl(trim(stepStr)))
        !if (cpuIdx==0) then
        !    write(*,*) "step",i
        !    write(*,*) MS(1)%v
        !    write(*,*) MS(1)%omega
        !    !write(*,*) MS(1)%rotMatGLo
        !endif
    endif

    enddo


    call PetscFinalize(ierr)
    call MPI_FINALIZE(ierr)
end program main
