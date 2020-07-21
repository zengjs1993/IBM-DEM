!-----------------------------------------------
subroutine InitializationDEM()
use vars
use ZjsVar
implicit none
    integer::req,status(MPI_STATUS_SIZE),info
    integer::cpuSend,cpuRecv,msgRecv

    call InitParameters()
    call CopyGridInformationDEM()
    call AllocateDEM()
    call SetCellType()
    call SetCellIndexMapDEM()
    call SetNeighCellIndexForCellsDEM()
    call SetBoundaryNormsAndPoints()
    call SetMirrorAndShiftCoords()
    !call InitParcelsDEM()
    call InitCylinderParcelsDEM()
    call SetWallProbe()

    call WriteDEM()
    call OutputCloudTecplot("0")
    
    call RebuildListMS()
    call RebuildParcelInCellListDEM()
    call OutputCloudTecplot("1") 
    call LogDEM(15,"Rebuild parcel collision list")
    call RebuildParcelCollisionListDEM()
    call LogDEM(17,"Update local cell weight list")
    call RebuildNeighborCellListForParcelsDEM()
    call UpdateLocalCellWeightListDEM()
    call LogDEM(19,"Update parcel info in ghost region") 
    call UpdateParcelInfoInGhostRegion() 
    call LogDEM(21,"Update ghost cell weight list") 
    call UpdateGhostCellWeightListDEM()
    !call AssembleForceMS()
    call WriteDEM()
return
endsubroutine
!-----------------------------------------------
subroutine SetCellType()
use vars
use ZjsVar
    cellType(:,:,:) = 0
    if (bType(1)==1)then
        cellType(LG1:SG1-1,:,:) = 1
    endif
    if (bType(2)==1)then
        cellType(EG1+1:UG1,:,:) = 1
    endif
    if (bType(3)==1)then
        cellType(:,LG2:SG2-1,:) = 1
    endif
    if (bType(4)==1)then
        cellType(:,EG2+1:UG2,:) = 1
    endif
    if (bType(5)==1)then
        cellType(:,:,LG3:SG3-1) = 1
    endif
    if (bType(6)==1)then
        cellType(:,:,EG3+1:UG3) = 1
    endif
    cellType(SG1:EG1,SG2:EG2,SG3:EG3) = -1
return
endsubroutine
!-----------------------------------------------
subroutine InitParameters()
use vars
use ZjsVar
    alphapInit = 0.65d0
    timeCntDEM = 0
    time = 0d0
    knCoeff = 1d4*4
    ResistCoeff = 0.0001d0
    coupleInterval = 10
    dtDEM = dt/coupleInterval
return
endsubroutine
!-----------------------------------------------
subroutine AllocateDEM()
use vars
use ZjsVar
implicit none

    integer::nGrid
    
    nGrid = ni*nj*nk
    
    Allocate(neighPartNo(1:maxPartNo))
    Allocate(neighPartIdx(1:maxCollisionNo,1:maxPartNo))
    Allocate(neighCellNo(1:maxPartNo))
    Allocate(neighCellIdx(1:maxNeighCellNo,1:maxPartNo))
    Allocate(neighCellWeight(1:maxNeighCellNo,1:maxPartNo,1:3))
    
    Allocate(npPartInCell(1:nGrid))
    Allocate(partIdxInCell(1:maxPartNoInCell,1:nGrid))
    
    Allocate(neighCellIdxForCell(1:maxNeighCellNo,1:nGrid))
    Allocate(map1Dto3D(1:3,1:nGrid))
    Allocate(map3Dto1D(LG1:UG1,LG2:UG2,LG3:UG3))
    Allocate(cellType(LG1:UG1,LG2:UG2,LG3:UG3))
    
    Allocate( srcMomX(LG1:UG1,LG2:UG2,LG3:UG3))
    Allocate( srcMomY(LG1:UG1,LG2:UG2,LG3:UG3))
    Allocate( srcMomZ(LG1:UG1,LG2:UG2,LG3:UG3))
    Allocate(  srcEng(LG1:UG1,LG2:UG2,LG3:UG3))
    Allocate(aveVelpX(LG1:UG1,LG2:UG2,LG3:UG3))
    Allocate(aveVelpY(LG1:UG1,LG2:UG2,LG3:UG3))
    Allocate(aveVelpZ(LG1:UG1,LG2:UG2,LG3:UG3))
    
    write(logFile,*) cpuLocal,"_dem.log"
    logFile = "./DEM/log/"//trim(adjustl(logFile))
return
endsubroutine
!-----------------------------------------------
subroutine SetBoundaryNormsAndPoints()
use vars
use ZjsVar
implicit none
    pointBound(1,1) = 0d0
    pointBound(2,1) = 0d0
    pointBound(3,1) = 0d0
    pointBound(1:3,3) = pointBound(1:3,1)
    pointBound(1:3,5) = pointBound(1:3,1)
    pointBound(1,2) = lxGlo
    pointBound(2,2) = lyGlo
    pointBound(3,2) = lzGlo
    pointBound(1:3,4) = pointBound(1:3,2)
    pointBound(1:3,6) = pointBound(1:3,2)
    
    normBound(1:3,1:6) = 0d0
    normBound(1,1) = 1d0
    normBound(1,2) = -1d0
    normBound(2,3) = 1d0
    normBound(2,4) = -1d0
    normBound(3,5) = 1d0
    normBound(3,6) = -1d0
return 
endsubroutine
!-----------------------------------------------
subroutine SetMirrorAndShiftCoords()
use vars
use ZjsVar
implicit none
    coordMirror(1) = 0    
    coordMirror(2) = lxGlo
    coordMirror(3) = 0    
    coordMirror(4) = lyGlo
    coordMirror(5) = 0    
    coordMirror(6) = lzGlo
    coordShift(1)  = lxGlo
    coordShift(2)  = -lxGlo
    coordShift(3)  = lyGlo 
    coordShift(4)  = -lyGlo
    coordShift(5)  = lzGlo 
    coordShift(6)  = -lzGlo
return
endsubroutine
!-----------------------------------------------
subroutine SetCellIndexMapDEM()
use vars
use ZjsVar
implicit none
   integer::i,j,k,idx
   !open(cpuLocal,file="./DEM/log/indexmap.log")
   do k=LG3,UG3
   do j=LG2,UG2
   do i=LG1,UG1
    idx = (k-LG3)*ni*nj+(j-LG2)*ni+i-LG1+1
    map3Dto1D(i,j,k) = idx
    map1Dto3D(1,idx) = i
    map1Dto3D(2,idx) = j
    map1Dto3D(3,idx) = k
    !write(cpuLocal,*) map3Dto1D(i,j,k), i,j,k,map1Dto3D(1,idx),map1Dto3D(2,idx),map1Dto3D(3,idx)
   enddo
   enddo
   enddo
return
endsubroutine
!-----------------------------------------------
subroutine SetNeighCellIndexForCellsDEM()
use vars
use ZjsVar
implicit none
    integer::i,j,k,ii,jj,kk,cellIdxLocal,idx1,idx2,idx3,neighIdx
    !open(cpuLocal,file=logFile,access="APPEND")
    !write(cpuLocal,*) LG1,UG1,LG2,UG2,LG3,UG3
    do k=LG3,UG3
    do j=LG2,UG2
    do i=LG1,UG1
        cellIdxLocal = map3Dto1D(i,j,k)
        !write(cpuLocal,*) cellIdxLocal
        do kk=-2,2
        do jj=-2,2
        do ii=-2,2
            idx1 = i+ii; idx2 = j+jj; idx3 = k+kk
            neighIdx = (kk+2)*25+(jj+2)*5+ii+2 + 1
            
            if (idx1<LG1.or.idx1>UG1.or.idx2<LG2.or.idx2>UG2.or.idx3<LG3.or.idx3>UG3)then
                neighCellIdxForCell(neighIdx,cellIdxLocal) = -1
            else
                neighCellIdxForCell(neighIdx,cellIdxLocal) = (map3Dto1D(idx1,idx2,idx3)) 
            endif
            !write(cpuLocal,*) idx1,idx2,idx3,cellIdxLocal,neighIdx,neighCellIdxForCell(neighIdx,cellIdxLocal)
        enddo
        enddo
        enddo
    enddo
    enddo
    enddo

    !close(cpuLocal)
return
endsubroutine
!-----------------------------------------------
subroutine CopyGridInformationDEM()
use vars
use ZjsVar
implicit none
    ni = UG1-LG1+1; nj = UG2-LG2+1; nk = UG3-LG3+1
    di = dx; dj = dy; dk = dz
    xStart = (SG1-1)*di; yStart = (SG2-1)*dj; zStart = (SG3-1)*dk
    XEnd = (EG1-1+1)*di; yEnd = (EG2-1+1)*dj; zEnd = (EG3-1+1)*dk
    cpuLocal = CPUIdx
    neighCpu(1) = CPUNeighIdx(1); neighCpu(2) = CPUNeighIdx(2)
    neighCpu(3) = CPUNeighIdx(3); neighCpu(4) = CPUNeighIdx(4)
    neighCpu(5) = CPUNeighIdx(5); neighCpu(6) = CPUNeighIdx(6)
    lxGlo = di*nTotX; lyGlo = dj*nTotY; lzGlo = dk*nTotZ
   
return
endsubroutine
!-----------------------------------------------
subroutine InitParcelsDEM()
use vars
use ZjsVar
implicit none
    character(len=100)::filename
    integer::i,j,k,idx
    real*8::xpc(1:3,1:64),x1,x2,x3,xc1,xc2,xc3
    real*8::pi,rhoMS,deltaL
    integer::nTplt
    
    !nMS = 0 
    !nMSLocal = 0
    !np = 0
    !npLocal = 0

    filename = "SPHERE_LP600.dat"
    nlpTplt = 600
    radSphere = 1.67d-3/2-0.3d0*dx
    !radSphere = 0.25d0-0.3d0*dx
    allocate(xtplt(1:3,1:nlpTplt))
    open(cpuIdx,file=adjustl(trim(filename)))
        do i=1,nlpTplt
            read(cpuIdx,*) xtplt(1,i), xtplt(2,i), xtplt(3,i)
        enddo
    close(cpuIdx)

    nMSGlo = 2
    xpc(1,1) = 5.05d-3; xpc(3,1) = 3.5d-2+1d-10; xpc(2,1) = 5.05d-3
    !xpc(1,1) = lx/2-1d-10; xpc(3,1) = lz/4-1d-10; xpc(2,1) = ly/2-1d-10
    xpc(1,2) = 4.95d-3; xpc(3,2) = 3.16d-2;      xpc(2,2) = 4.95d-3
    !xpc(1,2) = lx/2-1d-10; xpc(3,2) = lz/4*3-1d-10;      xpc(2,2) = ly/2-1d-10
    
    !deltaL = 2.5d-3
    !nMSGlo = 64
    !do k=1,4
    !do j=1,4
    !do i=1,4
    !    idx = 16*(k-1) + 4*(j-1) + i
    !    xpc(1,idx) = (i-1)*deltaL + deltaL/2d0 
    !    xpc(2,idx) = (j-1)*deltaL + deltaL/2d0 
    !    xpc(3,idx) = (k-1)*deltaL + deltaL/2d0 
    !enddo
    !enddo
    !enddo
    
    !deltaL = lx/2
    !nMSGlo = 8
    !do k=1,2
    !do j=1,2
    !do i=1,2
    !    idx = 4*(k-1) + 2*(j-1) + i
    !    xpc(1,idx) = (i-1)*deltaL + deltaL/2d0 
    !    xpc(2,idx) = (j-1)*deltaL + deltaL/2d0 
    !    xpc(3,idx) = (k-1)*deltaL + deltaL/2d0 
    !enddo
    !enddo
    !enddo
    
    rhoMS = 1.14d3

    np = 0; npLocal = 0
    nMS = 0; nMSLocal = 0
    pi = acos(-1d0)
    collisionDist = dx*0.75d0*2!+0.6d0*dx
    
    do j=1,nMSGlo
    xc1 = xpc(1,j); xc2 = xpc(2,j); xc3 = xpc(3,j)

    do i=1,nlpTplt
        x1 = radSphere*xtplt(1,i) + xpc(1,j)
        x2 = radSphere*xtplt(2,i) + xpc(2,j)
        x3 = radSphere*xtplt(3,i) + xpc(3,j)
        if (x1>xStart.and.x1<xEnd.and.x2>yStart.and.x2<yEnd.and.x3>zStart.and.x3<zEnd)then
            np = np + 1
            part(np)%x(1) = x1
            part(np)%x(2) = x2
            part(np)%x(3) = x3
            part(np)%v(:) = 0d0
            part(np)%radeff = dx*0.75d0
            part(np)%vol = 3.1415926536*dx/3/nlpTplt*(12*radSphere**2+dx**2)
            part(np)%rho = rhoMS
            part(np)%m = part(np)%rho*part(np)%vol
            part(np)%MSIdx = j
            part(np)%GlobalIdx = i

        endif
    enddo
    
    if (xc1>xStart.and.xc1<xEnd.and.xc2>yStart.and.xc2<yEnd.and.xc3>zStart.and.xc3<zEnd)then
        nMS = nMS + 1; nMSLocal = nMSLocal + 1
        MS(nMS)%x = xpc(:,j)
        MS(nMS)%v = 0d0
        !MS(nMS)%v(1) = 0.03d0*6*MS(nMS)%x(3)*(lz-MS(nMS)%x(3))/(lz**2)
        MS(nMS)%vold = 0d0
        MS(nMS)%rad = radSphere + 0.3d0*dx
        MS(nMS)%rho = rhoMS
        MS(nMS)%vol = 4d0/3*pi*MS(nMS)%rad**3
        MS(nMS)%m = MS(nMS)%vol*MS(nMS)%rho
        MS(nMS)%I = 0.4d0*MS(nMS)%m*MS(nMS)%rad**2
        MS(nMS)%globalIdx = j
        MS(nMS)%cpuIdx = cpuIdx
        MSGlo(j) = MS(nMS)
        MS(nMS)%omega = 0d0
        MS(nMS)%omega = 0d0
        write(*,*) "MS belongs to ",cpuIdx
        write(*,*) "vol ratio:",MS(nMS)%vol/part(np)%vol,part(np)%vol
    endif

    enddo

    
    npLocal = np   
    npPartInCell(1:ni*nj*nk) = 0
    call PutParcelsIntoCellsDEM(1,np,0)
return
endsubroutine
!-----------------------------------------------
subroutine InitCylinderParcelsDEM()
use vars
use ZjsVar
implicit none
    character(len=100)::filename
    integer::i,j,k,idx
    real*8::xpc(1:3,1:512),x1,x2,x3,xc1,xc2,xc3,xTmp(1:3)
    real*8::pi,rhoMS,deltaL,radratio,heightratio,height,rad
    real*8::AngInit(1:3),rotMatInv(1:9)
    integer::nTplt,flag,ntx,nty,ntz,nt,idxX,idxY,idxZ
    real*8::rnd(1:512)
    
    !filename = "CYLINDER_LP2295.dat"
    !nlpTplt = 2295
    filename = "CYLINDER_LP1125.dat"
    nlpTplt = 1125
    
    rad = 2.08e-4
    height = rad*4

    radratio = (rad-0.3d0*dx)*height/rad
    heightratio = (height-0.3d0*2*dx)
    
    allocate(xtplt(1:3,1:nlpTplt))
    open(cpuIdx,file=adjustl(trim(filename)))
        do i=1,nlpTplt
            read(cpuIdx,*) xtplt(1,i), xtplt(2,i), xtplt(3,i)
        enddo
    close(cpuIdx)

    !nMSGlo = 1 ! 2
    !!xpc(1,1) = 5.05d-3; xpc(3,1) = 3.5d-2+1d-10; xpc(2,1) = 5.05d-3
    !xpc(1,1) = lx/2-1d-10; xpc(3,1) = lz/2-1d-10; xpc(2,1) = ly/2-1d-10
    !!xpc(1,2) = 4.95d-3; xpc(3,2) = 3.16d-2;      xpc(2,2) = 4.95d-3
    !!xpc(1,2) = lx/2-1d-10; xpc(3,2) = lz/4*3-1d-10;      xpc(2,2) = ly/2-1d-10
    
    deltaL = 3d-3/3
    nMSGlo = 86
    ntx = 9; nty = 2; ntz = 9
    nt = ntx*nty*ntz
    do i=1,nMSGlo
        idx = mod(7*i-1,nt)+1
        idxZ = (idx-1)/(ntx*nty)
        idxY = (idx-1-idxZ*ntx*nty)/ntx
        idxX = mod(idx-1,ntx)
        !call random_number(rnd1)
        !call random_number(rnd2)
        !call random_number(rnd3)
        xpc(1,i) = (idxX)*deltaL + deltaL*0.5d0+1d-10 
        xpc(2,i) = (idxY)*deltaL + deltaL*0.5d0+1d-10
        xpc(2,i) = xpc(2,i)/deltaL*1.2d-3/2+1d-10
        xpc(3,i) = (idxZ)*deltaL + deltaL*0.5d0+1d-10
    enddo

    do i=1,nMSGlo
        call random_number(rnd(i))
    enddo

    rhoMS = 2.65d3

    np = 0; npLocal = 0
    nMS = 0; nMSLocal = 0
    pi = acos(-1d0)
    collisionDist = dx*0.75d0*2!+0.6d0*dx
    
    do j=1,nMSGlo
    xc1 = xpc(1,j); xc2 = xpc(2,j); xc3 = xpc(3,j)

    if (xc1>xStart.and.xc1<xEnd.and.xc2>yStart.and.xc2<yEnd.and.xc3>zStart.and.xc3<zEnd)then
        nMS = nMS + 1; nMSLocal = nMSLocal + 1
        MS(nMS)%x = xpc(:,j)
        MS(nMS)%v = 0d0
        MS(nMS)%vold = 0d0
        MS(nMS)%rad = rad
        MS(nMS)%rho = rhoMS
        MS(nMS)%vol = (pi*rad**2)*height
        MS(nMS)%m = MS(nMS)%vol*MS(nMS)%rho
        MS(nMS)%I = 0d0
        MS(nMS)%InertiaTensor(1:9) = 0
        MS(nMS)%InertiaTensor(1) = MS(nMS)%m*(3d0*rad**2+height**2)/12d0
        MS(nMS)%InertiaTensor(5) = MS(nMS)%InertiaTensor(1)
        MS(nMS)%InertiaTensor(9) = 0.5d0*MS(nMS)%m*rad**2
        AngInit(1:3) = 0d0; AngInit(2) = pi/2d0*rnd(j)
        call CalculateRotateMatrix(AngInit,MS(nMS)%rotMat)
        call InverseOrthogonalMatrix(MS(nMS)%rotMat,rotMatInv)
        
        MS(nMS)%globalIdx = j
        MS(nMS)%cpuIdx = cpuIdx
        MSGlo(j) = MS(nMS)
        MS(nMS)%omega = 0d0
        !MS(nMS)%omega(3) = 2*pi
        call TensorVectorMultiplication(rotMatInv,MS(nMS)%omega,MS(nMS)%AngMom)
        call TensorVectorMultiplication(MS(nMS)%InertiaTensor,MS(nMS)%AngMom,MS(nMS)%AngMom)
        call TensorVectorMultiplication(MS(nMS)%rotMat,MS(nMS)%AngMom,MS(nMS)%AngMom)
        MS(nMS)%AngMomOld = MS(nMS)%AngMom
        MS(nMS)%rotMatGlo = MS(nMS)%rotMat
        MS(nMS)%rotMatGloOld = MS(nMS)%rotMat

        write(*,*) "MS belongs to ",cpuIdx
    endif
    enddo

    call RebuildListMS()

    do j=1,nMSGlo

    ! Determine whether the current MS locates at neighbor procs
    flag = 1
    do i=1,nMS
        if (MS(i)%globalIdx==j) flag = 0
    enddo

    if (flag==0)then
    do i=1,nlpTplt
        x1 = radratio*xtplt(1,i)    + xpc(1,j)
        x2 = radratio*xtplt(2,i)    + xpc(2,j)
        x3 = heightratio*xtplt(3,i) + xpc(3,j)
        xTmp(1) = x1; xTmp(2) = x2; xTmp(3) = x3
        call RotateLagrangianPoints(xTmp,xpc(1:3,j),MSGlo(j)%rotMat)
        x1 = xTmp(1); x2 = xTmp(2); x3 = xTmp(3)
        if (x1>xStart.and.x1<xEnd.and.x2>yStart.and.x2<yEnd.and.x3>zStart.and.x3<zEnd)then
            np = np + 1
            part(np)%x(1) = x1
            part(np)%x(2) = x2
            part(np)%x(3) = x3
            part(np)%v(:) = 0d0
            part(np)%radeff = dx*0.75d0
            part(np)%vol = pi*dx/nlpTplt*(rad**2+2*rad*height+0.25d0*dx**2)
            part(np)%rho = rhoMS
            part(np)%m = part(np)%rho*part(np)%vol
            part(np)%MSIdx = j
            part(np)%GlobalIdx = i

        endif
    enddo 
    endif
    
    enddo
    
    npLocal = np   
    nMSLocal = nMS
    npPartInCell(1:ni*nj*nk) = 0
    call PutParcelsIntoCellsDEM(1,np,0)
return
endsubroutine
!-----------------------------------------------
subroutine SetWallProbe()
use ZjsVar
implicit none
    type(WallColProbe)::probe
    numWallProbe(1:6) = 0
    numWallProbe(2) = 1
    probe%points(1,1) = lxGlo; probe%points(2,1) = 4*dj; probe%points(3,1) = 4*dk
    probe%points(1,2) = lxGlo; probe%points(2,2) = 8*dj; probe%points(3,2) = 4*dk
    probe%points(1,3) = lxGlo; probe%points(2,3) = 8*dj; probe%points(3,3) = 8*dk
    probe%points(1,4) = lxGlo; probe%points(2,4) = 4*dj; probe%points(3,4) = 8*dk
    wallProbe(1,2) = probe
return
endsubroutine
!-----------------------------------------------
