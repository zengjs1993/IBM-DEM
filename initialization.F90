!----------------------------------------------------------
subroutine PartitionInitialization(cpux,cpuy,cpuz,nnx,nny,nnz,lenx,leny,lenz)
use vars
implicit none
    integer::cpux,cpuy,cpuz,nnx,nny,nnz
    integer::i,j,k,ii,jj,kk,ierr
    real*8::lenx,leny,lenz
    CPUnx = cpux; CPUny = cpuy; CPUnz = cpuz
    nTotX = nnx; nTotY = nny; nTotZ = nnz
    lx = lenx; ly = leny; lz = lenz
    dx = lx/nnx; dy = ly/nny; dz = lz/nnz
    nx = nnx/cpux; ny = nny/cpuy; nz = nnz/cpuz
    call MPI_COMM_RANK(MPI_COMM_WORLD,CPUIdx,ierr)
    CPUIdxZ = CPUIdx/(cpux*cpuy)
    CPUIdxY = mod(CPUIdx,cpux*cpuy)/cpux
    CPUIdxX = mod(CPUIdx,cpux)
    
    call SetBoundaryType()

    allocate(nxList(0:cpux),nyList(0:cpuy),nzList(0:cpuz))
    allocate(GlobalStartCellIndex(1:cpux,1:cpuy,1:cpuz))
    allocate(GlobalStartCellIndex1D(1:cpux*cpuy*cpuz))
    
    nxList(:) = nx; nyList(:) = ny; nzList(:) = nz
    nxList(0) = 0;  nyList(0) = 0;  nzList(0) = 0
    do i=1,mod(nnx,cpux)
        nxList(i) = nxList(i) + 1
    enddo
    do j=1,mod(nny,cpuy)
        nyList(j) = nyList(j) + 1
    enddo
    do k=1,mod(nnz,cpuz)
        nzList(k) = nzList(k) + 1
    enddo

    do k=1,cpuz
    do j=1,cpuy
    do i=1,cpux
        GlobalStartCellIndex(i,j,k) = nxList(i)*nyList(j)*nzList(k)
    enddo
    enddo
    enddo


    write(logFile,*) cpuIdx,".log"
    logFile = "./DEM/log/"//trim(adjustl(logFile))
    do k=1,cpuz
    do j=1,cpuy
    do i=1,cpux
        if (i==1)then
            ii = cpux
            if (j==1) then
                jj = cpuy; kk = k-1
            else
                jj = j-1; kk = k
            endif
        else
            ii = i-1; jj = j; kk=k
        endif
        if (i/=1.or.j/=1.or.k/=1)then
            GlobalStartCellIndex(i,j,k) = GlobalStartCellIndex(i,j,k) + GlobalStartCellIndex(ii,jj,kk)
        endif
        GlobalStartCellIndex1D((k-1)*cpux*cpuy+(j-1)*cpux+i) = &
        GlobalStartCellIndex(i,j,k)-nxList(i)*nyList(j)*nzList(k)
        open(cpuIdx,file=logFile,access="append")
        write(cpuIdx,*) GlobalStartCellIndex1D((k-1)*cpux*cpuy+(j-1)*cpux+i)
        close(cpuIdx)
    enddo
    enddo
    enddo

    nx = nxList(CPUIdxX+1)
    ny = nyList(CPUIdxY+1)
    nz = nzList(CPUIdxZ+1)

    call SetNeighborCPUIndex()

    do i=2,cpux
        nxList(i) = nxList(i) + nxList(i-1)
    enddo
    do j=2,cpuy
        nyList(j) = nyList(j) + nyList(j-1)
    enddo
    do k=2,cpuz
        nzList(k) = nzList(k) + nzList(k-1)
    enddo
    
    sg1 = nxList(CPUIdxX+1)-nx+1; eg1 = sg1+nx-1
    sg2 = nyList(CPUIdxY+1)-ny+1; eg2 = sg2+ny-1
    sg3 = nzList(CPUIdxZ+1)-nz+1; eg3 = sg3+nz-1
    lg1 = sg1-lap; ug1 = eg1+lap;
    lg2 = sg2-lap; ug2 = eg2+lap;
    lg3 = sg3-lap; ug3 = eg3+lap;

    sg(1) = sg1; sg(2) = sg2; sg(3) = sg3
    eg(1) = eg1; eg(2) = eg2; eg(3) = eg3
    lg(1) = lg1; lg(2) = lg2; lg(3) = lg3
    ug(1) = ug1; ug(2) = ug2; ug(3) = ug3
    
    call allocateMemory()
    call SetNeighCellIndex()
return
endsubroutine
!----------------------------------------------------------
subroutine SetNeighCellIndex()
use vars
implicit none
    integer::i,j,k,l,shift(1:3),is,js,ks,dim,flag1,flag2
    integer::cpuNeigh,startIdx
    integer::cpui,cpuj,cpuk,ii,jj,kk,ni,nj,nk

    GlobalNeighCellIndex(1:6,sg1:eg1,sg2:eg2,sg3:eg3) = -1
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1        
    do l=1,6
        cpuNeigh = CPUNeighIdx(l)
        shift(:) = 0
        dim = (l+1)/2
        shift(dim) = 2*mod(l+1,2) - 1
        is = i+shift(1)
        js = j+shift(2)
        ks = k+shift(3)

        flag1 = 0
        flag2 = 0
        if (dim==1.and.((is<1).or.(is>nTotX)))then
            if (bType(l)/=1)then
                is = mod(is+nTotX-1,nTotX)+1 
            else
                flag1 = 1
            endif
        endif
        if (dim==2.and.((js<1).or.(js>nTotY)))then
            if (bType(l)/=1)then
                js = mod(js+nTotY-1,nTotY)+1
            else
                flag1 = 1
            endif
        endif
        if (dim==3.and.((ks<1).or.(ks>nTotZ)))then
            if (bType(l)/=1)then
                ks = mod(ks+nTotZ-1,nTotZ)+1
            else
                flag1 = 1
            endif
        endif
        if (dim==1.and.((is<sg1).or.(is>eg1))) flag2 = 1
        if (dim==2.and.((js<sg2).or.(js>eg2))) flag2 = 1
        if (dim==3.and.((ks<sg3).or.(ks>eg3))) flag2 = 1
        
        if (flag1==0)then
            if (flag2==0) cpuNeigh = CPUIdx
            startIdx = GlobalStartCellIndex1D(cpuNeigh+1)
            cpuk = (cpuNeigh)/(CPUnx*CPUny) + 1
            cpuj = mod(cpuNeigh,CPUnx*CPUny)/CPUnx + 1
            cpui = mod(cpuNeigh,CPUnx) + 1
            ni = nxList(cpui)-nxList(cpui-1)
            nj = nyList(cpuj)-nyList(cpuj-1)
            nk = nzList(cpuk)-nzList(cpuk-1)
            ii = is-nxList(cpui-1)
            jj = js-nyList(cpuj-1)
            kk = ks-nzList(cpuk-1)
            GlobalNeighCellIndex(l,i,j,k) = startIdx + (kk-1)*ni*nj + (jj-1)*ni + ii
        else
            GlobalNeighCellIndex(l,i,j,k) = -1
        endif
    enddo
    enddo
    enddo
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine SetNeighborCPUIndex()
use vars
implicit none
    if (bw==0)then
        CPUNeighIdx(1) = CPUIdx+CPUnx-1
        CPUTransferIdx(1) = CPUIdx+CPUnx-1
    else if (bw==1)then
        CPUNeighIdx(1) = -1
        CPUTransferIdx(1) = CPUIdx+CPUnx-1
    else
        CPUNeighIdx(1) = CPUIdx-1 
        CPUTransferIdx(1) = CPUIdx-1 
    endif
    if (be==0)then
        CPUNeighIdx(2) = CPUIdx-(CPUnx-1)
        CPUTransferIdx(2) = CPUIdx-(CPUnx-1)
    else if (be==1)then
        CPUNeighIdx(2) = -1
        CPUTransferIdx(2) = CPUIdx-(CPUnx-1)
    else
        CPUNeighIdx(2) = CPUIdx+1
        CPUTransferIdx(2) = CPUIdx+1
    endif
    if (bb==0)then
        CPUNeighIdx(3) = CPUIdx+CPUnx*(CPUny-1)
        CPUTransferIdx(3) = CPUIdx+CPUnx*(CPUny-1)
    else if (bb==1)then
        CPUNeighIdx(3) = -1
        CPUTransferIdx(3) = CPUIdx+CPUnx*(CPUny-1)
    else
        CPUNeighIdx(3) = CPUIdx-CPUnx
        CPUTransferIdx(3) = CPUIdx-CPUnx
    endif
    if (bf==0)then
        CPUNeighIdx(4) = CPUIdx-CPUnx*(CPUny-1)
        CPUTransferIdx(4) = CPUIdx-CPUnx*(CPUny-1)
    else if (bf==1)then
        CPUNeighIdx(4) = -1
        CPUTransferIdx(4) = CPUIdx-CPUnx*(CPUny-1)
    else
        CPUNeighIdx(4) = CPUIdx+CPUnx
        CPUTransferIdx(4) = CPUIdx+CPUnx
    endif
    if (bs==0)then
        CPUNeighIdx(5) = CPUIdx+CPUnx*CPUny*(CPUnz-1)
        CPUTransferIdx(5) = CPUIdx+CPUnx*CPUny*(CPUnz-1)
    else if (bs==1)then
        CPUNeighIdx(5) = -1
        CPUTransferIdx(5) = CPUIdx+CPUnx*CPUny*(CPUnz-1)
    else
        CPUNeighIdx(5) = CPUIdx-CPUnx*CPUny
        CPUTransferIdx(5) = CPUIdx-CPUnx*CPUny
    endif
    if (bn==0)then
        CPUNeighIdx(6) = CPUIdx-CPUnx*CPUny*(CPUnz-1)
        CPUTransferIdx(6) = CPUIdx-CPUnx*CPUny*(CPUnz-1)
    else if (bn==1)then
        CPUNeighIdx(6) = -1
        CPUTransferIdx(6) = CPUIdx-CPUnx*CPUny*(CPUnz-1)
    else
        CPUNeighIdx(6) = CPUIdx+CPUnx*CPUny
        CPUTransferIdx(6) = CPUIdx+CPUnx*CPUny
    endif
return
endsubroutine
!----------------------------------------------------------
subroutine GridInitialization()
use vars
implicit none
    integer::i,j,k
    xoffset = (sg1-1)*dx
    yoffset = (sg2-1)*dy
    zoffset = (sg3-1)*dz
    do k=lg3,ug3
    do j=lg2,ug2
    do i=lg1,ug1
        xc(i,j,k) = (i-0.5d0)*dx
        yc(i,j,k) = (j-0.5d0)*dy
        zc(i,j,k) = (k-0.5d0)*dz
    enddo
    enddo
    enddo
    visc = 1d-6
    rhof = 1d3
    dt = 1d-4/2
    time = 0d0
    step = 0
    CPx = 0d0!0.03d0*visc*12/lz**2
    CPy = 0d0
    CPz = 0d0
    DefAveBulkUSwitch = 1
    DefAveBulkVSwitch = 0
    DefAveBulkWSwitch = 1
    bulkU = 0.1d0; bulkV = 0; bulkW = 0;
return
endsubroutine
!----------------------------------------------------------
subroutine SetCellIndexMap()
use vars
implicit none
    integer::i,j,k,idx
    do k=lg3,ug3
    do j=lg2,ug2
    do i=lg1,ug1
        ! 1D index starts from 1
        idx = (k-lg3)*(nx+2*lap)*(ny+2*lap)+(j-lg2)*(nx+2*lap)+i-lg1+1
        cellIdx3Dto1D(i,j,k) = idx
        cellIdx1Dto3D(1,idx) = i
        cellIdx1Dto3D(2,idx) = j
        cellIdx1Dto3D(3,idx) = k
    enddo
    enddo
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine FluidInitialization()
use vars
implicit none
    integer::i,j,k,ierr
    real*8::xTmp,yTmp,zTmp,Height
    do k=lg3,ug3
    do j=lg2,ug2
    do i=lg1,ug1
        !Taylor-Green vortex benchmark
        !xTmp = xc(i,j,k)
        !yTmp = yc(i,j,k)
        !zTmp = zc(i,j,k)
        !P(i,j,k) = -rhof/4d0*(cos(2d0*xTmp)+cos(2d0*yTmp))
        !omegaZ(i,j,k) = -2*cos(xTmp)*sin(yTmp)
        !xTmp = xc(i,j,k)-0.5d0*dx
        !u(i,j,k) = cos(xTmp)*sin(yTmp)
        !xTmp = xc(i,j,k)
        !yTmp = yc(i,j,k)-0.5d0*dy
        !v(i,j,k) = -sin(xTmp)*cos(yTmp)
        !w(i,j,k) = 0d0

        xTmp = xc(i,j,k)
        yTmp = yc(i,j,k)
        zTmp = zc(i,j,k)
        Height = 1d-2
        u(i,j,k) = -6d0*0.03d0/Height**2*zTmp*(zTmp-Height)
    enddo
    enddo
    enddo

    u(:,:,:) = 0
    v(:,:,:) = 0
    w(:,:,:) = 0
    P(:,:,:) = 0 
    fx(:,:,:) = 0
    fy(:,:,:) = 0
    fz(:,:,:) = 0
    RHSu(:,:,:) = 0
    RHSv(:,:,:) = 0
    RHSw(:,:,:) = 0
    RHSul(:,:,:) = 0
    RHSvl(:,:,:) = 0
    RHSwl(:,:,:) = 0
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call FillGhostRegionForFluid()
    call Boundary()
return
endsubroutine
!---------------------------------------------------------
subroutine SetBoundaryType()
use vars
implicit none
    bwg = 0 
    beg = 0 
    bbg = 1
    bfg = 1 
    bsg = 0 
    bng = 0 

    bw=-1; be=-1; bb=-1; bf=-1; bs=-1; bn=-1
    if (cpuIdxX==0)then
        bw = bwg
    endif
    if (cpuIdxX==CPUnx-1)then
        be = beg
    endif
    if (cpuIdxY==0)then
        bb = bbg
    endif
    if (cpuIdxY==CPUny-1)then
        bf = bfg
    endif
    if (cpuIdxZ==0)then
        bs = bsg
    endif
    if (cpuIdxZ==CPUnz-1)then
        bn = bng
    endif
    bType(1) = bw
    bType(2) = be
    bType(3) = bb
    bType(4) = bf
    bType(5) = bs
    bType(6) = bn
return
endsubroutine
!---------------------------------------------------------
