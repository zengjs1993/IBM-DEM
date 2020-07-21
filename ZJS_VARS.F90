module ZjsVar
    use vars
    implicit none
    integer,parameter::maxPartNo=100000,maxCollisionNo=200,max_neigh_cell=125,maxSendPartNo=5000
    integer,parameter::maxPartNoInCell=50,maxNeighCellNo=125
    type Parcel
        integer::globalIdx,localIdx,cellIdx,pType,MSIdx
        real*8::radeff,rad,rho,vol,m
        real*8,dimension(1:3)::x,v,fcol,fdrag,vgap
        real*8::wdrag
        real*8::wsum
    end type
    ! send and receive buffer
    integer,parameter::numSendInt=5,numSendReal=17
    integer,dimension(1:numSendInt*maxSendPartNo)::sendIntBuff,recvIntBuff
    real*8,dimension(1:numSendReal*maxSendPartNo)::sendRealBuff,recvRealBuff
    ! initial parameters
    real*8::alphapInit
    ! parcels
    type(Parcel),dimension(1:maxPartNo)::part
    ! total parcel number
    integer::np,npLocal
    ! sending numbers and start indice for six sides
    integer,dimension(1:6)::npSend,npRecv,startIdxRecv
    ! store sending parcel information 
    integer,dimension(1:maxSendPartNo,1:6)::partSendIdx
    type(Parcel),dimension(1:maxSendPartNo,1:6)::partSend
    
    ! neighbor parcel numbers and lists for collision 
    integer,dimension(:),allocatable::neighPartNo
    integer,dimension(:,:),allocatable::neighPartIdx
    ! neighbor cell numbers and lists for interpolation
    integer,dimension(:),allocatable::neighCellNo
    integer,dimension(:,:),allocatable::neighCellIdx
    real*8,dimension(:,:,:),allocatable::neighCellWeight
    
    ! parcel numbers and lists in cells 
    integer,dimension(:),allocatable::npPartInCell
    integer,dimension(:,:),allocatable::partIdxInCell
    ! neigh cell indice for cells
    integer,dimension(:,:),allocatable::neighCellIdxForCell
    ! mapping between 1D and 3D cell indice
    integer,dimension(:,:),allocatable::map1Dto3D
    integer,dimension(:,:,:),allocatable::map3Dto1D
    integer,dimension(:,:,:),allocatable::cellType
    
    ! Collision parameters
    real*8::collisionDist,collisionBufferRecord
    real*8,parameter::collisionBuffer=0.1d0
    real*8::knCoeff,ResistCoeff

    ! Time control
    integer::timeCntDEM
    integer::timeDEM
    real*8::dtDEM

    ! grid information
    ! local grid numbers
    integer::ni,nj,nk
    real*8::di,dj,dk
    real*8::xStart,yStart,zStart,xEnd,yEnd,zEnd
    integer::cpuLocal
    ! boundary type: -1-internal, 0-periodic, 1-solid
    integer,dimension(6)::neighCpu
    real*8::lxGlo,lyGlo,lzGlo
    real*8,dimension(6)::coordMirror,coordShift
    ! grid numbers of each halo region
    integer,parameter::haloX=2,haloY=2,haloZ=2
    real*8::errorX,errorY,errorZ
    ! solid boundary normals and points
    real*8,dimension(1:3,1:6)::normBound,pointBound

    ! Eulerian variables
    real*8,dimension(:,:,:),allocatable::srcMomX,srcMomY,srcMomZ,srcEng,aveVelpX,aveVelpY,aveVelpZ

    ! Fluid parameters
    real*8,parameter::viscf=1D-5
    
    ! templates of lagrangian points
    real*8::radSphere
    integer::nlpTplt 
    real*8,pointer,dimension(:,:)::xtplt

    ! Wall collision force probes
    type WallColProbe
        integer::type
        real*8,dimension(1:3,1:4)::points
        real*8,dimension(1:3)::norm
        real*8::area,pressure
        real*8,dimension(1:3)::force
    end type
    integer,parameter::maxNumWallProbe=10
    type(WallColProbe),dimension(1:maxNumWallProbe,1:6)::wallProbe
    integer,dimension(1:6)::numWallProbe

    ! Multisphere
    type MultiSphere
        integer::globalIdx,cpuIdx
        real*8,dimension(1:9)::rotMat,rotMatGlo,rotMatGloOld,InertiaTensor
        real*8,dimension(1:3)::x,xold,v,vold,omega,omegaold,fcol,fmom,fdrag,AngMom,AngMomOld,fhydro,fmomHydro,fmomCol
        real*8::I,rad,rho,vol,m
    end type
    
    integer::nMS,nMSlocal,nMSglo
    integer,parameter::maxMSNo=1000,maxGlobalMSNo=1000
    type(MultiSphere),dimension(1:maxMSNo)::MSglo,MS
    
    integer,dimension(1:6)::nMSSend,nMSRecv,startIdxMSRecv
    integer,dimension(1:maxMSNo,1:6)::MSSendIdx
    type(MultiSphere),dimension(1:maxMSNo,1:6)::MSSend

    integer,parameter::numSendMSInt=5,numSendMSReal=17
    integer,dimension(1:numSendInt*maxMSNo)::sendMSIntBuff,recvMSIntBuff
    real*8,dimension(1:numSendReal*maxMSNo)::sendMSRealBuff,recvMSRealBuff

    integer::coupleInterval
endmodule
