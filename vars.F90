module vars
!use petsc
!#include<petsc/finclude/petscsys.h>
!#include<petsc/finclude/petsc.h>
#include<petsc/finclude/petscksp.h>
use mpi
use petscksp
implicit none
    !include "mpif.h"
    !--Grid
    ! Co-located grid for pressure and velocity CVs
    ! Use Rhie-Chow's (1986) scheme to avoid checkboard pattern
    integer::nx,ny,nz
    integer,parameter::lap = 2
    real*8::dx,dy,dz
    real*8::xoffset,yoffset,zoffset
    real*8::lx,ly,lz
    real*8,pointer,dimension(:,:,:)::xc,yc,zc
    ! Boundary type: 0-period, 1-solid wall, -1-internal
    integer::bw,be,bb,bf,bs,bn
    ! Global boundary type
    integer::bwg,beg,bsg,bng,bfg,bbg
    integer,dimension(1:6)::bType
    integer,dimension(1:3)::sg,eg,lg,ug
    integer::sg1,sg2,sg3,eg1,eg2,eg3,lg1,lg2,lg3,ug1,ug2,ug3
    integer::nTotX,nTotY,nTotZ
    ! Store global grid partition information in each direction
    integer,pointer,dimension(:)::nxList,nyList,nzList
    ! Store global start cell index for each CPU
    integer,pointer,dimension(:,:,:)::GlobalStartCellIndex
    integer,pointer,dimension(:)::GlobalStartCellIndex1D
    integer,dimension(1:6)::CPUNeighIdx,CPUTransferIdx
    integer::CPUIdx,CPUIdxX,CPUIdxY,CPUIdxZ,CPUnx,CPUny,CPUnz
    integer,pointer,dimension(:,:,:)::GlobalCellIndex
    ! Store 6 global neighbor cell indice 
    integer,pointer,dimension(:,:,:,:)::GlobalNeighCellIndex


    !--Fluid
    real*8::visc,rhof
    real*8::dt,time
    integer::step
    ! Four variables in incompressible N-S system
    real*8,pointer,dimension(:,:,:)::u,v,w,P,Pcorr,velDiv,omegaX,omegaY,omegaZ
    real*8,pointer,dimension(:,:,:)::errorU,errorV,errorW,errorP
    real*8,pointer,dimension(:,:,:)::ustar,vstar,wstar
    real*8,pointer,dimension(:,:,:)::uold,vold,wold
    ! Right Hand Side temporary storage for two iteration steps in R-K loop
    ! following Breugem's work (2012, JCP), including viscous, constant pressure gradient
    ! and inertial forces 
    real*8,pointer,dimension(:,:,:)::RHSu,RHSv,RHSw,RHSul,RHSvl,RHSwl
    ! constant pressure gradient
    real*8::CPx,CPy,CPz
    ! Defined averaged bulk velocity
    integer::DefAveBulkUSwitch,DefAveBulkVSwitch,DefAveBulkWSwitch
    real*8::bulkU,bulkV,bulkW,modifyPx,modifyPy,modifyPz
    ! Inertial and viscous tensor gradient in N-S equation
    real*8,pointer,dimension(:,:,:)::TIux,TIuy,TIuz,TIvx,TIvy,TIvz,TIwx,TIwy,TIwz
    real*8,pointer,dimension(:,:,:)::TVux,TVuy,TVuz,TVvx,TVvy,TVvz,TVwx,TVwy,TVwz
    ! Pressure gradient in N-S equation
    real*8,pointer,dimension(:,:,:)::Px,Py,Pz
    ! Source due to existence of LPs
    real*8,pointer,dimension(:,:,:)::fx,fy,fz,fxOld,fyOld,fzOld
    
    ! Transform between 3D indice and 1D index for cells
    integer,pointer,dimension(:,:,:)::cellIdx3Dto1D
    integer,pointer,dimension(:,:)::cellIdx1Dto3D

    ! MPI transfer data buffer
    integer::maxSendNoFluid
    integer,pointer,dimension(:)::recvCntFluid,sendCntFluid
    integer,pointer,dimension(:,:)::recvIdxFluid,sendIdxFluid
    real*8,pointer,dimension(:)::recvBuffFluid,sendBuffFluid

    ! IO Infos
    character(len=100)::logFile

    ! Poisson equation
    KSP ksp
    Mat pMat
    Vec pCorrVec,rhsP,pVecLocal
    PetscInt,pointer,dimension(:)::toLocalIdx,fromGlobalIdx
    IS fromIS,toIS
    VecScatter scatter
    PetscScalar,pointer,dimension(:)::pCorr1D
    ! Const pressure grid index
    integer::pConstIdx
endmodule
