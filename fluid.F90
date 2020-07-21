!----------------------------------------------------------
subroutine ConvectionTensor()
use vars
implicit none
    integer::i,j,k
    real*8::uuE,uuW,uvB,uvF,uwN,uwS
    real*8::vuE,vuW,vvB,vvF,vwN,vwS
    real*8::wuE,wuW,wvB,wvF,wwN,wwS

    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        
            ! 2nd order center difference | normal velocity component on CV surface         
        uuE = 0.5d0*(u(i+1,j,k)+u(i,j,k)) * 0.5d0*(u(i+1,j  ,k  ) + u(i  ,j  ,k  ))
        uuW = 0.5d0*(u(i-1,j,k)+u(i,j,k)) * 0.5d0*(u(i-1,j  ,k  ) + u(i  ,j  ,k  ))
        uvB = 0.5d0*(u(i,j+1,k)+u(i,j,k)) * 0.5d0*(v(i-1,j+1,k  ) + v(i  ,j+1,k  ))
        uvF = 0.5d0*(u(i,j-1,k)+u(i,j,k)) * 0.5d0*(v(i-1,j  ,k  ) + v(i  ,j  ,k  ))
        uwN = 0.5d0*(u(i,j,k+1)+u(i,j,k)) * 0.5d0*(w(i-1,j  ,k+1) + w(i  ,j  ,k+1))
        uwS = 0.5d0*(u(i,j,k-1)+u(i,j,k)) * 0.5d0*(w(i-1,j  ,k  ) + w(i  ,j  ,k  ))

        vuE = 0.5d0*(v(i+1,j,k)+v(i,j,k)) * 0.5d0*(u(i+1,j-1,k  ) + u(i+1,j  ,k  ))
        vuW = 0.5d0*(v(i-1,j,k)+v(i,j,k)) * 0.5d0*(u(i  ,j-1,k  ) + u(i  ,j  ,k  ))
        vvB = 0.5d0*(v(i,j+1,k)+v(i,j,k)) * 0.5d0*(v(i  ,j+1,k  ) + v(i  ,j  ,k  ))
        vvF = 0.5d0*(v(i,j-1,k)+v(i,j,k)) * 0.5d0*(v(i  ,j-1,k  ) + v(i  ,j  ,k  ))
        vwN = 0.5d0*(v(i,j,k+1)+v(i,j,k)) * 0.5d0*(w(i  ,j-1,k+1) + w(i  ,j  ,k+1))
        vwS = 0.5d0*(v(i,j,k-1)+v(i,j,k)) * 0.5d0*(w(i  ,j-1,k  ) + w(i  ,j  ,k  ))

        wuE = 0.5d0*(w(i+1,j,k)+w(i,j,k)) * 0.5d0*(u(i+1,j  ,k-1) + u(i+1,j  ,k  ))
        wuW = 0.5d0*(w(i-1,j,k)+w(i,j,k)) * 0.5d0*(u(i  ,j  ,k-1) + u(i  ,j  ,k  ))
        wvB = 0.5d0*(w(i,j+1,k)+w(i,j,k)) * 0.5d0*(v(i  ,j+1,k-1) + v(i  ,j+1,k  ))
        wvF = 0.5d0*(w(i,j-1,k)+w(i,j,k)) * 0.5d0*(v(i  ,j  ,k-1) + v(i  ,j  ,k  ))
        wwN = 0.5d0*(w(i,j,k+1)+w(i,j,k)) * 0.5d0*(w(i  ,j  ,k+1) + w(i  ,j  ,k  ))
        wwS = 0.5d0*(w(i,j,k-1)+w(i,j,k)) * 0.5d0*(w(i  ,j  ,k-1) + w(i  ,j  ,k  ))

        TIux(i,j,k) = (uuE-uuW)/dx
        TIuy(i,j,k) = (uvB-uvF)/dy
        TIuz(i,j,k) = (uwN-uwS)/dz
        TIvx(i,j,k) = (vuE-vuW)/dx
        TIvy(i,j,k) = (vvB-vvF)/dy
        TIvz(i,j,k) = (vwN-vwS)/dz
        TIwx(i,j,k) = (wuE-wuW)/dx
        TIwy(i,j,k) = (wvB-wvF)/dy
        TIwz(i,j,k) = (wwN-wwS)/dz
    enddo
    enddo
    enddo
return
end subroutine
!----------------------------------------------------------
subroutine ViscousTensor()
use vars
implicit none
    integer::i,j,k
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        TVux(i,j,k) = (u(i+1,j,k)+u(i-1,j,k)-2*u(i,j,k))/(dx*dx)*visc
        TVuy(i,j,k) = (u(i,j+1,k)+u(i,j-1,k)-2*u(i,j,k))/(dy*dy)*visc
        TVuz(i,j,k) = (u(i,j,k+1)+u(i,j,k-1)-2*u(i,j,k))/(dz*dz)*visc
        TVvx(i,j,k) = (v(i+1,j,k)+v(i-1,j,k)-2*v(i,j,k))/(dx*dx)*visc
        TVvy(i,j,k) = (v(i,j+1,k)+v(i,j-1,k)-2*v(i,j,k))/(dy*dy)*visc
        TVvz(i,j,k) = (v(i,j,k+1)+v(i,j,k-1)-2*v(i,j,k))/(dz*dz)*visc
        TVwx(i,j,k) = (w(i+1,j,k)+w(i-1,j,k)-2*w(i,j,k))/(dx*dx)*visc
        TVwy(i,j,k) = (w(i,j+1,k)+w(i,j-1,k)-2*w(i,j,k))/(dy*dy)*visc
        TVwz(i,j,k) = (w(i,j,k+1)+w(i,j,k-1)-2*w(i,j,k))/(dz*dz)*visc
    enddo
    enddo
    enddo
return
end subroutine
!----------------------------------------------------------
subroutine PressureGradient()
use vars
implicit none
    integer::i,j,k
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
       Px(i,j,k) = (P(i,j,k)-P(i-1,j,k))/(dx)/rhof
       Py(i,j,k) = (P(i,j,k)-P(i,j-1,k))/(dy)/rhof
       Pz(i,j,k) = (P(i,j,k)-P(i,j,k-1))/(dz)/rhof
    enddo
    enddo
    enddo
return
end subroutine
!----------------------------------------------------------
subroutine VelocityDivergence(cRK)
use vars
implicit none
    integer::i,j,k
    real*8::cRK
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        velDiv(i,j,k) = (u(i+1,j,k)-u(i,j,k))/dx+ &
                        (v(i,j+1,k)-v(i,j,k))/dy+ &
                        (w(i,j,k+1)-w(i,j,k))/dz
        velDiv(i,j,k) = -velDiv(i,j,k)/dt/cRK*rhof      
    enddo
    enddo
    enddo
return
end subroutine
!----------------------------------------------------------
subroutine VelocityGuess(aRK,bRK,cRK)
use vars
implicit none
    integer::i,j,k
    real*8::aRK,bRK,cRK
    call ConvectionTensor()
    call ViscousTensor()
    call PressureGradient()
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        RHSu(i,j,k) = -(TIux(i,j,k)+TIuy(i,j,k)+TIuz(i,j,k))+ &
                       (TVux(i,j,k)+TVuy(i,j,k)+TVuz(i,j,k))+CPx+fx(i,j,k)
        RHSv(i,j,k) = -(TIvx(i,j,k)+TIvy(i,j,k)+TIvz(i,j,k))+ &
                       (TVvx(i,j,k)+TVvy(i,j,k)+TVvz(i,j,k))+CPy+fy(i,j,k)
        RHSw(i,j,k) = -(TIwx(i,j,k)+TIwy(i,j,k)+TIwz(i,j,k))+ &
                       (TVwx(i,j,k)+TVwy(i,j,k)+TVwz(i,j,k))+CPz+fz(i,j,k)
        u(i,j,k) = uold(i,j,k) + dt*(aRK*RHSu(i,j,k)+bRK*RHSul(i,j,k)-cRK*Px(i,j,k))
        v(i,j,k) = vold(i,j,k) + dt*(aRK*RHSv(i,j,k)+bRK*RHSvl(i,j,k)-cRK*Py(i,j,k))
        w(i,j,k) = wold(i,j,k) + dt*(aRK*RHSw(i,j,k)+bRK*RHSwl(i,j,k)-cRK*Pz(i,j,k))
    enddo
    enddo
    enddo
    RHSul(sg1:eg1,sg2:eg2,sg3:eg3) = RHSu(sg1:eg1,sg2:eg2,sg3:eg3)
    RHSvl(sg1:eg1,sg2:eg2,sg3:eg3) = RHSv(sg1:eg1,sg2:eg2,sg3:eg3)
    RHSwl(sg1:eg1,sg2:eg2,sg3:eg3) = RHSw(sg1:eg1,sg2:eg2,sg3:eg3)
endsubroutine
!----------------------------------------------------------
subroutine VelocityPostGuess()
use vars
implicit none
    integer::i,j,k,ierr
    !call ConvectionTensor()
    !call ViscousTensor()
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        u(i,j,k) = u(i,j,k) + dt*(fx(i,j,k)-fxOld(i,j,k))
        v(i,j,k) = v(i,j,k) + dt*(fy(i,j,k)-fyOld(i,j,k))
        w(i,j,k) = w(i,j,k) + dt*(fz(i,j,k)-fzOld(i,j,k))
        !u(i,j,k) = u(i,j,k) + dt*(0.5d0*fx(i,j,k)+0.5d0*fx(i-1,j,k))
        !v(i,j,k) = v(i,j,k) + dt*(0.5d0*fy(i,j,k)+0.5d0*fy(i,j-1,k))
        !w(i,j,k) = w(i,j,k) + dt*(0.5d0*fz(i,j,k)+0.5d0*fz(i,j,k-1))
    enddo
    enddo
    enddo
endsubroutine
!----------------------------------------------------------
subroutine ModifyBulkVelocity()
use vars
implicit none
    integer::i,j,k,ierr,nTot
    real*8::uSumLocal,vSumLocal,wSumLocal,uSumGlo,vSumGlo,wSumGlo
    if (DefAveBulkUSwitch == 1 .or. DefAveBulkVSwitch == 1 .or. DefAveBulkWSwitch == 1)then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        uSumLocal = sum(u(sg1:eg1,sg2:eg2,sg3:eg3))
        vSumLocal = sum(v(sg1:eg1,sg2:eg2,sg3:eg3))
        wSumLocal = sum(w(sg1:eg1,sg2:eg2,sg3:eg3))
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        call MPI_REDUCE(uSumLocal,uSumGlo,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(uSumGlo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call MPI_REDUCE(vSumLocal,vSumGlo,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vSumGlo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call MPI_REDUCE(wSumLocal,wSumGlo,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wSumGlo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       
        nTot = nTotX*nTotY*nTotZ

        modifyPx = 0d0; modifyPy = 0d0; modifyPz = 0d0

        if (DefAveBulkUSwitch==1) modifyPx = -(uSumGlo/nTot-bulkU)/dt
        if (DefAveBulkVSwitch==1) modifyPy = -(vSumGlo/nTot-bulkV)/dt
        if (DefAveBulkWSwitch==1) modifyPz = -(wSumGlo/nTot-bulkW)/dt

        do k=sg3,eg3
        do j=sg2,eg2
        do i=sg1,eg1
            u(i,j,k) = u(i,j,k) + dt*modifyPx
            v(i,j,k) = v(i,j,k) + dt*modifyPy
            w(i,j,k) = w(i,j,k) + dt*modifyPz
        enddo
        enddo
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    endif
return
endsubroutine
!----------------------------------------------------------
subroutine PressureCorrection()
use vars
implicit none
    integer::i,j,k
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        P(i,j,k) = P(i,j,k) + Pcorr(i,j,k)
    enddo
    enddo
    enddo
return
endsubroutine
!----------------------------------------------------------
subroutine VelocityCorrection(cRK)
use vars
implicit none
    integer::i,j,k
    real*8::cRK
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        u(i,j,k) = u(i,j,k) + cRK*Px(i,j,k)*dt - cRK*(P(i,j,k)-P(i-1,j,k))/(dx)*dt/rhof
        v(i,j,k) = v(i,j,k) + cRK*Py(i,j,k)*dt - cRK*(P(i,j,k)-P(i,j-1,k))/(dy)*dt/rhof
        w(i,j,k) = w(i,j,k) + cRK*Pz(i,j,k)*dt - cRK*(P(i,j,k)-P(i,j,k-1))/(dz)*dt/rhof
    enddo
    enddo
    enddo
    do k=sg3,eg3
    do j=sg2,eg2
    do i=sg1,eg1
        omegaX(i,j,k) = (w(i,j+1,k)-w(i,j,k))/dy-(v(i,j,k+1)-v(i,j,k))/dz
        omegaY(i,j,k) = (u(i,j,k+1)-u(i,j,k))/dz-(w(i+1,j,k)-w(i,j,k))/dx
        omegaZ(i,j,k) = (v(i+1,j,k)-v(i,j,k))/dx-(u(i,j+1,k)-u(i,j,k))/dy
    enddo
    enddo
    enddo

return
endsubroutine
!----------------------------------------------------------
subroutine Boundary()
use vars
implicit none
    if (bType(1)==1) u(lg1:sg1,:,:) = 0d0
    if (bType(2)==1) u(eg1+1:ug1,:,:) = 0d0
    if (bType(3)==1) v(:,lg2:sg2,:) = 0d0
    if (bType(4)==1) v(:,eg2+1:ug2,:) = 0d0
    if (bType(5)==1) w(:,:,lg3:sg3) = 0d0
    if (bType(6)==1) w(:,:,eg3+1:ug3) = 0d0
return
endsubroutine
!------------------------------------------------------------
subroutine FirstOrderForwardEuler()
use vars
use ZjsVar
implicit none
    integer::Ns,i,j,k,ierr
    Ns = 5
    uold = u; vold = v; wold = w

    do i=1,npLocal
        part(i)%fdrag(:) = 0d0
    enddo
    fx = 0d0; fy = 0d0; fz = 0d0
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i=1,Ns
        call UpdateFluidVelocityOnLagrangianPoints()
        call UpdateParcelInfoInGhostRegion()
        !open(cpuIdx,file=logFile,access="append")
        !do j=1,np
        !    write(cpuIdx,*)"part drag all:",part(j)%fdrag
        !enddo
        !close(cpuIdx)
        call UpdateMomentumSourceTerm()

        call VelocityPostGuess()
        call FillGhostRegionForFluid()
        call Boundary()

        !open(cpuIdx,file=logFile,access="append")
        !do j = 1,npLocal
        !    write(cpuIdx,*)"part drag", part(j)%GlobalIdx,part(j)%vgap, part(j)%fdrag
        !enddo
        !close(cpuIdx)
        
    enddo
    call VelocityGuess(1d0,0d0,1d0)
    call FillGhostRegionForFluid()
    call Boundary()

    call ModifyBulkVelocity()
    call FillGhostRegionForFluid()
    call Boundary()

    !open(cpuIdx,file=logFile,access="append")
    !do k=sg3,eg3
    !do j=sg2,eg2
    !do i=sg1,eg1
    !    if (i>8.and.i<25.and.j>8.and.j<25.and.k>8.and.k<25)then
    !        if (abs(fx(i,j,k))+abs(fy(i,j,k))+abs(fz(i,j,k))>1d-20)then
    !            write(cpuIdx,*) i,j,k,fx(i,j,k),fy(i,j,k),fz(i,j,k)
    !        endif
    !    endif
    !enddo
    !enddo
    !enddo
    !close(cpuidx)


    call VelocityDivergence(1d0)
    call FillGhostRegionForFluid()
    call Boundary()


    call SetPressureVector()
    call PoissonSolver()
    call PressureCorrection()
    call FillGhostRegionForFluid()
    call Boundary()


    call VelocityCorrection(1d0)
    call FillGhostRegionForFluid()
    call Boundary()
    
    call CloudEvolutionMS()
return
endsubroutine
!------------------------------------------------------------
subroutine ThirdOrderLowStorageRungeKutta()
use vars
implicit none
    real*8::aRK(1:3),bRK(1:3),cRK(1:3)
    integer::iRK
    aRK(1) = 32d0/60d0; bRK(1) = 0d0
    aRK(2) = 25d0/60d0; bRK(2) = -17d0/60d0
    aRK(3) = 45d0/60d0; bRK(3) = -25d0/60d0
    cRK = aRK+bRK
    if (step == 1)then
        open(cpuIdx,file=logFile,access="append")
        write(cpuIdx,*)aRK
        write(cpuIdx,*)bRK
        write(cpuIdx,*)cRK
        close(cpuIdx)
    endif

    do iRK=1,3
        call VelocityGuess(aRK(iRK),bRK(iRK),cRK(iRK))
        call FillGhostRegionForFluid()
        call Boundary()

        call VelocityDivergence(cRK(iRK))
        call FillGhostRegionForFluid()
        call Boundary()

        call SetPressureVector()
        call PoissonSolver()
        call PressureCorrection()
        call FillGhostRegionForFluid()
        call Boundary()

        call VelocityCorrection(cRK(iRK))
        call FillGhostRegionForFluid()
        call Boundary()
    enddo
return
endsubroutine
!------------------------------------------------------------
subroutine PostProcessing()
use vars
use ZjsVar
implicit none
   integer::i,j,k,ierr
   real*8::xTmp,yTmp,zTmp,pRefSim,pRefAna
   real*8::maxErrorUTmp,maxErrorVTmp,maxErrorPTmp,maxErrorU,maxErrorV,maxErrorP
   real*8::maxUTmp,maxVtmp,maxPTmp,maxU,maxV,maxP,ratio
   character(len=20)::msStr
   
   maxUTmp = sum(u(sg1:eg1,sg2:eg2,sg3:eg3))
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(maxUTmp,maxU,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   if (cpuIdx == 0 .and. (mod(step,50)==0 .or. step==0 .or. step==1))then
    write(*,*) "Darcy number:",maxU*1.67d-3/0.2336/visc/(nTotX*nTotY*nTotZ)
    write(*,*) "Reference ratio:",maxU/0.03d0/(nTotX*nTotY*nTotZ)
    write(*,*) "Modified pressure gradient:", modifyPx,modifyPy,modifyPz
    !write(*,*) "Ratio & Accelleration:",nTotX*nTotY*nTotZ*MS(1)%v(1)/maxU,MS(1)%fdrag/MS(1)%m*rhof
   endif
   if (nMSLocal>0 .and. (mod(step,50)==0 .or. step==0 .or. step==1))then
    write(*,*) "Angular velocity:",MS(1)%omega
    write(*,*) "Velocity:",MS(1)%v
    write(*,*)
   endif
   do i=1,nMSLocal
   ratio = (MS(i)%rho-rhof)/MS(i)%rho
    write(msStr,*) MS(i)%GlobalIdx
        open(cpuIdx,file="./DEM/"//trim(adjustl(msStr))//".probe",access="append")
        write(cpuIdx,*) step,MS(i)%x,MS(i)%v,MS(i)%omega,MS(i)%fdrag/MS(i)%m*rhof,MS(i)%fmom/MS(i)%I*rhof,MS(i)%fcol/MS(i)%m
        close(cpuIdx)
   enddo

   if (cpuIdx==0)then
        open(cpuIdx,file="./DEM/wallfriction.log",access="append")
        write(cpuIdx,*) step,modifyPx,modifyPy,modifyPz
        close(cpuIdx)
   endif

   !if (cpuIdx==0)then
   !    pRefAna = -rhof/4d0*(cos(2d0*0.5d0*dx)+cos(2d0*0.5d0*dy))*exp(-4d0*visc*time)
   !    pRefSim = P(1,1,1)
   !endif
   !
   !call MPI_BCAST(pRefSim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   !call MPI_BCAST(pRefAna,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

   !!if (cpuIdx==0)then
   !! write(*,*)"Analytical pressure reference value ",pRefAna
   !! write(*,*)"Simulated pressure reference value ",pRefSim
   !! write(*,*)"Error:",pRefSim-pRefAna
   !!endif

   !do k=sg3,eg3
   !do j=sg2,eg2
   !do i=sg1,eg1
   !     xTmp = xc(i,j,k)!
   !     yTmp = yc(i,j,k)!
   !     zTmp = zc(i,j,k)!
   !                     !
   !     errorP(i,j,k) = -rhof/4d0*(cos(2d0*xTmp)+cos(2d0*yTmp))*exp(-4d0*visc*time)-pRefAna
   !     
   !     xTmp = xc(i,j,k)-0.5d0*dx
   !     errorU(i,j,k) = cos(xTmp)*sin(yTmp)*exp(-2d0*visc*time)
   !     
   !     xTmp = xc(i,j,k)
   !     yTmp = yc(i,j,k)-0.5d0*dy
   !     errorV(i,j,k) = -sin(xTmp)*cos(yTmp)*exp(-2d0*visc*time)
   !       
   !enddo
   !enddo
   !enddo
   !errorU(sg1:eg1,sg2:eg2,sg3:eg3) = u(sg1:eg1,sg2:eg2,sg3:eg3) - errorU(sg1:eg1,sg2:eg2,sg3:eg3)
   !errorV(sg1:eg1,sg2:eg2,sg3:eg3) = v(sg1:eg1,sg2:eg2,sg3:eg3) - errorV(sg1:eg1,sg2:eg2,sg3:eg3)
   !errorP(sg1:eg1,sg2:eg2,sg3:eg3) = P(sg1:eg1,sg2:eg2,sg3:eg3) - errorP(sg1:eg1,sg2:eg2,sg3:eg3) - pRefSim 

   !maxErrorUTmp = maxval(errorU); maxErrorVTmp = maxval(errorV); maxErrorPTmp = maxval(abs(errorP))
   !call MPI_REDUCE(maxErrorUTmp,maxErrorU,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   !call MPI_REDUCE(maxErrorVTmp,maxErrorV,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   !call MPI_REDUCE(maxErrorPTmp,maxErrorP,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   !maxUTmp = maxval(U); maxVTmp = maxval(V); maxPTmp = maxval(P)
   !call MPI_REDUCE(maxUTmp,maxU,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   !call MPI_REDUCE(maxVTmp,maxV,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   !call MPI_REDUCE(maxPTmp,maxP,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

   !if (cpuIdx == 0 .and. (mod(step,50)==0 .or. step==0 .or. step==1))then
   !    write(*,*) "max U and Uerror",maxU,maxErrorU/maxU
   !    write(*,*) "max V and Verror",maxV,maxErrorV/maxV
   !    write(*,*) "max P and Perror",maxP,maxErrorP/(maxP-pRefSim)
   !endif
return
endsubroutine
!------------------------------------------------------------

