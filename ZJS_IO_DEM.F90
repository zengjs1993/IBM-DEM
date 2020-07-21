!-----------------------------------------------
subroutine WriteDEM()
use vars
use ZjsVar
implicit none
    character(len=100)::filename
    integer::i,j,idx
    integer::ierr
    write(filename,*) cpuLocal,".proc"
    open(cpuLocal,file="./DEM/result/"//trim(adjustl(filename)))
    write(cpuLocal,*),"ni, nj, nk = ",ni," ",nj," ",nk
    write(cpuLocal,*),"di, dj, dk = ",di," ",dj," ",dk
    write(cpuLocal,*),"xS, yS, zS = ",xStart," ",yStart," ",zStart 
    write(cpuLocal,*),"xE, yE, zE = ",xEnd,yEnd,zEnd
    write(cpuLocal,*),"boundary type = ",bType(1)," ",bType(2)," ",bType(3)," ",bType(4)," ",bType(5)," ",bType(6)
    write(cpuLocal,*),"SGi~EGi",SG1," ",EG1," ",SG2," ",EG2," ",SG3," ",EG3
    write(cpuLocal,*),"lxGlo, lyGlo, lzGlo = ",lxGlo," ",lyGlo," ",lzGlo
    write(cpuLocal,*),neighCpu(1),neighCpu(2),neighCpu(3),neighCpu(4),neighCpu(5),neighCpu(6)
    write(cpuLocal,*),np,npLocal
    write(cpuLocal,*),"KnCoeff, ResistCoeff, dtDEM:", knCoeff,resistCoeff,dtDEM
    do i=1,np
        write(cpuLocal,*),part(i)%globalIdx,part(i)%x(1),part(i)%x(2),part(i)%x(3)
    enddo
    do i=1,6
        write(cpuLocal,*)"Sending parcels from side ",i
        do j=1,npSend(i)
            idx = partSend(j,i)%cellIdx
            write(cpuLocal,*)partSend(j,i)%globalIdx,map1Dto3D(1,idx),map1Dto3D(2,idx),map1Dto3D(3,idx),&
                partSend(j,i)%x(1),partSend(j,i)%x(2),partSend(j,i)%x(3)
        enddo
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine LogDEM(line,msg)
use vars
use ZjsVar
implicit none
    integer::line
    character(len=*)::msg
    open(cpuLocal,file=logFile,access="APPEND")
    write(cpuLocal,*)line,msg
    close(cpuLocal)
return
endsubroutine
!-----------------------------------------------  
subroutine OutputCloudDEM(trec)
use vars
use ZjsVar
implicit none
    character(len=*)::trec
    character(len=100)::filename
    character(len=10)::cpuStr
    integer::i,j

    write(cpuStr,"(I10)"), cpuLocal
    write(filename,*), trim(adjustl("./DEM/proc"//trim(adjustl(cpuStr))//"/"//trec//".dat"))
    open(cpuLocal,file=trim(adjustl(filename)))
    write(cpuLocal,*)np,npLocal
    do i=1,np
        write(cpuLocal,"(1I5)",advance="no")part(i)%pType
        write(cpuLocal,"(7E15.7)"),part(i)%wsum+1d-20,part(i)%x(1)+1d-20,part(i)%x(2)+1d-20,part(i)%x(3)+1d-20,&
            part(i)%v(1)+1d-20,part(i)%v(2)+1d-20,part(i)%v(3)+1d-20
        !do j=1,neighCellNo(i)
        !    write(cpuLocal,"(I5)",advance="no") neighCellIdx(j,i)
        !enddo
        !write(cpuLocal,*) " " 
        !do j=1,neighCellNo(i)
        !    write(cpuLocal,"(E15.7)",advance="no") neighCellWeight(j,i)
        !enddo
        !write(cpuLocal,*) " "
    enddo
    close(cpuLocal)
    if (cpuLocal==0)then
        write(*,*) "STEP    ",trec
    endif
return
endsubroutine
!-----------------------------------------------

!------------MBQ-------------20181227-------------
subroutine OutputCloudTecplot(trec)
use vars
use ZjsVar
implicit none
integer::i,j,k
character(len=*)::trec
character(len=3)::cpuStr
real*8::lpx,lpy,lpz,lpvel
    if (npLocal>0)then
    write(cpuStr,"(I3)") CPUIdx
    open(88,file=adjustl("Results/Particle_"//trim(adjustl(trec))//"_"//trim(adjustl(cpuStr))//".dat"))
    write(88,*) "variables=x,y,z,u,v,w,index,ugap,vgap,wgap"
    write(88,*) "zone i=",npLocal, "j=",1
    do i=1,npLocal       
        lpx = part(i)%fdrag(1)*dt
        lpy = part(i)%fdrag(2)*dt
        lpz = part(i)%fdrag(3)*dt
        lpvel = sqrt(lpx**2+lpy**2+lpz**2)
        write(88,*)  part(i)%x(1),part(i)%x(2),part(i)%x(3),part(i)%v(1),part(i)%v(2),part(i)%v(3),part(i)%MSIdx,part(i)%vgap(1),part(i)%vgap(2),part(i)%vgap(3)
    end do
    close(88)
    else
    write(cpuStr,"(I3)") CPUIdx
    open(88,file=adjustl("Results/Particle_"//trim(adjustl(trec))//"_"//trim(adjustl(cpuStr))//".dat"))
    write(88,*) "variables=x,y,z,u,v,w,index,ugap,vgap,wgap"
    write(88,*) "zone i=",1, "j=",1
    write(88,*) 0.5d0*(xStart+xEnd),0.5d0*(yStart+yEnd),0.5d0*(zStart+zEnd),0d0,0d0,0d0,0,0d0,0d0,0d0
    do i=1,npLocal       
        lpx = part(i)%fdrag(1)*dt
        lpy = part(i)%fdrag(2)*dt
        lpz = part(i)%fdrag(3)*dt
        lpvel = sqrt(lpx**2+lpy**2+lpz**2)
        write(88,*)  part(i)%x(1),part(i)%x(2),part(i)%x(3),part(i)%v(1),part(i)%v(2),part(i)%v(3),part(i)%MSIdx,part(i)%vgap(1),part(i)%vgap(2),part(i)%vgap(3)
    end do
    close(88)

    endif
    
return
endsubroutine
!---------------------------------------------------
