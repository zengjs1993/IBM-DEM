!----------------------------------------------------------
subroutine OutputGridAndPartitionInformation()
use vars
implicit none
    character(len=100)::filename
    write(filename,*) cpuIdx,".partition"
    open(cpuIdx,file="./DEM/log/"//trim(adjustl(filename)),access="append")
    write(cpuIdx,*) "Global CPU regions: ",CPUnx,CPUny,CPUnz
    write(cpuIdx,*) "Local CPU Index: ",CPUIdx,CPUIdxX,CPUIdxY,CPUIdxZ
    write(cpuIdx,*) "Neighbor CPU Index: ",CPUNeighIdx(1),CPUNeighIdx(2),CPUNeighIdx(3),CPUNeighIdx(4),CPUNeighIdx(5),CPUNeighIdx(6)
    write(cpuIdx,*) "Local Grid numbers: ",nx,ny,nz
    write(cpuIdx,*) "Local Grid intervals: ",sg1,eg1,sg2,eg2,sg3,eg3
    write(cpuIdx,*) "Global Grid intervals: ",lg1,ug1,lg2,ug2,lg3,ug3
    write(cpuIdx,*) "Boundary type:", bw,be,bb,bf,bs,bn
    close(cpuIdx)
return
end subroutine
!----------------------------------------------------------
subroutine OutputFluidVariables(trec)
use vars
implicit none
    character(len=100)::filename
    character(len=*)::trec
    character(len=10)::cpuStr
    integer::i,j,k
    
    write(cpuStr,"(I3)"), cpuIdx
    write(filename,*), trim(adjustl("./CFD/proc"//trim(adjustl(cpuStr))//"/"//trim(adjustl(trec))//".dat"))
    
    open(cpuIdx,file=trim(adjustl(filename)))
    do k=lg3,ug3
    do j=lg2,ug2
    do i=lg1,ug1
        write(cpuIdx,"(3I3)",advance="no") i,j,k
        write(cpuIdx,"(5E16.7)",advance="no")  1d-20+u(i,j,k),1d-20+v(i,j,k),1d-20+w(i,j,k),1d-20+p(i,j,k),1d-40+velDiv(i,j,k)
        write(cpuIdx,*) " "
    enddo
    enddo
    enddo
return
end subroutine
!----------------------------------------------------------
subroutine OutputFluidTecplot(cnt)
use vars
implicit none
    character(len=100)::filename
    character(len=10)::trec
    character(len=10)::cpuStr
    integer::i,j,k
    integer::cnt
    real*8::vel,omega
    if(cpuIdx==0) write(*,*) "time step = ",cnt
    write(trec,"(I10)"),cnt
    write(cpuStr,"(I3)"), cpuIdx
    write(filename,*), trim(adjustl("./CFD/tecplot/"//trim(adjustl(trec))//"_"//trim(adjustl(cpuStr))//".dat"))
    
    open(cpuIdx,file=trim(adjustl(filename)))
    !write(cpuIdx,*) "variables=x,y,z,u,v,w,omega,vel,P,fx,fy,fz,velDiv,Pcorr,errorU,errorV,errorP"
    write(cpuIdx,*) "variables=x,y,z,u,v,w,vel,P"
    write(cpuIdx,*) "zone i=",ug1-lg1+1,"j=",ug2-lg2+1,"k=",ug3-lg3+1
    !write(cpuIdx,*) "zone i=",eg1-sg1+1,"j=",eg2-sg2+1,"k=",eg3-sg3+1
    do k=lg3,ug3
    do j=lg2,ug2
    do i=lg1,ug1

    !do k=sg3,eg3
    !do j=sg2,eg2
    !do i=sg1,eg1
        vel = sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
        !write(cpuIdx,*)xc(i,j,k),yc(i,j,k),zc(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),omegaZ(i,j,k),vel,P(i,j,k), &
        !    fx(i,j,k),fy(i,j,k),fz(i,j,k),velDiv(i,j,k),Pcorr(i,j,k),errorU(i,j,k),errorV(i,j,k),errorP(i,j,k)
        write(cpuIdx,*)xc(i,j,k),yc(i,j,k),zc(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),vel,P(i,j,k)
    enddo
    enddo
    enddo
return
end subroutine
!----------------------------------------------------------
