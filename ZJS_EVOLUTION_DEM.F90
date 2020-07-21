!-----------------------------------------------
subroutine CloudEvolutionDEM()
use vars
use ZjsVar
    integer::ierr

    call RebuildAllListOrNot()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    !call ShiftDEM()

    !call Integration()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    call AdjustParcelPosition()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    call UpdateLocalCellWeightListDEM()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    call UpdateParcelInfoInGhostRegion()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    call UpdateGhostCellWeightListDEM()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
    call AssembleForceMS()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   

return
endsubroutine
!-----------------------------------------------
subroutine LogVolumeFraction(timeInt)
use vars
use ZjsVar

return
endsubroutine
!-----------------------------------------------
subroutine RebuildAllListOrNot()
use vars
use ZjsVar
    call RebuildParcelInCellListDEM()
    call RebuildParcelCollisionListDEM()
    call RebuildNeighborCellListForParcelsDEM()
return
endsubroutine
!-----------------------------------------------
subroutine Integration()
use vars
use ZjsVar
implicit none
    integer::i,j,idx1,idx2,idx3,cellNoTmp,cellIdxTmp
    real*8::mi,vel1,vel2,vel3,f1,f2,f3
    real*8::weight,rhofLocal,alphapLocal,rhofTmp,rhopi
    real*8,dimension(1:3)::ff,beta,cellVelLocal,vpCellLocal,dragi

    dtDEM=dt
    !dtDEM = 2d-5
    call CalculateCollisionForce()
    call CalculateSolidBoundaryCollisionForce()
    do i=1,npLocal
        f1 = part(i)%fcol(1)
        f2 = part(i)%fcol(2)
        f3 = part(i)%fcol(3)
        mi = part(i)%m
        vel1 = dtDEM*(f1/mi)!
        vel2 = dtDEM*(f2/mi)!
        vel3 = dtDEM*(f3/mi)!

        part(i)%v(1) = vel1
        part(i)%v(2) = vel2
        part(i)%v(3) = vel3

        part(i)%x(1) = part(i)%x(1) + dtDEM*vel1
        part(i)%x(2) = part(i)%x(2) + dtDEM*vel2
        part(i)%x(3) = part(i)%x(3) + dtDEM*vel3

    enddo
endsubroutine
!-----------------------------------------------
subroutine ShiftDEM()
use vars
use ZjsVar
implicit none
  integer::i
  real*8::rnd1,rnd2,rnd3
  do i=1,npLocal
    call random_number(rnd1)
    call random_number(rnd2)
    call random_number(rnd3)
    !part(i)%v(1) =  (part(i)%x(2)-0.5d0)*0.5d0/0.25d0
    !part(i)%v(2) = -(part(i)%x(1)-0.5d0)*0.5d0/0.25d0
    part(i)%v(1) = 0.5d0
    part(i)%x(1) = part(i)%x(1) + 0.5d0*dt!dtDEM*2.5d-1*(rnd1-0.5)!2.5d-1*di!
  enddo
return
endsubroutine
!-----------------------------------------------
subroutine CloudEvolution()
use vars
use ZjsVar
implicit none
    integer::i,j,ierr
    character(len=10)::idx
    dtDEM = 1d-5
    do i=1,100
        timeCntDEM = timeCntDEM+1
        timeDEM = timeDEM+dtDEM
        call RebuildParcelInCellListDEM()
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call RebuildParcelCollisionListDEM()
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call ShiftDEM()
        !call Integration()
        call AdjustParcelPosition()
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call UpdateParcelInfoInGhostRegion()
        if (mod(i,1)==0)then
            write(idx,"(I10)") i
            call OutputCloudDEM(trim(adjustl(idx)))
        endif
    enddo
return
endsubroutine
!-----------------------------------------------
subroutine AdjustParcelPosition()
use vars
use ZjsVar
implicit none
    integer::i

    errorX = di*1d-20; errorY = dj*1d-20; errorZ = dk*1d-20
    do i=1,npLocal
        if (abs(part(i)%x(1)-xStart)<errorX)then
            part(i)%x(1) = xStart+errorX
        endif
        if (abs(part(i)%x(1)-xEnd)<errorX)then
            part(i)%x(1) = xEnd+errorX
        endif
        if (abs(part(i)%x(2)-yStart)<errorY)then
            part(i)%x(2) = yStart+errorY
        endif
        if (abs(part(i)%x(2)-yEnd)<errorY)then
            part(i)%x(2) = yEnd+errorY
        endif
        if (abs(part(i)%x(3)-zStart)<errorZ)then
            part(i)%x(3) = zStart+errorZ
        endif
        if (abs(part(i)%x(3)-zEnd)<errorZ)then
            part(i)%x(3) = zEnd+errorZ
        endif

    enddo
return
endsubroutine
!-----------------------------------------------
subroutine ProcessDEM()
use vars
use ZjsVar
implicit none
    character(len=10)::idx
    call InitializationDEM()
!    write(idx,"(I10)")0
!    call OutputCloudDEM(trim(adjustl(idx)))
    call CloudEvolution()
!    call WriteDEM()
return
endsubroutine
!-----------------------------------------------  
