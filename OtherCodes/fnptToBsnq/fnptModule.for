      module fnptModule
      implicit none

        integer,protected::fnptnz,fnptnt,fnptsti,fnptLast
        double precision,protected::fnptdt,fnptDep
        double precision,allocatable,protected::fnptData(:,:)
        logical::fnptEOF
        character(len=256)::fnptFile

      contains
        
        subroutine readFNPT
        !use mpi
        implicit none

          character(len=256)::text1,text2
          integer::i,tmpi,tmpi2
          integer::myrank,ierr
          double precision::tmpr(1:8),tmpr2
          logical::ex

          !call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
          myrank=0

          open(601+myrank,file='fnptSetup.txt')
          read(601+myrank,*)text1
          read(601+myrank,*)fnptFile
          read(601+myrank,*)text1
          read(601+myrank,*)fnptnt
          read(601+myrank,*)text1
          read(601+myrank,*)fnptdt
          read(601+myrank,*)text1
          read(601+myrank,*)tmpi
          read(601+myrank,*)text1
          read(601+myrank,*)fnptnz
          read(601+myrank,*)text1
          read(601+myrank,*)fnptDep
          close(601+myrank)          

          allocate(fnptData(fnptnz*(fnptnt+1),5))
          !! fnptData(:,1) = Z Coordinate
          !! fnptData(:,2) = Pressure
          !! fnptData(:,3) = U
          !! fnptData(:,4) = W
          !! fnptData(:,5) = Time (s)

          inquire(file=trim(fnptFile),exist=ex)
          if(ex) then
            open((601+myrank),file=trim(fnptFile))
          else
            write(*,*)"[ERR] Missing FNPT data file"
            stop
          endif          
          do i=1,(fnptnz*fnptnt)
            !read((601+myrank),*)tmpr,tmpr,fnptData(i,1:4),tmpr,tmpr
            read((601+myrank),*)tmpr(1:8)
            !fnptData(i,1:4)=tmpr(3:6)
            !fnptData(i,1)=tmpr(3)+fnptDep+myrank
            fnptData(i,1)=tmpr(3)+fnptDep
            fnptData(i,2)=tmpr(4)
            fnptData(i,3)=tmpr(5)
            fnptData(i,4)=tmpr(6)
            fnptData(i,5)=tmpr(1)
          enddo
          close((601+myrank))          

          fnptsti=1
          tmpi=max(1,tmpi)
          tmpi=min(fnptnt-1,tmpi)
          fnptsti=(tmpi-1)*fnptnz+1
          call etaFNPT(tmpr2)
          call findStartFNPT(tmpi2)

          if(myrank.eq.0)then
          !if(.true.)then
            write(*,*)
            write(*,'("[INF] ------- FNPT Info -------")')
            write(*,'("[---] ",a30,I15)')"Num of Time-Steps ",fnptnt
            write(*,'("[---] ",a30,F15.6)')"Time step (s) ",fnptdt
            write(*,'("[---] ",a30,I15)')"Num of vert nodes ",fnptnz
            write(*,'("[---] ",a30,F15.6)')"Water depth (m) ",fnptDep
            write(*,'("[---] ",a30,I15)')"Suggested Start time-step ",
     &        tmpi2
            write(*,'("[---] ",a30,I15)')"Start time-step ",tmpi
            write(*,'("[---] ",a30,F15.6)')"Start time (s) ",
     &        (tmpi-1)*fnptdt
            write(*,'("[---] ",a30,F15.6)')"Elevation at t0 (m) ",
     &        tmpr2
            if(tmpr2.gt.fnptDep)then
              write(*,'("[ERR] Note that elevation at t0 is > 0")')
              write(*,'("[---] Suggest starting at earlier timestep")')
            endif
            write(*,*)
            !write(*,*)myrank,fnptData(fnptnz,1:4)
            !stop
          endif

          fnptLast=fnptnt*fnptnz
          fnptEOF=.false.

        end subroutine readFNPT


        subroutine findStartFNPT(startI)
        implicit none

          integer::i,j
          integer,intent(out)::startI

          do i=1,fnptnt
            j=(i-1)*fnptnz+1
            if(abs(fnptData(j+fnptnz-1,1)-fnptDep).gt.0) exit
          enddo          
          startI=max(i-1,1)          

        end subroutine findStartFNPT


        subroutine updateFNPT
        !use mpi
        implicit none

          integer::myrank,ierr
          real(kind=8)::tmpr1,tmpr2,tmpr3
          
          !call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
          myrank=0
          
          if(.not.fnptEOF)then
            if(myrank.eq.0)then
              call etaFNPT(tmpr1)
              call velFNPT(0.15d0,tmpr2,1)
              call velFNPT(0.15d0,tmpr3,3)
              write(601,'(4E20.10)')(fnptsti+fnptnz-1)/fnptnz*fnptdt,
     &          tmpr1,tmpr2,tmpr3
            endif          
            
            fnptsti=fnptsti+fnptnz
            !fnptsti=fnptsti+myrank
C           if(myrank.eq.0) then
C             write(*,'("[INF] FNPT Start I = ",i10)')fnptsti
C           endif
            
            if(fnptsti.gt.fnptLast) fnptEOF=.true.          
          endif

        end subroutine updateFNPT


        subroutine etaFNPT(zMax)
        implicit none

          double precision,intent(out)::zMax

          if(fnptEOF)then
            zMax=fnptDep
            return
          endif
          
          zMax=fnptData(fnptsti+fnptnz-1,1)
            
        end subroutine etaFNPT


        subroutine velFNPT(zIn,velOut,velIndex)
        implicit none

          integer,intent(in)::velIndex
          double precision,intent(in)::zIn
          double precision,intent(out)::velOut

          integer::i,i2
          double precision:: dr,driav,tmpr

          if(fnptEOF)then
            velOut=0d0
            return
          endif

          do i2=0,fnptnz-1
            i=fnptsti+i2

            if(i2.eq.0)then
              driav=(fnptData(i+1,1)-fnptData(i,1))/2d0
            elseif(i2.eq.(fnptnz-1))then
              driav=(fnptData(i,1)-fnptData(i-1,1))/2d0
            else
              driav=(fnptData(i+1,1)-fnptData(i-1,1))/4d0
            endif

            dr=dabs(zIn-fnptData(i,1))
            if(dr.le.driav) exit
          enddo

          if(i.eq.fnptsti+fnptnz-1)then
            if((zIn-fnptData(i,1)).gt.0d0)then
              velOut=0d0
              return
            elseif(zIn.lt.fnptData(1,1))then
              i=fnptsti+1
            else
              i=i-1
            endif
          endif

          if(i.eq.fnptsti) i=i+1

          if(velIndex.eq.1) i2=3
          if(velIndex.eq.3) i2=4

          if(zIn.lt.fnptData(fnptsti,1))then
            tmpr=fnptData(fnptsti,1)+abs(zIn-fnptData(fnptsti,1))
            i=fnptsti+1
            velOut=fnptData(i,i2)
     &        +(fnptData(i+1,i2)-fnptData(i-1,i2))
     &        *(tmpr-fnptData(i,1))
     &        /(fnptData(i+1,1)-fnptData(i-1,1))
          else
            velOut=fnptData(i,i2)
     &        +(fnptData(i+1,i2)-fnptData(i-1,i2))
     &        *(zIn-fnptData(i,1))
     &        /(fnptData(i+1,1)-fnptData(i-1,1))
          endif

        end subroutine velFNPT

        subroutine depIntVelFNPT(tOut,etaOut,velOut)
        implicit none
                    
          double precision,intent(out)::velOut,tOut,etaOut

          integer::n(3),i,i2
          double precision:: dz,h,w(3)

          tOut=0d0          
          etaOut=0d0
          velOut=0d0          

          if(fnptEOF)then
            return
          endif

          call etaFNPT(h)
          etaOut=h-fnptDep
          dz=h/(fnptnz-1)
          
          w=(/ 1d0, 4d0, 1d0 /)                            
          do i2=0,fnptnz-3,2
            n(1)=fnptsti+i2
            n(2)=n(1)+1
            n(3)=n(1)+2

            do i=1,3
              velOut=velOut+(w(i)*fnptData(n(i),3))
            enddo            
          enddo  
          velOut=dz/3d0*velOut

          tOut=fnptData(fnptsti,5)

        end subroutine depIntVelFNPT

      end module fnptModule