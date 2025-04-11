      subroutine druku(y,i0,j0,ly,ldeg,difmat,rr1,rr2,nrr0)
      implicit none
      integer*4 i0,j0,ly,ldeg,nrr0
      real*8 rr1,rr2
      real*8 y(i0,j0)
      external difmat
c
      integer*4 nrrmax
      parameter(nrrmax=1000000)
c
      integer*4 i,j,k,irr,nrr,nn
      real*8 rr,drr,epsilon
      real*8 yabs(6,3)
      real*8 y1(6,3),y2(6,3),yk(6,4),mat(6,6,3)
      save nn
c
      real*8 eps
      data eps/1.0d-04/
c
      do j=1,j0
        do i=1,i0
          y2(i,j)=y(i,j)
        enddo
      enddo
      nrr=max0(2,nrr0/4)
100   continue
      drr=(rr2-rr1)/dble(nrr)
      do j=1,j0
        do i=1,i0
          y1(i,j)=y2(i,j)
          y2(i,j)=y(i,j)
        enddo
      enddo
      do irr=1,nrr
        rr=rr1+dble(irr-1)*drr
        call difmat(ly,ldeg,rr,mat(1,1,1))
        call difmat(ly,ldeg,rr+0.5d0*drr,mat(1,1,2))
        call difmat(ly,ldeg,rr+drr,mat(1,1,3))
        do j=1,j0
          do i=1,i0
            yk(i,1)=0.d0
            do k=1,i0
              yk(i,1)=yk(i,1)+drr*mat(i,k,1)*y2(k,j)
            enddo
          enddo
c
          do i=1,i0
            yk(i,2)=0.d0
            do k=1,i0
              yk(i,2)=yk(i,2)+drr*mat(i,k,2)*(y2(k,j)+yk(k,1)/2.d0)
            enddo
          enddo
c
          do i=1,i0
            yk(i,3)=0.d0
            do k=1,i0
              yk(i,3)=yk(i,3)+drr*mat(i,k,2)*(y2(k,j)+yk(k,2)/2.d0)
            enddo
          enddo
c
          do i=1,i0
            yk(i,4)=0.d0
            do k=1,i0
              yk(i,4)=yk(i,4)+drr*mat(i,k,3)*(y2(k,j)+yk(k,3))
            enddo
          enddo
c
          do i=1,i0
            y2(i,j)=y2(i,j)+(yk(i,1)+2.d0*(yk(i,2)+yk(i,3))+yk(i,4))
     &                     /6.d0
          enddo
        enddo
      enddo
c
      do j=1,j0
        do i=1,i0
          yabs(i,j)=dmax1(dabs(y(i,j)),dabs(y2(i,j)))
        enddo
      enddo
c
      if(nrr.lt.nrrmax)then
        do j=1,j0
          do i=1,i0
            epsilon=eps*(dabs(y2(i,j)-y(i,j))+eps*yabs(i,j))
            if(dabs(y2(i,j)-y1(i,j)).gt.epsilon)then
              nrr=nrr*2
              goto 100
            endif
          enddo
        enddo
      else
        print *, ' Warning in druku: Convergence problem!'
      endif
c
      do j=1,j0
        do i=1,i0
          y(i,j)=y2(i,j)
        enddo
      enddo
      if(nn.lt.nrr)then
        nn=nrr
c        write(*,'(a,i4,a,i3,a,i6)')'  ruku: deg = ',ldeg,
c     &                             ', layer = ',ly,', steps = ',nn
      endif
c
      nrr0=nrr
c
      return
      end