      subroutine sruku(y,i0,j0,ly,ldeg,difmat,fac,rr0,gdds,rr1,rr2)
      implicit none
      integer*4 i0,j0,ly,ldeg
      real*8 fac,rr0,gdds,rr1,rr2
      real*8 y(i0,j0)
      external difmat
c
      integer*4 nrrmax
      parameter(nrrmax=131072)
c
      integer*4 i,j,k,irr,nrr
      real*8 drr,epsilon
      real*8 x0,x1,x2,s0,s1,s2,rr0up,rr0lw
      real*8 yabs(6,4)
      real*8 y1(6,4),y2(6,4),yk(6,4),mat(6,6,4)
c
      real*8 eps
      data eps/1.0d-04/
c
      rr0up=rr0+0.5d0*gdds
      rr0lw=rr0-0.5d0*gdds
c
      do j=1,j0
        do i=1,i0
          y2(i,j)=y(i,j)
        enddo
      enddo
      nrr=16
100   continue
      drr=(rr2-rr1)/dble(nrr)
      do j=1,j0
        do i=1,i0
          y1(i,j)=y2(i,j)
          y2(i,j)=y(i,j)
        enddo
      enddo
      do irr=1,nrr
        x0=rr1+dble(irr-1)*drr
        x1=x0+0.5d0*drr
        x2=x0+drr
        call difmat(ly,ldeg,x0,mat(1,1,1))
        call difmat(ly,ldeg,x1,mat(1,1,2))
        call difmat(ly,ldeg,x2,mat(1,1,3))
        if(x0.ge.rr0)then
          s0=fac*((rr0up-x0)/x0)**2
          s1=fac*((rr0up-x1)/x1)**2
          s2=fac*((rr0up-x2)/x2)**2
        else
          s0=fac*((x0-rr0lw)/x0)**2
          s1=fac*((x1-rr0lw)/x1)**2
          s2=fac*((x2-rr0lw)/x2)**2
        endif
        do j=1,j0
          do i=1,i0
            yk(i,1)=0.d0
            do k=1,i0
              yk(i,1)=yk(i,1)+mat(i,k,1)*y2(k,j)
            enddo
          enddo
          yk(j,1)=yk(j,1)+s0
c
          do i=1,i0
            yk(i,2)=0.d0
            do k=1,i0
              yk(i,2)=yk(i,2)+mat(i,k,2)*(y2(k,j)+yk(k,1)*drr/2.d0)
            enddo
          enddo
          yk(j,2)=yk(j,2)+s1
c
          do i=1,i0
            yk(i,3)=0.d0
            do k=1,i0
              yk(i,3)=yk(i,3)+mat(i,k,2)*(y2(k,j)+yk(k,2)*drr/2.d0)
            enddo
          enddo
          yk(j,3)=yk(j,3)+s1
c
          do i=1,i0
            yk(i,4)=0.d0
            do k=1,i0
              yk(i,4)=yk(i,4)+mat(i,k,3)*(y2(k,j)+drr*yk(k,3))
            enddo
          enddo
          yk(j,4)=yk(j,4)+s2
c
          do i=1,i0
            y2(i,j)=y2(i,j)
     &             +(yk(i,1)+2.d0*(yk(i,2)+yk(i,3))+yk(i,4))*drr/6.d0
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
        do j=1,j0
          do i=1,i0
            epsilon=dsqrt(eps)*(dabs(y2(i,j)-y(i,j))
     &                         +dsqrt(eps)*yabs(i,j))
            if(dabs(y2(i,j)-y1(i,j)).gt.epsilon)then
              print *, ' Warning in sruku: Convergence problem!'
              goto 200
            endif
          enddo
        enddo
      endif
c
200   continue
c
      do j=1,j0
        do i=1,i0
          y(i,j)=y2(i,j)
        enddo
      enddo
c
      return
      end