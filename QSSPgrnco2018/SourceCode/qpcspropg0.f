      subroutine qpcspropg0(lyup,lylw)
      use qpcalloc
      implicit none
c
c     calculation of spheroidal response for ldeg = 0
c
      integer*4 lyup,lylw
c
c     work space
c
      integer*4 i,j,istp,ly,iy,iys,iz,ily,nly,key
      real*8 f,rr1,rr2,dlnr,h,fac,ca,cb
      real*8 maiup(2,2),mailw(2,2),swap6x6(6,6)
      real*8 wave(2),coef2(2,2),b2(2,2),swap(3)
      external qpcdifmat0,qpcmaig0
c
      do iz=1,nz
        do istp=1,4
          do i=1,6
            ypsvg(i,istp,iz)=0.d0
          enddo
        enddo
      enddo
c
      if(depr(1).lt.0.d0.and.depr(nz).lt.0.d0)return
c
      do ly=1,ly0
        do i=1,3
          yup3(i,ly)=0.d0
          ylw3(i,ly)=0.d0
        enddo
      enddo
c
c===============================================================================
c
c     propagation from surface to source
c
      if(lyup.eq.1)then
        yup3(1,lyup)=1.d0
        yup3(2,lyup)=0.d0
        yup3(3,lyup)=0.d0
      else
        call qpcsmat0(lyup,lyup)
        do i=1,3
          yup3(i,lyup)=mas3x3(i,2,lyup)
        enddo
      endif
c
      do ly=lyup,lys2-1
        do i=1,3
          swap(i)=yup3(i,ly)
        enddo
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h/rrlw(ly))
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
          call druku(swap,3,1,ly,0,qpcdifmat0,rr1,rr2,ndruku(0,ly))
        enddo
c
        do i=1,3
          yup3(i,ly+1)=swap(i)
        enddo
      enddo
c
c===============================================================================
c
c     propagation from bottom to source
c
      call qpcsmat0(lylw,lylw)
      do i=1,3
        ylw3(i,lylw)=mas3x3(i,1,lylw)
      enddo
c
      do ly=lylw-1,lys1,-1
        do i=1,3
          swap(i)=ylw3(i,ly+1)
        enddo
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h/rrlw(ly))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
          call druku(swap,3,1,ly,0,qpcdifmat0,rr1,rr2,ndruku(0,ly))
        enddo
c
        do i=1,3
          ylw3(i,ly)=swap(i)
        enddo
      enddo
c
      do ly=lyup,lys2
        yup3(1,ly)=yup3(1,ly)/rrup(ly)
        yup3(2,ly)=yup3(2,ly)/rrup(ly)**2
c       yup3(3,ly)=yup3(3,ly)
      enddo
c
      do ly=lys1,lylw
        ylw3(1,ly)=ylw3(1,ly)/rrup(ly)
        ylw3(2,ly)=ylw3(2,ly)/rrup(ly)**2
c       ylw3(3,ly)=ylw3(3,ly)
      enddo
c
      if(lys1.eq.lys2)then
c
c       for point source
c
        b2(1,1)=1.d0
        b2(2,1)=0.d0
        b2(1,2)=0.d0
        b2(2,2)=1.d0
        do i=1,2
          coef2(i,1)=yup3(i,lys1)
          coef2(i,2)=-ylw3(i,lys1)
        enddo
        key=0
        call svd500(coef2,b2,2,2,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpcprop0: anormal exit from svd500!'
          return
        endif
        do ly=lyup,lys1
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,2
              do i=1,2
                ypsvg(i,istp,iz)=ypsvg(i,istp,iz)+b2(1,istp)*yup3(i,ly)
              enddo
c              ypsvg(5,istp,iz)=b2(1,istp)*yup3(3,ly)
c              ypsvg(6,istp,iz)=ypsvg(5,istp,iz)/rrup(ly)
            enddo
          endif
        enddo
        do ly=lys1+1,lylw
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,2
              do i=1,2
                ypsvg(i,istp,iz)=ypsvg(i,istp,iz)+b2(2,istp)*ylw3(i,ly)
              enddo
c              ypsvg(5,istp,iz)=b2(2,istp)*ylw3(3,ly)
c              ypsvg(6,istp,iz)=ypsvg(5,istp,iz)/rrup(ly)
            enddo
          endif
        enddo
      else
c
c       for line source
c
        fac=-rrup(lys)**2/gddsnorm
c
        do iys=lys1,lys2-1
          do istp=1,2
            do i=1,2
              b2(i,istp)=0.d0
            enddo
          enddo
          call sruku(b2,2,2,iys,0,qpcmaig0,
     &               fac,rrup(lys),gdds0,rrup(iys),rrlw(iys))
c
          do i=1,2
            coef2(i,1)=yup3(i,iys+1)
            coef2(i,2)=-ylw3(i,iys+1)
          enddo
          key=0
          call svd500(coef2,b2,2,2,0.d0,key)
          if(key.eq.0)then
            print *,' Warning in qpcspropg0: anormal exit from svd500!'
            return
          endif
c
          do ly=lyup,iys
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,2
                do i=1,2
                  ypsvg(i,istp,iz)=ypsvg(i,istp,iz)
     &                            +b2(1,istp)*yup3(i,ly)
                enddo
c                ypsvg(5,istp,iz)=ypsvg(5,istp,iz)+b2(1,istp)*yup3(3,ly)
c                ypsvg(6,istp,iz)=ypsvg(5,istp,iz)/rrup(ly)
              enddo
            endif
          enddo
          do ly=iys+1,lylw
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,2
                do i=1,2
                  ypsvg(i,istp,iz)=ypsvg(i,istp,iz)
     &                            +b2(2,istp)*ylw3(i,ly)
                enddo
c                ypsvg(5,istp,iz)=ypsvg(5,istp,iz)+b2(2,istp)*ylw3(3,ly)
c                ypsvg(6,istp,iz)=ypsvg(5,istp,iz)/rrup(ly)
              enddo
            endif
          enddo
        enddo
      endif
      return
      end