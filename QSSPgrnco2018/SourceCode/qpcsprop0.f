      subroutine qpcsprop0(lyup,lylw)
      use qpcalloc
      implicit none
c
c     calculation of spheroidal response for ldeg = 0
c
      integer*4 lyup,lylw
c
c     work space
c
      integer*4 i,j,istp,ly,iy,iys,iz,key
      real*8 fac,ca,cb
      real*8 maiup(2,2),mailw(2,2)
      real*8 wave(2),coef2(2,2),b2(2,2),c(3)
      external qpcmai0
c
      do iz=1,nz
        do istp=1,4
          do i=1,6
            ypsv(i,istp,iz)=0.d0
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
        do i=1,3
          yup3(i,lyup)=mas3x3(i,2,lyup)
        enddo
      endif
c
      do ly=lyup,lys2-1
        call axb(mas3x3inv(1,1,ly),yup3(1,ly),3,3,1,c)
        wave(1)=(rrlw(ly)/rrup(ly))**2
        wave(2)=rrlw(ly)/rrup(ly)
        do iy=lyup,ly
          do i=1,3
            yup3(i,iy)=yup3(i,iy)*wave(2)
          enddo
        enddo
c
	  c(1)=c(1)*wave(2)*wave(1)
c       c(2)=c(2)
        c(3)=c(3)*wave(2)
c
        call axb(mas3x3(1,1,ly),c,3,3,1,yup3(1,ly+1))
      enddo
c
c===============================================================================
c
c     propagation from bottom to source
c
      do i=1,3
        ylw3(i,lylw)=mas3x3(i,1,lylw)
      enddo
c
      do ly=lylw-1,lys1,-1
        call axb(mas3x3inv(1,1,ly),ylw3(1,ly+1),3,3,1,c)
        wave(1)=(rrlw(ly)/rrup(ly))**2
        wave(2)=rrlw(ly)/rrup(ly)
c
        do iy=lylw,ly+1,-1
          do i=1,3
            ylw3(i,iy)=ylw3(i,iy)*wave(1)
          enddo
        enddo
c
c       c(1)=c(1)
        c(2)=c(2)*wave(1)*wave(2)
        c(3)=c(3)*wave(1)
c
        call axb(mas3x3(1,1,ly),c,3,3,1,ylw3(1,ly))
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
                ypsv(i,istp,iz)=ypsv(i,istp,iz)+b2(1,istp)*yup3(i,ly)
              enddo
c              ypsv(5,istp,iz)=b2(1,istp)*yup3(3,ly)
c              ypsv(6,istp,iz)=ypsv(5,istp,iz)/rrup(ly)
            enddo
          endif
        enddo
        do ly=lys1+1,lylw
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,2
              do i=1,2
                ypsv(i,istp,iz)=ypsv(i,istp,iz)+b2(2,istp)*ylw3(i,ly)
              enddo
c              ypsv(5,istp,iz)=b2(2,istp)*ylw3(3,ly)
c              ypsv(6,istp,iz)=ypsv(5,istp,iz)/rrup(ly)
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
          call sruku(b2,2,2,iys,0,qpcmai0,
     &               fac,rrup(lys),gdds0,rrup(iys),rrlw(iys))
c
          do i=1,2
            coef2(i,1)=yup3(i,iys+1)
            coef2(i,2)=-ylw3(i,iys+1)
          enddo
          key=0
          call svd500(coef2,b2,2,2,0.d0,key)
          if(key.eq.0)then
            print *,' Warning in qpcsprop0: anormal exit from svd500!'
            return
          endif
c
          do ly=lyup,iys
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,2
                do i=1,2
                  ypsv(i,istp,iz)=ypsv(i,istp,iz)+b2(1,istp)*yup3(i,ly)
                enddo
c                ypsv(5,istp,iz)=ypsv(5,istp,iz)+b2(1,istp)*yup3(3,ly)
c                ypsv(6,istp,iz)=ypsv(5,istp,iz)/rrup(ly)
              enddo
            endif
          enddo
          do ly=iys+1,lylw
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,2
                do i=1,2
                  ypsv(i,istp,iz)=ypsv(i,istp,iz)+b2(2,istp)*ylw3(i,ly)
                enddo
c                ypsv(5,istp,iz)=ypsv(5,istp,iz)+b2(2,istp)*ylw3(3,ly)
c                ypsv(6,istp,iz)=ypsv(5,istp,iz)/rrup(ly)
              enddo
            endif
          enddo
        enddo
      endif
c
      return
      end
