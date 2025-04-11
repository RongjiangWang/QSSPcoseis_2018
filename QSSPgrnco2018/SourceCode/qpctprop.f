      subroutine qpctprop(ldeg,lyup,lylw)
      use qpcalloc
      implicit none
c
c
c     calculation of toroidal response
c
      integer*4 ldeg,lyup,lylw
c
c     work space
c
      integer*4 i,j,istp,iy,ly,iys,iz,key
      real*8 fac,ca,cb
      real*8 maiup(2,2),mailw(2,2),swap6x6(6,6)
      real*8 wave(2),coef2(2,2),b2(2,2),c(2)
      external qpcmaish
c
      do iz=1,nz
        do istp=1,2
          do i=1,2
            ysh(i,istp,iz)=0.d0
          enddo
        enddo
      enddo
c
      if(depr(1).lt.0.d0.and.depr(nz).lt.0.d0)return
c
      do ly=1,ly0
        do i=1,2
          yup2(i,ly)=0.d0
          ylw2(i,ly)=0.d0
        enddo
      enddo
c===============================================================================
c
c     propagation from surface to source
c
      if(lyup.eq.1)then
        yup2(1,lyup)=1.d0
        yup2(2,lyup)=0.d0
      else
        yup2(1,lyup)=mat2x2(1,2,lyup)
        yup2(2,lyup)=mat2x2(2,2,lyup)
      endif
c
      do ly=lyup,lys2-1
        call axb(mat2x2inv(1,1,ly),yup2(1,ly),2,2,1,c)
        wave(1)=(rrlw(ly)/rrup(ly))**ldeg
        wave(2)=(rrlw(ly)/rrup(ly))**(ldeg+1)
        do iy=lyup,ly
          do i=1,2
            yup2(i,iy)=yup2(i,iy)*wave(2)
          enddo
        enddo
c
	  c(1)=c(1)*wave(1)*wave(2)
c       c(2)=c(2)
c
        call axb(mat2x2(1,1,ly),c,2,2,1,yup2(1,ly+1))
      enddo
c
c===============================================================================
c
c     propagation from bottom to source
c
      if(lylw.eq.lycm)then
        ylw2(1,lylw)=1.d0
        ylw2(2,lylw)=0.d0
      else
        ylw2(1,lylw)=mat2x2(1,1,lylw)
        ylw2(2,lylw)=mat2x2(2,1,lylw)
      endif
c
      do ly=lylw-1,lys1,-1
        call axb(mat2x2inv(1,1,ly),ylw2(1,ly+1),2,2,1,c)
        wave(1)=(rrlw(ly)/rrup(ly))**ldeg
        wave(2)=(rrlw(ly)/rrup(ly))**(ldeg+1)
c
        do iy=lylw,ly+1,-1
          do i=1,2
            ylw2(i,iy)=ylw2(i,iy)*wave(1)
          enddo
        enddo
c
c       c(1)=c(1)
        c(2)=c(2)*wave(1)*wave(2)
c
        call axb(mat2x2(1,1,ly),c,2,2,1,ylw2(1,ly))
      enddo
c
      do ly=lyup,lys2
        yup2(1,ly)=yup2(1,ly)/rrup(ly)
        yup2(2,ly)=yup2(2,ly)/rrup(ly)**2
      enddo
c
      do ly=lys1,lylw
        ylw2(1,ly)=ylw2(1,ly)/rrup(ly)
        ylw2(2,ly)=ylw2(2,ly)/rrup(ly)**2
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
          coef2(i,1)=yup2(i,lys1)
          coef2(i,2)=-ylw2(i,lys1)
        enddo
        key=0
        call svd500(coef2,b2,2,2,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpctprop: anormal exit from svd500!'
          return
        endif
        do ly=lyup,lys1
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,2
              do i=1,2
                ysh(i,istp,iz)=ysh(i,istp,iz)+b2(1,istp)*yup2(i,ly)
              enddo
            enddo
          endif
        enddo
        do ly=lys1+1,lylw
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,2
              do i=1,2
                ysh(i,istp,iz)=ysh(i,istp,iz)+b2(2,istp)*ylw2(i,ly)
              enddo
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
          call sruku(b2,2,2,iys,ldeg,qpcmaish,
     &               fac,rrup(lys),gdds0,rrup(iys),rrlw(iys))
          do i=1,2
            coef2(i,1)=yup2(i,iys+1)
            coef2(i,2)=-ylw2(i,iys+1)
          enddo
          key=0
          call svd500(coef2,b2,2,2,0.d0,key)
          if(key.eq.0)then
            print *,' Warning in qpctprop: anormal exit from svd500!'
            return
          endif
c
          do ly=lyup,iys
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,2
                do i=1,2
                  ysh(i,istp,iz)=ysh(i,istp,iz)+b2(1,istp)*yup2(i,ly)
                enddo
              enddo
            endif
          enddo
          do ly=iys+1,lylw
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,2
                do i=1,2
                  ysh(i,istp,iz)=ysh(i,istp,iz)+b2(2,istp)*ylw2(i,ly)
                enddo
              enddo
            endif
          enddo
        enddo
      endif
c
      return
      end