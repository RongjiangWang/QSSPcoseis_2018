      subroutine qpcsprop(ldeg,lyup,lylw)
      use qpcalloc
      implicit none
c
c     calculation of spheroidal response
c
      integer*4 ldeg,lyup,lylw
c
c     work space
c
      integer*4 i,j,j0,istp,ly,iy,iys,iz,key
      real*8 y4max
      real*8 cdet,cyswap,ca,cb,fac,r0
      real*8 c0(6,3),c1(6,3),cc0(4,2),cc1(4,2)
      real*8 y1(6,3),yc(6,3),ypot(4)
      real*8 yupc(4,2),ylwc(4,2)
      real*8 wave(6),orth(3,3),orthc(2,2)
      real*8 coef6(6,6),b6(6,4)
      external qpcmaipsv
c
      do iz=1,nz
        do istp=1,4
          do i=1,6
            ypsv(i,istp,iz)=0.d0
          enddo
        enddo
      enddo
c
      do ly=1,ly0
        do j=1,3
          do i=1,6
            yup6(i,j,ly)=0.d0
            ylw6(i,j,ly)=0.d0
          enddo
        enddo
      enddo
c
c===============================================================================
c
c     propagation from free surface to source
c
      if(ldeg.le.1)return
c
      if(lyup.eq.1)then
        do j=1,3
          do i=1,6
            yup6(i,j,lyup)=0.d0
          enddo
        enddo
        yup6(1,1,lyup)=1.d0
        yup6(3,2,lyup)=1.d0
        if(ldeg.eq.1)then
          yup6(1,3,lyup)=1.d0
          yup6(1,3,lyup)=1.d0
        endif
        yup6(5,3,lyup)=-grup(1)
      else if(lyup.lt.lycm)then
        do j=1,3
          do i=1,6
            yup6(i,j,lyup)=mas6x6(i,2*j,lyup)
          enddo
        enddo
      else
        stop ' Errot in qpssprop: wrong source depth!'
      endif
c
      do ly=lyup,lys2-1
        wave(1)=(rrlw(ly)/rrup(ly))**ldeg
        wave(2)=(rrlw(ly)/rrup(ly))**(ldeg+1)
        wave(3)=(rrlw(ly)/rrup(ly))**(ldeg+2)
        wave(4)=(rrlw(ly)/rrup(ly))**(ldeg-1)
        wave(5)=(rrlw(ly)/rrup(ly))**ldeg
        wave(6)=(rrlw(ly)/rrup(ly))**(ldeg+1)
c
        call axb(mas6x6inv(1,1,ly),yup6(1,1,ly),6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(2,1)*c0(4,2)*c0(6,3)
     &      +c0(4,1)*c0(6,2)*c0(2,3)
     &      +c0(6,1)*c0(2,2)*c0(4,3)
     &      -c0(6,1)*c0(4,2)*c0(2,3)
     &      -c0(4,1)*c0(2,2)*c0(6,3)
     &      -c0(2,1)*c0(6,2)*c0(4,3)
        orth(1,1)=(c0(4,2)*c0(6,3)-c0(4,3)*c0(6,2))/cdet
        orth(2,1)=(c0(4,3)*c0(6,1)-c0(4,1)*c0(6,3))/cdet
        orth(3,1)=(c0(4,1)*c0(6,2)-c0(4,2)*c0(6,1))/cdet
        orth(1,2)=(c0(2,3)*c0(6,2)-c0(2,2)*c0(6,3))/cdet
        orth(2,2)=(c0(2,1)*c0(6,3)-c0(2,3)*c0(6,1))/cdet
        orth(3,2)=(c0(2,2)*c0(6,1)-c0(2,1)*c0(6,2))/cdet
        orth(1,3)=(c0(2,2)*c0(4,3)-c0(2,3)*c0(4,2))/cdet
        orth(2,3)=(c0(2,3)*c0(4,1)-c0(2,1)*c0(4,3))/cdet
        orth(3,3)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))/cdet
c
        call axb(c0,orth,6,3,3,c1)
c
        do iy=lyup,ly
c
c         orthonormalization of the receiver vectors
c
          call axb(yup6(1,1,iy),orth,6,3,3,y1)
          call memcpy(y1,yup6(1,1,iy),18)
          do j=1,3
            do i=1,6
              yup6(i,j,iy)=yup6(i,j,iy)*wave(2*j)
            enddo
          enddo
        enddo
c
        c1(1,1)=c1(1,1)*wave(1)*wave(2)
        c1(2,1)=1.d0
        c1(3,1)=c1(3,1)*wave(3)*wave(2)
        c1(4,1)=0.d0
        c1(5,1)=c1(5,1)*wave(5)*wave(2)
        c1(6,1)=0.d0
c
        c1(1,2)=c1(1,2)*wave(1)*wave(4)
        c1(2,2)=0.d0
        c1(3,2)=c1(3,2)*wave(3)*wave(4)
        c1(4,2)=1.d0
        c1(5,2)=c1(5,2)*wave(5)*wave(4)
        c1(6,2)=0.d0
c
        c1(1,3)=c1(1,3)*wave(1)*wave(6)
        c1(2,3)=0.d0
        c1(3,3)=c1(3,3)*wave(3)*wave(6)
        c1(4,3)=0.d0
        c1(5,3)=c1(5,3)*wave(5)*wave(6)
        c1(6,3)=1.d0
c
        call axb(mas6x6(1,1,ly),c1,6,6,3,yup6(1,1,ly+1))
      enddo
c
c===============================================================================
c
c     propagation within inner core
c
      if(lylw.ge.lycc)then
c
c       lowest layer is within inner core
c
        do j=1,3
          do i=1,6
            ylw6(i,j,lylw)=mas6x6(i,2*j-1,lylw)
          enddo
        enddo
      endif
c
      do ly=lylw-1,lycc,-1
        wave(1)=(rrlw(ly)/rrup(ly))**ldeg
        wave(2)=(rrlw(ly)/rrup(ly))**(ldeg+1)
        wave(3)=(rrlw(ly)/rrup(ly))**(ldeg+2)
        wave(4)=(rrlw(ly)/rrup(ly))**(ldeg-1)
        wave(5)=(rrlw(ly)/rrup(ly))**ldeg
        wave(6)=(rrlw(ly)/rrup(ly))**(ldeg+1)
c
        call axb(mas6x6inv(1,1,ly),ylw6(1,1,ly+1),6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call axb(c0,orth,6,3,3,c1)
c
        do iy=lylw,ly+1,-1
c
c         orthonormalization of the receiver vectors
c
          call axb(ylw6(1,1,iy),orth,6,3,3,y1)
          call memcpy(y1,ylw6(1,1,iy),18)
          do j=1,3
            do i=1,6
              ylw6(i,j,iy)=ylw6(i,j,iy)*wave(2*j-1)
            enddo
          enddo
        enddo
c
        c1(1,1)=1.d0
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=0.d0
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=0.d0
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=0.d0
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=1.d0
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=0.d0
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=0.d0
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=0.d0
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=1.d0
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call axb(mas6x6(1,1,ly),c1,6,6,3,ylw6(1,1,ly))
      enddo
c
c===============================================================================
c
c     propagation within outer core
c
      if(lylw.ge.lycc)then
c
c       interface conditions: solid to liquid
c
        y4max=dabs(ylw6(4,3,lycc))
        j0=3
        do j=1,2
          if(y4max.lt.dabs(ylw6(4,j,lycc)))then
            y4max=dabs(ylw6(4,j,lycc))
            j0=j
          endif
        enddo
        do i=1,6
          cyswap=ylw6(i,j0,lycc)
          ylw6(i,j0,lycc)=ylw6(i,3,lycc)
          ylw6(i,3,lycc)=cyswap
        enddo
        do j=1,2
          ylwc(1,j)=ylw6(1,j,lycc)
     &             -ylw6(4,j,lycc)*ylw6(1,3,lycc)/ylw6(4,3,lycc)
          ylwc(2,j)=ylw6(2,j,lycc)
     &             -ylw6(4,j,lycc)*ylw6(2,3,lycc)/ylw6(4,3,lycc)
          ylwc(3,j)=ylw6(5,j,lycc)
     &             -ylw6(4,j,lycc)*ylw6(5,3,lycc)/ylw6(4,3,lycc)
          ylwc(4,j)=ylw6(6,j,lycc)
     &             -ylw6(4,j,lycc)*ylw6(6,3,lycc)/ylw6(4,3,lycc)
        enddo
c
        do ly=lylw,lycc
          do i=1,6
            cyswap=ylw6(i,j0,ly)
            ylw6(i,j0,ly)=ylw6(i,3,ly)
            ylw6(i,3,ly)=cyswap
          enddo
          do j=1,2
            do i=1,6
              ylw6(i,j,ly)=ylw6(i,j,ly)
     &                   -ylw6(4,j,ly)*ylw6(i,3,ly)/ylw6(4,3,ly)
            enddo
          enddo
          do i=1,6
            ylw6(i,3,ly)=0.d0
          enddo
        enddo
c
        ca=ylwc(2,2)
        cb=ylwc(2,1)
        do i=1,4
          ylwc(i,1)=ca*ylwc(i,1)-cb*ylwc(i,2)
          ylwc(i,2)=0.d0
        enddo
        ylwc(2,1)=0.d0
        ylwc(2,2)=1.d0
c
        do ly=lylw,lycc
          do i=1,6
            ylw6(i,1,ly)=ca*ylw6(i,j,ly)-cb*ylw6(i,2,ly)
            ylw6(i,2,ly)=0.d0
          enddo
        enddo
      else if(lylw.ge.lycm)then
        do j=1,2
          do i=1,4
            ylwc(i,j)=mas4x4(i,2*j-1,lylw)
          enddo
        enddo
c
        do j=1,2
          ylw6(1,j,lylw)=ylwc(1,j)
          ylw6(2,j,lylw)=0.d0
          ylw6(3,j,lylw)=ylwc(2,j)
          ylw6(4,j,lylw)=0.d0
          ylw6(5,j,lylw)=ylwc(3,j)
          ylw6(6,j,lylw)=ylwc(4,j)
        enddo
        do i=1,6
          ylw6(i,3,lylw)=0.d0
        enddo
      endif
c
      do ly=min0(lylw,lycc)-1,lycm,-1
        wave(1)=(rrlw(ly)/rrup(ly))**ldeg
        wave(2)=(rrlw(ly)/rrup(ly))**(ldeg+1)
        wave(3)=(rrlw(ly)/rrup(ly))**ldeg
        wave(4)=(rrlw(ly)/rrup(ly))**(ldeg+1)
c
        call axb(mas4x4inv(1,1,ly),ylwc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(1,1)*cc0(3,2)-cc0(3,1)*cc0(1,2)
        orthc(1,1)=cc0(3,2)/cdet
        orthc(1,2)=-cc0(1,2)/cdet
        orthc(2,1)=-cc0(3,1)/cdet
        orthc(2,2)=cc0(1,1)/cdet
c
        call axb(cc0,orthc,4,2,2,cc1)
c
        do iy=lylw,ly+1,-1
c
c         orthonormalization of the receiver vectors
c
          call axb(ylw6(1,1,iy),orthc,6,2,2,y1)
          call memcpy(y1,ylw6(1,1,iy),12)
          do j=1,2
            do i=1,6
              ylw6(i,j,iy)=ylw6(i,j,iy)*wave(2*j-1)
            enddo
          enddo
        enddo
c
        cc1(1,1)=1.d0
        cc1(2,1)=cc1(2,1)*wave(2)*wave(1)
        cc1(3,1)=0.d0
        cc1(4,1)=cc1(4,1)*wave(4)*wave(1)
c
        cc1(1,2)=0.d0
        cc1(2,2)=cc1(2,2)*wave(2)*wave(3)
        cc1(3,2)=1.d0
        cc1(4,2)=cc1(4,2)*wave(4)*wave(3)
c
        call axb(mas4x4(1,1,ly),cc1,4,4,2,ylwc)
        do j=1,2
          ylw6(1,j,ly)=ylwc(1,j)
          ylw6(2,j,ly)=0.d0
          ylw6(3,j,ly)=ylwc(2,j)
          ylw6(4,j,ly)=0.d0
          ylw6(5,j,ly)=ylwc(3,j)
          ylw6(6,j,ly)=ylwc(4,j)
        enddo
        do i=1,6
          ylw6(i,3,ly)=0.d0
        enddo
      enddo
c
c===============================================================================
c
c     propagation from core-mantle boundary to source
c
      if(lylw.ge.lycm)then
c
c       interface conditions: liquid to solid
c
        do j=1,2
          ylw6(1,j,lycm)=ylwc(1,j)
          ylw6(2,j,lycm)=0.d0
          ylw6(3,j,lycm)=0.d0
          ylw6(4,j,lycm)=0.d0
          ylw6(5,j,lycm)=ylwc(3,j)
          ylw6(6,j,lycm)=ylwc(4,j)
        enddo
        ylw6(1,3,lycm)=0.d0
        ylw6(2,3,lycm)=0.d0
        ylw6(3,3,lycm)=1.d0
        ylw6(4,3,lycm)=0.d0
        ylw6(5,3,lycm)=0.d0
        ylw6(6,3,lycm)=0.d0
      else
        do j=1,3
          do i=1,6
            ylw6(i,j,lylw)=mas6x6(i,2*j-1,lylw)
          enddo
        enddo
      endif
c
      do ly=min0(lylw,lycm)-1,lys1,-1
        wave(1)=(rrlw(ly)/rrup(ly))**ldeg
        wave(2)=(rrlw(ly)/rrup(ly))**(ldeg+1)
        wave(3)=(rrlw(ly)/rrup(ly))**(ldeg+2)
        wave(4)=(rrlw(ly)/rrup(ly))**(ldeg-1)
        wave(5)=(rrlw(ly)/rrup(ly))**ldeg
        wave(6)=(rrlw(ly)/rrup(ly))**(ldeg+1)
c
        call axb(mas6x6inv(1,1,ly),ylw6(1,1,ly+1),6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call axb(c0,orth,6,3,3,c1)
c
        do iy=lylw,ly+1,-1
c
c         orthonormalization of the receiver vectors
c
          call axb(ylw6(1,1,iy),orth,6,3,3,y1)
          call memcpy(y1,ylw6(1,1,iy),18)
          do j=1,3
            do i=1,6
              ylw6(i,j,iy)=ylw6(i,j,iy)*wave(2*j-1)
            enddo
          enddo
        enddo
c
        c1(1,1)=1.d0
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=0.d0
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=0.d0
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=0.d0
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=1.d0
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=0.d0
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=0.d0
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=0.d0
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=1.d0
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call axb(mas6x6(1,1,ly),c1,6,6,3,ylw6(1,1,ly))
      enddo
c
      do ly=lyup,lys2
        do j=1,3
          yup6(1,j,ly)=yup6(1,j,ly)/rrup(ly)
          yup6(2,j,ly)=yup6(2,j,ly)/rrup(ly)**2
          yup6(3,j,ly)=yup6(3,j,ly)/rrup(ly)
          yup6(4,j,ly)=yup6(4,j,ly)/rrup(ly)**2
c         yup6(5,j,ly)=yup6(5,j,ly)
          yup6(6,j,ly)=yup6(6,j,ly)/rrup(ly)
        enddo
      enddo
c
      do ly=lys1,lylw
        do j=1,3
          ylw6(1,j,ly)=ylw6(1,j,ly)/rrup(ly)
          ylw6(2,j,ly)=ylw6(2,j,ly)/rrup(ly)**2
          ylw6(3,j,ly)=ylw6(3,j,ly)/rrup(ly)
          ylw6(4,j,ly)=ylw6(4,j,ly)/rrup(ly)**2
c         ylw6(5,j,ly)=ylw6(5,j,ly)
          ylw6(6,j,ly)=ylw6(6,j,ly)/rrup(ly)
        enddo
      enddo
c
      if(lys1.eq.lys2)then
c
c       for point source
c
        do istp=1,4
          do i=1,6
            b6(i,istp)=0.d0
          enddo
          b6(istp,istp)=1.d0
        enddo
        do j=1,3
          do i=1,6
            coef6(i,j)=yup6(i,j,lys1)
            coef6(i,j+3)=-ylw6(i,j,lys1)
          enddo
        enddo
        key=0
        call svd500(coef6,b6,6,4,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpcsprop: anormal exit from svd500!'
          return
        endif
        do ly=lyup,lys1
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,4
              do i=1,6
                do j=1,3
                  ypsv(i,istp,iz)=ypsv(i,istp,iz)
     &                           +b6(j,istp)*yup6(i,j,ly)
                enddo
              enddo
            enddo
          endif
        enddo
c
        do ly=lys1+1,lylw
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,4
              do i=1,6
                do j=1,3
                  ypsv(i,istp,iz)=ypsv(i,istp,iz)
     &                           +b6(j+3,istp)*ylw6(i,j,ly)
                enddo
              enddo
            enddo
          endif
        enddo
c
        if(depr(1).lt.0.d0.or.depr(nz).lt.0.d0)then
c
c         for potential at receivers outside the earth
c
          do istp=1,4
            ypot(istp)=0.d0
            do j=1,3
              ypot(istp)=ypot(istp)+b6(j,istp)*yup6(5,j,1)
            enddo
          enddo
          do iz=1,nz
            if(depr(iz).lt.0.d0)then
              ca=(1.d0/(1.d0-depr(iz)/rearth))**(ldeg+1)
              do istp=1,4
                ypsv(5,istp,iz)=ypsv(5,istp,iz)+ypot(istp)*ca
              enddo
            endif
          enddo
        endif
      else
c
c       for line source
c
        fac=-rrup(lys)**2/gddsnorm
c
        do iys=lys1,lys2-1
          do istp=1,4
            do i=1,6
              b6(i,istp)=0.d0
            enddo
          enddo
          call sruku(b6,6,4,iys,ldeg,qpcmaipsv,
     &               fac,rrup(lys),gdds0,rrup(iys),rrlw(iys))
          do j=1,3
            do i=1,6
              coef6(i,j)=yup6(i,j,iys+1)
              coef6(i,j+3)=-ylw6(i,j,iys+1)
            enddo
          enddo
c
          key=0
          call svd500(coef6,b6,6,4,0.d0,key)
          if(key.eq.0)then
            print *,' Warning in qpcsprop: anormal exit from svd500!'
            return
          endif
c
          do ly=lyup,iys
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,4
                do i=1,6
                  do j=1,3
                    ypsv(i,istp,iz)=ypsv(i,istp,iz)
     &                             +b6(j,istp)*yup6(i,j,ly)
                  enddo
                enddo
              enddo
            endif
          enddo
c
          do ly=iys+1,lylw
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,4
                do i=1,6
                  do j=1,3
                    ypsv(i,istp,iz)=ypsv(i,istp,iz)
     &                             +b6(j+3,istp)*ylw6(i,j,ly)
                  enddo
                enddo
              enddo
            endif
          enddo
c
          if(depr(1).lt.0.d0.or.depr(nz).lt.0.d0)then
c
c           for potential at receivers outside the earth
c
            do istp=1,4
              ypot(istp)=0.d0
              do j=1,3
                ypot(istp)=ypot(istp)+b6(j,istp)*yup6(5,j,1)
              enddo
            enddo
            do iz=1,nz
              if(depr(iz).lt.0.d0)then
                ca=(1.d0/(1.d0-depr(iz)/rearth))**(ldeg+1)
                do istp=1,4
                  ypsv(5,istp,iz)=ypsv(5,istp,iz)+ypot(istp)*ca
                enddo
              endif
            enddo
          endif
        enddo
      endif
      return
      end