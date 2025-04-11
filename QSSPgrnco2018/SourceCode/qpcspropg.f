      subroutine qpcspropg(ldeg,lyup,lylw)
      use qpcalloc
      implicit none
c
c     calculation of spheroidal response
c
      integer*4 ldeg,lyup,lylw
c
c     work space
c
      integer*4 i,j,j0,istp,ly,iy,iys,iz,ily,nly,key
      real*8 y4max,rr1,rr2,dlnr,h,f
      real*8 cdet,alf,bet,cyabs,cyswap,ca,cb,fac
      real*8 yupc(4,2),ylwc(4,2),ypot(4)
      real*8 coef6(6,6),b6(6,4),c(2),swap(6,3)
      external qpcdifmatl,qpcdifmats,qpcmaipsvg
c
      do iz=1,nz
        do istp=1,4
          do i=1,6
            ypsvg(i,istp,iz)=0.d0
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
      if(lyup.eq.1)then
        do j=1,3
          do i=1,6
            yup6(i,j,lyup)=0.d0
          enddo
        enddo
        yup6(1,1,lyup)=1.d0
        yup6(3,2,lyup)=1.d0
        yup6(5,3,lyup)=1.d0
        if(ldeg.eq.1)then
          yup6(6,3,lyup)=3.d0
        endif
      else
        call qpcsmat(ldeg,lyup,lyup)
        do j=1,3
          do i=1,6
            yup6(i,j,lyup)=mas6x6(i,2*j,lyup)
          enddo
        enddo
      endif
c
c===============================================================================
c
c     propagation from atmosphere/ocean bottom to source
c
      do ly=lyup,lys2-1
        call memcpy(yup6(1,1,ly),swap,18)
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h*dble(ldeg)/rrlw(ly))
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,1)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,1)=swap(i,1)*cyabs
          enddo
          do iy=lyup,ly
            do i=1,6
              yup6(i,1,iy)=yup6(i,1,iy)*cyabs
            enddo
          enddo
c
          alf=0.d0
          do i=1,6
            alf=alf+swap(i,2)*swap(i,1)/cypnorm(i,ly)**2
          enddo
          do i=1,6
            swap(i,2)=swap(i,2)-alf*swap(i,1)
          enddo
          do iy=lyup,ly
            do i=1,6
              yup6(i,2,iy)=yup6(i,2,iy)-alf*yup6(i,1,iy)
            enddo
          enddo
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,2)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,2)=swap(i,2)*cyabs
          enddo
          do iy=lyup,ly
            do i=1,6
              yup6(i,2,iy)=yup6(i,2,iy)*cyabs
            enddo
          enddo
c
          alf=0.d0
          bet=0.d0
          do i=1,6
            alf=alf+swap(i,3)*swap(i,1)/cypnorm(i,ly)**2
            bet=bet+swap(i,3)*swap(i,2)/cypnorm(i,ly)**2
          enddo
          do i=1,6
            swap(i,3)=swap(i,3)-alf*swap(i,1)-bet*swap(i,2)
          enddo
          do iy=lyup,ly
            do i=1,6
              yup6(i,3,iy)=yup6(i,3,iy)-alf*yup6(i,1,iy)
     &                                 -bet*yup6(i,2,iy)
            enddo
          enddo
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,3)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,3)=swap(i,3)*cyabs
          enddo
          do iy=lyup,ly
            do i=1,6
              yup6(i,3,iy)=yup6(i,3,iy)*cyabs
            enddo
          enddo
c
          call druku(swap,6,3,ly,ldeg,
     &               qpcdifmats,rr1,rr2,ndruku(ldeg,ly))
        enddo
        call memcpy(swap,yup6(1,1,ly+1),18)
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
        call qpcsmat(ldeg,lylw,lylw)
        do j=1,3
          do i=1,6
            ylw6(i,j,lylw)=mas6x6(i,2*j-1,lylw)
          enddo
        enddo
      endif
c
      do ly=lylw-1,lycc,-1
        call memcpy(ylw6(1,1,ly+1),swap,18)
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h*dble(ldeg)/rrlw(ly))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,1)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,1)=swap(i,1)*cyabs
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,1,iy)=ylw6(i,1,iy)*cyabs
            enddo
          enddo
c
          alf=0.d0
          do i=1,6
            alf=alf+swap(i,2)*swap(i,1)/cypnorm(i,ly)**2
          enddo
          do i=1,6
            swap(i,2)=swap(i,2)-alf*swap(i,1)
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,2,iy)=ylw6(i,2,iy)-alf*ylw6(i,1,iy)
            enddo
          enddo
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,2)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,2)=swap(i,2)*cyabs
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,2,iy)=ylw6(i,2,iy)*cyabs
            enddo
          enddo
c
          alf=0.d0
          bet=0.d0
          do i=1,6
            alf=alf+swap(i,3)*swap(i,1)/cypnorm(i,ly)**2
            bet=bet+swap(i,3)*swap(i,2)/cypnorm(i,ly)**2
          enddo
          do i=1,6
            swap(i,3)=swap(i,3)-alf*swap(i,1)-bet*swap(i,2)
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,3,iy)=ylw6(i,3,iy)-alf*ylw6(i,1,iy)
     &                                 -bet*ylw6(i,2,iy)
            enddo
          enddo
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,3)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,3)=swap(i,3)*cyabs
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,3,iy)=ylw6(i,3,iy)*cyabs
            enddo
          enddo
c
          call druku(swap,6,3,ly,ldeg,
     &               qpcdifmats,rr1,rr2,ndruku(ldeg,ly))
        enddo
        call memcpy(swap,ylw6(1,1,ly),18)
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
c       y2 = Ut
c
        ca=rholw(lycc-1)*grlw(lycc-1)*rrlw(lycc-1)
        cb=rholw(lycc-1)*rrlw(lycc-1)**2
        do j=1,2
          c(j)=ca*ylwc(1,j)-ylwc(2,j)-cb*ylwc(3,j)
        enddo
        do i=1,4
          ylwc(i,1)=c(2)*ylwc(i,1)-c(1)*ylwc(i,2)
        enddo
        ylwc(2,1)=0.d0
        do i=1,4
          ylwc(i,2)=0.d0
        enddo
        ylwc(2,2)=c(2)
c
        do ly=lylw,lycc+1
          do i=1,6
            ylw6(i,1,ly)=c(2)*ylw6(i,1,ly)-c(1)*ylw6(i,2,ly)
            ylw6(i,2,ly)=0.d0
          enddo
        enddo
      else if(lylw.ge.lycm)then
        call qpcsmatc(ldeg,lylw,lylw)
        do j=1,2
          do i=1,4
            ylwc(i,j)=mas4x4(i,2*j-1,lylw)
          enddo
        enddo
        ca=rhoup(lylw)*grup(lylw)*rrup(lylw)
        cb=rhoup(lylw)*rrup(lylw)**2
        do j=1,2
          ylw6(1,j,lylw)=ylwc(1,j)
          ylw6(2,j,lylw)=ca*ylwc(1,j)-cb*ylwc(3,j)
          ylw6(3,j,lylw)=ylwc(2,j)
          ylw6(5,j,lylw)=ylwc(3,j)
          ylw6(6,j,lylw)=ylwc(4,j)
        enddo
        do i=1,6
          ylw6(i,3,lylw)=0.d0
        enddo
      endif
c
      do ly=min0(lylw-1,lycc-1),lycm,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h*dble(ldeg)/rrlw(ly))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=0.d0
          do i=1,4
            cyabs=cyabs+ylwc(i,1)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,4
            ylwc(i,1)=ylwc(i,1)*cyabs
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,1,iy)=ylw6(i,1,iy)*cyabs
            enddo
          enddo
c
          alf=0.d0
          do i=1,4
            alf=alf+ylwc(i,2)*ylwc(i,1)/cypnorm(i,ly)**2
          enddo
          do i=1,4
            ylwc(i,2)=ylwc(i,2)-alf*ylwc(i,1)
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,2,iy)=ylw6(i,2,iy)-alf*ylw6(i,1,iy)
            enddo
          enddo
c
          cyabs=0.d0
          do i=1,4
            cyabs=cyabs+ylwc(i,2)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,4
            ylwc(i,2)=ylwc(i,2)*cyabs
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,2,iy)=ylw6(i,2,iy)*cyabs
            enddo
          enddo
c
          call druku(ylwc,4,2,ly,ldeg,
     &               qpcdifmatl,rr1,rr2,ndruku(ldeg,ly))
        enddo
        ca=rhoup(ly)*grup(ly)*rrup(ly)
        cb=rhoup(ly)*rrup(ly)**2
        do j=1,2
          ylw6(1,j,ly)=ylwc(1,j)
          ylw6(2,j,ly)=ca*ylwc(1,j)-cb*ylwc(3,j)
          ylw6(3,j,ly)=ylwc(2,j)
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
c     propagation from core-mantle boundary to source or ocean bottom
c
      if(lylw.ge.lycm)then
c
c       interface conditions: liquid to solid
c
        ca=rhoup(lycm)*grup(lycm)*rrup(lycm)
        cb=rhoup(lycm)*rrup(lycm)**2
        do j=1,2
          ylw6(1,j,lycm)=ylwc(1,j)
          ylw6(2,j,lycm)=ca*ylwc(1,j)-cb*ylwc(3,j)
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
        call qpcsmat(ldeg,lylw,lylw)
        do j=1,3
          do i=1,6
            ylw6(i,j,lylw)=mas6x6(i,2*j-1,lylw)
          enddo
        enddo
      endif
c
      do ly=min0(lylw-1,lycm-1),lys1,-1
        call memcpy(ylw6(1,1,ly+1),swap,18)
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h*dble(ldeg)/rrlw(ly))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,1)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,1)=swap(i,1)*cyabs
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,1,iy)=ylw6(i,1,iy)*cyabs
            enddo
          enddo
c
          alf=0.d0
          do i=1,6
            alf=alf+swap(i,2)*swap(i,1)/cypnorm(i,ly)**2
          enddo
          do i=1,6
            swap(i,2)=swap(i,2)-alf*swap(i,1)
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,2,iy)=ylw6(i,2,iy)-alf*ylw6(i,1,iy)
            enddo
          enddo
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,2)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,2)=swap(i,2)*cyabs
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,2,iy)=ylw6(i,2,iy)*cyabs
            enddo
          enddo
c
          alf=0.d0
          bet=0.d0
          do i=1,6
            alf=alf+swap(i,3)*swap(i,1)/cypnorm(i,ly)**2
            bet=bet+swap(i,3)*swap(i,2)/cypnorm(i,ly)**2
          enddo
          do i=1,6
            swap(i,3)=swap(i,3)-alf*swap(i,1)-bet*swap(i,2)
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,3,iy)=ylw6(i,3,iy)-alf*ylw6(i,1,iy)
     &                                 -bet*ylw6(i,2,iy)
            enddo
          enddo
c
          cyabs=0.d0
          do i=1,6
            cyabs=cyabs+swap(i,3)**2/cypnorm(i,ly)**2
          enddo
          cyabs=1.d0/dsqrt(cyabs)
          do i=1,6
            swap(i,3)=swap(i,3)*cyabs
          enddo
          do iy=lylw,ly+1,-1
            do i=1,6
              ylw6(i,3,iy)=ylw6(i,3,iy)*cyabs
            enddo
          enddo
c
          call druku(swap,6,3,ly,ldeg,
     &               qpcdifmats,rr1,rr2,ndruku(ldeg,ly))
        enddo
        call memcpy(swap,ylw6(1,1,ly),18)
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
          print *,' Warning in qpcspropg: anormal exit from svd500!'
          return
        endif
        do ly=lyup,lys1
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,4
              do i=1,6
                do j=1,3
                  ypsvg(i,istp,iz)=ypsvg(i,istp,iz)
     &                           +b6(j,istp)*yup6(i,j,ly)
                enddo
              enddo
            enddo
          endif
        enddo
        do ly=lys1+1,lylw
          iz=izrly(ly)
          if(iz.gt.0)then
            do istp=1,4
              do i=1,6
                do j=1,3
                  ypsvg(i,istp,iz)=ypsvg(i,istp,iz)
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
                ypsvg(5,istp,iz)=ypsvg(5,istp,iz)+ypot(istp)*ca
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
          call sruku(b6,6,4,iys,ldeg,qpcmaipsvg,
     &               fac,rrup(lys),gdds0,rrup(iys),rrlw(iys))
c
          do j=1,3
            do i=1,6
              coef6(i,j)=yup6(i,j,iys+1)
              coef6(i,j+3)=-ylw6(i,j,iys+1)
            enddo
          enddo
          key=0
          call svd500(coef6,b6,6,4,0.d0,key)
          if(key.eq.0)then
            print *,' Warning in qpcspropg: anormal exit from svd500!'
            return
          endif
c
          do ly=lyup,iys
            iz=izrly(ly)
            if(iz.gt.0)then
              do istp=1,4
                do i=1,6
                  do j=1,3
                    ypsvg(i,istp,iz)=ypsvg(i,istp,iz)
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
                    ypsvg(i,istp,iz)=ypsvg(i,istp,iz)
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
                  ypsvg(5,istp,iz)=ypsvg(5,istp,iz)+ypot(istp)*ca
                enddo
              endif
            enddo
          endif
        enddo
      endif
      return
      end