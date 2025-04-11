      subroutine qpcsmatc(ldeg,ly,lylw)
      use qpcalloc
      implicit none
c
c     calculate 4x4 spheroidal layer matrix for a liquid shell
c
      integer*4 ldeg,ly,lylw
c
      integer*4 i,j,key
      real*8 ga,dldeg,dlp1,d2lp1,d3lp2
      real*8 mas(4,4)
c
      dldeg=dble(ldeg)
      dlp1=dldeg+1.d0
      d2lp1=2.d0*dldeg+1.d0
      d3lp2=3.d0*dldeg+2.d0
      ga=2.d0*PI*BIGG*(rhoup(ly)+rholw(ly))
c
c     y1 <- y1 (normal displacement)
c     y2 <- y3 (horizontal displacement)
c     y3 <- y5 (potential)
c     y4 <- y6 (gravity)
c
      mas4x4(1,1,ly)=dldeg
      mas4x4(2,1,ly)=1.d0
      mas4x4(3,1,ly)=0.d0
      mas4x4(4,1,ly)=-ga*dldeg
c
      mas4x4(1,2,ly)=-dlp1
      mas4x4(2,2,ly)=1.d0
      mas4x4(3,2,ly)=ga
      mas4x4(4,2,ly)=0.d0
c
      mas4x4(1,3,ly)=0.d0
      mas4x4(2,3,ly)=0.d0
      mas4x4(3,3,ly)=1.d0
      mas4x4(4,3,ly)=d2lp1
c
      mas4x4(1,4,ly)=0.d0
      mas4x4(2,4,ly)=0.d0
      mas4x4(3,4,ly)=0.d0
      mas4x4(4,4,ly)=1.d0
c
      if(ly.eq.lylw)return
c
      do j=1,4
        do i=1,4
          mas(i,j)=mas4x4(i,j,ly)
          mas4x4inv(i,j,ly)=0.d0
        enddo
        mas4x4inv(j,j,ly)=1.d0
      enddo
      key=0
      call svd500(mas,mas4x4inv(1,1,ly),4,4,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpcsmatc: anormal exit from svd500!'
        return
      endif
c
      return
      end