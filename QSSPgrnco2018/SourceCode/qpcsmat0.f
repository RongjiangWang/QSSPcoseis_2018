      subroutine qpcsmat0(ly,lylw)
      use qpcalloc
      implicit none
c
c     calculate 3x3 spheroidal layer matrix for a solid shell in case of degree l = 0
c
      integer*4 ly,lylw
c
      integer*4 i,j,key
      real*8 d3ksi,d4mue,ga
      real*8 mas(3,3)
c
      d4mue=2.d0*(mueup(ly)+muelw(ly))
      d3ksi=1.5d0*(kapup(ly)+kaplw(ly))
      ga=2.d0*PI*BIGG*(rhoup(ly)+rholw(ly))
c
      mas3x3(1,1,ly)=1.d0
      mas3x3(2,1,ly)=d3ksi
      mas3x3(3,1,ly)=ga/2.d0
c
      mas3x3(1,2,ly)=1.d0
      mas3x3(2,2,ly)=-d4mue
      mas3x3(3,2,ly)=-ga
c
      mas3x3(1,3,ly)=0.d0
      mas3x3(2,3,ly)=0.d0
      mas3x3(3,3,ly)=1.d0
c
      if(ly.eq.lylw)return
c
c     calculate inverse matrix
c
      do j=1,3
        do i=1,3
          mas(i,j)=mas3x3(i,j,ly)
          mas3x3inv(i,j,ly)=0.d0
        enddo
        mas3x3inv(j,j,ly)=1.d0
      enddo
      key=0
      call svd500(mas,mas3x3inv(1,1,ly),3,3,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpcsmat0: anormal exit from svd500!'
        return
      endif
c
      return
      end