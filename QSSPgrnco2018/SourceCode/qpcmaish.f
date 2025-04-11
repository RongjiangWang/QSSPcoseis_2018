      subroutine qpcmaish(ly,ldeg,rr,mat)
      use qpcalloc
      implicit none
      integer*4 ly,ldeg,idiff
      real*8 rr
      real*8 mat(6,6)
c
c     3x3 coefficient matrix for spheroidal mode l = 0
c
      integer*4 i,j,key
      real*8 up,lw,muerr
c
      up=(rr-rrlw(ly))/(rrup(ly)-rrlw(ly))
      lw=1.d0-up
c
      muerr=up*mueup(ly)+lw*muelw(ly)
c
      mat(1,1)=1.d0/rr
      mat(1,2)=1.d0/muerr
c
      mat(2,1)=dble(ldeg-1)*dble(ldeg+2)*muerr/rr**2
      mat(2,2)=-3.d0/rr
      return
      end