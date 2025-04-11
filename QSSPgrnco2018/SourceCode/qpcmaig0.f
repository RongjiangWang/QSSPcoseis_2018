      subroutine qpcmaig0(ly,ldeg,rr,mat)
      use qpcalloc
      implicit none
      integer*4 ly,ldeg,idiff
      real*8 rr
      real*8 mat(6,6)
c
c     3x3 coefficient matrix for spheroidal mode l = 0
c
      integer*4 i,j,key
      real*8 up,lw,xx,x1,rr1
      real*8 rhorr,kaprr,ksirr,lamrr,muerr
      real*8 mass,drho,rho1,grrr
c
      up=(rr-rrlw(ly))/(rrup(ly)-rrlw(ly))
      lw=1.d0-up
      rr1=rrlw(ly)
c
      drho=(rhoup(ly)-rholw(ly))/(rrup(ly)-rrlw(ly))
      rho1=rholw(ly)-drho*rrlw(ly) 
      mass=PI*(rr-rr1)*((4.d0/3.d0)*rho1*(rr**2+rr*rr1+rr1**2)
     &    +drho*(rr**3+rr**2*rr1+rr*rr1**2+rr1**3))
      grrr=(grlw(ly)*rr1**2+BIGG*mass)/rr**2
      rhorr=up*rhoup(ly)+lw*rholw(ly)
      kaprr=up*kapup(ly)+lw*kaplw(ly)
      muerr=up*mueup(ly)+lw*muelw(ly)
c
      lamrr=kaprr-muerr*2.d0/3.d0
      ksirr=lamrr+2.d0*muerr
c
      mat(1,1)=-2.d0*lamrr/ksirr/rr
      mat(1,2)=1.d0/ksirr
c
      mat(2,1)=4.d0*(muerr*(1.d0+2.d0*lamrr/ksirr)/rr-rhorr*grrr)/rr
      mat(2,2)=-4.d0*muerr/ksirr/rr
      return
      end