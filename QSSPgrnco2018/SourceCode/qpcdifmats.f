      subroutine qpcdifmats(ly,ldeg,rr,mat)
      use qpcalloc
      implicit none
      integer*4 ly,ldeg
      real*8 rr
      real*8 mat(6,6)
c
c     6x6 coefficient matrix for spheroidal mode l > 0 in solid media
c
      real*8 up,lw,xx,x1,rr1
      real*8 grrr,garr,rhorr,kaprr,ksirr,lamrr,muerr
      real*8 mass,drho,rho1
      real*8 dldeg,dlp1,dll1
c
      dldeg=dble(ldeg)
      dlp1=dldeg+1.d0
      dll1=dldeg*dlp1
c
      up=(rr-rrlw(ly))/(rrup(ly)-rrlw(ly))
      lw=1.d0-up
c
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
      lamrr=kaprr-muerr*2.d0/3.d0
      ksirr=lamrr+2.d0*muerr
      garr=4.d0*PI*BIGG*rhorr
c
      mat(1,1)=(1.d0-2.d0*lamrr/ksirr)/rr
      mat(1,2)=1.d0/ksirr/rr
      mat(1,3)=dll1*lamrr/ksirr/rr
      mat(1,4)=0.d0
      mat(1,5)=0.d0
      mat(1,6)=0.d0
c
      mat(2,1)=4.d0*(muerr*(1.d0+2.d0*lamrr/ksirr)/rr-rhorr*grrr)
      mat(2,2)=2.d0*lamrr/ksirr/rr
      mat(2,3)=dll1*(rhorr*grrr-2.d0*muerr*(1.d0+2.d0*lamrr/ksirr)/rr)
      mat(2,4)=dll1/rr
      mat(2,5)=rhorr*dlp1*rr
      mat(2,6)=-rhorr*rr
c
      mat(3,1)=-1.d0/rr
      mat(3,2)=0.d0
      mat(3,3)=2.d0/rr
      mat(3,4)=1.d0/muerr/rr
      mat(3,5)=0.d0
      mat(3,6)=0.d0
c
      mat(4,1)=rhorr*grrr-2.d0*muerr*(1.d0+2.d0*lamrr/ksirr)/rr
      mat(4,2)=-lamrr/ksirr/rr
      mat(4,3)=2.d0*muerr*(2.d0*dll1*(1.d0-muerr/ksirr)-1.d0)/rr
      mat(4,4)=-1.d0/rr
      mat(4,5)=-rhorr*rr
      mat(4,6)=0.d0
c
      mat(5,1)=garr/rr
      mat(5,2)=0.d0
      mat(5,3)=0.d0
      mat(5,4)=0.d0
      mat(5,5)=-dlp1/rr
      mat(5,6)=1.d0/rr
c
      mat(6,1)=garr*dlp1/rr
      mat(6,2)=0.d0
      mat(6,3)=-garr*dll1/rr
      mat(6,4)=0.d0
      mat(6,5)=0.d0
      mat(6,6)=dldeg/rr
      return
      end