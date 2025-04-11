      subroutine qpcdifmatl(ly,ldeg,rr,mat)
      use qpcalloc
      implicit none
      integer*4 ly,ldeg
      real*8 rr
      real*8 mat(6,6)
c
c     4x4 coefficient matrix for spheroidal mode l > 0 in liquid media
c
      real*8 up,lw,xx,x1,rr1
      real*8 grrr,garr,rhorr,kaprr,mass,drho,rho1
      real*8 dldeg,dlp1,dlll1
c
      dldeg=dble(ldeg)
      dlp1=dldeg+1.d0
      dlll1=dldeg*dlp1
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
c
      garr=4.d0*PI*BIGG*rhorr
c
      mat(1,1)=rhorr*grrr/kaprr-1.d0/rr
      mat(1,2)=dlll1/rr
      mat(1,3)=-rr*rhorr/kaprr
      mat(1,4)=0.d0
c
      mat(2,1)=1.d0/rr
      mat(2,2)=0.d0
      mat(2,3)=0.d0
      mat(2,4)=0.d0
c
      mat(3,1)=garr/rr
      mat(3,2)=0.d0
      mat(3,3)=-dlp1/rr
      mat(3,4)=1.d0/rr
c
      mat(4,1)=garr*dlp1/rr
      mat(4,2)=-garr*dlll1/rr
      mat(4,3)=0.d0
      mat(4,4)=dldeg/rr
      return
      end