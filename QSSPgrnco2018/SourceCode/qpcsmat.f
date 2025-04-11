      subroutine qpcsmat(ldeg,ly,lylw)
      use qpcalloc
      implicit none
c
c     calculate 6x6 spheroidal layer matrix for a solid shell
c
      integer*4 ldeg,ly,lylw
c
      integer*4 i,j,key
      real*8 dldeg,dlp1,dlp2,dlm1,dllm1,dllp1,d2lp1,dllp2
      real*8 da,d2mue,dksi,deta,dga
      real*8 mas(6,6)
c
      dldeg=dble(ldeg)
      d2mue=0.5d0*(mueup(ly)+muelw(ly))
      dga=2.d0*PI*BIGG*(rhoup(ly)+rholw(ly))
c
      dksi=d2mue/(kapup(ly)+kaplw(ly)+d2mue*4.d0/3.d0)
      deta=1.d0-dksi
c
      dlp1=dldeg+1.d0
      dlp2=dldeg+2.d0
      dlm1=dldeg-1.d0
      dllm1=dldeg*dlm1
      dllp1=dldeg*dlp1
      dllp2=dldeg*dlp2
      d2lp1=2.d0*dldeg+1.d0
c
      mas6x6(1,1,ly)=dldeg
      mas6x6(2,1,ly)=d2mue*dllm1
      mas6x6(3,1,ly)=1.d0
      mas6x6(4,1,ly)=d2mue*dlm1
      mas6x6(5,1,ly)=dga
      mas6x6(6,1,ly)=dga*dlp1
c
      mas6x6(1,2,ly)=-dlp1
      mas6x6(2,2,ly)=d2mue*dlp1*dlp2
      mas6x6(3,2,ly)=1.d0
      mas6x6(4,2,ly)=-d2mue*dlp2
      mas6x6(5,2,ly)=dga
      mas6x6(6,2,ly)=dga*dlp1
c
      mas6x6(1,3,ly)=dldeg*deta-2.d0*dksi
      mas6x6(2,3,ly)=d2mue*(dllm1*deta+4.d0*dksi-3.d0)
      mas6x6(3,3,ly)=deta+2.d0/dlp1
      mas6x6(4,3,ly)=d2mue*(dlm1*deta-2.d0*dksi+d2lp1/dlp1)
      mas6x6(5,3,ly)=-dga*dksi
      mas6x6(6,3,ly)=dga*(dlp1*deta-d2lp1)
c
      mas6x6(1,4,ly)=dllp1*deta+2.d0*dldeg*dksi
      mas6x6(2,4,ly)=d2mue*dldeg*(-dlp1*dlp2*deta-4.d0*dksi+3.d0)
      mas6x6(3,4,ly)=-dldeg*deta+2.d0
      mas6x6(4,4,ly)=d2mue*(dllp2*deta+2.d0*dldeg*dksi-d2lp1)
      mas6x6(5,4,ly)=dga*dldeg*dksi
      mas6x6(6,4,ly)=-dga*dllp1*deta
c
      mas6x6(1,5,ly)=0.d0
      mas6x6(2,5,ly)=0.d0
      mas6x6(3,5,ly)=0.d0
      mas6x6(4,5,ly)=0.d0
      mas6x6(5,5,ly)=1.d0
      mas6x6(6,5,ly)=d2lp1
c
      mas6x6(1,6,ly)=0.d0
      mas6x6(2,6,ly)=0.d0
      mas6x6(3,6,ly)=0.d0
      mas6x6(4,6,ly)=0.d0
      mas6x6(5,6,ly)=1.d0
      mas6x6(6,6,ly)=0.d0
c
      if(ly.eq.lylw)return
c
      do j=1,6
        do i=1,6
          mas(i,j)=mas6x6(i,j,ly)
          mas6x6inv(i,j,ly)=0.d0
        enddo
        mas6x6inv(j,j,ly)=1.d0
      enddo
      key=0
      call svd500(mas,mas6x6inv(1,1,ly),6,6,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpcsmat: anormal exit from svd500!'
        return
      endif
c
      return
      end