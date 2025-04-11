      subroutine qpctmat(ldeg,ly,lylw)
      use qpcalloc
      implicit none
c
      integer*4 ldeg,ly,lylw
c
      real*8 dlm1,dlp2,d2lp1
      real*8 mue
c
      dlm1=dble(ldeg-1)
      dlp2=dble(ldeg+2)
      d2lp1=dble(2*ldeg+1)
      mue=0.5d0*(mueup(ly)+muelw(ly))
c
      mat2x2(1,1,ly)=1.d0
      mat2x2(2,1,ly)=mue*dlm1
      mat2x2(1,2,ly)=-1.d0
      mat2x2(2,2,ly)=mue*dlp2
c
      if(ly.eq.lylw)return
c
c     calculate inverse matrix at upper radius
c
      mat2x2inv(1,1,ly)=dlp2/d2lp1
      mat2x2inv(1,2,ly)=1.d0/(mue*d2lp1)
      mat2x2inv(2,1,ly)=-dlm1/d2lp1
      mat2x2inv(2,2,ly)=1.d0/(mue*d2lp1)
c
      return
      end