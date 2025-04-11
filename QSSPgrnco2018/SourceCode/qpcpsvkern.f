      subroutine qpcpsvkern(ldeg)
      use qpcalloc
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c
      integer*4 ldeg
c
      integer*4 i,iz,ly,istp,lyup,lylw
      
c
      do iz=1,nz
        do istp=1,4
          do i=1,6
            ypsv(i,istp,iz)=0.d0
          enddo
        enddo
      enddo
c
      lyup=lyuppsv(ldeg)
      lylw=lylwpsv(ldeg)
c
      if(ldeg.eq.1)return
c
      if(ldeg.eq.0)then
        do ly=lyup,lylw
          call qpcsmat0(ly,lylw)
        enddo
        call qpcsprop0(lyup,lylw)
      else
        do ly=lyup,min0(lycm-1,lylw)
          call qpcsmat(ldeg,ly,lylw)
        enddo
        do ly=max0(lyup,lycm),min0(lycc-1,lylw)
          call qpcsmatc(ldeg,ly,lylw)
        enddo
        do ly=max0(lyup,lycc),min0(ly0,lylw)
          call qpcsmat(ldeg,ly,lylw)
        enddo
        call qpcsprop(ldeg,lyup,lylw)
      endif
      return
      end
