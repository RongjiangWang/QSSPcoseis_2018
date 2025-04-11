      subroutine qpcwvint(ierr)
      use qpcalloc
      implicit none
      integer*4 ierr
c
      integer*4 ir,iz,ig,id,nd,istp,ldeg,ldegup,ldeg3,ldeg4
      real*8 cp0,cp1,cp2,fac,ca,cb,drs,rrr,muer,lamr,ksir,srt,srp
c
      do ir=1,nr
        ssd(ir)=dsin(disrad(ir))
        csd(ir)=dcos(disrad(ir))
        ssf(ir)=2.d0*dsin(0.5d0*disrad(ir))**2
      enddo
c
      call qpclegendre(ierr)
c
      do ig=1,ngrn
c
        write(*,'(a)')' '
        write(*,'(a)')'  open Green function data base: '
     &              //grnfile(ig)(1:40)
        write(*,'(a)')'  ... please wait ...'
c
        open(21,file=uspecfile(ig),form='unformatted',status='old')
        open(22,file=vspecfile(ig),form='unformatted',status='old')
        open(23,file=wspecfile(ig),form='unformatted',status='old')
        open(24,file=especfile(ig),form='unformatted',status='old')
        open(25,file=fspecfile(ig),form='unformatted',status='old')
        open(26,file=gspecfile(ig),form='unformatted',status='old')
        open(27,file=pspecfile(ig),form='unformatted',status='old')
        open(28,file=qspecfile(ig),form='unformatted',status='old')
c
        open(30,file=spatgrnfile(ig),form='unformatted',
     &       status='unknown')
c
        read(21)ldegup
        read(22)ldegup
        read(23)ldegup
        read(24)ldegup
        read(25)ldegup
        read(26)ldegup
        read(27)ldegup
        read(28)ldegup
c
        ldeg4=ldegup-1
        ldeg3=(ldegup-ndmax)*4/5
        do ldeg=0,ldeg3
          tap(ldeg)=1.d0
        enddo
        do ldeg=1+ldeg3,ldegup-ndmax
          tap(ldeg)=dcos(0.5d0*PI*dble(ldeg-ldeg3)
     &                           /dble(ldegup-ndmax-ldeg3))**2
        enddo
        do ldeg=1+ldegup-ndmax,ldegup
          tap(ldeg)=0.d0
        enddo
c
        do istp=1,4
          do ir=1,nr
            do iz=1,nz
              spatur(iz,ir,istp)=0.d0
              spatut(iz,ir,istp)=0.d0
              spatup(iz,ir,istp)=0.d0
              spaterr(iz,ir,istp)=0.d0
              spatert(iz,ir,istp)=0.d0
              spaterp(iz,ir,istp)=0.d0
              spatetr(iz,ir,istp)=0.d0
              spatett(iz,ir,istp)=0.d0
              spatetp(iz,ir,istp)=0.d0
              spatepr(iz,ir,istp)=0.d0
              spatept(iz,ir,istp)=0.d0
              spatepp(iz,ir,istp)=0.d0
              spatgr(iz,ir,istp)=0.d0
              spatgt(iz,ir,istp)=0.d0
              spatgp(iz,ir,istp)=0.d0
              spatpo(iz,ir,istp)=0.d0
            enddo
          enddo
        enddo
c
        do iz=1,nz
          rrr=rearth-depr(iz)
          if(depr(iz).ge.0.d0)then
            drs=PI*rearth/dble(1+ldeg4)
     &         +5.d0*dabs(rrup(lygrn(ig))-rrr)
          else
            drs=PI*rearth/dble(1+ldeg4)
     &         +5.d0*dabs(rrup(lygrn(ig))-rearth)
          endif
c
          nd=0
          do ir=1,nr
            idr(ir)=min0(ndmax,
     &                   max0(idnint(dlog(dism(ir)/drs)/dlog(5.d0)),0))
            nd=max0(nd,idr(ir))
          enddo
c
          if(depr(iz).ge.0.d0)then
            muer=mueup(lyr(iz))
            lamr=kapup(lyr(iz))-muer*2.d0/3.d0
            ksir=lamr+2.d0*muer
          else
            muer=mueup(1)
            lamr=kapup(1)-muer*2.d0/3.d0
            ksir=lamr+2.d0*muer
          endif
c
          read(21)((ul0(ldeg,istp),ldeg=0,ldegup),istp=1,4)
          read(22)((vl0(ldeg,istp),ldeg=0,ldegup),istp=1,4)
          read(23)((wl0(ldeg,istp),ldeg=0,ldegup),istp=2,3)
          read(24)((el0(ldeg,istp),ldeg=0,ldegup),istp=1,4)
          read(25)((fl0(ldeg,istp),ldeg=0,ldegup),istp=1,4)
          read(26)((gl0(ldeg,istp),ldeg=0,ldegup),istp=2,3)
          read(27)((pl0(ldeg,istp),ldeg=0,ldegup),istp=1,4)
          read(28)((ql0(ldeg,istp),ldeg=0,ldegup),istp=1,4)
c
c         explosion and clvd sources
c         ur,ut,up(=0) -- displacement vector
c         err,ert,erp(=0),etr,ett,etp(=0),epr(=0),ept(=0),epp  --strain tensor
c         gr,gt,gp(=0)  -- gravity vector
c         po -- geopotential
c
          do istp=1,4,3
            do ldeg=0,ldeg4
              urlm(ldeg,istp,0)=ul0(ldeg,istp)                          !m=0,  1
              utlm(ldeg,istp,0)=-vl0(ldeg,istp)                         !m=1,  1
c
              errlm(ldeg,istp,0)=(el0(ldeg,istp)                        !m=0,  1
     &                     +(lamr/rrr)*(-2.d0*ul0(ldeg,istp)
     &                     +dble(ldeg*(ldeg+1))*vl0(ldeg,istp)))/ksir
              ertlm(ldeg,istp,0)=(-ul0(ldeg,istp)                       !m=1,  1
     &                            +vl0(ldeg,istp))/rrr
c
              if(muer.gt.0.d0)then
                srt=-fl0(ldeg,istp)
                etrlm(ldeg,istp,0)=srt/muer-ertlm(ldeg,istp,0)          !m=1,  1
              else
                etrlm(ldeg,istp,0)=ertlm(ldeg,istp,0)                   !m=1,  1
              endif
              ett0lm(ldeg,istp,0)=(ul0(ldeg,istp)
     &                        -dble(ldeg*(ldeg+1))*vl0(ldeg,istp))/rrr  !m=0,  1
              ettalm(ldeg,istp,0)=vl0(ldeg,istp)/rrr                    !m=1,  cos(t)/sin(t)
c
              epp0lm(ldeg,istp,0)= ul0(ldeg,istp)/rrr                   !m=0,  1
              eppalm(ldeg,istp,0)=-vl0(ldeg,istp)/rrr                   !m=1,  cos(t)/sin(t)
c
              grlm(ldeg,istp,0)=ql0(ldeg,istp)                          !m=0,  1
              gtlm(ldeg,istp,0)=-pl0(ldeg,istp)/rrr                     !m=1,  1
              polm(ldeg,istp,0)=pl0(ldeg,istp)                          !m=0,  1
            enddo
          enddo
c
          do istp=2,3
            do ldeg=0,3-istp
              urlm(ldeg,istp,0)=0.d0
              utlm(ldeg,istp,0)=0.d0
              uplm(ldeg,istp,0)=0.d0
              errlm(ldeg,istp,0)=0.d0
              ertlm(ldeg,istp,0)=0.d0
              erplm(ldeg,istp,0)=0.d0
              etrlm(ldeg,istp,0)=0.d0
              ett0lm(ldeg,istp,0)=0.d0
              ettalm(ldeg,istp,0)=0.d0
              ettblm(ldeg,istp,0)=0.d0
              etp0lm(ldeg,istp,0)=0.d0
              etpalm(ldeg,istp,0)=0.d0
              etpblm(ldeg,istp,0)=0.d0
              eprlm(ldeg,istp,0)=0.d0
              ept0lm(ldeg,istp,0)=0.d0
              eptalm(ldeg,istp,0)=0.d0
              eptblm(ldeg,istp,0)=0.d0
              epp0lm(ldeg,istp,0)=0.d0
              eppalm(ldeg,istp,0)=0.d0
              eppblm(ldeg,istp,0)=0.d0
              grlm(ldeg,istp,0)=0.d0
              gtlm(ldeg,istp,0)=0.d0
              gplm(ldeg,istp,0)=0.d0
              polm(ldeg,istp,0)=0.d0
            enddo
          enddo
c
c         dip-slip
c
          do ldeg=1,ldeg4
            ca=dble(ldeg-1)**2/dble(2*ldeg-1)
            cb=dble(ldeg+2)**2/dble(2*ldeg+3)
            urlm(ldeg,3,0)=ul0(ldeg,3)                                  !m=1,  1
            utlm(ldeg,3,0)= ca*vl0(ldeg-1,3)-cb*vl0(ldeg+1,3)           !m=1,  1/sin(t)
     &                     +wl0(ldeg,3)
            uplm(ldeg,3,0)=-ca*wl0(ldeg-1,3)+cb*wl0(ldeg+1,3)           !m=1,  1/sin(t)
     &                     -vl0(ldeg,3)
c
            errlm(ldeg,3,0)=(el0(ldeg,3)                                !m=1,  1
     &                     +(lamr/rrr)*(-2.d0*ul0(ldeg,3)
     &                     +dble(ldeg*(ldeg+1))*vl0(ldeg,3)))/ksir
c           ertlm(ldeg,3,0) s. below
            erplm(ldeg,3,0)=(-ul0(ldeg,3)-uplm(ldeg,3,0))/rrr           !m=1,  1/sin(t)
c
            if(muer.gt.0.d0)then
              etrlm(ldeg,3,0)=ca*(fl0(ldeg-1,3)/muer                    !m=1,  1/sin(t)
     &                            +(vl0(ldeg-1,3)-ul0(ldeg-1,3))/rrr)
     &                       -cb*(fl0(ldeg+1,3)/muer
     &                            +(vl0(ldeg+1,3)-ul0(ldeg+1,3))/rrr)
     &                       +(wl0(ldeg,3)/rrr+gl0(ldeg,3)/muer)
              srt=ca*fl0(ldeg-1,3)-cb*fl0(ldeg+1,3)+gl0(ldeg,3)
              ertlm(ldeg,3,0)=srt/muer-etrlm(ldeg,3,0)                  !m=1,  1/sin(t)
            else
              etrlm(ldeg,3,0)=ca*(ul0(ldeg-1,3)-vl0(ldeg-1,3))/rrr      !m=1,  1/sin(t)
     &                         -cb*(ul0(ldeg+1,3)-vl0(ldeg+1,3))/rrr
              ertlm(ldeg,3,0)=etrlm(ldeg,3,0)                           !m=1,  1/sin(t)
            endif
c
            ett0lm(ldeg,3,0)=(ul0(ldeg,3)                               !m=1,  1
     &                +(1.d0-dble(ldeg*(ldeg+1)))*vl0(ldeg,3))/rrr
            ettalm(ldeg,3,0)= vl0(ldeg,3)/rrr                           !m=2,  cos(t)/sin(t)
            ettblm(ldeg,3,0)=-wl0(ldeg,3)/rrr                           !m=2,  1/sin(t)
c
            etp0lm(ldeg,3,0)=-wl0(ldeg,3)/rrr                           !m=1,  1
            etpalm(ldeg,3,0)= vl0(ldeg,3)/rrr                           !m=2,  1/sin(t)
            etpblm(ldeg,3,0)=-wl0(ldeg,3)/rrr                           !m=2,  cos(t)/sin(t)
c
            if(muer.gt.0.d0)then
              srp=-ca*gl0(ldeg-1,3)+cb*gl0(ldeg+1,3)-fl0(ldeg,3)
              eprlm(ldeg,3,0)=srp/muer-erplm(ldeg,3,0)                  !m=1,  1/sin(t)
            else
              eprlm(ldeg,3,0)=erplm(ldeg,3,0)                           !m=1,  1/sin(t)
            endif
c
            ept0lm(ldeg,3,0)=-(1.d0-dble(ldeg*(ldeg+1)))                !m=1,  1
     &                       *wl0(ldeg,3)/rrr
            eptalm(ldeg,3,0)=-wl0(ldeg,3)/rrr                           !m=2,  cos(t)/sin(t)
            eptblm(ldeg,3,0)= vl0(ldeg,3)/rrr                           !m=2,  1/sin(t)
            
c
            epp0lm(ldeg,3,0)=(ul0(ldeg,3)-vl0(ldeg,3))/rrr              !m=1,  1
            eppalm(ldeg,3,0)=-vl0(ldeg,3)/rrr                           !m=2,  cos(t)/sin(t)
            eppblm(ldeg,3,0)= wl0(ldeg,3)/rrr                           !m=2,  1/sin(t)
c
            grlm(ldeg,3,0)=ql0(ldeg,3)                                  !m=1,  1
            gtlm(ldeg,3,0)=(ca*pl0(ldeg-1,3)-cb*pl0(ldeg+1,3))/rrr      !m=1,  1/sin(t)
            gplm(ldeg,3,0)=-pl0(ldeg,3)/rrr                             !m=1,  1/sin(t)
            polm(ldeg,3,0)=pl0(ldeg,3)                                  !m=1,  1
          enddo
c
c         strike-slip
c
          do ldeg=2,ldeg4
            ca=dble((ldeg-1)*(ldeg-2))/dble(2*ldeg-1)
            cb=dble((ldeg+2)*(ldeg+3))/dble(2*ldeg+3)
            urlm(ldeg,2,0)=ul0(ldeg,2)                                  !m=2,  1
            utlm(ldeg,2,0)= ca*vl0(ldeg-1,2)-cb*vl0(ldeg+1,2)           !m=2,  1/sin(t)
     &                     -2.d0*wl0(ldeg,2)
            uplm(ldeg,2,0)=-ca*wl0(ldeg-1,2)+cb*wl0(ldeg+1,2)           !m=2,  1/sin(t)
     &                     +2.d0*vl0(ldeg,2)
c
            errlm(ldeg,2,0)=(el0(ldeg,2)                                !m=2,  1
     &                     +(lamr/rrr)*(-2.d0*ul0(ldeg,2)
     &                     +dble(ldeg*(ldeg+1))*vl0(ldeg,2)))/ksir
c           ertlm(ldeg,2,0) s. below
            erplm(ldeg,2,0)=(2.d0*ul0(ldeg,2)-uplm(ldeg,2,0))/rrr       !m=2,  1/sin(t)
c
            if(muer.gt.0.d0)then
              etrlm(ldeg,2,0)=ca*(fl0(ldeg-1,2)/muer                    !m=2,  1/sin(t)
     &                          +(vl0(ldeg-1,2)-ul0(ldeg-1,2))/rrr)
     &                       -cb*(fl0(ldeg+1,2)/muer
     &                          +(vl0(ldeg+1,2)-ul0(ldeg+1,2))/rrr)
     &                     -2.d0*(wl0(ldeg,2)/rrr+gl0(ldeg,2)/muer)
              srt=ca*fl0(ldeg-1,2)-cb*fl0(ldeg+1,2)-2.d0*gl0(ldeg,2)
              ertlm(ldeg,2,0)=srt/muer-etrlm(ldeg,2,0)                  !m=2,  1/sin(t)
            else
              etrlm(ldeg,2,0)=ca*(ul0(ldeg-1,2)-vl0(ldeg-1,2))/rrr      !m=2,  1/sin(t)
     &                       -cb*(ul0(ldeg+1,2)-vl0(ldeg+1,2))/rrr
              ertlm(ldeg,2,0)=etrlm(ldeg,2,0)                           !m=2  1/sin(t)
            endif
c
            ett0lm(ldeg,2,0)=(ul0(ldeg,2)                               !m=2,  1
     &                 -dble(ldeg*(ldeg+1))*vl0(ldeg,2))/rrr
            ettalm(ldeg,2,0)=-utlm(ldeg,2,0)/rrr                        !m=2,  cos(t)/sin^2(t)
            ettblm(ldeg,2,0)=2.d0*uplm(ldeg,2,0)/rrr                    !m=2,  1/sin^2(t)
c
            etp0lm(ldeg,2,0)=2.d0*utlm(ldeg,2,0) /rrr                   !m=2,  1/sin^2(t)
            etpalm(ldeg,2,0)=-uplm(ldeg,2,0)/rrr                        !m=2,  cos(t)/sin^2(t)
c
            if(muer.gt.0.d0)then
              srp=-ca*gl0(ldeg-1,2)+cb*gl0(ldeg+1,2)+2.d0*fl0(ldeg,2)
              eprlm(ldeg,2,0)=srp/muer-erplm(ldeg,2,0)                  !m=2,  1/sin(t)
            else
              eprlm(ldeg,2,0)=erplm(ldeg,2,0)                           !m=2,  1/sin(t)
            endif
c
            ept0lm(ldeg,2,0)=-dble(ldeg+2)*uplm(ldeg,2,0)/rrr           !m=2,  cos(t)/sin^2(t)
            eptalm(ldeg,2,0)=dble(ldeg-2)*uplm(ldeg-1,2,0)/rrr          !m=2,  1/sin^2(t)
c
            epp0lm(ldeg,2,0)=ul0(ldeg,2)/rrr                            !m=2,  1
            eppalm(ldeg,2,0)=utlm(ldeg,2,0)/rrr                         !m=2,  cos(t)/sin^2(t)
            eppblm(ldeg,2,0)=-2.d0*uplm(ldeg,2,0)/rrr                   !m=2,  1/sin^2(t)
c
            grlm(ldeg,2,0)=ql0(ldeg,2)                                  !m=2,  1
            gtlm(ldeg,2,0)=(ca*pl0(ldeg-1,2)-cb*pl0(ldeg+1,2))/rrr      !m=2,  1/sin(t)
            gplm(ldeg,2,0)=2.d0*pl0(ldeg,2)/rrr                         !m=2,  1/sin(t)
            polm(ldeg,2,0)=pl0(ldeg,2)                                  !m=2,  1
          enddo
c
c         use differential transform to suppress spatial aliasing
c
          do id=1,nd
c
c           m = 0
c
            do istp=1,4,3
              urlm(0,istp,id)=urlm(0,istp,id-1)
     &                       -urlm(1,istp,id-1)/3.d0
              errlm(0,istp,id)=errlm(0,istp,id-1)
     &                        -errlm(1,istp,id-1)/3.d0
              ett0lm(0,istp,id)=ett0lm(0,istp,id-1)
     &                         -ett0lm(1,istp,id-1)/3.d0
              epp0lm(0,istp,id)=epp0lm(0,istp,id-1)
     &                         -epp0lm(1,istp,id-1)/3.d0
              grlm(0,istp,id)=grlm(0,istp,id-1)
     &                       -grlm(1,istp,id-1)/3.d0
              polm(0,istp,id)=polm(0,istp,id-1)
     &                       -polm(1,istp,id-1)/3.d0
            enddo
c
            do ldeg=1,ldeg4-id
              ca=dble(ldeg+1)/dble(2*ldeg+3)
              cb=dble(ldeg)/dble(2*ldeg-1)
              do istp=1,4,3
                urlm(ldeg,istp,id)=urlm(ldeg,istp,id-1)
     &                         -ca*urlm(ldeg+1,istp,id-1)
     &                         -cb*urlm(ldeg-1,istp,id-1)
                errlm(ldeg,istp,id)=errlm(ldeg,istp,id-1)
     &                         -ca*errlm(ldeg+1,istp,id-1)
     &                         -cb*errlm(ldeg-1,istp,id-1)
                ett0lm(ldeg,istp,id)=ett0lm(ldeg,istp,id-1)
     &                         -ca*ett0lm(ldeg+1,istp,id-1)
     &                         -cb*ett0lm(ldeg-1,istp,id-1)
                epp0lm(ldeg,istp,id)=epp0lm(ldeg,istp,id-1)
     &                         -ca*epp0lm(ldeg+1,istp,id-1)
     &                         -cb*epp0lm(ldeg-1,istp,id-1)
                grlm(ldeg,istp,id)=grlm(ldeg,istp,id-1)
     &                         -ca*grlm(ldeg+1,istp,id-1)
     &                         -cb*grlm(ldeg-1,istp,id-1)
                polm(ldeg,istp,id)=polm(ldeg,istp,id-1)
     &                         -ca*polm(ldeg+1,istp,id-1)
     &                         -cb*polm(ldeg-1,istp,id-1)
              enddo
            enddo
c
c           m = 1
c
            do istp=1,4,3
              utlm(0,istp,id)=0.d0
              ertlm(0,istp,id)=0.d0
              etrlm(0,istp,id)=0.d0
              ettalm(0,istp,id)=0.d0
              eppalm(0,istp,id)=0.d0
              gtlm(0,istp,id)=0.d0
            enddo
            urlm(0,3,id)=0.d0
            utlm(0,3,id)=0.d0
            uplm(0,3,id)=0.d0
            errlm(0,3,id)=0.d0
            ertlm(0,3,id)=0.d0
            erplm(0,3,id)=0.d0
            etrlm(0,3,id)=0.d0
            ett0lm(0,3,id)=0.d0
            etp0lm(0,3,id)=0.d0
            eprlm(0,3,id)=0.d0
            ept0lm(0,3,id)=0.d0
            epp0lm(0,3,id)=0.d0
            grlm(0,3,id)=0.d0
            gtlm(0,3,id)=0.d0
            gplm(0,3,id)=0.d0
            polm(0,3,id)=0.d0
c
            do ldeg=1,ldegup-id
              ca=dble(ldeg+2)/dble(2*ldeg+3)
              cb=dble(ldeg-1)/dble(2*ldeg-1)
              do istp=1,4,3
                utlm(ldeg,istp,id)=utlm(ldeg,istp,id-1)
     &                         -ca*utlm(ldeg+1,istp,id-1)
     &                         -cb*utlm(ldeg-1,istp,id-1)
                ertlm(ldeg,istp,id)=ertlm(ldeg,istp,id-1)
     &                           -ca*ertlm(ldeg+1,istp,id-1)
     &                           -cb*ertlm(ldeg-1,istp,id-1)
                etrlm(ldeg,istp,id)=etrlm(ldeg,istp,id-1)
     &                           -ca*etrlm(ldeg+1,istp,id-1)
     &                           -cb*etrlm(ldeg-1,istp,id-1)
                ettalm(ldeg,istp,id)=ettalm(ldeg,istp,id-1)
     &                           -ca*ettalm(ldeg+1,istp,id-1)
     &                           -cb*ettalm(ldeg-1,istp,id-1)
                eppalm(ldeg,istp,id)=eppalm(ldeg,istp,id-1)
     &                           -ca*eppalm(ldeg+1,istp,id-1)
     &                           -cb*eppalm(ldeg-1,istp,id-1)
                gtlm(ldeg,istp,id)=gtlm(ldeg,istp,id-1)
     &                         -ca*gtlm(ldeg+1,istp,id-1)
     &                         -cb*gtlm(ldeg-1,istp,id-1)
              enddo
              urlm(ldeg,3,id)=urlm(ldeg,3,id-1)
     &                    -ca*urlm(ldeg+1,3,id-1)
     &                    -cb*urlm(ldeg-1,3,id-1)
              utlm(ldeg,3,id)=utlm(ldeg,3,id-1)
     &                    -ca*utlm(ldeg+1,3,id-1)
     &                    -cb*utlm(ldeg-1,3,id-1)
              uplm(ldeg,3,id)=uplm(ldeg,3,id-1)
     &                    -ca*uplm(ldeg+1,3,id-1)
     &                    -cb*uplm(ldeg-1,3,id-1)
              errlm(ldeg,3,id)=errlm(ldeg,3,id-1)
     &                     -ca*errlm(ldeg+1,3,id-1)
     &                     -cb*errlm(ldeg-1,3,id-1)
              ertlm(ldeg,3,id)=ertlm(ldeg,3,id-1)
     &                     -ca*ertlm(ldeg+1,3,id-1)
     &                     -cb*ertlm(ldeg-1,3,id-1)
              erplm(ldeg,3,id)=erplm(ldeg,3,id-1)
     &                     -ca*erplm(ldeg+1,3,id-1)
     &                     -cb*erplm(ldeg-1,3,id-1)
              etrlm(ldeg,3,id)=etrlm(ldeg,3,id-1)
     &                     -ca*etrlm(ldeg+1,3,id-1)
     &                     -cb*etrlm(ldeg-1,3,id-1)
              ett0lm(ldeg,3,id)=ett0lm(ldeg,3,id-1)
     &                      -ca*ett0lm(ldeg+1,3,id-1)
     &                      -cb*ett0lm(ldeg-1,3,id-1)
              etp0lm(ldeg,3,id)=etp0lm(ldeg,3,id-1)
     &                      -ca*etp0lm(ldeg+1,3,id-1)
     &                      -cb*etp0lm(ldeg-1,3,id-1)
              eprlm(ldeg,3,id)=eprlm(ldeg,3,id-1)
     &                     -ca*eprlm(ldeg+1,3,id-1)
     &                     -cb*eprlm(ldeg-1,3,id-1)
              ept0lm(ldeg,3,id)=ept0lm(ldeg,3,id-1)
     &                      -ca*ept0lm(ldeg+1,3,id-1)
     &                      -cb*ept0lm(ldeg-1,3,id-1)
              epp0lm(ldeg,3,id)=epp0lm(ldeg,3,id-1)
     &                      -ca*epp0lm(ldeg+1,3,id-1)
     &                      -cb*epp0lm(ldeg-1,3,id-1)
              grlm(ldeg,3,id)=grlm(ldeg,3,id-1)
     &                    -ca*grlm(ldeg+1,3,id-1)
     &                    -cb*grlm(ldeg-1,3,id-1)
              gtlm(ldeg,3,id)=gtlm(ldeg,3,id-1)
     &                    -ca*gtlm(ldeg+1,3,id-1)
     &                    -cb*gtlm(ldeg-1,3,id-1)
              gplm(ldeg,3,id)=gplm(ldeg,3,id-1)
     &                    -ca*gplm(ldeg+1,3,id-1)
     &                    -cb*gplm(ldeg-1,3,id-1)
              polm(ldeg,3,id)=polm(ldeg,3,id-1)
     &                    -ca*polm(ldeg+1,3,id-1)
     &                    -cb*polm(ldeg-1,3,id-1)
            enddo
c
c           m = 2
c
            do ldeg=0,1
              ettalm(ldeg,3,id)=0.d0
              ettblm(ldeg,3,id)=0.d0
              etpalm(ldeg,3,id)=0.d0
              etpblm(ldeg,3,id)=0.d0
              eptalm(ldeg,3,id)=0.d0
              eptblm(ldeg,3,id)=0.d0
              eppalm(ldeg,3,id)=0.d0
              eppblm(ldeg,3,id)=0.d0
c
              urlm(ldeg,2,id)=0.d0
              utlm(ldeg,2,id)=0.d0
              uplm(ldeg,2,id)=0.d0
              errlm(ldeg,2,id)=0.d0
              ertlm(ldeg,2,id)=0.d0
              erplm(ldeg,2,id)=0.d0
              etrlm(ldeg,2,id)=0.d0
              ett0lm(ldeg,2,id)=0.d0
              ettalm(ldeg,2,id)=0.d0
              ettblm(ldeg,2,id)=0.d0
              etp0lm(ldeg,2,id)=0.d0
              etpalm(ldeg,2,id)=0.d0
              eprlm(ldeg,2,id)=0.d0
              ept0lm(ldeg,2,id)=0.d0
              eptalm(ldeg,2,id)=0.d0
              epp0lm(ldeg,2,id)=0.d0
              eppalm(ldeg,2,id)=0.d0
              eppblm(ldeg,2,id)=0.d0
              grlm(ldeg,2,id)=0.d0
              gtlm(ldeg,2,id)=0.d0
              gplm(ldeg,2,id)=0.d0
              polm(ldeg,2,id)=0.d0
            enddo
c
            do ldeg=2,ldegup-id
              ca=dble(ldeg+3)/dble(2*ldeg+3)
              cb=dble(ldeg-2)/dble(2*ldeg-1)
              ettalm(ldeg,3,id)=ettalm(ldeg,3,id-1)
     &                      -ca*ettalm(ldeg+1,3,id-1)
     &                      -cb*ettalm(ldeg-1,3,id-1)
              ettblm(ldeg,3,id)=ettblm(ldeg,3,id-1)
     &                      -ca*ettblm(ldeg+1,3,id-1)
     &                      -cb*ettblm(ldeg-1,3,id-1)
              etpalm(ldeg,3,id)=etpalm(ldeg,3,id-1)
     &                      -ca*etpalm(ldeg+1,3,id-1)
     &                      -cb*etpalm(ldeg-1,3,id-1)
              etpblm(ldeg,3,id)=etpblm(ldeg,3,id-1)
     &                      -ca*etpblm(ldeg+1,3,id-1)
     &                      -cb*etpblm(ldeg-1,3,id-1)
              eptalm(ldeg,3,id)=eptalm(ldeg,3,id-1)
     &                      -ca*eptalm(ldeg+1,3,id-1)
     &                      -cb*eptalm(ldeg-1,3,id-1)
              eptblm(ldeg,3,id)=eptblm(ldeg,3,id-1)
     &                      -ca*eptblm(ldeg+1,3,id-1)
     &                      -cb*eptblm(ldeg-1,3,id-1)
              eppalm(ldeg,3,id)=eppalm(ldeg,3,id-1)
     &                      -ca*eppalm(ldeg+1,3,id-1)
     &                      -cb*eppalm(ldeg-1,3,id-1)
              eppblm(ldeg,3,id)=eppblm(ldeg,3,id-1)
     &                      -ca*eppblm(ldeg+1,3,id-1)
     &                      -cb*eppblm(ldeg-1,3,id-1)
c
              urlm(ldeg,2,id)=urlm(ldeg,2,id-1)
     &                    -ca*urlm(ldeg+1,2,id-1)
     &                    -cb*urlm(ldeg-1,2,id-1)
              utlm(ldeg,2,id)=utlm(ldeg,2,id-1)
     &                    -ca*utlm(ldeg+1,2,id-1)
     &                    -cb*utlm(ldeg-1,2,id-1)
              uplm(ldeg,2,id)=uplm(ldeg,2,id-1)
     &                    -ca*uplm(ldeg+1,2,id-1)
     &                    -cb*uplm(ldeg-1,2,id-1)
              errlm(ldeg,2,id)=errlm(ldeg,2,id-1)
     &                     -ca*errlm(ldeg+1,2,id-1)
     &                     -cb*errlm(ldeg-1,2,id-1)
              ertlm(ldeg,2,id)=ertlm(ldeg,2,id-1)
     &                     -ca*ertlm(ldeg+1,2,id-1)
     &                     -cb*ertlm(ldeg-1,2,id-1)
              erplm(ldeg,2,id)=erplm(ldeg,2,id-1)
     &                     -ca*erplm(ldeg+1,2,id-1)
     &                     -cb*erplm(ldeg-1,2,id-1)
              etrlm(ldeg,2,id)=etrlm(ldeg,2,id-1)
     &                     -ca*etrlm(ldeg+1,2,id-1)
     &                     -cb*etrlm(ldeg-1,2,id-1)
              ett0lm(ldeg,2,id)=ett0lm(ldeg,2,id-1)
     &                      -ca*ett0lm(ldeg+1,2,id-1)
     &                      -cb*ett0lm(ldeg-1,2,id-1)
              ettalm(ldeg,2,id)=ettalm(ldeg,2,id-1)
     &                      -ca*ettalm(ldeg+1,2,id-1)
     &                      -cb*ettalm(ldeg-1,2,id-1)
              ettblm(ldeg,2,id)=ettblm(ldeg,2,id-1)
     &                      -ca*ettblm(ldeg+1,2,id-1)
     &                      -cb*ettblm(ldeg-1,2,id-1)
              etp0lm(ldeg,2,id)=etp0lm(ldeg,2,id-1)
     &                      -ca*etp0lm(ldeg+1,2,id-1)
     &                      -cb*etp0lm(ldeg-1,2,id-1)
              etpalm(ldeg,2,id)=etpalm(ldeg,2,id-1)
     &                      -ca*etpalm(ldeg+1,2,id-1)
     &                      -cb*etpalm(ldeg-1,2,id-1)
              eprlm(ldeg,2,id)=eprlm(ldeg,2,id-1)
     &                     -ca*eprlm(ldeg+1,2,id-1)
     &                     -cb*eprlm(ldeg-1,2,id-1)
              ept0lm(ldeg,2,id)=ept0lm(ldeg,2,id-1)
     &                      -ca*ept0lm(ldeg+1,2,id-1)
     &                      -cb*ept0lm(ldeg-1,2,id-1)
              eptalm(ldeg,2,id)=eptalm(ldeg,2,id-1)
     &                      -ca*eptalm(ldeg+1,2,id-1)
     &                      -cb*eptalm(ldeg-1,2,id-1)
              epp0lm(ldeg,2,id)=epp0lm(ldeg,2,id-1)
     &                      -ca*epp0lm(ldeg+1,2,id-1)
     &                      -cb*epp0lm(ldeg-1,2,id-1)
              eppalm(ldeg,2,id)=eppalm(ldeg,2,id-1)
     &                      -ca*eppalm(ldeg+1,2,id-1)
     &                      -cb*eppalm(ldeg-1,2,id-1)
              eppblm(ldeg,2,id)=eppblm(ldeg,2,id-1)
     &                      -ca*eppblm(ldeg+1,2,id-1)
     &                      -cb*eppblm(ldeg-1,2,id-1)
              grlm(ldeg,2,id)=grlm(ldeg,2,id-1)
     &                    -ca*grlm(ldeg+1,2,id-1)
     &                    -cb*grlm(ldeg-1,2,id-1)
              gtlm(ldeg,2,id)=gtlm(ldeg,2,id-1)
     &                    -ca*gtlm(ldeg+1,2,id-1)
     &                    -cb*gtlm(ldeg-1,2,id-1)
              gplm(ldeg,2,id)=gplm(ldeg,2,id-1)
     &                    -ca*gplm(ldeg+1,2,id-1)
     &                    -cb*gplm(ldeg-1,2,id-1)
              polm(ldeg,2,id)=polm(ldeg,2,id-1)
     &                    -ca*polm(ldeg+1,2,id-1)
     &                    -cb*polm(ldeg-1,2,id-1)
            enddo
          enddo
c
          do ir=1,nr
            id=idr(ir)
            do ldeg=0,ldeg4-id
              cp0=tap(ldeg)*plm(ldeg,0,ir)
              cp1=tap(ldeg)*plm(ldeg,1,ir)
              cp2=tap(ldeg)*plm(ldeg,2,ir)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c             explosion & clvd
c
              do istp=1,4,3
                spatur(iz,ir,istp)=spatur(iz,ir,istp)
     &                            +urlm(ldeg,istp,id)*cp0
                spatut(iz,ir,istp)=spatut(iz,ir,istp)
     &                            +utlm(ldeg,istp,id)*cp1*ssd(ir)
c
                spaterr(iz,ir,istp)=spaterr(iz,ir,istp)
     &                            +errlm(ldeg,istp,id)*cp0
                spatert(iz,ir,istp)=spatert(iz,ir,istp)
     &                            +ertlm(ldeg,istp,id)*cp1*ssd(ir)
                spatetr(iz,ir,istp)=spatetr(iz,ir,istp)
     &                            +etrlm(ldeg,istp,id)*cp1*ssd(ir)
                spatett(iz,ir,istp)=spatett(iz,ir,istp)
     &                            +ett0lm(ldeg,istp,id)*cp0
     &                            +ettalm(ldeg,istp,id)*cp1*csd(ir)
                spatepp(iz,ir,istp)=spatepp(iz,ir,istp)
     &                            +epp0lm(ldeg,istp,id)*cp0
     &                            +eppalm(ldeg,istp,id)*cp1*csd(ir)
c
                spatgr(iz,ir,istp)=spatgr(iz,ir,istp)
     &                            +grlm(ldeg,istp,id)*cp0
                spatgt(iz,ir,istp)=spatgt(iz,ir,istp)
     &                            +gtlm(ldeg,istp,id)*cp1*ssd(ir)
                spatpo(iz,ir,istp)=spatpo(iz,ir,istp)
     &                            +polm(ldeg,istp,id)*cp0
              enddo
c
c             dip-slip
c
              spatur(iz,ir,3)=spatur(iz,ir,3)
     &                       +urlm(ldeg,3,id)*cp1*ssd(ir)
              spatut(iz,ir,3)=spatut(iz,ir,3)
     &                       +utlm(ldeg,3,id)*cp1
              spatup(iz,ir,3)=spatup(iz,ir,3)
     &                       +uplm(ldeg,3,id)*cp1
c
              spaterr(iz,ir,3)=spaterr(iz,ir,3)
     &                        +errlm(ldeg,3,id)*cp1*ssd(ir)
              
              spatert(iz,ir,3)=spatert(iz,ir,3)
     &                        +ertlm(ldeg,3,id)*cp1
              spaterp(iz,ir,3)=spaterp(iz,ir,3)
     &                        +erplm(ldeg,3,id)*cp1
              spatetr(iz,ir,3)=spatetr(iz,ir,3)
     &                        +etrlm(ldeg,3,id)*cp1
              spatett(iz,ir,3)=spatett(iz,ir,3)
     &                        +ett0lm(ldeg,3,id)*cp1*ssd(ir)
     &                        +ettalm(ldeg,3,id)*cp2*ssd(ir)*csd(ir)
     &                        +ettblm(ldeg,3,id)*cp2*ssd(ir)
              spatetp(iz,ir,3)=spatetp(iz,ir,3)
     &                        +etp0lm(ldeg,3,id)*cp1*ssd(ir)
     &                        +etpalm(ldeg,3,id)*cp2*ssd(ir)
     &                        +etpblm(ldeg,3,id)*cp2*ssd(ir)*csd(ir)
              spatepr(iz,ir,3)=spatepr(iz,ir,3)
     &                        +eprlm(ldeg,3,id)*cp1
              spatept(iz,ir,3)=spatept(iz,ir,3)
     &                        +ept0lm(ldeg,3,id)*cp1*ssd(ir)
     &                        +eptalm(ldeg,3,id)*cp2*ssd(ir)*csd(ir)
     &                        +eptblm(ldeg,3,id)*cp2*ssd(ir)
              spatepp(iz,ir,3)=spatepp(iz,ir,3)
     &                        +epp0lm(ldeg,3,id)*cp1*ssd(ir)
     &                        +eppalm(ldeg,3,id)*cp2*ssd(ir)*csd(ir)
     &                        +eppblm(ldeg,3,id)*cp2*ssd(ir)
c
              spatgr(iz,ir,3)=spatgr(iz,ir,3)
     &                       +grlm(ldeg,3,id)*cp1*ssd(ir)
              spatgt(iz,ir,3)=spatgt(iz,ir,3)
     &                       +gtlm(ldeg,3,id)*cp1
              spatgp(iz,ir,3)=spatgp(iz,ir,3)
     &                       +gplm(ldeg,3,id)*cp1
              spatpo(iz,ir,3)=spatpo(iz,ir,3)
     &                       +polm(ldeg,3,id)*cp1*ssd(ir)
c
c             strike-slip
c
              spatur(iz,ir,2)=spatur(iz,ir,2)
     &                       +urlm(ldeg,2,id)*cp2*ssd(ir)**2
              spatut(iz,ir,2)=spatut(iz,ir,2)
     &                       +utlm(ldeg,2,id)*cp2*ssd(ir)
              spatup(iz,ir,2)=spatup(iz,ir,2)
     &                       +uplm(ldeg,2,id)*cp2*ssd(ir)
c
              spaterr(iz,ir,2)=spaterr(iz,ir,2)
     &                       +errlm(ldeg,2,id)*cp2*ssd(ir)**2
              spatert(iz,ir,2)=spatert(iz,ir,2)
     &                       +ertlm(ldeg,2,id)*cp2*ssd(ir)
              spaterp(iz,ir,2)=spaterp(iz,ir,2)
     &                       +erplm(ldeg,2,id)*cp2*ssd(ir)

              spatetr(iz,ir,2)=spatetr(iz,ir,2)
     &                       +etrlm(ldeg,2,id)*cp2*ssd(ir)
              spatett(iz,ir,2)=spatett(iz,ir,2)
     &                       +ett0lm(ldeg,2,id)*cp2*ssd(ir)**2
     &                       +ettalm(ldeg,2,id)*cp2*csd(ir)
     &                       +ettblm(ldeg,2,id)*cp2
              spatetp(iz,ir,2)=spatetp(iz,ir,2)
     &                       +etp0lm(ldeg,2,id)*cp2
     &                       +etpalm(ldeg,2,id)*cp2*csd(ir)
              spatepr(iz,ir,2)=spatepr(iz,ir,2)
     &                       +eprlm(ldeg,2,id)*cp2*ssd(ir)
              spatept(iz,ir,2)=spatept(iz,ir,2)
     &                       +ept0lm(ldeg,2,id)*cp2*csd(ir)
     &                       +eptalm(ldeg,2,id)*cp2
              spatepp(iz,ir,2)=spatepp(iz,ir,2)
     &                       +epp0lm(ldeg,2,id)*cp2*ssd(ir)**2
     &                       +eppalm(ldeg,2,id)*cp2*csd(ir)
     &                       +eppblm(ldeg,2,id)*cp2
c
              spatgr(iz,ir,2)=spatgr(iz,ir,2)
     &                       +grlm(ldeg,2,id)*cp2*ssd(ir)**2
              spatgt(iz,ir,2)=spatgt(iz,ir,2)
     &                       +gtlm(ldeg,2,id)*cp2*ssd(ir)
              spatgp(iz,ir,2)=spatgp(iz,ir,2)
     &                       +gplm(ldeg,2,id)*cp2*ssd(ir)
              spatpo(iz,ir,2)=spatpo(iz,ir,2)
     &                       +polm(ldeg,2,id)*cp2*ssd(ir)**2
            enddo
c
            if(id.gt.0)then
              fac=1.d0/ssf(ir)**id
              do istp=1,4
                spatur(iz,ir,istp)=spatur(iz,ir,istp)*fac
                spatut(iz,ir,istp)=spatut(iz,ir,istp)*fac
c
                spaterr(iz,ir,istp)=spaterr(iz,ir,istp)*fac
                spatert(iz,ir,istp)=spatert(iz,ir,istp)*fac
                spatetr(iz,ir,istp)=spatetr(iz,ir,istp)*fac
                spatett(iz,ir,istp)=spatett(iz,ir,istp)*fac
                spatepp(iz,ir,istp)=spatepp(iz,ir,istp)*fac
c
                spatgr(iz,ir,istp)=spatgr(iz,ir,istp)*fac
                spatgt(iz,ir,istp)=spatgt(iz,ir,istp)*fac
                spatpo(iz,ir,istp)=spatpo(iz,ir,istp)*fac
              enddo
              do istp=2,3
                spatup(iz,ir,istp)=spatup(iz,ir,istp)*fac
c
                spaterp(iz,ir,istp)=spaterp(iz,ir,istp)*fac
                spatetp(iz,ir,istp)=spatetp(iz,ir,istp)*fac
                spatepr(iz,ir,istp)=spatepr(iz,ir,istp)*fac
                spatept(iz,ir,istp)=spatept(iz,ir,istp)*fac
c
                spatgp(iz,ir,istp)=spatgp(iz,ir,istp)*fac
              enddo
            endif
          enddo
        enddo
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)
        close(26)
        close(27)
        close(28)
        write(*,'(i6,a)')nz,' depth dependent spectra read from '
     &                    //grnfile(ig)(1:40)
        do ir=1,nr
c
c---------explosion source
c
          write(30)(spatur(iz,ir,1),iz=1,nz)
          write(30)(spatut(iz,ir,1),iz=1,nz)
          write(30)(spaterr(iz,ir,1),iz=1,nz)
          write(30)(spatert(iz,ir,1),iz=1,nz)
          write(30)(spatetr(iz,ir,1),iz=1,nz)
          write(30)(spatett(iz,ir,1),iz=1,nz)
          write(30)(spatepp(iz,ir,1),iz=1,nz)
          write(30)(spatgr(iz,ir,1),iz=1,nz)
          write(30)(spatgt(iz,ir,1),iz=1,nz)
          write(30)(spatpo(iz,ir,1),iz=1,nz)
c
c---------strike-slip source
c
          write(30)(spatur(iz,ir,2),iz=1,nz)
          write(30)(spatut(iz,ir,2),iz=1,nz)
          write(30)(spatup(iz,ir,2),iz=1,nz)
          write(30)(spaterr(iz,ir,2),iz=1,nz)
          write(30)(spatert(iz,ir,2),iz=1,nz)
          write(30)(spaterp(iz,ir,2),iz=1,nz)
          write(30)(spatetr(iz,ir,2),iz=1,nz)
          write(30)(spatett(iz,ir,2),iz=1,nz)
          write(30)(spatetp(iz,ir,2),iz=1,nz)
          write(30)(spatepr(iz,ir,2),iz=1,nz)
          write(30)(spatept(iz,ir,2),iz=1,nz)
          write(30)(spatepp(iz,ir,2),iz=1,nz)
          write(30)(spatgr(iz,ir,2),iz=1,nz)
          write(30)(spatgt(iz,ir,2),iz=1,nz)
          write(30)(spatgp(iz,ir,2),iz=1,nz)
          write(30)(spatpo(iz,ir,2),iz=1,nz)
c
c---------dip-slip source
c
          write(30)(spatur(iz,ir,3),iz=1,nz)
          write(30)(spatut(iz,ir,3),iz=1,nz)
          write(30)(spatup(iz,ir,3),iz=1,nz)
          write(30)(spaterr(iz,ir,3),iz=1,nz)
          write(30)(spatert(iz,ir,3),iz=1,nz)
          write(30)(spaterp(iz,ir,3),iz=1,nz)
          write(30)(spatetr(iz,ir,3),iz=1,nz)
          write(30)(spatett(iz,ir,3),iz=1,nz)
          write(30)(spatetp(iz,ir,3),iz=1,nz)
          write(30)(spatepr(iz,ir,3),iz=1,nz)
          write(30)(spatept(iz,ir,3),iz=1,nz)
          write(30)(spatepp(iz,ir,3),iz=1,nz)
          write(30)(spatgr(iz,ir,3),iz=1,nz)
          write(30)(spatgt(iz,ir,3),iz=1,nz)
          write(30)(spatgp(iz,ir,3),iz=1,nz)
          write(30)(spatpo(iz,ir,3),iz=1,nz)
c
c---------cldv source
c
          write(30)(spatur(iz,ir,4),iz=1,nz)
          write(30)(spatut(iz,ir,4),iz=1,nz)
          write(30)(spaterr(iz,ir,4),iz=1,nz)
          write(30)(spatert(iz,ir,4),iz=1,nz)
          write(30)(spatetr(iz,ir,4),iz=1,nz)
          write(30)(spatett(iz,ir,4),iz=1,nz)
          write(30)(spatepp(iz,ir,4),iz=1,nz)
          write(30)(spatgr(iz,ir,4),iz=1,nz)
          write(30)(spatgt(iz,ir,4),iz=1,nz)
          write(30)(spatpo(iz,ir,4),iz=1,nz)
        enddo
        close(30)
      enddo
      return
      end