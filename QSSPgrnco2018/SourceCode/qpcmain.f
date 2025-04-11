      program qpcmain
      use qpcalloc
      implicit none
c
c     work space
c
      integer*4 ig,ierr,runtime
      integer*4 time
c
c     read input file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#      QQQ    SSSS   SSSS  PPPP    GGGG   CCCC       #'
      print *,'#     Q   Q  S      S      P   P  G      C           #'
      print *,'#     Q   Q   SSS    SSS   PPPP   G GGG  C           #'
      print *,'#     Q  QQ      S      S  P      G   G  C           #'
      print *,'#      QQQQ  SSSS   SSSS   P       GGGG   CCCC       #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                   to generate                      #'
      print *,'#             Green Function database for            #'
      print *,'#        coseismic displacement, strain/stress,      #'
      print *,'#  ground rotation, and gravity/geopotential changes #'
      print *,'#                      based on                      #'
      print *,'#   a spherically symmetric, elastic earth model     #'
      print *,'#                                                    #'
      print *,'#                  (Version 2018)                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#             Last modified: May 2018                #'
      print *,'######################################################'
      print *,'                                                      '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      open(10,file=inputfile,status='old')
      call qpcgetinp(10)
      close(10)
c
      call qpcsublayer(ierr)
      call qpcdocu(ierr)
c
      do ig=1,ngrn
        if(grnsel(ig).eq.1)then
          lys=lygrn(ig)
          lys1=lygrn1(ig)
          lys2=lygrn2(ig)
          call qpcgrnspec(ig)
        endif
      enddo
      call qpcwvint(ierr)
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #    End of computations with qsspgrnco     #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
      stop
      end
