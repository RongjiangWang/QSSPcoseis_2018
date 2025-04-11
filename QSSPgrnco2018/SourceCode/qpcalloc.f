      module qpcalloc
c
c     CONSTANTS
c     =========
      real*8 PI
      parameter(PI=3.14159265358979d0)
      real*8 DEG2RAD,KM2M
      parameter(DEG2RAD=1.745329251994328d-02)
      parameter(KM2M=1.0d+03)
      real*8 BIGG
      parameter(BIGG=6.6732d-11)
      real*8 RESOLUT
      parameter(RESOLUT=0.05d0)
      real*8 FLTAPER
      parameter(FLTAPER=0.2d0)
      real*8 THICKMAX
      parameter(THICKMAX=5.0d+05)
      integer*4 ndmax
      parameter(ndmax=2)
c
c     GREEN FUNCTION PARAMETERS
c     =========================
c
      logical*2 nogravity
      integer*4 ngrn,nz,nr,lyadd,ldeggr,ldegcut,ldegmax,lymax
      character*80 inputfile,infofile,specgrntmp,specgrndir,
     &             spatgrndir
c
      integer*4, allocatable:: lygrn(:),lygrn1(:),lygrn2(:),grnsel(:)
      real*8, allocatable:: grndep(:),depr(:)
      character*80, allocatable:: spatgrnfile(:),grnfile(:),
     &             uspecfile(:),vspecfile(:),wspecfile(:),especfile(:),
     &             fspecfile(:),gspecfile(:),pspecfile(:),qspecfile(:)
c
      real*8 rearth,gravity,zr1,zr2,dzr,dr1,dr2,ddr1,ddr2
      real*8 gdds0,gddsnorm,mues,ksis,kaps
c
      integer*4, allocatable:: lyr(:)
      real*8, allocatable:: gdds(:),grrs(:)
      real*8, allocatable:: ul0(:,:),vl0(:,:),wl0(:,:),
     &       el0(:,:),fl0(:,:),gl0(:,:),pl0(:,:),ql0(:,:)
c
      real*8, allocatable:: urlm(:,:,:),utlm(:,:,:),uplm(:,:,:),
     &       grlm(:,:,:),gtlm(:,:,:),gplm(:,:,:),polm(:,:,:),
     &       errlm(:,:,:),ertlm(:,:,:),erplm(:,:,:),etrlm(:,:,:),
     &       ett0lm(:,:,:),ettalm(:,:,:),ettblm(:,:,:),etp0lm(:,:,:),
     &       etpalm(:,:,:),etpblm(:,:,:),eprlm(:,:,:),ept0lm(:,:,:),
     &       eptalm(:,:,:),eptblm(:,:,:),epp0lm(:,:,:),eppalm(:,:,:),
     &       eppblm(:,:,:)
      real*8, allocatable:: spatur(:,:,:),spatut(:,:,:),spatup(:,:,:),
     &       spatgr(:,:,:),spatgt(:,:,:),spatgp(:,:,:),spatpo(:,:,:),
     &       spaterr(:,:,:),spatert(:,:,:),spaterp(:,:,:),
     &       spatetr(:,:,:),spatett(:,:,:),spatetp(:,:,:),
     &       spatepr(:,:,:),spatept(:,:,:),spatepp(:,:,:)
c
c     SOURCE-RECEIVER CONFIGURATION PARAMETERS
c     ========================================
c
      integer*4, allocatable:: idr(:)
      real*8, allocatable:: dism(:),disrad(:),ssd(:),csd(:),ssf(:)
c
c     LEGENDRE POLYNOMIALS TABLES
c     ===========================
c
c     plm = associated Legendre polynomials divided by sin(x)**m
c
      real*8, allocatable:: plm(:,:,:),tap(:),disk(:)
c
c     ORIGINAL MODEL PARAMETERS
c     =========================
c
c     lys = layer index of source
c     lycm = layer index of core-mantle boundary
c     lycc = layer index of inner and outer core boundary
c     ly0 = max. layer number for integration
c     lyuppsv = upper starting layer index for p-sv solution
c     lyupt = upper starting layer index for sh solution
c     lylwpsv = lower starting layer index for p-sv solution
c     lylwt = lower starting layer index for sh solution
c
      integer*4 lys,lys1,lys2,lycm,lycc,ly0
      integer*4, allocatable:: lylwpsv(:),lylwt(:),lyuppsv(:),lyupt(:)
c
c     dp = depth (up = top of layer, lw = bottom of layer)
c     vp, vs = p and s wave velocity
c     rho = density
c
      integer*4 l0
c
      real*8, allocatable:: dp0(:),kap0(:),mue0(:),rho0(:)
      real*8, allocatable:: dp0up(:),kap0up(:),mue0up(:),rho0up(:)
      real*8, allocatable:: dp0lw(:),kap0lw(:),mue0lw(:),rho0lw(:)
c
      integer*4, allocatable:: izrly(:),ldegpsv(:),ldegsh(:),
     &                         ndruku(:,:)
      real*8, allocatable:: rrup(:),rrlw(:),kapup(:),kaplw(:),
     &        mueup(:),muelw(:),rhoup(:),rholw(:),grup(:),grlw(:)
c
c     LAYER MATRICES
c     ==============
c
c     mat2x2 = 2x2 toroidal solution matrix
c     mas3x3 = 3x3 spheroidal solution matrix (l = 0)
c     mas4x4 = 4x4 spheroidal solution matrix (l > 0, in liquid)
c     mas6x6 = 6x6 spheroidal solution matrix (l > 0, in solid)
c     mas(t)inv = inverse solution matrix
c
      real*8, allocatable:: mat2x2(:,:,:),mat2x2inv(:,:,:),
     &                      mas3x3(:,:,:),mas3x3inv(:,:,:),
     &                      mas4x4(:,:,:),mas4x4inv(:,:,:),
     &                      mas6x6(:,:,:),mas6x6inv(:,:,:),
     &                      cypnorm(:,:),ypsv(:,:,:),ypsvg(:,:,:),
     &                      ysh(:,:,:)
c
      real*8, allocatable:: yup2(:,:),ylw2(:,:),yup3(:,:),ylw3(:,:),
     &                      yup6(:,:,:),ylw6(:,:,:)
c
      end module
