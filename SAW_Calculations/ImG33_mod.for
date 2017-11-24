c     last change:  a     1 sep 2015   11:36 am
c           pangim
c*********************************************
c this coding has been developed by a.g. every
c*********************************************
c     program for calculating half space green's function
c     for anisotropic solid
c********************************************
      implicit real*8 (a-h,o-z)
      real*8 atr(3,3),del(3,3),eps(3,3,3)
      real*8 cl(3,3,3,3)
      real*8 bl(0:6,0:6,0:6),cml(6,6)
      integer ipar(2)
      common /cipar/ ipar
      common /cdep/ del,eps
c===============================
      rt13=dsqrt(1.0d0/3.0d0)
      rt16=dsqrt(1.0d0/6.0d0)
      pi=3.141592653589793d0
      twpi=2.0d0*pi
      rt12=dsqrt(1.0d0/2.0d0)
c===============================
      do 543 i=1,6
      do 543 j=1,6
      cml(i,j)=0.0d0
543   continue
c ***************************
c     orientation of lab frame wrt crystallographic axes.
c  3 axis normal to crystal surface
c**********************************
c     111 crystal
      atr(1,1)=rt16
      atr(1,2)=rt16
      atr(1,3)=-2.0d0*rt16
      atr(2,1)=-rt12
      atr(2,2)=rt12
      atr(2,3)=0.0d0
      atr(3,1)=rt13
      atr(3,2)=rt13
      atr(3,3)=rt13
c ***********************************************
c     100 crystal, x1 axis at angle pph (in degrees) to 100 direction
c     x2 axis at angle ppp to 010 direction
c     x3 axis in 001 direction
c      pph=0.00001d0*pi/180.0d0
c      cpph=dcos(pph)
c      spph=dsin(pph)
c      atr(1,1)=cpph
c      atr(1,2)=spph
c      atr(1,3)=0.0d0
c      atr(2,1)=-spph
c      atr(2,2)=cpph
c      atr(2,3)=0.0d0
c      atr(3,1)=0.0d0
c      atr(3,2)=0.0d0
c      atr(3,3)=1.0d0
c     crystal with axes rotated by phi about x3 and then eta about x2
c      eta=0.00001d0*pi/180.0d0
c     phi=45.0001d0*pi/180.0d0
c      phi=0.00001d0*pi/180.0d0
c      ceta=dcos(eta)
c      seta=dsin(eta)
c      cphi=dcos(phi)
c      sphi=dsin(phi)
c      atr(1,1)=cphi*ceta
c      atr(1,2)=sphi*ceta
c      atr(1,3)=seta
c      atr(2,1)=-sphi
c      atr(2,2)=cphi
c      atr(2,3)=0.0d0
c      atr(3,1)=-cphi*seta
c      atr(3,2)=-sphi*seta
c      atr(3,3)=ceta
c ***************************
c     110 crystal, x1 axis in-110 direction
c     x2 axis 001 direction
c     x3 axis in 110 direction
c      atr(1,1)=-rt12
c      atr(1,2)=rt12
c      atr(1,3)=0.0d0
c      atr(2,1)=0.0d0
c      atr(2,2)=0.0d0
c      atr(2,3)=1.0d0
c      atr(3,1)=rt12
c      atr(3,2)=rt12
c      atr(3,3)=0.0d0
c ***************************
c     materials constants
c ***************************
c     elastic constants and density of substrate
c      silicon cubic symmetry.
c      rhol=2.332d0
c      cml(1,1)=165.0d0
c      cml(1,2)=63.0d0
c      cml(4,4)=79.1d0
c======================
c     elastic constants and density of al w/ cubic symmetry
c     sourced from vallin (1964) - room temp.
c      rhol=2.7d0
c      cml(1,1)=116.3d0
c      cml(1,2)=64.8d0
c      cml(4,4)=30.9d0
C     New set of elastic constants from a website that says they're from a book by Simmons and Wang (1970)
c      cml(1,1)=107.30d0
c      cml(1,2)=60.90d0
c      cml(4,4)=28.30d0
c======================
C     elastic constants and density of Nb with cubic symmetry
c     sourced from Carroll, JAP (1965)
c      rhol=8.5605d0;
c      cml(1,1)=245.6d0;
c      cml(1,2)=138.7d0; 
c     cml(4,4)=29.30d0;
c======================
c     elastic constants from MD Al with cubic symmetry. Taken from Winey (2009) 325 K
c      rhol=2.7d0
c      cml(1,1)=107.7d0
c      cml(1,2)=57.9d0
c      cml(4,4)=26.2d0
c======================
c     elastic constants and density of cu w/ cubic symmetry
c     sourced from ledbetter (1974) - room temp.
c      rhol=8.96d0
c      cml(1,1)=169.1d0
c      cml(1,2)=122.2d0
c      cml(4,4)=75.42d0
c density with a 12.8% reduction from nominal based on image processing over an area a grating wavelength deep to match 90 dpa data
c BASED ON OVERTON starting point, not LEDBETTER
      rhol=7.81312d0
      cml(1,1)=142.6740d0
      cml(1,2)=102.8560d0
      cml(4,4)=52.9272d0
c density with a 0.79% reduction to match the 5 dpa case
c      rhol=8.889216d0
c      cml(1,1)=165.9d0
c      cml(1,2)=119.6d0
c      cml(4,4)=73.51d0
c     second set from overton (1954) for 300 K corresponds to 'new' saved .txt files
c     cml(1,1)=168.4d0
c      cml(1,2)=121.4d0
c      cml(4,4)=75.39d0
c Base density, optimizing for unirradiated case starting from Overton data set
c      cml(1,1)=165.8740d0
c      cml(1,2)=119.5790d0
c      cml(4,4)=73.5053d0
c======================
c     elastic constants from MD Cu with cubic symmetry. Taken from Mishin (2001) 300 K
c      rhol=8.96d0
c      cml(1,1)=169.9d0
c      cml(1,2)=122.6d0
c      cml(4,4)=76.2d0
c======================
c     elastic constants and density of W w/ cubic symmetry
c     sourced from my Mathematica notebook, which includes references
c      rhol=19.25d0
c      cml(1,1)=522.796d0
c      cml(1,2)=203.468d0
c      cml(4,4)=160.668d0
c======================
c     elastic constants and density of Mo w/ cubic symmetry
c     see source in elastic_constants.xlsx
c      rhol=19.3d0
c      cml(1,1)=440.77d0
c      cml(1,2)=172.43d0
c      cml(4,4)=121.65d0
c======================
c     elastic constants and density of \alpha Fe w/ cubic symmetry
c     see source in elastic_constants.xlsx
c      rhol=7.874d0
c      cml(1,1)=233.1d0
c      cml(1,2)=135.44d0
c      cml(4,4)=117.83d0
c======================
c     elastic constants and density of Ni w/ cubic symmetry
c     see source in elastic_constants.xlsx
c      rhol=8.9d0
c      cml(1,1)=252d0
c      cml(1,2)=152d0
c      cml(4,4)=123d0
c======================
c     elastic constants and density of Au w/ cubic symmetry
c     see source in elastic_constants.xlsx
c      rhol=19.3d0
c      cml(1,1)=192.34d0
c      cml(1,2)=163.14d0
c      cml(4,4)=41.95d0
c======================
c     elastic constants and density of Ag w/ cubic symmetry
c     see sources in elastic_constants.xlsx
c      rhol=10.49d0
c      cml(1,1)=123.99d0
c      cml(1,2)=93.67d0
c      cml(4,4)=46.12d0
c======================
c     elastic constants and density of Pb w/ cubic symmetry
c     see sources in elastic_constants.xlsx
c      rhol=11.344d0
c      cml(1,1)=49.66d0
c      cml(1,2)=42.31d0
c      cml(4,4)=14.977d0
c======================
c     elastic constants of SrTiO3, strontium titanate from Beattie (1971) at 295 K
c      rhol=4.81d0
c      cml(1,1)=317.6d0
c      cml(1,2)=102.5d0
c      cml(4,4)=123.5d0
c======================
c======================
c     Constants from different cases from Felix's first paper in 2015. Values taken from Table 2
c     - Density is constant for all cases
c      rhol=19.26d0
c     Elastic constant values for base case
c      cml(1,1)=522.8d0
c      cml(1,2)=203.5d0
c      cml(4,4)=160.7d0
c     Elastic constant values for SIAs
c      cml(1,1)=498.6d0
c      cml(1,2)=229.8d0
c      cml(4,4)=147.8d0
c     Elastic constant values for HeV clusters
c      cml(1,1)=504.8d0
c      cml(1,2)=203.4d0
c      cml(4,4)=147.6d0
c     Elastic constant values for both defect types
c      cml(1,1)=514.4d0
c      cml(1,2)=208.7d0
c      cml(4,4)=155.5d0
c======================
c======================
c for an isotropic solid, avoid perfectly isotropic elastic constants
c otherwise run into numerical instabilities
c      polymer
c      rhol=0.92d0
c      cml(1,1)=3.50d0
c      cml(4,4)=0.27d0
c      cml(1,2)=(cml(1,1)-2.0d0*cml(4,4))*0.999d0
c=======================
      cml(1,3)=cml(1,2)
      cml(2,2)=cml(1,1)
      cml(2,3)=cml(1,2)
      cml(3,3)=cml(1,1)
      cml(5,5)=cml(4,4)
      cml(6,6)=cml(4,4)
      cml(3,1)=cml(1,2)
      cml(2,1)=cml(1,2)
      cml(3,2)=cml(1,2)
c==================
c   this section easily modified for lower symmetry crystal. e.g. for gan
c      rhol=6.0950d0
c      cml(1,1)=384.54d0
c      cml(2,2)=384.54d0
c      cml(3,3)=398.0d0
c      cml(4,4)=94.660d0
c      cml(5,5)=94.660d0      
c      cml(6,6)=119.61d0       
c      cml(6,6)=119.88d0     
c      cml(1,2)=144.78d0
c      cml(1,3)=114.0d0
c      cml(2,3)=114.0d0
c      cml(3,2)=cml(2,3)
c      cml(2,1)=cml(1,2)
c      cml(3,1)=cml(1,3)
c     ==============================
      do 987 i=1,6
      do 987 j=1,6
      cml(i,j)=cml(i,j)/rhol
 987  continue
c     ===================================
c     elastic constants originally in gpa and rho in g/cm**3
c     c(i,j)are now elastic constants divided by density
c     velocities come out in mm/microsec. thick in microns
c    =======================================================
c     i and j components of green's function gij
      ipar(1)=3
      ipar(2)=3
c     ===============================
c     initialising various arrays
      call sinit(atr,bl,cml,cl)
c     ==============================
c     main calculation.
      call modeconv(bl,cl,rhol)
c     ==================================
      stop
      end


      subroutine modeconv(bl,cl,rhol)
c     ***************************************************************
c     a sequence of values of s// are generated and the equation
c     of the slowness surface is then solved to determine the
c     corresponding values of sz.
c     ***************************************************************
      implicit real*8 (a-h,o-z)
      real*8 cl(3,3,3,3),del(3,3),eps(3,3,3)
      complex*16 s1p(0:6),s2p(0:6)   
      real*8 bl(0:6,0:6,0:6)
      complex*16 xii,s1,s2,hnw,gdet,amp1
      integer ipar(2),gray(900)
      integer*2 intns(900,900)
      common /cipar/ ipar
      common /cdep/ del,eps
c  ==============================
      xii=(0.0d0,1.0d0)
      s1p(0)=(1.0d0,0.0d0)
      s2p(0)=(1.0d0,0.0d0)
      pi=3.141592653589793d0
      twpi=2.0d0*pi
c      open(unit=10,file="g1im.dat",status="new")
      open(unit=10,file="g1im.txt",status="new")
c
c the rayleigh wave is a delta function which is too narrow to show up in plots
c  of img33, so a small amount of damping is introduced below. too small and  you
c won't see rw, too large and the rw will be too broadened
      damping=0.001d0
c the range of s1 and s2 set suitably. usually broad enough initially to capture
c whole plot of img33, and then narrower range to get accurate rayleigh slowness
c      s1min=-0.25d0
      s1min=0.0d0
c      s1max=0.25d0
      s1max=0.75d0
      s1del=(s1max-s1min)/900.0d0
c      s2min=-0.25d0
      s2min=0.00d0      
c      s2max=0.25d0
      s2max=0.75d0
      s2del=(s2max-s2min)/900.0d0
c      smax=0.5d0
ccc      rrr=rrand()

c      bsimax=0.0000001d0
c img33 is now calculated for a 900x900 array of s1,s2 values and output array in file g1im.dat
c i use the plotting package origin to generate the image.
      do 750 jkd=1,900
      write (6,*)jkd
      do 75 i=1,900
      s1=s1min+(dfloat(jkd)-0.5d0)*s1del
      s2=s2min+(dfloat(i)-0.5d0)*s2del
c      write (6,*)jkd,cdabs(s1)
c     ==================
      call hval(hnw,bl,cl,rhol,s1,s2,s1p,s2p,gdet,amp1)
c     ==================
c    re(hnw)=im(g), im(hnw)=re(g). so hnw=ig* for fixed frequency omega!!!
      bsi=dimag(xii/(1.0d0/hnw+damping))
      if (bsi.gt.bsimax) then
      bsimax=bsi
      end if
c  you will have to try different values of numerical factor below for best results
      bsi=bsi*1000.0d0
      nadd=int(bsi)
      frac=bsi-nadd
c      rrr=rnd()
      rrr = 0.5
      if(rrr.lt.frac)nadd=nadd+1
      if(nadd.gt.98)nadd=99
      if(nadd.lt.1)nadd=1
      intns(i,jkd)=nadd
75    continue
750   continue
       do 100 j=1,900
       do 400 i=1,900
       gray(i)=100-intns(j,i)
 400   continue
      write (10,60) gray
100   continue
60    format(900i3)
999   format(1x,f12.4,f12.4)
      return
      end

      subroutine hval(hnw,bl,cl,rhol,s1,s2,s1p,s2p,gdet,amp1)
c*************************      
c     calculation of hnw
c*************************
      implicit real*8 (a-h,o-z)
      real*8 cu(3,3,3,3),cl(3,3,3,3),del(3,3),eps(3,3,3)
      real*8 bu(0:6,0:6,0:6)       
      real*8 bl(0:6,0:6,0:6) 
      complex*16 s1p(0:6),s2p(0:6)   
      complex*16 e3u(0:6),e3l(0:6)
      complex*16 au(7),acu(4)       
      complex*16 al(7),acl(4)
      complex*16 aau(7),aal(7)
      complex*16 s3u(6),ucu(6,3),s3cu(3),ucmu(6,3)       
      complex*16 s3l(6),ucl(6,3),s3cl(3),ucml(3,3)  
      complex*16 h(6),smu(6,3),sml(3,3),hnw,gdet
      complex*16 s1,s2,amp1
      integer ipar(2)
      integer ivnl(3)
      logical polish
      common /cipar/ ipar
      common /cdep/ del,eps
c  ==============================
      polish=.true.
      xii=(0.0d0,1.0d0)
      pi=3.141592653589793d0
      twpi=2.0d0*pi
c  =================================================
c     calculation of powers of s1 and s2 and of array a
      do 570 k=1,6
      s1p(k)=s1p(k-1)*s1
      s2p(k)=s2p(k-1)*s2
570   continue
      do 571 k=0,6
      e3l(k)=0.0d0
      do 572 l=0,6,2
      if(k.gt.l)go to 573
      do 575 m=0,l-k
      n=l-m-k
      bmnk=bl(m,n,k)
      if(bmnk.eq.0.0d0)go to 575
      e3l(k)=e3l(k)+bl(m,n,k)*s1p(m)*s2p(n)            
575   continue
573   continue
572   continue
      al(7-k)=e3l(k)/bl(0,0,6)
      aal(k+1)=e3l(k)/bl(0,0,6)
571   continue
c      write (6,*)aal(2)
c     for sz perpendicular to a mirror
c     plane, equation of degree 3 in sz**2
c      acu(4)=au(1)
c      acu(3)=au(3)
c      acu(2)=au(5)
c      acu(1)=au(7)
c      acl(4)=al(1)
c      acl(3)=al(3)
c      acl(2)=al(5)
c      acl(1)=al(7)
c     solving for sz
      call zroots(aal,6,s3l,polish)
c      write (6,*)s3l
c      write (6,*)s3(1)
c      call cubrt(acu,s3cu)
c      s3u(1)=cdsqrt(s3cu(1))
c      s3u(2)=-s3u(1)
c      s3u(3)=cdsqrt(s3cu(2))
c      s3u(4)=-s3u(3)
c      s3u(5)=cdsqrt(s3cu(3))
c      s3u(6)=-s3u(5)
c      call cubrt(acl,s3cl)
c      s3l(1)=cdsqrt(s3cl(1))
c      s3l(2)=-s3l(1)
c      s3l(3)=cdsqrt(s3cl(2))
c      s3l(4)=-s3l(3)
c      s3l(5)=cdsqrt(s3cl(3))
c     s3l(6)=-s3l(5)
123   continue
c     ===============================================================
c     calculation of eigenvectors uc(ibr,j) for 6 solutions
c      call sub2(cu,ucu,s1,s2,s3u)
      call sub2(cl,ucl,s1,s2,s3l)   
c     ===============================================================
c     calculation of number nrw of outgoing real waves and ncw of
c     outgoing complex waves. total number of outgoing waves ivn(1-3)
      call sub3(1,bl,ivnl,s3l,e3l,s1p,s2p,s1,s2)      
c     ===============================================================
c      do 774 ibr=1,6
c      smu(ibr,1)=s1
c      smu(ibr,2)=s2
c      smu(ibr,3)=s3u(ibr)
c      do 7741 jbr=1,3
c      ucmu(ibr,jbr)=ucu(ibr,jbr)
c7741  continue
c774   continue
      do 773 ibr=1,3
      sml(ibr,1)=s1
      sml(ibr,2)=s2
      sml(ibr,3)=s3l(ivnl(ibr))
      do 7731 jbr=1,3
      ucml(ibr,jbr)=ucl(ivnl(ibr),jbr)
7731  continue
773   continue
c     ============================================================
c     calculation of g(ix,iy),determinants,amplitudes
c     of outgoing waves, h(i)
      call sub4(cl,rhol,ucml,sml,h,gdet,amp1)
c     ============================================================
c      write (6,*)h
      hnw=(0.0d0,0.0d0)      
      do 44477 nw=1,3
       hnw=hnw+h(nw)
44477  continue
      return
      end


c     ************************************************************
c     calculation of g(ix,iy),determinants,amplitudes 
c     of outgoing waves h(i)
      subroutine sub4(cl,rhol,ucml,sml,h,gdet,amp1)
c     ************************************************************
      implicit real*8 (a-h,o-z)
      real*8 cu(3,3,3,3),cl(3,3,3,3),gr(6,6),deltafr(6,3)
      integer ipar(2)
      complex*16 g(3,3),h(3),smu(6,3),sml(3,3),xii,gdet
      complex*16 adj31,adj32,adj33,amp1
      complex*16 ucmu(6,3),ucml(3,3)
      common /cipar/ ipar
c
      igf=ipar(1)
      jgf=ipar(2)
      xii=(0.0d0,1.0d0)
      do 210 i=1,6
      do 210 j=1,3
      deltafr(i,j)=0.0d0
210   continue
c*****************************************************
c 1,2,3 for outer surface of film, 4,5,6 for interface
      deltafr(1,1)=1.0d0        
      deltafr(2,2)=1.0d0        
      deltafr(3,3)=1.0d0
      do 212 i=1,3
      do 212 j=1,3
      g(i,j)=0.0d0
212   continue
c     ============================================
      do 777 ix=1,3
      do 777 iy=1,3
      do 776 ip=1,3
      do 776 iq=1,3
      c3xpq=cl(3,ix,ip,iq)
      if(c3xpq.eq.0.0d0)go to 776
      g(ix,iy)=g(ix,iy)
     1 -cl(3,ix,ip,iq)*ucml(iy,ip)*sml(iy,iq)*rhol
776   continue
777   continue
      adj31=g(1,2)*g(2,3)-g(2,2)*g(1,3)
      adj32=g(1,3)*g(2,1)-g(1,1)*g(2,3)
      adj33=g(1,1)*g(2,2)-g(1,2)*g(2,1)
      gdet=g(3,1)*adj31+g(3,2)*adj32+g(3,3)*adj33
c     =============================================
c     g is entered and its inverse is returned as g
c      call gaussj2(g,9,9,deltaf,3,3)      
      do 310 i=1,3
      do 310 j=1,3
      gr(i,j)=dimag(xii*g(i,j))
      gr(i,j+3)=-dimag(g(i,j))
      gr(i+3,j)=dimag(g(i,j))
      gr(i+3,j+3)=dimag(xii*g(i,j))
310   continue
      call gaussj(gr,6,6,deltafr,3,3)
      amp1=deltafr(1,jgf)+xii*deltafr(1+3,jgf)
c     =============================================      
      do 870 ii=1,3
c****************************************************************
c for outer surface
      h(ii)=-ucml(ii,igf)*(deltafr(ii,jgf)+xii*deltafr(ii+3,jgf))
870   continue
      return
      end





c     *************************************************************
c     calculation of gam and eigenvectors uc(ibr,j) for 6 solutions
      subroutine sub2(c,uc,s1,s2,s3)
c     *************************************************************
      implicit real*8 (a-h,o-z)
      real*8 c(3,3,3,3)
      complex*16 s3(6),gam(3,3),s(3),uc(6,3),s1,s2
      s(1)=s1   
      s(2)=s2
      do 670 ibr=1,6
      s(3)=s3(ibr)
      do 671 k=1,3
      do 671 l=1,3
      gam(k,l)=0.0d0
      do 672 m=1,3
      do 672 n=1,3
      ckmln=c(k,m,l,n)
      if(ckmln.eq.0.0d0)go to 672
      gam(k,l)=gam(k,l)+c(k,m,l,n)*s(m)*s(n)
672   continue
671   continue
      uc(ibr,1)=gam(1,2)*gam(2,3)-gam(1,3)*(gam(2,2)-1.0d0)
      uc(ibr,2)=gam(1,3)*gam(2,1)-gam(2,3)*(gam(1,1)-1.0d0)
      uc(ibr,3)=(gam(1,1)-1.0d0)*(gam(2,2)-1.0d0)-gam(1,2)*gam(2,1)
670   continue
      return
      end


c     ****************************************************************
c     calculation of number nrw of outgoing real waves and ncw of
c     outgoing complex waves. total number of outgoing waves=3
      subroutine sub3(isp,b,ivn,s3,e3,s1p,s2p,s1,s2)
c     ****************************************************************
      implicit real*8 (a-h,o-z)
      complex*16 s1p(0:6),s2p(0:6),s3p(0:6)
      complex*16 e1(0:6),e2(0:6),e3(0:6)
      real*8 b(0:6,0:6,0:6)
      complex*16 s3(6),xii,s1,s2,sdvg,vgmibr1,vgmibr2,vgmibr3
      integer ivn(3),ir(3),jr(3)
      nrw=0
      ncw=0
      xii=dcmplx(0.0d0,1.0d0)
      s3p(0)=dcmplx(1.0d0,0.0d0)
      do 772 ibr=1,6
      if(dabs(dimag(s3(ibr))).gt.0.00000001d0)goto 770
c
      s3(ibr)=s3(ibr)-xii*dimag(s3(ibr))
      do 575 k=1,6
      s3p(k)=s3p(k-1)*s3(ibr)
575   continue
      do 576 k=0,6
      e1(k)=0.0d0
      e2(k)=0.0d0
      do 577 l=0,6,2
      if(k.gt.l) go to 578
      do 579 m=0,l-k
      n=l-m-k
      bkmn=b(k,m,n)
      if(bkmn.eq.0.0d0)go to 5791
      e1(k)=e1(k)+bkmn*s2p(m)*s3p(n)
5791  bmkn=b(m,k,n)
      if(bmkn.eq.0.0d0)go to 579
      e2(k)=e2(k)+bmkn*s1p(m)*s3p(n)
579   continue
578   continue
577   continue
576   continue
c
      vgmibr3=6.d0*e3(6)*s3p(5)+5.d0*e3(5)*s3p(4)+4.d0*e3(4)*s3p(3)
     1          +3.d0*e3(3)*s3p(2)+2.d0*e3(2)*s3p(1)+1.d0*e3(1)*s3p(0)
      vgmibr2=6.0d0*e2(6)*s2p(5)+5.0d0*e2(5)*s2p(4)+4.0d0*e2(4)*s2p(3)
     1       +3.0d0*e2(3)*s2p(2)+2.0d0*e2(2)*s2p(1)+1.0d0*e2(1)*s2p(0)
      vgmibr1=6.0d0*e1(6)*s1p(5)+5.0d0*e1(5)*s1p(4)+4.0d0*e1(4)*s1p(3)
     1       +3.0d0*e1(3)*s1p(2)+2.0d0*e1(2)*s1p(1)+1.0d0*e1(1)*s1p(0)
      sdvg=s1*vgmibr1+s2*vgmibr2+s3(ibr)*vgmibr3
      vgmibr3=vgmibr3/sdvg
c
      if(isp.eq.1.and.dimag(xii*vgmibr3).gt.0.0d0)nrw=nrw+1
      if(isp.eq.1.and.dimag(xii*vgmibr3).gt.0.0d0)ir(nrw)=ibr
      if(isp.eq.-1.and.dimag(xii*vgmibr3).lt.0.0d0)nrw=nrw+1
      if(isp.eq.-1.and.dimag(xii*vgmibr3).lt.0.0d0)ir(nrw)=ibr
c
      go to 771
770   continue
      if(isp.eq.1.and.dimag(s3(ibr)).gt.0.0d0)ncw=ncw+1
      if(isp.eq.1.and.dimag(s3(ibr)).gt.0.0d0)jr(ncw)=ibr
      if(isp.eq.-1.and.dimag(s3(ibr)).lt.0.0d0)ncw=ncw+1
      if(isp.eq.-1.and.dimag(s3(ibr)).lt.0.0d0)jr(ncw)=ibr
771   continue
772   continue
      if(ncw.eq.0)then
      ivn(1)=ir(1)
      ivn(2)=ir(2)
      ivn(3)=ir(3)
      end if
      if(ncw.eq.1)then
      ivn(1)=ir(1)
      ivn(2)=ir(2)
      ivn(3)=jr(1)
      end if
      if(ncw.eq.2)then
      ivn(1)=ir(1)
      ivn(2)=jr(1)
      ivn(3)=jr(2)
      end if
      if(ncw.eq.3)then
      ivn(1)=jr(1)
      ivn(2)=jr(2)
      ivn(3)=jr(3)
      end if
      return
      end


      subroutine sinit(atr,b,cm,c)
c     *******************************
c
      implicit real*8 (a-h,o-z)
      real*8 c(3,3,3,3),del(3,3),eps(3,3,3)
      real*8 atr(3,3),ccc(3,3,3,3)
      real*8 a4(3,3,3,3),a6(3,3,3,3,3,3)
      real*8 b(0:6,0:6,0:6),cm(6,6)
      integer ii(9)
      common /cdep/ del,eps
c     ********************************************************
c     tensor elastic constants in crystal frame
      do 63 i=1,3
      do 63 j=1,3
      do 63 k=1,3
      do 63 l=1,3
      m=9-i-j
      if(i.eq.j)m=i
      n=9-k-l
      if(k.eq.l)n=k
      ccc(i,j,k,l)=cm(m,n)
 63   continue
c ==========================================
c     tensor elastic constants in lab frame
      do 663 i=1,3
      do 663 j=1,3
      do 663 k=1,3
      do 663 l=1,3
      c(i,j,k,l)=0.0d0
      do 664 ia=1,3
      do 664 ja=1,3
      do 664 ka=1,3
      do 664 la=1,3
      zy=atr(i,ia)*atr(j,ja)*atr(k,ka)*atr(l,la)*ccc(ia,ja,ka,la)
      c(i,j,k,l)=c(i,j,k,l)+zy
 664  continue
 663  continue
c  ==========================================================
c     kronecker delta and alternating tensors
      do 3 i=1,3
      do 3 j=1,3
      del(i,j)=0.0d0
      do 3 k=1,3
      eps(i,j,k)=0.0d0
 3    continue
      del(1,1)=1.0d0
      del(2,2)=1.0d0
      del(3,3)=1.0d0
      eps(1,2,3)=1.0d0
      eps(2,3,1)=1.0d0
      eps(3,1,2)=1.0d0
      eps(3,2,1)=-1.0d0
      eps(2,1,3)=-1.0d0
      eps(1,3,2)=-1.0d0

c     ==========================================================
c     procedure for calculating b(i,j,k), the matrix of coefficients
c     of the equation of the slowness surface
      do 31 i=1,3
      do 31 j=1,3
      do 31 k=1,3
      do 31 l=1,3
      a4(i,j,k,l)=0.0d0
      do 34 ip=1,3
      do 34 iq=1,3
      do 34 irr=1,3
      do 34 is=1,3
      do 34 it=1,3
      a4(i,j,k,l)=a4(i,j,k,l)+eps(ip,iq,it)*eps(irr,is,it)*
     1c(ip,i,irr,j)*c(iq,k,is,l)/2.0d0
34    continue
      do 32 m=1,3
      do 32 n=1,3
      a6(i,j,k,l,m,n)=0.0d0
      do 33 ip=1,3
      do 33 iq=1,3
      do 33 irr=1,3
      a6(i,j,k,l,m,n)=a6(i,j,k,l,m,n)+eps(ip,iq,irr)*
     1c(1,i,ip,j)*c(2,k,iq,l)*c(3,m,irr,n)
33    continue
32    continue
31    continue
c
      do 30 i=0,6
      do 30 j=0,6
      do 30 k=0,6
      b(i,j,k)=0.0d0
30    continue
c
      do 40 i=0,6
      do 40 j=0,6-i
      k=6-i-j
      do 41 il=1,3
      do 41 im=1,3
      do 41 ip=1,3
      do 41 iq=1,3
      do 41 ir=1,3
      do 41 is=1,3
      ii(1)=0
      ii(2)=0
      ii(3)=0
      ii(il)=ii(il)+1
      ii(im)=ii(im)+1
      ii(ip)=ii(ip)+1
      ii(iq)=ii(iq)+1
      ii(ir)=ii(ir)+1
      ii(is)=ii(is)+1
      if(i.eq.ii(1).and.j.eq.ii(2).and.k.eq.ii(3))
     1b(i,j,k)=b(i,j,k)+a6(il,im,ip,iq,ir,is)
41    continue
40    continue
      do 43 i=0,4
      do 43 j=0,4-i
      k=4-i-j
      do 44 il=1,3
      do 44 im=1,3
      do 44 ir=1,3
      do 44 is=1,3
      ii(1)=0
      ii(2)=0
      ii(3)=0
      ii(il)=ii(il)+1
      ii(im)=ii(im)+1
      ii(ir)=ii(ir)+1
      ii(is)=ii(is)+1
      if(i.eq.ii(1).and.j.eq.ii(2).and.k.eq.ii(3))
     1b(i,j,k)=b(i,j,k)-a4(il,im,ir,is)
44    continue
43    continue
      b(2,0,0)=c(1,1,1,1)+c(2,1,2,1)+c(3,1,3,1)
      b(0,2,0)=c(1,2,1,2)+c(2,2,2,2)+c(3,2,3,2)
      b(0,0,2)=c(1,3,1,3)+c(2,3,2,3)+c(3,3,3,3)
      b(1,1,0)=c(1,1,1,2)+c(2,1,2,2)+c(3,1,3,2)
     1        +c(1,2,1,1)+c(2,2,2,1)+c(3,2,3,1)
      b(1,0,1)=c(1,1,1,3)+c(2,1,2,3)+c(3,1,3,3)
     1        +c(1,3,1,1)+c(2,3,2,1)+c(3,3,3,1)
      b(0,1,1)=c(1,2,1,3)+c(2,2,2,3)+c(3,2,3,3)
     1        +c(1,3,1,2)+c(2,3,2,2)+c(3,3,3,2)
c
      b(0,0,0)=-1.0d0
c
c     =================================
      return
      end


      subroutine cubrt(a,cz)
c
c     roots of a cubic equation a(i)z**(i-1)
c     z**3+a2*z**2+a1*z+a0 as in abramowitz and stegun
c
      real*8 third,sn60,twn7
      complex*16 a(4),cz(3)
      complex*16 a0,a1,a2,g,h,a23,sud,q,omeg,omeg2,g2,h4,hq
      third=1.0d0/3.0d0
      twn7=1.0d0/27.0d0
      sn60=0.5d0*dsqrt(3.0d0)
c      write(6,*)third,twn7,sn60
      a2=a(3)/a(4)
      a1=a(2)/a(4)
      a0=a(1)/a(4)
      a23=a2*third
      h=(-a2*a2*third+a1)*third
      g=2.0d0*a2*a2*a2*twn7-a1*a2*third+a0
      g2=g*g
      h4=4.0d0*h*h*h
      if(cdabs(h4).gt.0.000000001d0*cdabs(g2))then
      sud=(-g+cdsqrt(g2+h4))/2.0d0
      else
      sud=h4/(4.0d0*g)
      end if
      q=sud**third
      omeg=dcmplx(-0.5d0,sn60)
      omeg2=dcmplx(-0.5d0,-sn60)
      hq=h/q
      cz(1)=q-hq-a23
      cz(2)=omeg*q-omeg2*hq-a23
      cz(3)=omeg2*q-omeg*hq-a23
c      write(6,*)omeg,omeg*omeg*omeg
c      write(6,*)omeg2,omeg2*omeg2*omeg2
      return
      end

      subroutine gaussj(a,n,np,b,m,mp)
      parameter (nmax=18)
      implicit real*8 (a-h,o-z)
      real*8 a(np,np),b(np,mp) 
      integer ipiv(nmax),indxr(nmax),indxc(nmax)
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.0d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (dabs(a(j,k)).ge.big)then
                  big=dabs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
c                pause 'singular matrix'
                write(*,'("singular matrix")') !bp
                read(*,*)
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.0d0) then
c           pause 'singular matrix.'
        write(*,'("singluar matrix.")') !bp
        read(*,*)
      endif
        pivinv=1.0d0/a(icol,icol)
        a(icol,icol)=1.0d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.0d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      end


      subroutine gaussj2(a,n,np,b,m,mp)
      parameter (nmax=9)
      implicit real*8 (a-h,o-z)
      complex*16 dum,pivinv
      complex*16 a(np,np),b(np,mp) 
      integer ipiv(nmax),indxr(nmax),indxc(nmax)
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.0d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (cdabs(a(j,k)).ge.big)then
                  big=cdabs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
c                pause 'singular matrix'
             write(*,'("singular matrix")') !bp
             read(*,*)
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.(0.0d0,0.0d0)) then
c       pause 'singular matrix.'
        write(*,'("singular matrix.")') !bp
        read(*,*)
      endif
        pivinv=1.0d0/a(icol,icol)
        a(icol,icol)=1.0d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.0d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      end


      subroutine zroots(a,m,roots,polish)
      parameter (eps=1.d-6,maxm=101)
      complex*16 a(*),roots(m),ad(maxm),x,b,c,xi
      logical polish
      xi=dcmplx(0.0d0,1.0d0)
      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      do 13 j=m,1,-1
        x=dcmplx(0.,0.)
        call laguer(ad,j,x,eps,.false.)
        if(dabs(dimag(x)).le.2.*eps**2*dabs(dimag(xi*x))) 
     1        x=dcmplx(dimag(xi*x),0.)
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
13    continue
      if (polish) then
        do 14 j=1,m
          call laguer(a,m,roots(j),eps,.true.)
14      continue
      endif
      do 16 j=2,m
        x=roots(j)
        do 15 i=j-1,1,-1
          if(dimag(xi*roots(i)).le.dimag(xi*x))go to 10
          roots(i+1)=roots(i)
15      continue
        i=0
10      roots(i+1)=x
16    continue
      return
      end

      subroutine laguer(a,m,x,eps,polish)
      complex*16 a(*),x,dx,x1,b,d,f,g,h,sq,gp,gm,g2,zero,xi
      logical polish
      parameter (zero=(0.0d0,0.0d0),epss=6.d-12,maxit=100)
      dxold=cdabs(x)
      do 12 iter=1,maxit
        b=a(m+1)
        err=cdabs(b)
        d=zero
        f=zero
        abx=cdabs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=cdabs(b)+abx*err
11      continue
        err=epss*err
        if(cdabs(b).le.err) then
          dx=zero
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.*f/b
          sq=cdsqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          if(cdabs(gp).lt.cdabs(gm)) gp=gm
          dx=m/gp
        endif
        x1=x-dx
        if(x.eq.x1)return
        x=x1 
        cdx=cdabs(dx)
        if(iter.gt.6.and.cdx.ge.dxold)return
        dxold=cdx
        if(.not.polish)then
          if(cdabs(dx).le.eps*cdabs(x))return
        endif
12    continue
c      pause 'too many iterations'
      write(*,'("too many iterations")') !bp
      read(*,*)
      return
      end

