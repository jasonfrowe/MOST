C Contents
C
C spdetrend
C phasept
C jmfour
C plotphbin
C phase
C nyquest

c23456789012345678901234567890123456789012345678901234567890123456789012
C        1         2         3         4         5         6         7

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function nyquest(npt,time)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ndt,i
      parameter(nmax=600000)
      real time(npt),dt(nmax),mindt,maxdt,mode,modecal,mean,std,stdev,
     .     sigcut

      sigcut=3.0

      ndt=0
      mindt= 99.9e30
      maxdt=-99.9e30
      mean=0.
      do 10 i=2,npt
         ndt=ndt+1
         dt(ndt)=time(i)-time(i-1)
         mindt=min(dt(ndt),mindt)
         maxdt=max(dt(ndt),maxdt)
         mean=mean+dt(ndt)
 10   continue
      mean=mean/real(ndt)
      std=stdev(ndt,dt,mean)
      mindt=max(mindt,mean-sigcut*std)
      maxdt=min(maxdt,mean+sigcut*std)

      mode=modecal(ndt,dt,mindt,maxdt)
      
c      write(6,*) "mean,mode,mindt,maxdt:",mean,mode,mindt,maxdt

      nyquest=1.0/(2.0*mode)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine spdetrend(npt1,time1,mag1,merr1,npt2,time2,mag2,merr2,
     .     tbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt1,npt2,nmax,i,bnpt
      parameter(nmax=600000)
      real time1(npt1),mag1(npt1),merr1(npt1),bdata(nmax),bderr(nmax),
     .     yp1,yp2,y2(nmax),btime(nmax),x,y,tbin,time2(npt2),mag2(npt2),
     .     merr2(npt2)

      bnpt=npt2
      do 10 i=1,npt2
         btime(i)=time2(i)
         bdata(i)=mag2(i)
         bderr(i)=merr2(i)
 10   continue

      call bindt(bnpt,btime,bdata,bderr,tbin,2,2)

      yp1=1.0e30
      yp2=1.0e30
      call spline(btime,bdata,bnpt,yp1,yp2,y2)
      
      do 20 i=1,npt1
         x=time1(i)
         call splint(btime,bdata,y2,bnpt,x,y)
         mag1(i)=mag1(i)-y
 20   continue

      do 21 i=1,npt2
         x=time2(i)
         call splint(btime,bdata,y2,bnpt,x,y)
         mag2(i)=mag2(i)-y
 21   continue

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phasept(npt,time,mag,phase,period)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     phase - returned phase of observation
C     period - fixed period for data
      implicit none
      integer npt
      real time(npt),mag(npt),phase(npt),period

      integer i
      real temp

      do 10 i=1,npt
c         temp=time(i)/period
c         phase(i)=temp-int(temp)
         temp=time(i)!-time(1)
         phase(i)=temp/period-int(temp/period)
c         write(6,*) time(i),phase(i)

c         phase(i)=phase(i)+0.5
c         if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
 10   continue
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jmfourw(npt,time,mag,merr,per1,per2,steps,bper,btheta,
     .     panx,pany,plot,avgamp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,steps,panx,pany,nmax,stepmax,k,j,plot
      parameter(stepmax=2000)
      real time(npt),mag(npt),per1,per2,bper,btheta,percold,perc,by(4),
     .     nx(stepmax),ny1(stepmax),ny2(stepmax),px(3),py(3),bx(4),tbin,
     .     orbper,dum,weight,merr(npt),avgamp
      parameter(nmax=600000)
      character*80 tmpc

      INTEGER I,NUMPTS,NMDENU
      DOUBLE PRECISION JD (nmax),B(nmax),LONU,DELTNU,TWOPI
      DOUBLE PRECISION SIDEL(nmax),CODEL(nmax),SI(nmax),CO(nmax)
      DOUBLE PRECISION SISUM,COSUM,SN,NU,AMP,PHI
      CHARACTER *20 INFILE,SPFILE
      DOUBLE PRECISION JDTMP,LOJD,HIJD,BTMP,JDZRO,Davgamp

      if((steps.lt.stepmax).and.(plot.eq.1)) then
         write(6,*) "Increasing number of steps to 2000"
         steps=stepmax
      endif


C     for calculating the average amplitude
      Davgamp=0.

      tbin=0.0

      btheta=-99.9e10
      do 5 i=1,stepmax
         ny1(i)=-99.9e10
         ny2(i)= 99.9e10
 5    continue

c      open(unit=21,file="junk.dat")

      LOJD=0.
      HIJD=100.
      JDZRO=0.

      LONU=dble(per1)
      DELTNU=dble(per2-per1)/dble(steps)
      NMDENU=dble(steps)

      do 3 i=1,npt
         jd(i)=dble(time(i))
         b(i)=dble(mag(i))
 3    continue
      NUMPTS = npt

      TWOPI = 6.2831853
      NU = LONU
      DO 61 I=1,NUMPTS
        SI(I) = SIN(TWOPI*LONU*JD(I))
        CO(I) = COS(TWOPI*LONU*JD(I))
        SIDEL(I) = SIN(TWOPI*DELTNU*JD(I))
        CODEL(I) = COS(TWOPI*DELTNU*JD(I))
   61 CONTINUE

      percold=0.0      
      DO 51 J=1,NMDENU
         perc=100.0*real(j)/real(steps)
         if(perc-percold.gt.0.1) then
c            write(tmpc,500) "Percent done ", perc
            percold=perc
c            call ovrwrt(tmpc,2)
         endif
 500     format(A13,F6.1)
         SISUM = 0.
         COSUM = 0.
         weight=0.0
         DO 41 I=1,NUMPTS
            SISUM = SISUM + B(I)*SI(I)/merr(i)
            COSUM = COSUM + B(I)*CO(I)/merr(i)
            weight=weight+1.0/merr(i)
 41      CONTINUE
C     WRITE (50,*) NU,SISUM,COSUM
         AMP = (2./weight)*SQRT(SISUM**2. + COSUM**2.)
         PHI = ATAN(-SISUM/COSUM)
C     calculate average of amplitudes
         Davgamp=Davgamp+AMP
         IF (COSUM.LT.0) PHI = PHI + TWOPI/2.
         IF ((SISUM.GT.0).AND.(COSUM.GE.0)) PHI = PHI + TWOPI
c     
         if(plot.eq.1) then
            k=1+real(stepmax)/real(steps)*(j-1)
            nx(k)=real(nu)
            if(real(amp).gt.ny1(k)) ny1(k)=real(amp)
            if(real(amp).lt.ny2(k)) ny2(k)=real(amp)
         endif
         if(tbin.gt.0.) then
         if((real(amp).gt.btheta).and.(real(nu).gt.1440.0/tbin)) then
            btheta=real(amp)
            bper=real(nu)
         endif
         else
         if(real(amp).gt.btheta) then
            btheta=real(amp)
            bper=real(nu)
         endif
         endif 
c     call pgpt1(NU,AMP,1)
c     WRITE (21,*) NU,AMP,PHI
c     
         DO 31 I=1,NUMPTS
            SN = SI(I)
            SI(I) = SN*CODEL(I) + CO(I)*SIDEL(I)
            CO(I) = CO(I)*CODEL(I) - SN*SIDEL(I)
 31      CONTINUE
         NU = NU + DELTNU
 51   CONTINUE
      
      avgamp=Davgamp/real(steps)
      
      if(plot.eq.1) then
c         call pgpage()
         call pgpanl(panx,pany)
         call pgwindow(per1,per2,0.0,btheta)
         call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
c         call pglabel("Frequency","AMP"," ")
         call pgbbuf()
      
c         close(21)
         do 30 i=2,stepmax
            px(1)=nx(i-1)
            py(1)=ny2(i-1)
            px(2)=nx(i)
            px(3)=nx(i)
            py(2)=ny1(i)
            py(3)=ny2(i)
            call pgline(3,px,py)
 30      continue
         call pgebuf()


         call pgsls(2)
         call pgsci(2)
         orbper=14.199363
         i=per1/orbper+1
         j=per2/orbper
      
         do 100 k=i,j
            px(1)=k*orbper
            px(2)=px(1)
            py(1)=0.0
            py(2)=btheta
            call pgline(2,px,py)
 100     continue

         call pgsls(1)
         call pgsci(1)

         dum=1000.0/86400.0
         call pgwindow(dum*per1,dum*per2,0.0,btheta)
         call pgbox('CMTS1',0.0,0,'',0.0,0)

c      px(1)=360.0/0.997268
c      px(2)=px(1)
c      py(1)=0.0
c      py(2)=btheta
c      call pgline(2,px,py)
c      call pgsls(1)
c      if(tbin.gt.0.0) then
c         bx(1)=per1
c         by(1)=0.0
c         bx(2)=(60.0*24.0)/tbin
c         by(2)=0.0
c         bx(3)=bx(2)
c         by(3)=btheta
c         bx(4)=per1
c         by(4)=btheta
c         call pgsfs(3)
c         call pgpoly(4,bx,by)
c         call pgsfs(1)
c      endif
c      call pgsci(1)

      endif
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jmfour(npt,time,mag,per1,per2,steps,bper,btheta,
     .     panx,pany,plot)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,steps,panx,pany,nmax,stepmax,k,j,plot
      parameter(stepmax=2000)
      real time(npt),mag(npt),per1,per2,bper,btheta,percold,perc,by(4),
     .     nx(stepmax),ny1(stepmax),ny2(stepmax),px(2),py(2),bx(4),tbin,
     .     orbper,dum
      parameter(nmax=600000)
      character*80 tmpc

      INTEGER I,NUMPTS,NMDENU
      DOUBLE PRECISION JD (nmax),B(nmax),LONU,DELTNU,TWOPI
      DOUBLE PRECISION SIDEL(nmax),CODEL(nmax),SI(nmax),CO(nmax)
      DOUBLE PRECISION SISUM,COSUM,SN,NU,AMP,PHI
      CHARACTER *20 INFILE,SPFILE
      DOUBLE PRECISION JDTMP,LOJD,HIJD,BTMP,JDZRO

      if(steps.lt.stepmax) then
         write(6,*) "Increasing number of steps to 2000"
         steps=stepmax
      endif


      tbin=0.0

      btheta=-99.9e10
      do 5 i=1,stepmax
         ny1(i)=-99.9e10
         ny2(i)= 99.9e10
 5    continue

c      open(unit=21,file="ft.dat")

      LOJD=0.
      HIJD=100.
      JDZRO=0.

      LONU=dble(per1)
      DELTNU=dble(per2-per1)/dble(steps)
      NMDENU=dble(steps)

      do 3 i=1,npt
         jd(i)=dble(time(i))
         b(i)=dble(mag(i))
 3    continue
      NUMPTS = npt

      TWOPI = 6.2831853
      NU = LONU
      DO 61 I=1,NUMPTS
        SI(I) = SIN(TWOPI*LONU*JD(I))
        CO(I) = COS(TWOPI*LONU*JD(I))
        SIDEL(I) = SIN(TWOPI*DELTNU*JD(I))
        CODEL(I) = COS(TWOPI*DELTNU*JD(I))
   61 CONTINUE

      percold=0.0      
      DO 51 J=1,NMDENU
         perc=100.0*real(j)/real(steps)
         if(perc-percold.gt.0.1) then
            write(tmpc,500) "Percent done ", perc
            percold=perc
            call ovrwrt(tmpc,2)
         endif
 500     format(A13,F6.1)
         SISUM = 0.
         COSUM = 0.
         DO 41 I=1,NUMPTS
            SISUM = SISUM + B(I)*SI(I)
            COSUM = COSUM + B(I)*CO(I)
 41      CONTINUE
C     WRITE (50,*) NU,SISUM,COSUM
         AMP = (2./NUMPTS)*SQRT(SISUM**2. + COSUM**2.)
         PHI = ATAN(-SISUM/COSUM)
         IF (COSUM.LT.0) PHI = PHI + TWOPI/2.
         IF ((SISUM.GT.0).AND.(COSUM.GE.0)) PHI = PHI + TWOPI
c     
         k=1+real(stepmax)/real(steps)*(j-1)
         nx(k)=real(nu)
         if(real(amp).gt.ny1(k)) ny1(k)=real(amp)
         if(real(amp).lt.ny2(k)) ny2(k)=real(amp)
         if(tbin.gt.0.) then
         if((real(amp).gt.btheta).and.(real(nu).gt.1440.0/tbin)) then
            btheta=real(amp)
            bper=real(nu)
         endif
         else
         if(real(amp).gt.btheta) then
            btheta=real(amp)
            bper=real(nu)
         endif
         endif 
c     call pgpt1(NU,AMP,1)
c     WRITE (21,*) NU,AMP,PHI
c     
         DO 31 I=1,NUMPTS
            SN = SI(I)
            SI(I) = SN*CODEL(I) + CO(I)*SIDEL(I)
            CO(I) = CO(I)*CODEL(I) - SN*SIDEL(I)
 31      CONTINUE
         NU = NU + DELTNU
 51   CONTINUE

      if(plot.eq.0) goto 901
      
c      call pgpage()
      call pgpanl(panx,pany)
      call pgwindow(per1,per2,0.0,btheta)
      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Frequency","AMP"," ")
      call pgbbuf()
      
c     close(21)
      do 30 i=1,stepmax
         px(1)=nx(i)
         px(2)=nx(i)
         py(1)=ny1(i)
         py(2)=ny2(i)
         call pgline(2,px,py)
 30   continue
      call pgebuf()


      call pgsls(2)
      call pgsci(2)
      orbper=14.199363
      i=per1/orbper+1
      j=per2/orbper
      
      do 100 k=i,j
         px(1)=k*orbper
         px(2)=px(1)
         py(1)=0.0
         py(2)=btheta
         call pgline(2,px,py)
 100  continue

      call pgsls(1)
      call pgsci(1)

      dum=1000.0/86400.0
      call pgwindow(dum*per1,dum*per2,0.0,btheta)
      call pgbox('CMTS1',0.0,0,'',0.0,0)

c      px(1)=360.0/0.997268
c      px(2)=px(1)
c      py(1)=0.0
c      py(2)=btheta
c      call pgline(2,px,py)
c      call pgsls(1)
c      if(tbin.gt.0.0) then
c         bx(1)=per1
c         by(1)=0.0
c         bx(2)=(60.0*24.0)/tbin
c         by(2)=0.0
c         bx(3)=bx(2)
c         by(3)=btheta
c         bx(4)=per1
c         by(4)=btheta
c         call pgsfs(3)
c         call pgpoly(4,bx,by)
c         call pgsfs(1)
c      endif
c      call pgsci(1)
      
 901  continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phase(npt,time,pht,per)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This routine will change the a time into a phase using the 
C     supplied period.  I do not remove the integer part, use the
C     routine phasept to do that!
C
C     npt - number of points
C     time - the date of the observation
C     pht - the phase of the data (returned)
C     per - the period for phasing the data

      integer npt,i,nmax
      parameter(nmax=600000)
      real time(npt),pht(nmax),per,mint

C     find the start time of the data set.
C     setting some really low value
      mint= 99.9e30
      do 10 i=1,npt
         mint=min(time(i),mint)
 10   continue
      
      do 20 i=1,npt
         pht(i)=(time(i)-mint)/per
 20   continue
     

      return
      end
