c23456789012345678901234567890123456789012345678901234567890123456789012
C        1         2         3         4         5         6         7
      program mostfit3
      implicit none
      integer nmax,npt,nphas,nunit,i,nph,j,nfitmax,nfit,nfitc,dumi,
     .     ntest,bins,xcoo(2),ycoo(2),ntime,k,nzpts,nt
      parameter(nmax=600000,nfitmax=20)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .     fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .     ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),orbfreq,per,
     .     pht(nmax),xt(nmax),yt(nmax),gain,ans(nfitmax),minph(nmax),
     .     minphr(nmax),modecal,phoffset,aa(nmax,nfitmax),yerr(nmax),
     .     atime(nmax),at(nmax),rmavg,bx(nmax),by(nmax),skystd(nmax),
     .     x1,x2,y1,y2,px(nmax),py(nmax),etimes(nmax),etimeW,zpts,
     .     xcoom,ycoom,pmin,pmax
      double precision dtime(nmax)
      character*80 filename,cline
      common /funct1/ npt,etimeW,itime,gain,time,mag,merr,sky

      pmin=-1.0
      pmax=-1.0
      if(iargc().lt.1) goto 903 !need at least one input arguement
      call getarg(1,filename)
      if(iargc().ge.3) then
      	 call getarg(2,cline)
      	 read(cline,*) pmin
         call getarg(3,cline)
         read(cline,*) pmax
      endif

      rmavg=0.0

C     setting initial zeros
      do 5 i=1,nmax
         minph(i)= 99.9e30
         minphr(i)= 99.9
 5    continue

C     The gain of the detector
      gain=6.1

c      write(6,*) "Enter Photometry filename"
c      read(5,500) filename
 500  format(A80)
      
C     set the unit number to 10 and open it.
      nunit=10
      open(unit=nunit,file=filename,status='old',err=901)
      
C     the orbital frequency
c      orbfreq=14.199363
      orbfreq=14.18139

C     first lets read in the raw photometry database
      write(6,*) "Reading in data set"
      call readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,fy,fxy,
     .   tboard,mfield,ftotal,aflux,nap,itime,skystd)

C     bias correction stuff
c      do 4 i=1,npt
c         sky(i)=sky(i)+tboard(i)
c 4    continue

C     get all the exposure times
      call getetimes(ntime,etimes,npt,itime)

C     change the time into a phase
      
C     the period of orbit (in days)
      per=1.0/orbfreq
      
C     save the phase info in pht
      write(6,*) "Phasing data"
      call phase(npt,time,pht,per)

C     Now we can find the sky/flux depedence for each occurance of the
C     straylight.  So first lets see how many phases we have
      nphas=0
      do 10 i=1,npt
         nphas=max(int(pht(i)),nphas)
 10   continue
      write(6,*) "Number of stray light events:",nphas

      nph=0
      xcoom=0.
      ycoom=0.
C     Now we calculate the average minimum of the phase bins.
      do 40 i=1,npt
         xcoom=xcoom+xc(i)
         ycoom=ycoom+yc(i)
         if(int(pht(i)).gt.0) then
            if(minph(int(pht(i))).gt.mag(i)) then
               minph(int(pht(i)))=mag(i)
               minphr(int(pht(i)))=pht(i)-int(pht(i))
            endif
         endif
 40   continue
      xcoo(1)=xcoom/real(npt)-10
      ycoo(1)=ycoom/real(npt)-10
      xcoo(2)=xcoom/real(npt)+10
      ycoo(2)=xcoom/real(npt)+10

C     Pixel Detrending!!!
      if(pmin.ne.pmax) then !see if user wants pixel detrending
      	call pixeldetrend2(npt,mag,xc,yc,pmin,pmax)
      	goto 900
      endif

C	  I use the mode of all the calculated phases to find lowest point
C	  in the phase.
      phoffset=modecal(nphas,minphr,0.0,1.0)
c      write(6,*) "phoffset:",phoffset

C     now we shift it 0.5 in phase and check for any phases greater than
C     one
      nph=0
      do 50 i=1,npt
         pht(i)=pht(i)+1.5-phoffset
         if(int(pht(i)).eq.1) then
            nph=nph+1
            xt(nph)=pht(i)
            yt(nph)=mag(i)
         endif
 50   continue

      call pgopen('?')
      call pgpage()

C     we will calculate the average zero point
      rmavg=0.
C     for counting number of valid zero points
      nzpts=0
C     loop over the results for each exposure time.
      write(6,*) "ntime:",ntime
      do 19 k=1,ntime
         etimeW=etimes(k)
         write(6,*) "exposure time:",etimeW
         nfitc=0
C        why was this 0->nphas?
         do 20 i=1,nphas
C     nph counts the number of points in each phase
            nph=0
C     atime is used to calculate the average time in each phase.
            atime(i)=0.0
            do 25 j=1,npt
               if((int(pht(j)).eq.i).and.(itime(j).eq.etimeW)) then
C     inserted for HD 209458
                  if((sky(j).gt.1.0).and.(sky(j).lt.1500.0))then
                     nph=nph+1
                     atime(i)=atime(i)+time(j)
                     xt(nph)=sky(j)
                     yt(nph)=10**((mag(j)-25.0)/(-2.5))/gain
                  endif
               endif
 25         continue
 
C        if you get division by zero here, try increasing the 
C        minimum value of nph for the fit to procede.
            if(nph.gt.20) then
c            write(6,*) "i",i,nph
               atime(i)=atime(i)/real(nph)
               nfit=2
               call fitskydep(nph,xt,yt,nfit,ans)
               dumi=0
               do 21 j=1,nfit               
                  dumi=dumi+ntest(ans(j))
 21            continue
               if(dumi.eq.0) then
                  nfitc=nfitc+1
                  do 22 j=1,nfit
                     aa(nfitc,j)=ans(j)
 22               continue
                  at(nfitc)=atime(i)
c                  px(nfitc)=atime(i)
c                  py(nfitc)=ans(2)
               endif
            endif
 20      continue
 
         x1=0.0
         x2=0.0
         y1=0.2
         y2=0.2
c      call pgpage()
c      call plotpoints(nfitc,px,py,17,0,x1,x2,y1,y2)

 
C     now we remove the overall gradient.
c         write(6,*) "Remove evolving linear fit"
c         if(nfitc.ge.3) then
c            call removesky(npt,time,mag,sky,itime,nfitc,nfit,at,
c     .         aa,gain,etimeW)
c         endif

C     now we remove a spline of the relationship.
         write(6,*) "Remove stationary spline"
         call removespline(npt,time,mag,sky,itime,mfield,gain,0,etimeW)
c         rmavg=0.
c         goto 900
C     Now we fit the relationship between CCD position and magnitude
c      call posfit(npt,time,mag,merr,xc,yc)

         nph=0
         do 26 i=1,npt
            nt=0
            nt=nt+ntest(time(i))
            nt=nt+ntest(mag(i))
            nt=nt+ntest(merr(i))
            if((itime(i).eq.etimeW).and.(mag(i).lt.90.0).and.(nt.eq.0)
     .            )then
               nph=nph+1
               xt(nph)=time(i)
               yt(nph)=mag(i)
               yerr(nph)=merr(i)
            endif
 26      continue
         if(nph.gt.20) then
            nzpts=nzpts+1
            call rejhilow(nph,xt,yt,yerr,2,2)
            call sigclip(nph,xt,yt,yerr,2.5)
            call avgrm(nph,yt,zpts)
            write(6,*) "zpts",zpts
            rmavg=rmavg+zpts
         else
C        mark errors
            zpts=-99.9
         endif
         do 27 i=1,npt
           if(itime(i).eq.etimeW) mag(i)=mag(i)-zpts
 27      continue 
 19   continue

c      call removespline(npt,time,mag,sky,itime,mfield,gain,0,0.)

      rmavg=rmavg/real(nzpts)
      write(6,*) "rmavg:",rmavg
 900  filename="test.dat"
c      do 3 i=1,npt
c         sky(i)=sky(i)-tboard(i)
c 3    continue
      call exportdata(npt,dtime,mag,merr,filename,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,rmavg)

      call pgclos()
      
      close(10)
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 903  write(6,*) "Usage: mostfit3 <filename> [mmin] [mmax]"
      write(6,*) "<filename> : MOST Photometry file"
      write(6,*) "[mmin] [mmax] : (optional) if defined pixel detrending
     . is performed instead"
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getetimes(ntime,etimes,npt,itime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer ntime,npt,i,j,mm
      real etimes(npt),itime(npt)
      
      ntime=1
      etimes(ntime)=itime(1)
      
      do 10 i=2,npt
         mm=0
         do 20 j=1,ntime
            if(etimes(j).eq.itime(i))then
               mm=1
            endif
 20      continue
         if(mm.eq.0) then
            ntime=ntime+1
            etimes(ntime)=itime(i)
            write(6,*) ntime,etimes(ntime)
         endif 
 10   continue
  
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pixeldetrend2(npt,mag,xc,yc,maglo,maghi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j
      real mag(npt),xc(npt),yc(npt),x1,x2,ymag,dy,nm,temp,maglo,maghi
      
      
c      maglo=-0.005
c      maghi= 0.005
      write(6,*) "maglo,maghi",maglo,maghi
      nm=4.0
      do 10 i=1,npt
         x1=xc(i)
         x2=yc(i)
         ymag=0.
         dy=0.
         do 48 j=1,npt
            if((mag(j).gt.maglo).and.(mag(j).lt.maghi))then
               temp=(x1-xc(j))**nm + (x2-yc(j))**nm
               if(temp.gt.0)then
                  ymag=ymag+mag(j)/temp
                  dy=dy+1.0/temp
               endif
c               write(6,*) ymag,dy
            endif
 48      continue
         ymag=ymag/dy
         mag(i)=mag(i)-ymag
c         write(6,*) x1,x2,ymag,mag(i)
c         read(5,*)
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pixeldetrend(npt,mag,xc,yc,gain,xcoo,ycoo)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,fx,fy,i,j,yn,k,mm,nn,xcoo(2),ycoo(2),nct,tgrdn
      parameter(fx=60,fy=60)
      integer fgridn(fx,fy),tgrdx(fx*fy),tgrdy(fx*fy)
      real mag(npt),xc(npt),yc(npt),fgrid(fx,fy),ymean,gain,
     .     tgrdv(fx*fy),x1,x2,dy,nm,temp,x1a(fx),x2a(fy),ymag

      do 41 i=1,fx
         do 42 j=1,fy
            fgrid(i,j)=0.0
            fgridn(i,j)=0
 42      continue
 41   continue

      ymean=0.0
      yn=0
      do 40 i=1,npt
         j=int(xc(i)+0.5)-xcoo(1)
         k=int(yc(i)+0.5)-ycoo(1)
c         write(6,*) j,k,xc(i),yc(i),xcoo(1),ycoo(1)
         if((j.gt.0).and.(j.lt.fx).and.(k.gt.0).and.(k.lt.fy).and.
     .        (mag(i).lt.14.95).and.(mag(i).gt.14.89))then
            fgrid(j,k)=fgrid(j,k)+mag(i)
            fgridn(j,k)=fgridn(j,k)+1
            ymean=ymean+mag(i)
            yn=yn+1
         endif
c         read(5,*)
 40   continue
      ymean=ymean/real(yn)
c      write(6,*) "ymean",ymean
      mm=xcoo(2)-xcoo(1)+1
      nn=ycoo(2)-ycoo(1)+1

      do 47 i=1,mm
         x1a(i)=real(i+xcoo(1)-1)
 47   continue

      do 44 i=1,nn
         x2a(i)=real(i+ycoo(1)-1)
 44   continue

      nct=20
      tgrdn=0
      do 43 i=1,mm
         do 46 j=1,nn
            write(6,*) fgridn(i,j)
            if(fgridn(i,j).gt.nct) then
               fgrid(i,j)=fgrid(i,j)/real(fgridn(i,j))
               tgrdn=tgrdn+1
               tgrdx(tgrdn)=x1a(i)
               tgrdy(tgrdn)=x2a(j)
               tgrdv(tgrdn)=fgrid(i,j)
               write(6,*) tgrdn,tgrdx(tgrdn),tgrdy(tgrdn),tgrdv(tgrdn)
            else
               fgrid(i,j)=ymean
            endif
c           if(fgridn(i,j).gt.nct) write(6,*) i,j,fgrid(i,j),fgridn(i,j)
 46      continue
 43   continue

      do 10 i=1,npt

         x1=xc(i)
         x2=yc(i)
c     call polin2(x1a,x2a,fgrid,mm,nn,x1,x2,ymag,dy)
         ymag=0.
         dy=0.
         nm=1.0
         do 48 j=1,tgrdn
            temp=(x1-real(tgrdx(j)))**nm + (x2-real(tgrdy(j)))**nm
            ymag=ymag+tgrdv(j)/temp
            dy=dy+1.0/temp
 48      continue
         ymag=ymag/dy
         ymag=ymag-ymean
c         write(6,*) ymag,dy,mag(i)
         mag(i)=mag(i)-ymag
         if((mag(i).gt.0.0).and.(mag(i).lt.90.0)) then
            mag(i)=mag(i)
         else
            mag(i)=99.9
         endif
 10   continue

      return
      end
