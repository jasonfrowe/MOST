      program comparecurve
      implicit none
      integer nunit,npt,nmax,nunit2,npt2,pnpt,i,j,i2,j2,nfit,nfitmax
      parameter(nmax=120000,nfitmax=2)
      integer ia(nfitmax)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),time2(nmax),
     .   mag2(nmax),merr2(nmax),sky2(nmax),xc2(nmax),yc2(nmax),
     .   fx2(nmax),fy2(nmax),fxy2(nmax),tboard2(nmax),mfield2(nmax),
     .   ftotal2(nmax),aflux2(nmax),nap2(nmax),itime2(nmax),px(nmax),
     .   py(nmax),perr(nmax),tx(nmax),ty(nmax),x1,x2,y1,y2,
     .   covar(nfitmax,nfitmax),ans(nfitmax),chisq,rmavg,rmavg2,sig
      character*80 filename
     
      nunit=10
      write(6,*) "Enter filename1:"
      read(5,500) filename
      open(unit=nunit,file=filename,status='old',err=901)
 500  format(A80)
 
      nunit2=11
      write(6,*) "Enter filename2:"
      read(5,500) filename
      open(unit=nunit2,file=filename,status='old',err=901)
      
C     okay files are open.  Lets read in all that data

C     first file
      call readdata(nunit,npt,time,mag,merr,
     .     sky,xc,yc,fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime)
C     second file
      call readdata(nunit2,npt2,time2,mag2,merr2,sky2,xc2,yc2,fx2,fy2,
     .   fxy2,tboard2,mfield2,ftotal2,aflux2,nap2,itime2)

C     remove means
      call avgrm(npt,mag,rmavg)
      call avgrm(npt2,mag2,rmavg2)

C     sig clip both data sets
      sig=3.0
      call sigclip2(npt,time,mag,merr,sig,sky,xc,yc,fx,fy,
     .     fxy,tboard,mfield,ftotal,aflux,nap,itime,rmavg)
      call sigclip2(npt,time,mag,merr,sig,sky,xc,yc,fx,fy,
     .     fxy,tboard,mfield,ftotal,aflux,nap,itime,rmavg)
      call sigclip2(npt2,time2,mag2,merr2,sig,sky2,xc2,yc2,fx2,fy2,
     .     fxy2,tboard2,mfield2,ftotal2,aflux2,nap2,itime2,rmavg2)
      call sigclip2(npt2,time2,mag2,merr2,sig,sky2,xc2,yc2,fx2,fy2,
     .     fxy2,tboard2,mfield2,ftotal2,aflux2,nap2,itime2,rmavg2)

     
C     okay let plot some stuff

C     do not want to destroy time2, so we will save a copy of it.
      do 5 i=1,npt2
         tx(i)=time2(i)
 5    continue
      j2=npt2
      pnpt=0
      do 10 i=1,npt
         i2=0
         do 20 j=1,j2   
            if(time(i).eq.tx(j)) then
               pnpt=pnpt+1
               px(pnpt)=mag(i)
               py(pnpt)=mag2(j)
               perr(pnpt)=sqrt(merr(i)**2.0 + merr2(i)**2.0)
            else
               i2=i2+1
               ty(i2)=tx(j)
            endif
 20      continue
 
         do 25 j=1,i2
            tx(i2)=ty(i2)
 25      continue
         j2=i2
        
 10   continue      

      nfit=2
      do 30 i=1,nfit
         ia(i)=1
 30   continue

      call lfit(px,py,perr,pnpt,ans,ia,nfit,covar,nfitmax,chisq)
         
      call pgopen('?')
      x1=0.0
      x2=0.0
      y1=0.0
      y2=0.0
      write(6,*) "Number of points to plot:",pnpt
      call plotpoints(pnpt,px,py,-1,0,x1,x2,y1,y2)
      call plotfit(nfit,ans,x1,x2)
      call pgclos()
      
C     End of program.. clean up time
      close(10)
C     no errors so skip messages and end
      goto 999
C     cannot file the requestede file
 901  write(6,*) "File not found:",filename
C     error goto exit
      goto 999
C     EXIT!
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sigclip2(npt,time,mag,merr,sig,sky,xc,yc,fx,fy,
     .     fxy,tboard,peak,ftotal,aflux,nap,itime,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ntmp,i,nclip,niter,nitmax
      parameter (nmax=120000,nitmax=50)
      real time(npt),mag(npt),merr(npt),sig,std,stdev,tmp1(nmax),
     .     tmp2(nmax),tmp3(nmax),mean,dum,tmp4(nmax),tmp5(nmax),
     .     tmp6(nmax),tmp7(nmax),tmp8(nmax),tmp9(nmax),tmp10(nmax),
     .     tmp11(nmax),tmp12(nmax),sky(npt),xc(npt),omean,
     .     yc(npt),fx(npt),fy(npt),fxy(npt),tboard(npt),peak(npt),rmavg,
     .     tmp13(nmax),tmp14(nmax),tmp15(nmax),tmp16(nmax),ftotal(nmax),
     .     aflux(nmax),nap(nmax),itime(nmax)


C     watch out for infinite loops
      niter=0

C     find mean of data set
      ntmp=0
      mean=0.
      do 4 i=1,npt
         if(mag(i).lt.90.0) then
            ntmp=ntmp+1
            mean=mean+mag(i)
         endif
 4    continue
      mean=mean/real(ntmp)
C     find standard dev. of data set.
      std=stdev(npt,mag,mean)

      write(6,*) "sigclip:",npt,mean,std

C     count number of clipped points
 6    nclip=0
C     count number of new points
      ntmp=0
      do 10 i=1,npt
         dum=abs(mag(i)-mean)
         if(dum.lt.sig*std) then
            ntmp=ntmp+1
            tmp1(ntmp)=time(i)
            tmp2(ntmp)=mag(i)
            tmp3(ntmp)=merr(i)
            tmp5(ntmp)=sky(i)
            tmp6(ntmp)=xc(i)
            tmp7(ntmp)=yc(i)
            tmp8(ntmp)=fx(i)
            tmp9(ntmp)=fy(i)
            tmp10(ntmp)=fxy(i)
            tmp11(ntmp)=tboard(i)
            tmp12(ntmp)=peak(i)
            tmp13(ntmp)=ftotal(i)
            tmp14(ntmp)=aflux(i)
            tmp15(ntmp)=nap(i)
            tmp16(ntmp)=itime(i)
         else
            nclip=nclip+1
         endif
 10   continue
C     if nothing is clipped, no point in continuing
      if(nclip.eq.0) goto 15
      
C     save copy of old mean
      omean=mean

C     find mean of new data set
      mean=0.
      do 5 i=1,ntmp
         mean=mean+tmp2(i)
 5    continue
      mean=mean/real(ntmp)
C     if mean doesn't change, we're done.
      if(mean.eq.omean) goto 15

C     if mean doesn't change, we're done
      if(nclip.gt.0) goto 15 
C     find st. dev. of new data set
      std=stdev(ntmp,tmp2,mean)
C     now restart clipping
      niter=niter+1
C     if we loop too much.. get out. (probably an infinite loop)
      if(niter.gt.nitmax) goto 15
      goto 6

 15   do 20 i=1,ntmp
         time(i)=tmp1(i)
         mag(i)=tmp2(i)-mean
         merr(i)=tmp3(i)
         sky(i)=tmp5(i)
         xc(i)=tmp6(i)
         yc(i)=tmp7(i)
         fx(i)=tmp8(i)
         fy(i)=tmp9(i)
         fxy(i)=tmp10(i)
         tboard(i)=tmp11(i)
         peak(i)=tmp12(i)
         ftotal(i)=tmp13(i)
         aflux(i)=tmp14(i)
         nap(i)=tmp15(i)
         itime(i)=tmp16(i)
 20   continue
      npt=ntmp
C     update zero point
      rmavg=rmavg+mean

      return
      end
     
     
     