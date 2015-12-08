      program autocut
      implicit none
      integer iargc,nunit,nmax,npt,nbin,nbinmax,i,ncut,j,niter,i1,i2
      parameter(nmax=600000,nbinmax=500)
      integer nd(nmax),flag(nmax)
      real rp(nmax),bdatax(nbinmax),bdatay(nbinmax),errs(6),ave,std
      double precision work(nmax),sigcut
      double precision time(nmax),mag(nmax),merr(nmax),sky(nmax),
     .  xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),
     .  mfield(nmax),ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),
     .  skystd(nmax)
      character*80 filename,title,cline,outname

      if(iargc().lt.2) goto 901 !get command line arguments. 
      call getarg(1,filename)
      call getarg(2,outname)

      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      call readdata(nunit,npt,time,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
      close(nunit)

      do 10 i=1,npt
        flag(i)=0
 10   continue

      nbin=30
  
      sigcut=3.0
      if(iargc().ge.3)then
        call getarg(3,cline)
        read(cline,*) sigcut
      endif
      if(sigcut.le.1.5)then
        i1=1
        i2=2
      elseif(sigcut.le.2.5)then
        i1=3
        i2=4
      else
        i1=5
        i2=6
      endif

      niter=5
      if(iargc().ge.4)then
        call getarg(4,cline)
        read(cline,*) niter
        if(niter.lt.1) goto 901
        if(niter.gt.5) niter=5
      endif

      call pgopen('?') !open PGPlot device
      call pgask(.true.) !don't ask for new page.. just do it.
      call PGPAP ( 7.0 ,1.0) !paper size
      call pgsubp(5,5)  !break up plot into grid

C     Eliminate passages through the SAA
      do 16 i=1,npt
        if(mfield(i).lt.20000) flag(i)=1
 16   continue

      do 20 j=1,niter
      write(0,*) "Iter:",j

      title="Xcoo"
      call histogram(npt,rp,xc,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,flag)
      write(6,500) ave,std,(errs(i),i=1,6)
 500  format(9(F11.6,1X))
      if((errs(i1).ne.0).and.(errs(i2).ne.0))then
        do 11 i=1,npt
            if(xc(i).gt.ave+errs(i1))then
                flag(i)=1
            elseif(xc(i).lt.ave+errs(i2))then
                flag(i)=1
c        else
c            flag(i)=0
            endif
 11     continue
      endif
 
      title="Ycoo"
      call histogram(npt,rp,yc,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,flag)
      write(6,500) ave,std,(errs(i),i=1,6)
      ncut=0
      if((errs(i1).ne.0).and.(errs(i2).ne.0))then
        do 12 i=1,npt
            if(yc(i).gt.ave+errs(i1))then
                flag(i)=1
            elseif(yc(i).lt.ave+errs(i2))then
                flag(i)=1
c        else
c            flag(i)=0
            endif
            ncut=ncut+flag(i)
 12     continue
      endif
c      write(6,*) "ncut:",ncut
 
      title="sig_x"
      call histogram(npt,rp,fx,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,flag)
      write(6,500) ave,std,(errs(i),i=1,6)
      if((errs(i1).ne.0).and.(errs(i2).ne.0))then
        do 13 i=1,npt
            if(fx(i).gt.ave+errs(i1))then
                flag(i)=1
            elseif(fx(i).lt.ave+errs(i2))then
                flag(i)=1
c        else
c            flag(i)=0
            endif
 13     continue
      endif

      title="sig_y"
      call histogram(npt,rp,fy,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,flag)
      write(6,500) ave,std,(errs(i),i=1,6)
      if((errs(i1).ne.0).and.(errs(i2).ne.0))then
        do 14 i=1,npt
            if(fy(i).gt.ave+errs(i1))then
                flag(i)=1
            elseif(fy(i).lt.ave+errs(i2))then
                flag(i)=1
c        else
c            flag(i)=0
            endif
 14     continue
      endif
      
      title="sig_xy"
      call histogram(npt,rp,fxy,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,flag)
      write(6,500) ave,std,(errs(i),i=1,6)
      if((errs(i1).ne.0).and.(errs(i2).ne.0))then
        do 15 i=1,npt
            if(fxy(i).gt.ave+errs(i1))then
                flag(i)=1
            elseif(fxy(i).lt.ave+errs(i2))then
                flag(i)=1
c        else
c            flag(i)=0
            endif
 15     continue
      endif
      
 20   continue
      
      open(unit=nunit,file=outname) !export detrended data
c      write(0,*) "hello1"
      call exportdata(nunit,npt,time,mag,merr,sky,xc,yc,fx,fy,fxy,
     .  tboard,mfield,ftotal,aflux,nap,itime,skystd,flag)
      close(nunit)

      call pgclos()

      goto 999
 901  write(0,*) "Usage: autocut photometry.dat output.dat [sigcut] [nit
     .er]"
      write(0,*) "  photometry.dat - MOST photometry"
      write(0,*) "  output.dat     - output filename"
      write(0,*) "  sigcut         - sigma level for cutting (1,2 or 3)"
      write(0,*) "  niter          - number of iterations (1-5)"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportdata(nunit,npt,time,mag,merr,sky,xc,yc,fx,fy,fxy,
     .  tboard,mfield,ftotal,aflux,nap,itime,skystd,flag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,npt,i,flag(npt)
      double precision time(npt),mag(npt),merr(npt),sky(npt),xc(npt),
     .  yc(npt),fx(npt),fy(npt),fxy(npt),tboard(npt),mfield(npt),
     .  ftotal(npt),aflux(npt),nap(npt),itime(npt),skystd(npt)
     
c      write(0,*) "hello2",npt
      do 10 i=1,npt
        if(flag(i).eq.0)write(nunit,510) time(i),mag(i),merr(i),sky(i),
     .      xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),
     .      ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
 10   continue
 
 510  format(F13.8,1X,2(F9.6,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F10.2),1X,F8.3,1X,F6.2,1X,F8.4)
     
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL*8 ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL*8 s,ep
      ave=0.0
      m=min(1000,n)
      do 11 j=1,m
         ave=ave+data(j)
 11   continue
      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      return
      END 
      
c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real*8      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real*8      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine histogram(npt,rp,dp,dp2,np,nbin,nbinmax,bdatax,bdatay,
     .  title,rmed,std,errs,flag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j,npt2,nbin,nbinmax
      integer np(npt),flag(npt)
      real rp(npt),datamin,datamax,bdatax(nbinmax),bdatay(nbinmax),bmax,
     .  rave,errs(6),std,rmed
      double precision dp(npt),dp2(npt),ave,var,med
      character*80 title
      
C     First we remove the average (helps a lot with double -> real)
      j=0
      do 9 i=1,npt
        if(flag(i).eq.0)then
            j=j+1
            dp2(j)=dp(i)
        endif
 9    continue
      npt2=j

      call avevar(dp2,npt2,ave,var)
      rave=real(ave) !convert to real*4
      std=real(sqrt(var))
      
C     Now we convert dp to rp
      j=0 !because we have sigma-clipping, we need a counter.
      do 10 i=1,npt
        if((abs(dp(i)-ave).lt.4.0*sqrt(var)).and.(flag(i).eq.0))then
            j=j+1
            dp2(j)=dp(i)-ave
            rp(j)=real(dp2(j))
        endif
 10   continue
      npt2=j

C     Find median          
      call rqsort(npt2,dp2,np) !changed npt to k
      i=npt2/2
      if(i.le.0) i=1
      med=dp2(np(i))
      rmed=real(med)

C     Find datarange
      datamin=rp(1)
      datamax=rp(1)
      do 12 i=2,npt2
        datamin=min(rp(i),datamin)
        datamax=max(rp(i),datamax)
 12   continue
 
      call bindata(nbin,npt2,rp,bdatax,bdatay,datamin,datamax,bmax)
      
      call pgpage() !fresh plotting surface
      call pgslw(1)
      call pgsch(2.0)
c         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
c         call pgvport(xb1,xb2,yb1,yb2)
      call pgvport(0.2,1.0,0.2,0.9)
c      fc=2.*(datamax-datamin)/real(nbin) !center histograms
      call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size
         
C        Add axis labels
      call pglabel(title,"Relative Probability","")
      call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
      
C     Shift axis scale to account for average removal
      call pgwindow(datamin+rave,datamax+rave,0.,1.0+0.1*1.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
      
      call errorest(npt2,rp,nbin,bdatax,bdatay,bmax,rmed,errs)
      rmed=rmed+rave !correct for average removal
C     Need to recalulate standard deviation at this point!      

      return
      end      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest(npt,dd,nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,i,j,mdpnt,k
      real dd(npt),bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),shgt,
     .  intl,inth,per,perold,e1,e2,inthold,intlold
      logical loop
      
C     find out which bin contains the middle.      
      do 5 i=1,nbin-1
c        write(0,*) frsol,bdatax(i),bdatax(i+1)
        if((bdatax(i).le.frsol).and.(bdatax(i+1).gt.frsol))then
            mdpnt=i
        endif
 5    continue
C     case when "bestfit" is outside histogram (i.e. eccentricity)
      if(frsol.le.bdatax(1)) mdpnt=1
      if(frsol.ge.bdatax(nbin)) mdpnt=nbin 
c      write(0,*) "mdpnt:",mdpnt
      
      perold=0 !initalization of variables
      intlold=0
      inthold=0
      do 14 i=1,6  !initialize variable to zero.
        errs(i)=0.0
 14   continue
C     bmax is the value of the largest histogram bin
      do 13 k=1,bmax-1
        shgt=real(bmax-k)
C       first we move in the forward direction      
        loop=.true. !loop until loop is false
        i=mdpnt  !removed -1
        do 10 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i+1).lt.shgt))then
                loop=.false.
c               inth=bdatax(i)
                inth=(bdatay(i+1)-shgt)/(bdatay(i+1)-bdatay(i))
     .              *(bdatax(i+1)-bdatax(i))+
     .              bdatax(i)
            endif
            if(i.ge.nbin-1) then
                loop=.false.
                inth=bdatax(nbin)
            endif
            i=i+1
c           write(6,*) "i:",i
 10     enddo
       
C       now we move in the reverse direction
        loop=.true. !loop until loop is false
        i=mdpnt+1
        do 11 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i-1).lt.shgt))then
                loop=.false.
c               intl=bdatax(i)
                intl=bdatax(i)-
     .              (bdatay(i)-shgt)/(bdatay(i)-bdatay(i-1))
     .              *(bdatax(i)-bdatax(i-1))
            endif
            if(i.le.2) then
                loop=.false.
                intl=bdatax(1)
            endif
            i=i-1
 11     enddo

        j=0
        do 12 i=1,npt
            if((dd(i).gt.intl).and.(dd(i).lt.inth))then
                j=j+1    
            endif
 12     continue
 
        per=real(j)/real(npt)
        if((per.ge.0.683).and.(perold.lt.0.683))then
            e1=inth-(per-0.683)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.683)/(per-perold)*(intlold-intl)
            errs(1)=e1-frsol
            errs(2)=e2-frsol
c           write(6,*) "1 sigma:",frsol,errs(1),errs(2)
        endif
        if((per.ge.0.954).and.(perold.lt.0.954))then
            e1=inth-(per-0.954)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.954)/(per-perold)*(intlold-intl)
            errs(3)=e1-frsol
            errs(4)=e2-frsol
c           write(6,*) "2 sigma:",frsol,errs(3),errs(4)
        endif 
        if((per.ge.0.9973).and.(perold.lt.0.9973))then
            e1=inth-(per-0.9973)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.9973)/(per-perold)*(intlold-intl)
            errs(5)=e1-frsol
            errs(6)=e2-frsol
c           write(6,*) "3 sigma:",frsol,errs(5),errs(6)
        endif       
        perold=per
        intlold=intl
        inthold=inth
c       write(6,*) per,intl,inth,shgt
 13   continue
      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum

C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)
      
C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i-1)*binsize !added -1
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum
 
      return
      end
 
 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdata(nunit,npt,time,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in data
C     npt - number of data points (returned)
C     id - id of star (returned)
C     xcoo,ycoo - co-ordinates of star (returned)
C     time,mag,merr - the data and associated error.
      implicit none
      integer nmax,nunit
      parameter(nmax=600000)
      integer npt,i,ntest,rn
      double precision time(nmax),mag(nmax),merr(nmax),
     .     sky(nmax),xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),
     .     tboard(nmax),mfield(nmax),ftotal(nmax),aflux(nmax),nap(nmax),
     .     itime(nmax),skystd(nmax)

      i=1

   9  continue
      if(i.gt.1) write(6,*) "error reading",i
  10  read(nunit,*,err=9,end=20) time(i),mag(i),merr(i),sky(i),
     . xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),
     . ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
        merr(i)=abs(merr(i))
        if(merr(i).le.0.0) goto 10
        rn=ntest(mag(i))
        if(rn.eq.1) goto 10
        rn=ntest(sky(i))
        if(rn.eq.1) goto 10
C        
C       This is where the outlier rejection routine will be!! 
C

c         if((yc(i).gt.11.0).or.(yc(i).lt.9.0)) goto 10

c        if((fx(i).lt.1.2d0).and.(fy(i).lt.1.3d0).and.
c     .      (fxy(i).gt.-0.03d0).and.(fxy(i).lt.0.06d0).and.
c     .      (tboard(i).gt.292.0d0).and.(mfield(i).gt.20000.0d0).and.
c     .      (yc(i).gt.39.0d0).and.(yc(i).lt.45.0d0).and.
c     .      (xc(i).gt.10.0d0).and.(xc(i).lt.15.0d0))then
            i=i+1
c        endif
      goto 10
 20   continue

 510  format(F13.8,1X,2(F9.6,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F10.2),1X,F8.3,1X,F6.2,1X,F8.4)
      npt=i-1
      write(6,*) "Points read: ",npt

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function ntest(rr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this routine is used to test for NANs and other nasty numbers.
      implicit none
      double precision rr,thi,tlow

C     this defines the numerical boundary of anything that might be 
C     exciting.

      thi=  99.9d30
      tlow=-99.9d30

      if((rr.gt.tlow).and.(rr.lt.thi))then
         ntest=0
      else
         ntest=1
      endif

      return
      end