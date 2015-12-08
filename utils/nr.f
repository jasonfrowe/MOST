C Contents
C
C mnbrak
C brent 
C mrqmin
C mrqcof
C funcs
C medfit
C rofunc
C select
C moment
C spline
C splint
C funcsl
C lfit
C gaussj
C covsrt

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc)
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      implicit none
c      REAL ax,bx,cx,fa,fb,fc,f,GOLD,GLIMIT,TINY
c      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
c      REAL dum,fu,q,r,u,ulim
c      fa=f(ax)
c      fb=f(bx)
c      if(fb.gt.fa)then
c         dum=ax
c         ax=bx
c         bx=dum
c         dum=fb
c         fb=fa
c         fa=dum
c      endif
c      cx=bx+GOLD*(bx-ax)
c      fc=f(cx)
c 1    if(fb.ge.fc)then 
c         r=(bx-ax)*(fb-fc)
c         q=(bx-cx)*(fb-fa)
c         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
c         ulim=bx+GLIMIT*(cx-bx)
c         if((bx-u)*(u-cx).gt.0.)then 
c            fu=f(u)
c            if(fu.lt.fc)then
c               ax=bx
c               fa=fb
c               bx=u
c               fb=fu
c               return
c            else if(fu.gt.fb)then
c               cx=u
c               fc=fu
c               return
c            endif
c            u=cx+GOLD*(cx-bx)
c            fu=f(u)
c         elseif((cx-u)*(u-ulim).gt.0.)then
c            fu=f(u)
c            if(fu.lt.fc)then
c               bx=cx
c               cx=u
c               u=cx+GOLD*(cx-bx)
c               fb=fc
c               fc=fu
c               fu=f(u)
c            endif
c         elseif((u-ulim)*(ulim-cx).ge.0.)then 
c            u=ulim
c            fu=f(u)
c         else
c            u=cx+GOLD*(cx-bx)
c            fu=f(u)
c         endif
c         ax=bx
c         bx=cx
c         cx=u
c         fa=fb
c         fb=fc
c         fc=fu
c         goto 1
c      endif
c      return
c      END


cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      FUNCTION brent(ax,bx,cx,tol,xmin)
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      implicit none
c      INTEGER ITMAX
c      REAL brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
c      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
c      INTEGER iter
c      REAL a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
c      a=min(ax,cx)
c      b=max(ax,cx)
c      v=bx 
c      w=v
c      x=v
c      e=0.
c      fx=f(x)
c      fv=fx
c      fw=fx
c      do 11 iter=1,ITMAX
c         xm=0.5*(a+b)
c         tol1=tol*abs(x)+ZEPS
c         tol2=2.*tol1
c         if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
c         if(abs(e).gt.tol1) then
c            r=(x-w)*(fx-fv)
c            q=(x-v)*(fx-fw)
c            p=(x-v)*q-(x-w)*r
c            q=2.*(q-r)
c            if(q.gt.0.) p=-p
c            q=abs(q)
c            etemp=e
c            e=d
c            if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.
c     .         p.ge.q*(b-x)) goto 1
c            d=p/q
c            u=x+d
c            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
c            goto 2
c         endif
c 1       if(x.ge.xm) then 
c            e=a-x
c         else
c            e=b-x
c         endif
c         d=CGOLD*e
c 2       if(abs(d).ge.tol1) then 
c            u=x+d
c         else
c            u=x+sign(tol1,d)
c         endif
c         fu=f(u)
c         if(fu.le.fx) then
c            if(u.ge.x) then
c               a=x
c            else
c               b=x
c            endif
c            v=w
c            fv=fw
c            w=x
c            fw=fx
c            x=u
c            fx=fu
c         else
c            if(u.lt.x) then
c               a=u
c            else
c               b=u
c            endif
c            if(fu.le.fw .or. w.eq.x) then
c               v=w
c               fv=fw
c               w=u
c               fw=fu
c            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
c               v=u
c               fv=fu
c            endif
c         endif
c 11   continue
c      write(6,*) "Whoops.. too many iterations"
cc      pause
c 3    xmin=x
c      brent=fx
c      return
c      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca, 
     *     chisq,alamda) 
ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,nca,ndata,ia(ma),MMAX 
      REAL alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca), 
     *     sig(ndata),x(ndata),y(ndata) 
      PARAMETER (MMAX=1000) 
      INTEGER j,k,l,mfit 
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX) 
      SAVE ochisq,atry,beta,da,mfit 
      if(alamda.lt.0.)then 
         mfit=0 
         do 11 j=1,ma 
            if (ia(j).ne.0) mfit=mfit+1 
 11      enddo  
         alamda=0.001 
         call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq) 
         ochisq=chisq 
         do 12 j=1,ma 
            atry(j)=a(j) 
 12      enddo 
      endif 
      do 14 j=1,mfit 
         do 13 k=1,mfit
            covar(j,k)=alpha(j,k) 
 13      enddo 
         covar(j,j)=alpha(j,j)*(1.+alamda) 
         da(j)=beta(j) 
 14   enddo 
      call gaussj(covar,mfit,nca,da,1,1)  

      if(alamda.eq.0.)then 
         call covsrt(covar,nca,ma,ia,mfit) 
         call covsrt(alpha,nca,ma,ia,mfit) 
         return 
      endif 
      j=0 
      do 15 l=1,ma 
         if(ia(l).ne.0) then 
            j=j+1 
            atry(l)=a(l)+da(j) 
         endif 
 15   enddo 
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq) 
      if(chisq.lt.ochisq)then 
         alamda=0.1*alamda 
         ochisq=chisq 
         do 17 j=1,mfit 
            do 16 k=1,mfit 
               alpha(j,k)=covar(j,k) 
 16         enddo 
            beta(j)=da(j) 
 17      enddo 
         do 18 l=1,ma 
            a(l)=atry(l) 
 18      enddo 
      else 
         alamda=10.*alamda 
         chisq=ochisq 
      endif 
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp, 
     *     chisq) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,nalp,ndata,ia(ma),MMAX 
      REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata), 
     *     y(ndata) 
      PARAMETER (MMAX=1000) 
      INTEGER mfit,i,j,k,l,m 
      REAL dy,sig2i,wt,ymod,dyda(MMAX) 
      mfit=0 
      do 11 j=1,ma 
         if (ia(j).ne.0) mfit=mfit+1 
 11   enddo 
      do 13 j=1,mfit 
         do 12 k=1,j
            alpha(j,k)=0. 
 12      enddo 
         beta(j)=0. 
 13   enddo 
      chisq=0. 
      do 16 i=1,ndata 
         call funcs(x(i),a,ymod,dyda,ma)
         sig2i=1./(sig(i)*sig(i)) 
         dy=y(i)-ymod 
         j=0 
         do 15 l=1,ma 
            if(ia(l).ne.0) then 
               j=j+1 
               wt=dyda(l)*sig2i 
               k=0 
               do 14 m=1,l 
                  if(ia(m).ne.0) then 
                     k=k+1 
                     alpha(j,k)=alpha(j,k)+wt*dyda(m) 
                  endif 
 14            enddo 
               beta(j)=beta(j)+dy*wt 
            endif 
 15      enddo 
         chisq=chisq+dy*dy*sig2i  
 16   enddo
      do 18 j=2,mfit 
         do 17 k=1,j-1 
            alpha(k,j)=alpha(j,k) 
 17      enddo
 18   enddo
      return 
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE funcs(x,a,y,dyda,na) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER na,nstarmax
      parameter(nstarmax=20)
      REAL x,y,a(na),dyda(na) 
      real arg, cosarg,w(nstarmax),sinarg
      integer i,j,nstar,k,cc
      common /period/ nstar


      y=0
      j=(na-nstar)/nstar
      do 20 k=1,nstar
         cc=(j+1)*(k-1)+2
         dyda(cc-1)=0.
         do 21 i=cc,cc+j-2,2
            arg=real((i-cc+2)/2)*a(cc-1)*x+a(i+1)
            cosarg=cos(arg)
            sinarg=sin(arg)
            y=y+cosarg*a(i)
            dyda(i)=cosarg
            dyda(i+1)=-a(i)*sinarg
            dyda(cc-1)=dyda(cc-1)-a(i)*real((i-cc+2)/2)*x*sinarg
 21      continue
 20   continue

c      j=(na-nstar)/nstar
c      do 42 k=1,nstar
c         cc=(j+1)*(k-1)+2
c         write(6,500) "Amp", (a(i)  ,i=cc,cc+j-2,2)
c         write(6,500) "Psi", (a(i+1),i=cc,cc+j-2,2)
c 42   continue

 500  format(A3,1X,20(F9.6))

      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE medfit(x,y,ndata,a,b,abdev)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ndata,NMAX,ndatat 
      PARAMETER (NMAX=600000) 
      REAL a,abdev,b,x(ndata),y(ndata), 
     *     arr(NMAX),xt(NMAX),yt(NMAX),aa,abdevt 
      COMMON /arrays/ xt,yt,arr,aa,abdevt,ndatat 
C USES rofunc 
      INTEGER j 
      REAL b1,b2,bb,chisq,del,f,f1,f2,sigb,sx,sxx,sxy,sy,rofunc 
      sx=0. 
      sy=0. 
      sxy=0. 
      sxx=0. 
      do 11 j=1,ndata 
         xt(j)=x(j) 
         yt(j)=y(j) 
         sx=sx+x(j) 
         sy=sy+y(j) 
         sxy=sxy+x(j)*y(j) 
         sxx=sxx+x(j)**2 
 11   enddo 
      ndatat=ndata 
      del=ndata*sxx-sx**2 
      aa=(sxx*sy-sx*sxy)/del 
      bb=(ndata*sxy-sx*sy)/del 
      chisq=0.
      do 12 j=1,ndata 
         chisq=chisq+(y(j)-(aa+bb*x(j)))**2 
 12   enddo 
      sigb=sqrt(chisq/del) 
      b1=bb 
      f1=rofunc(b1) 
      if(sigb.gt.0.)then 
         b2=bb+sign(3.*sigb,f1) 
         f2=rofunc(b2) 
         if(b2.eq.b1)then 
            a=aa 
            b=bb 
            abdev=abdevt/ndata 
            return 
         endif 
 1       if(f1*f2.gt.0.)then 
            bb=b2+1.6*(b2-b1) 
            b1=b2 
            f1=f2 
            b2=bb 
            f2=rofunc(b2) 
            goto 1 
         endif 
         sigb=0.01*sigb 
 2       if(abs(b2-b1).gt.sigb)then 
            bb=b1+0.5*(b2-b1) 
            if(bb.eq.b1.or.bb.eq.b2)goto 3 
            f=rofunc(bb) 
            if(f*f1.ge.0.)then 
               f1=f 
               b1=bb 
            else 
               f2=f 
               b2=bb 
            endif 
            goto 2 
         endif
      endif 
 3    a=aa 
      b=bb 
      abdev=abdevt/ndata 
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION rofunc(b) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER NMAX 
      REAL rofunc,b,EPS 
      PARAMETER (NMAX=600000,EPS=1.e-7) 
C USES select 
      INTEGER j,ndata 
      REAL aa,abdev,d,sum,arr(NMAX),x(NMAX),y(NMAX),select 
      COMMON /arrays/ x,y,arr,aa,abdev,ndata 
      do 11 j=1,ndata 
         arr(j)=y(j)-b*x(j) 
 11   enddo 
      if (mod(ndata,2).eq.0) then 
         j=ndata/2 
         aa=0.5*(select(j,ndata,arr)+select(j+1,ndata,arr)) 
      else 
         aa=select((ndata+1)/2,ndata,arr) 
      endif 
      sum=0. 
      abdev=0. 
      do 12 j=1,ndata 
         d=y(j)-(b*x(j)+aa) 
         abdev=abdev+abs(d) 
         if (y(j).ne.0.) d=d/abs(y(j)) 
         if (abs(d).gt.EPS) sum=sum+x(j)*sign(1.0,d) 
 12   enddo 
      rofunc=sum 
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION select(k,n,arr) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER k,n 
      REAL select,arr(n) 
      INTEGER i,ir,j,l,mid 
      REAL a,temp 
      l=1 
      ir=n 
 1    if(ir-l.le.1)then 
         if(ir-l.eq.1)then 
            if(arr(ir).lt.arr(l))then 
               temp=arr(l) 
               arr(l)=arr(ir) 
               arr(ir)=temp 
            endif 
         endif 
         select=arr(k) 
         return 
      else 
         mid=(l+ir)/2 
         temp=arr(mid) 
         arr(mid)=arr(l+1) 
         arr(l+1)=temp 
         if(arr(l).gt.arr(ir))then 
            temp=arr(l) 
            arr(l)=arr(ir) 
            arr(ir)=temp 
         endif 
         if(arr(l+1).gt.arr(ir))then 
            temp=arr(l+1) 
            arr(l+1)=arr(ir) 
            arr(ir)=temp 
         endif 
         if(arr(l).gt.arr(l+1))then 
            temp=arr(l) 
            arr(l)=arr(l+1) 
            arr(l+1)=temp 
         endif 
         i=l+1 
         j=ir 
         a=arr(l+1) 
 3       continue
         i=i+1 
         if(arr(i).lt.a)goto 3 
 4       continue 
         j=j-1
         if(arr(j).gt.a)goto 4 
         if(j.lt.i)goto 5 
         temp=arr(i) 
         arr(i)=arr(j) 
         arr(j)=temp 
         goto 3 
 5       arr(l+1)=arr(j) 
         arr(j)=a 
         if(j.ge.k)ir=j-1 
         if(j.le.k)l=i 
      endif 
      goto 1 
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n
      REAL adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      REAL p,s,ep
      if(n.le.1)pause 'n must be at least 2 in moment'
      s=0.
      do 11 j=1,n
         s=s+data(j)
 11   enddo
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         adev=adev+abs(s)
         p=s*s
         var=var+p
         p=p*s
         skew=skew+p
         p=p*s
         curt=curt+p
 12   enddo
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
         skew=skew/(n*sdev**3)
         curt=curt/(n*var**2)-3.
      else
         pause 'no skew or kurtosis when zero variance in moment'
      endif
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE spline(x,y,n,yp1,ypn,y2) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n,NMAX 
      REAL yp1,ypn,x(n),y(n),y2(n) 
      PARAMETER (NMAX=600000) 
      INTEGER i,k 
      REAL p,qn,sig,un,u(NMAX) 
      if (yp1.gt..99e30) then 
         y2(1)=0. 
         u(1)=0. 
      else 
         y2(1)=-0.5 
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      endif 
      do 11 i=2,n-1 
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1)) 
         p=sig*y2(i-1)+2. 
         y2(i)=(sig-1.)/p 
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) 
     .        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p 
 11   continue
      if (ypn.gt..99e30) then 
         qn=0. 
         un=0. 
      else 
         qn=0.5 
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif 
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.) 
      do 12 k=n-1,1,-1 
         y2(k)=y2(k)*y2(k+1)+u(k) 
 12   continue
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE splint(xa,ya,y2a,n,x,y) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n 
      REAL x,y,xa(n),y2a(n),ya(n) 
      INTEGER k,khi,klo 
      REAL a,b,h 
      klo=1 
      khi=n 
 1    if (khi-klo.gt.1) then 
         k=(khi+klo)/2 
         if(xa(k).gt.x)then 
            khi=k 
         else 
            klo=k 
         endif 
      goto 1 
      endif 
      h=xa(khi)-xa(klo) 
      if (h.eq.0.) pause 'bad xa input in splint'  
      a=(xa(khi)-x)/h 
      b=(x-xa(klo))/h 
      y=a*ya(klo)+b*ya(khi)+ 
     .     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6. 
      return 
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FUNCSl(X,P,NP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
      DIMENSION P(NP)
      P(1)=1.
      DO 11 J=2,NP
         P(J)=P(J-1)*(X)
c         P(J)=P(J-1)/X
 11   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,ia(ma),npc,ndat,MMAX 
      REAL chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat) 
      EXTERNAL funcsl
      PARAMETER (MMAX=50)

      INTEGER i,j,k,l,m,mfit 
      REAL sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX) 
      mfit=0 
      do 11 j=1,ma 
         if(ia(j).ne.0) mfit=mfit+1 
 11   continue
      if(mfit.eq.0) pause 'lfit: no parameters to be fitted'
      do 13 j=1,mfit 
         do 12 k=1,mfit 
            covar(j,k)=0. 
 12      continue
         beta(j)=0. 
 13   continue
      do 17 i=1,ndat
         call funcsl(x(i),afunc,ma) 
         ym=y(i) 
         if(mfit.lt.ma) then 
            do 14 j=1,ma 
               if(ia(j).eq.0) ym=ym-a(j)*afunc(j) 
 14         continue
         endif 
         sig2i=1./sig(i)**2 
         j=0 
         do 16 l=1,ma 
            if (ia(l).ne.0) then 
               j=j+1 
               wt=afunc(l)*sig2i 
               k=0 
               do 15 m=1,l 
                  if (ia(m).ne.0) then 
                     k=k+1 
                     covar(j,k)=covar(j,k)+wt*afunc(m) 
                  endif 
 15            continue
               beta(j)=beta(j)+ym*wt 
            endif 
 16      continue 
 17   continue 
      do 19 j=2,mfit 
         do 18 k=1,j-1 
            covar(k,j)=covar(j,k) 
 18      continue
 19   continue 
      call gaussj(covar,mfit,npc,beta,1,1) 
      j=0 
      do 21 l=1,ma 
         if(ia(l).ne.0) then 
            j=j+1 
            a(l)=beta(j) 
         endif 
 21   continue
      chisq=0.
      do 23 i=1,ndat 
         call funcsl(x(i),afunc,ma) 
         sum=0. 
         do 22 j=1,ma 
            sum=sum+a(j)*afunc(j) 
 22      continue 
         chisq=chisq+((y(i)-sum)/sig(i))**2 
 23   continue 
      call covsrt(covar,npc,ma,ia,mfit)
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE gaussj(a,n,np,b,m,mp) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER m,mp,n,np,NMAX 
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=50) 
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),
     *     ipiv(NMAX) 
      REAL big,dum,pivinv
      do 11 j=1,n 
         ipiv(j)=0
 11   enddo 
      do 22 i=1,n
         big=0. 
         do 13 j=1,n
            if(ipiv(j).ne.1)then
               do 12 k=1,n
                  if (ipiv(k).eq.0) then
                     if (abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  else if (ipiv(k).gt.1) then
                     pause 'singular matrix in gaussj'
                  endif
 12            enddo
            endif
 13      enddo
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then
            do 14 l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
 14         enddo
            do 15 l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
 15         enddo
         endif
         indxr(i)=irow
         indxc(i)=icol
         if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
         pivinv=1./a(icol,icol)
         a(icol,icol)=1.
         do 16 l=1,n
            a(icol,l)=a(icol,l)*pivinv
 16      enddo
         do 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
 17      enddo
         do 21 ll=1,n
            if(ll.ne.icol)then
               dum=a(ll,icol)
               a(ll,icol)=0.
               do 18 l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
 18            enddo
               do 19 l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
 19            enddo
            endif
 21      enddo
 22   enddo
      do 24 l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            do 23 k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
 23         enddo
         endif
 24   enddo
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,mfit,npc,ia(ma)
      REAL covar(npc,npc)
      INTEGER i,j,k
      REAL swap
      do 12 i=mfit+1,ma
         do 11 j=1,i
            covar(i,j)=0.
            covar(j,i)=0.
 11      enddo
 12   enddo
      k=mfit
      do 15 j=ma,1,-1
         if(ia(j).ne.0)then
            do 13 i=1,ma
               swap=covar(i,k)
               covar(i,k)=covar(i,j)
               covar(i,j)=swap
 13         enddo
            do 14 i=1,ma
               swap=covar(k,i)
               covar(k,i)=covar(j,i)
               covar(j,i)=swap
 14         enddo
            k=k-1
         endif
 15   enddo
      return
      END
