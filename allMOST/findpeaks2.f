CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findpeaks2(fitsdata,naxes,datamin,datamax,
     .     stn,stx,sty,pcoo,hasstar)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer naxes(2),xmax,ymax,stmax,nsmax,i,j,k,mx,my,maxiter,ii,
     .  iter,ns,xc,yc
      parameter(xmax=600,ymax=600,stmax=20,nsmax=100,ns=5)
      integer stn,stx(stmax,2),sty(stmax,2),hasstar(stmax),nx,ny
      real fitsdata(xmax,ymax),datamin,datamax,pcoo(stmax,nsmax,2),
     .  mpeak,sumx,sumy,fsum,tempdata(xmax,ymax),pval(stn),oldcoo(2),
     .  meandata,ptemp(2,2),crit
      
C     Avoid infinite loops
      maxiter=50 !maximum number of iterations for find centroid
      

      nx=naxes(1)  !number of x-pixels
      ny=naxes(2)  !number of y-pixels

      do 5 i=1,naxes(1)
        do 6 j=1,naxes(2)
            tempdata(i,j)=fitsdata(i,j)
 6      continue
 5    continue

      do 10 i=1,stn
      do 11 ii=1,hasstar(i)  
      
         mpeak=-1.0e30 !some low number
         sumx=0.0 !summation in x-direction
         sumy=0.0 !summation in y-direction
         fsum=0.0 !sum over all pixels
         meandata=0.0
         do 20 j=1,naxes(1) !loop over all pixels
            do 30 k=1,naxes(2)
C              Now we check to see that the pixels are in the subraster
               if((j.ge.stx(i,1)).and.(j.le.stx(i,2)).and.
     .              (k.ge.sty(i,1)).and.(k.le.sty(i,2))) then
C              Check that fits values are inside datamin/datamax bounds
                  if((tempdata(j,k).gt.datamin).and.
     .                 (tempdata(j,k).lt.datamax))then
                     if(tempdata(j,k).gt.mpeak) then !find max-pixel
                        mpeak=tempdata(j,k)  !store peak pixel value
                        mx=j
                        my=k
                     endif
C                    All of this is pretty much useless.
                     meandata=meandata+tempdata(j,k)
                     fsum=fsum+1.0
c                     sumx=sumx+tempdata(j,k)*real(j)
c                     sumy=sumy+tempdata(j,k)*real(k)
c                     fsum=fsum+tempdata(j,k)
                  endif
               endif
 30         continue
 20      continue
         pval(i)=mpeak
         pcoo(i,ii,1)=real(mx)
         pcoo(i,ii,2)=real(my)
         meandata=meandata/fsum

         iter=0
 42      oldcoo(1)=pcoo(i,ii,1)
         oldcoo(2)=pcoo(i,ii,2)
         iter=iter+1

         sumx=0
         sumy=0
         fsum=0
         do 40 j=1,ns
            do 41 k=1,ns
               xc=int(pcoo(i,ii,1))-3+j
               yc=int(pcoo(i,ii,2))-3+k
               if((xc.ge.stx(i,1)).and.(xc.le.stx(i,2)).and.
     .              (yc.ge.sty(i,1)).and.(yc.le.sty(i,2))) then
                  if((tempdata(xc,yc).gt.datamin).and.
     .                 (tempdata(xc,yc).lt.datamax))then
                     sumx=sumx+tempdata(xc,yc)*real(xc)
                     sumy=sumy+tempdata(xc,yc)*real(yc)
                     fsum=fsum+tempdata(xc,yc)
                     tempdata(xc,yc)=meandata
                  endif
               endif
 41         continue
 40      continue
         
         pcoo(i,ii,1)=sumx/fsum
         pcoo(i,ii,2)=sumy/fsum
c         write(6,*)i,ii,pcoo(i,ii,1),pcoo(i,ii,2)
c         read(5,*)

         if(iter.gt.maxiter) goto 43
         crit=((pcoo(i,ii,1)-oldcoo(1))**2.0 + 
     .      (pcoo(i,ii,2)-oldcoo(2))**2.0)**0.5
         if(crit.gt.1.0) goto 42
 43      continue
         pcoo(i,ii,1)=pcoo(i,ii,1)-real(stx(i,1))+1.0
         pcoo(i,ii,2)=pcoo(i,ii,2)-real(sty(i,1))+1.0

c         write(6,*) i, "peak =",pcoo(i,1),pcoo(i,2),pval(i)

 11   continue
 10   continue
 
      if((pcoo(1,1,1).lt.pcoo(1,2,1)).and.
     .  (pcoo(1,1,2).lt.pcoo(1,2,2)))then
        ptemp(1,1)=pcoo(1,1,1)
        ptemp(1,2)=pcoo(1,1,2)
        pcoo(1,1,1)=pcoo(1,2,1)
        pcoo(1,1,2)=pcoo(1,2,2)
        pcoo(1,2,1)=ptemp(1,1)
        pcoo(1,2,2)=ptemp(1,2)
      endif
      
c      write(6,*) pcoo(1,1,1),pcoo(1,1,2)
c      write(6,*) pcoo(1,2,1),pcoo(1,2,2)
c      read(5,*)

      return
      end
      
      