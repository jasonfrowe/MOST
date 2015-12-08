C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      program allmost
C     Writen by Jason Rowe (rowe@astro.ubc.ca)
C     Last Updated: 20110622
C     Performs photometry on MOST SDS2 fits files.

C     Version 0.95c

C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THIS PROGRAM IS CURRENTLY, "USE AT YOUR OWN RISK!!!"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Log:
C-----------------------------------------------------------------------
C     Jun 22   Fabry output name is now properly named
C     Jun 04   Added GuideStar Photometry Output
C     Jan 28   This program has been written with the idea of 
C            a nice toolbox of subroutines for use in any 
C            program that deal with MOST STSD2 data.  This is v0.5
C     Feb 03   Fixed gpsf.  Bug where some profiles removed sky and 
C            others did not (it should not!)
C              PSF fitting routines are sorta working.  Not stable to 
C            make it though entire dataset, but getting there. 
C     Feb 10   Adding multi-object fits to correct for guiding errors
C            This is currently #3 on the to do list.  As this is a 
C            major add on, upping to v0.6a
C     Feb 16   Fixed labeling of pixel numbers in plots.  Added
C            NGauss and NMoffat functions for fitting.  Working on
C            image deconvolution with median PSF
C     Feb 17 Added new centroid and sky estimation routines
C     Feb 18 Deconvolution is now working.
C     Aug 13 Moved to allMOST version.  Too many changes to note!
C     Mar 01 Updated to handle new FITS format
C     Mar 30 Fixed code to properly handle hr1217
C     Mar 31 Fixed pcool segmentation fault. (for 51 Peg)
C            Added exposure time correction
C            Added check for when there are zero pixels to fit.
C              this is important when image is saturated. See NSTAR
C     Apr 02 Added TIMECORR
C     Feb 16 dates are now HJD, program also does Fabry photometry now
C            made improvement to starident for centroiding.
C              it is not senitive to the users clicking ability. 
C     Feb 17 Added sky standard deviation output, useful for tracking
C              cross talk between CCDs  
C     Feb 23 Not having a flatfield file does not result in an error
C     Jun 18 Bug fixes in centroids routine.
C     Sep 18 Added start of PGPLOT control so CPU is not wasted with
C              graphics output.
C
C-----------------------------------------------------------------------
C
C     To Do List:
C     no         2. Add PSF lookup tables
C     no         5. Make program interactive
C     no         7. Remove "Numerical Recipes" routines
C
C CODE BEGINS HERE -------- CODE BEGINS HERE -------- CODE BEGINS HERE C
      
C     NSMAX controls the maximum number of stars that can be 
C      simultaniously fit. 
C     STMAX controls the maximum number of secondary subrasters
C     NRECMAX controls the number of fits header entries.
C     XMAX,YMAX control the maximum dimensions of the fits file.
      implicit none
      integer stmax,readwrite,unitfits,status,nkeys,nrecmax,np,
     .     naxes(2),tix(2),tiy(2),xmax,ymax,stn,nsmax,ns,nfitm,nfitpsf,
     .     fitsky,i,stfrm,j,nbegin,id,ifun1,datasty,coott,nstack,gns,
     .     gsfields,jj
      parameter(stmax=20,nrecmax=700,xmax=600,ymax=600,nsmax=100,
     .     nfitm=200)
      integer hasstar(stmax),stx(stmax,2),sty(stmax,2),rstx(stmax,2),
     .     rsty(stmax,2),ifun,bpeak(stmax),munit(nsmax),sint,
     .     sflag(stmax),iplot,iplotted,gsunit(nsmax)
      real gain,rnoise,sat(stmax),fitsdata(xmax,ymax),datamin,datamax,
     .     pcoo(stmax,nsmax,2),pcool(stmax,nsmax,2),pval(stmax),fstar,
     .     ans(nfitm),rad,skyrad,rflux,zsky,stdsky,fitrad,ftotal,
     .     pmag(stmax,nsmax),peaks(stmax,2),date,tboard,
     .     boardtemp,magfield,mfield,ans1(nfitm),mag,merr,xcoo,
     .     ycoo,flat(xmax,ymax),dpix(xmax,ymax),skyprint,magprint,nfl,
     .     merprint,pans(5),skymap(xmax,ymax),prad,itime,exptime,
     .     expcor,fabmag,faberr,fabsky,fabans(nfitm),fabflat(xmax*ymax),
     .     fabdpix(xmax*ymax),pstdsky,gsmag,gsmerr
      double precision readdate,ddate
      character*80 listname,oname,header(nrecmax),filename,flatfile,
     .     basename,command,datestamp
      integer tt
      character*80 pnum,pfile

      open(unit=17,file='hdcor.dat')
 
C     iplot controls the display.  It's a waste of CPU to display
C     every frame.
C     So: iplot=0 : waste CPU and plot away!
C         iplot=-1 : only plot the first frame
C         iplot= n : plot every nth frame
      iplot=-1
      iplotted=0 !counts the number of frames between plots

C     This variable is for logging errors.  As long as it is zero, then
C     everything is fine.
      status=0

C     just marking how many frames we have gone through
      nbegin=0

C     Reads in parameter file and/or sets defaults
      call setpars(gain,rnoise,listname,sat,oname,flatfile,ifun,rad,
     .   prad,skyrad,fitsky,fitrad,basename,sint,datasty,coott)

C     This reads in the flatfield map
      call flatmap(flat,dpix,flatfile)
      call fabrymap(fabflat,fabdpix,"fabcor.dat")

C     Flat fielding data goes into ftest.dat
      open(unit=14,file="ftest.dat")
C     Photometry goes into "oname"
      open(unit=13,file=oname)
C     Listname contains the list of files to extract photometry from
      open(unit=10,file=listname,status='old',err=901)

CP    lines commented with cp are added if you want to produce a
CP    series of postscript images.  I used this individual frames
CP    to produce movies of data reduction "in action".  It is neat.
CP    If you want to use it, uncomment all lines starting with "cp"
CP    and you have to comment two PGPLOT subroutines.
CP    You can identify PGPLOT subroutines as they all start 
CP    with "pg".  You will have to do this twice.  Comment out the first
CP    and only occurance of pgopen below.  It is replaced below, so that
CP    a new postscipt image is produced for each frame.  The second line
CP    to comment out is "pgclos" occuring at the end of the MAIN 
CP    program.  It is moved into the main loop to close out the 
CP    current postscript file before we start a few one.
CP    I might add precompiler flags to have this done with a compiler 
CP    flag.  Bug me if you want this and I will do it.
      call pgopen('/xwindow')
cc      call pgopen('ccdlayout.ps/cps')
      call pgask(.false.)

C     depreciated
C     tt is a counter used by the subroutine MAKEFLAT.  It is increased
C     by one everytime MAKEFLAT is run from the subroutine APPHOT
C     Just lets me know how many frames have been averaged together.
c      tt=0

C     np is a counter used for naming the postscipts files if want to
C     make a movie (see CP comments above)
      np=0

C     Mainloop
 400  read(10,500,end=401) filename
 500  format(A80)
        nbegin=nbegin+1

CP    MOVIESTUFF-START  - do not worry about this stuff.  For JR only
cp      tt=tt+1
cp      np=np+1
cp      if(np.lt.10) then
cp         write(pnum,501) "0000",np,".png/png"
cp      elseif((np.ge.10).and.(np.lt.100))then
cp         write(pnum,502) "000",np,".png/png"
cp      elseif((np.ge.100).and.(np.lt.1000))then
cp         write(pnum,503) "00",np,".png/png"
cp      elseif((np.ge.1000).and.(np.lt.10000))then
cp         write(pnum,504) "0",np,".png/png"
cp      else
cp         write(pnum,513) np,".png/png"
cp      endif
cp 501  format(A4,I1,A8)
cp 502  format(A3,I2,A8)
cp 503  format(A2,I3,A8)
cp 504  format(A1,I4,A8)
cp 513  format(I5,A8)
cp      pfile="/media/Etobicoke/MOST/163_55Cnc_2011/"//pnum
cp      call pgopen(pfile)
cp      call pgpage()
cp      call pgsvp(0.0,1.0,0.0,1.0)
cp      call pgsci(0)
cp      call pgsfs(1)
cp      call pgrect(0.0,1.0,0.0,1.0)
cp      call pgsci(1)
cp
CP    MOVIESTUFF-END

C     makes lines thicker when printing
c       call pgslw(5)

CA    IN THE OLD FITS FORMAT::
CA    Some MOST fits files are missing NAXISn so be careful 
CA    I am using wfits in IRAF to prep the fits for usage.
CA    it will have status=223  (use hedit in IRAF)
        if(datasty.eq.1) then

C     This routine opens the fits file and assigns a unit number
         call openfits(readwrite,unitfits,filename,status)
         if(status.ne.0) write(6,506) "openfits error :",status
 506     format(A16,I3)
C     This routine reads in all of the header information
         call readheader(unitfits,status,header,nkeys)
         if(status.ne.0) write(6,507) "readheader error :",status
 507     format(A18,I3)
C     This routine reads in all the pixel information from the fits file
         call readfits(unitfits,status,fitsdata,naxes,datamin,datamax)
         if(status.ne.0) write(6,508) "readfits error :",status
 508     format(A16,I3)
C     This routine closes access to the fits file and releases the unit
C     number
         call closfits(unitfits,status)
         if(status.ne.0) write(6,509) "closfits error :",status
 509     format(A16,I3)

C     Read info about location of primary image
         call readti(header,nrecmax,nkeys,tix,tiy)
C     Read info about location of secondary images
         call readst(header,nrecmax,nkeys,stn,stx,sty,rstx,rsty)

        elseif(datasty.eq.0) then
         call fitsreader(filename,stn,stx,sty,rstx,rsty,tix,tiy,nkeys,
     .      header,fitsdata,datamin,datamax,naxes)
        endif
        
      if(stn.eq.0) then
         write(6,*) "filename:",filename
         pause
      endif         
          

C     read in the date of observation from the header
        ddate=readdate(header,nrecmax,nkeys,datasty)
c        if(ddate.eq.0.) then
c           command="/bin/rm "//filename
c           call system(command)
c           goto 400
c        endif
C     get the CCD board temperature
        tboard=boardtemp(header,nrecmax,nkeys)
C     get the measurement of the earths magnetic field
        mfield=magfield(header,nrecmax,nkeys)
C     get the exposure time
        itime=exptime(header,nrecmax,nkeys,nstack)
C     correct exposure time for shutterless system
C     expcor will probably depend of the position of the subraster on
C     the CCD, but keep it simple for now
c        expcor=0.0192*real*nstack
c        itime=itime+expcor

C     Set saturation level
        call satset(sat,nstack)

C     Correct Bias and Dark current
        call bdcor(header,nrecmax,nkeys,fitsdata,naxes,datamin,datamax,
     .       sat,stmax,stn,stx,sty,rstx,rsty,flat,dpix,sflag,nstack)

C     Display the image
        iplotted=iplotted+1 !increase counter for frames between plots
        if((nbegin.eq.1).or.(iplot.eq.iplotted))then 
            call displayfits(fitsdata,naxes,datamin,datamax,tix,tiy,stn,
     .          stx,sty,nstack)
        endif

C     pick the stars that you are interested in.
        if (nbegin.eq.1) then
        
C     gsfields counts the number of guidestar fields
           gns=gsfields(header,nrecmax,nkeys)
           write(0,*) "Guide Stars: ",gns
           call gsfiles(gns,gsunit,basename) !open gsfiles for output
        
           call starident(fitsdata,naxes,datamin,datamax,
     .          stn,stx,sty,pcoo,hasstar,rstx,rsty,ns)

C     TEMP-this line only!
C           call tempfunction(pcoo,hasstar)

C     ns contains the total number of stars.
           ns=0
           do 5 i=1,stn
              ns=ns+hasstar(i)
 5         continue

C     magfiles opens up the photometry files for writing.
C     if the file already exists it WILL be overwritten.
           call magfiles(ns,munit,basename)

        endif

C     Find location of stars in image
        if (nbegin.gt.1) then
           call findpeaks(fitsdata,naxes,datamin,datamax,
     .       stn,stx,sty,peaks,hasstar,pval)
           call coooffset(stn,pcoo,pcool,peaks,bpeak,hasstar,coott)
C          This routine was special for HD80606
c           call findpeaks2(fitsdata,naxes,datamin,datamax,
c     .      stn,stx,sty,pcoo,hasstar)
        endif

C     Mark the positions of the stars.
        if((nbegin.eq.1).or.(iplot.eq.iplotted))then
            call displaypeaks(pcoo,stn,naxes,hasstar,stx,sty)
            write(datestamp,514) ddate
 514        format(F13.8)
            call pgtext(23.0,57.0,datestamp)
			if(iplot.gt.0) iplotted=0 !reset counter
        endif

C     Skyfitting routine
c        call skyfit(fitsdata,naxes,datamin,datamax,stn,stx,sty,
c     .       pcoo,stfrm,skymap,ans,nfitpsf,ifun)        

C     extract PSF fits from all the identified sources.

        ns=0
        do 10 i=1,stn
           if(hasstar(i).gt.0) then
              stfrm=i
c              call skyfit(fitsdata,naxes,datamin,datamax,stn,stx,sty,
c     .          pcoo,stfrm,skymap,ans1,nfitpsf,ifun1)
              call nstar(fitsdata,naxes,datamin,datamax,stn,stx,sty,
     .             hasstar,ns,ans,ifun,nfitpsf,zsky,stdsky,sat,
     .             fitsky,fitrad,skyrad,pcoo,rstx,rsty,stfrm,sflag,sint)
              if((sflag(i).eq.0).or.((sint.eq.1).and.(sflag(i).eq.1)))
     .             then
                 call psfint(pmag,prad,nfitpsf,ans,zsky,rflux,gain,
     .                rnoise,ifun,hasstar,stfrm,itime)
                 ans1(1)=ans(1)
                 ans1(2)=ans(2)
                 ans1(3)=ans(3)
                 ans1(4)=ans(4)
                 do 11 j=1,hasstar(stfrm)
                    ns=ns+1
                    id=3*j
                    ans1(5)=ans(id+2)
                    ans1(6)=ans(id+3)
                    ans1(7)=ans(id+4)
                    ifun1=ifun-4
                    xcoo=pcoo(i,j,1)
                    ycoo=pcoo(i,j,2)
                   call apphot2(mag,merr,stfrm,rad,skyrad,ans1,fitsdata,
     .                   naxes,datamin,datamax,ifun1,nfitpsf,gain,
     .                   rnoise,zsky,stx,sty,xcoo,ycoo,rflux,rstx,rsty,
     .                   itime,nfl,ftotal,fstar,nstack)
                    if((i.eq.1).and.(j.eq.1)) then
                       skyprint=zsky
                       pstdsky=stdsky
                       magprint=mag
                       merprint=merr
                       pans(1)=ans(3*(j-1)+6)
                       pans(2)=ans(3*(j-1)+7)
                       pans(3)=ans(2)
                       pans(4)=ans(3)
                       pans(5)=ans(4)
                    endif
                    rflux=0.
                    pmag(i,j)=25.0-2.5*log10((pmag(i,j)+
     .                   rflux)/itime*gain)
c                    mfield=fitsdata(12,72)-zsky
c                    ftotal=pval(stfrm)
                    if(ftotal.gt.9999999.0) ftotal=9999999.0
c                    write(munit(ns),510) ddate,pmag(i,j),abs(merr),zsky,
                    jj=j
                    write(munit(ns),510) ddate,mag,abs(merr),zsky,
     .                   ans(3*(jj-1)+6),ans(3*(jj-1)+7),ans(2),ans(3),
     .                   ans(4),tboard,mfield,ftotal,fstar,nfl,itime,
     .                   stdsky
 11              continue
              else
                 do 12 j=1,hasstar(stfrm)
                    ns=ns+1
 12              continue
              endif
           endif
 10     continue
 
        write(6,512) filename(8:20),magprint,pmag(1,1),merprint,
     .   skyprint,pans(1),pans(2),pans(3),pans(4),pans(5),mfield,nfl,
     .   nstack
c        read(5,*)        

C     fabry photometry
c         write(6,*) "fabtime"
        call fabphot(naxes,tix,tiy,fitsdata,itime,gain,rnoise,fabmag,
     .   faberr,fabsky,fabans,fabflat,fabdpix,nstack,nbegin,nfl)
c         write(6,*) "fabdone"
c        write(13,511) ddate,fabmag,faberr,fabsky,(fabans(i),i=1,4),itime
        write(13,510) ddate,fabmag,faberr,fabsky,(pans(j),j=1,5),tboard,
     .      mfield,0.0,fstar,nfl,itime,0.0

C       Export guide-stars
        do 13 i=1,gns
            call getgsmag(gsmag,gsmerr,fstar,nfl,itime,header,nrecmax,
     .          nkeys,i,gain,rnoise,skyprint,nstack)
            write(gsunit(i),510) ddate,gsmag,abs(gsmerr),skyprint,
     .          (pans(j),j=1,5),tboard,mfield,0.0,fstar,nfl,itime,0.0
 13     continue
   

 505  format(F10.5,1X,3(F7.4,1X),F8.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2)
 512  format(A13,1X,3(F7.4,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F9.2,1X,F6.2,1X,I3)
 510  format(F13.8,1X,2(F9.6,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F10.2),1X,F8.3,1X,F6.2,1X,F8.4)
 511  format(F13.8,1X,2(F9.6,1X),F9.2,5(1X,F10.3),1X,F6.2)
        if(nbegin.eq.1) call fbpeak(stn,pmag,hasstar,bpeak,pcoo,pcool)

CP    MOVIESTUFF-START
cp        call pgclos()
CP    MOVIESTUFF-END
      goto 400
CCCCCCCCCCCCCCCCCCCCCC
C        end of main loop
CCCCCCCCCCCCCCCCCCCCCC
 401  continue

      call pgclos()
      close(13)
      close(14)
      close(10)
      close(17)
      call cmagfiles(ns,munit)

      goto 999
 901  write(6,*) "Cannot open ",listname
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine satset(sat,nstack)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,stmax,nstack
      parameter(stmax=20)
      real sat(stmax)
      
      
      do 30 i=1,stmax
         sat(i)=5000000.0!real(nstack)*2.0**14
 30   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine setpars(gain,rnoise,listname,sat,oname,flatfile,ifun,
     .   rad,prad,skyrad,fitsky,fitrad,basename,sint,datasty,coott)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer stmax,ifun,datasty,nline,coott,i,ln
      parameter(stmax=20)
      integer fitsky,sint,match
      real gain,rnoise,sat(stmax),rad,skyrad,fitrad,prad
      character*80 listname,oname,flatfile,basename,line,command,argc

C     first I set all the default parameter values.  Then I read 
C     everything in from allmost.opt

C     Setting the data format. 
C     0 - NEW MEF fits format
C     1 - old fits format (pre-archive)
      datasty=0

C     Setting ifun sets general analytic function
C     ifun = 5 NGauss
C     ifun = 6 NMoffat
      ifun=6
C     prad: set the radius to sum up the PSF photometry
      prad=15.0
C     rad set the radius to sum up photometry
      rad=6.0
C     skyrad sets the skyradius to determine the sky value
      skyrad=6.0
C     set fitsky=0 to not fit sky in NSTAR
C     set fitsky=1 to fit sky in NSTAR
      fitsky=0
C     radius for PSF fit.
      fitrad=3.0

C     sint decides on what to do with satuation.
C     set sint=0 to avoid frame when saturation present
C     set sint=1 to attemp PSF fit using the wings of the star
      sint=1


C     base name for photometry output.
      basename="mostlc"

C     CCD characteristics
C     gain in e/ADU
      gain=6.1
C     readnoise in ADU
      rnoise=1.1 

C     Not in header. 
C     These variables set the saturation level for each subraster
C     MOVE INTO PARAMETER FILE
c      sat(1)=16384.
c      sat(2)=16384.
c      sat(3)=16384.
c      sat(4)=16384.
c      sat(5)=16384.
c      sat(6)=16384.
c      sat(7)=16384.
c      sat(8)=16384.
c      sat(9)=16384.     

      do 30 i=1,stmax
         sat(i)=2.0**14
 30   continue


C     file that contains the list of files for photometry
      listname="files.list"

cC     Output filename
cc      oname="fab.dat"
c      ln=len_trim(basename)
c      oname=basename(1:ln)//"_fab.dat" 

C     file containing flatfield information.
      flatfile="flat.dat"

C     do we want to refine centroids? 0=yes, 1=no
C     (no is conservative and always works)
      coott=1
      
C     Default parameters have been set, so lets read in the opt file
C     count number of lines
      nline=0
      open(unit=16,file="allmost.pars",status='old',err=901)
 10   read(16,500,end=20,err=902) line
      match=0
      nline=nline+1
C     comments have a # symbol in the first column.
      if(line(1:1).eq."#") goto 10
C     if it is not a comment, read in the command and arguement      
      read(line,*) command,argc
C     if the command is ifun, we parse it
      if(command.eq."ifun") then
         match=1
         read(argc,*,err=902) ifun
         if((ifun.lt.5).or.(ifun.gt.6)) then
            write(6,*) "ifun value is invalid"
            write(6,*) "Setting to Moffat"
            ifun=6
         endif
      endif
C     if the command is fitrad we parse it
      if(command.eq."fitrad") then
         match=1
         read(argc,*,err=902) fitrad
      endif
C     if the command is prad we parse it
      if(command.eq."prad") then
         match=1
         read(argc,*,err=902) prad
      endif
C     if the command is rad we parse it
      if(command.eq."rad") then
         match=1
         read(argc,*,err=902) rad
      endif
C     if the command is skyrad we parse it
      if(command.eq."skyrad") then
         match=1
         read(argc,*,err=902) skyrad
      endif
C     if the command is fitsky we parse it
      if(command.eq."fitsky") then
         match=1
         read(argc,*,err=902) fitsky
         if((fitsky.lt.0).or.(fitsky.gt.1)) then
            write(6,*) "fitsky value is invalid"
            write(6,*) "Setting to no fit (0)"
            fitsky=0
         endif
      endif
C     if command is "sint" we parse it.
      if(command.eq."sint") then
         match=1
         read(argc,*,err=902) sint
         if((sint.lt.0).or.(sint.gt.1)) then
            write(6,*) "sint value is invalid"
            write(6,*) "Setting to fit wings (1)"
            sint=0
         endif
      endif
C     if command is "basename" we parse it
      if(command.eq."basename") then
         match=1
         read(argc,*,err=902) basename
      endif
C     if command is "gain" we parse it
      if(command.eq."gain") then
         match=1
         read(argc,*,err=902) gain
      endif
C     if command is "rnoise" we parse it.
      if(command.eq."rnoise") then
         match=1
         read(argc,*,err=902) rnoise
      endif
C     if command is "listname" we parse it
      if(command.eq."listname") then
         match=1
         read(argc,*,err=902) listname
      endif
C     if command is "flatfile" we parse it
      if(command.eq."flatfile") then
         match=1
         read(argc,*,err=902) flatfile
      endif
C     if command is "datasty" we parse it
      if(command.eq."datasty") then
         match=1
         read(argc,*,err=902) datasty
         if((datasty.lt.0).or.(datasty.gt.1)) then
            write(6,*) "datasty value is invalid"
            write(6,*) "Setting to MEF format (0)"
            datasty=0
         endif
      endif
      
C     if command is "coott" we parse it
      if(command.eq."coott") then
         match=1
         read(argc,*,err=902) coott
         if((coott.lt.0).or.(coott.gt.1)) then
            write(6,*) "coott value is invalid"
            write(6,*) "Setting Co-ordinate tracker to simple (1)"
            coott=1
         endif
      endif

      if(match.eq.0) write(6,*) "No match for command:",command
    
C     read in next line
      goto 10
    
 500  format(a80)
      
 20   close(16)
C     done reading file, now we can exit

C     Output filename
c      oname="fab.dat"
      ln=len_trim(basename)
      oname=basename(1:ln)//"_fab.dat" 


      goto 999

C     if the option file is not found issue a warning and continue
 901  write(6,*) "allmost.pars not found, using defaults"
      goto 999
C     if we encounter an error reading in the parameter file
 902  write(6,501) "Error on line:",nline,"in parameter file"
      write(6,500) line
      goto 999
 501  format(A14,1X,I3,1X,A17)
 999  return 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine skyfit(fitsdata,naxes,datamin,datamax,stn,stx,sty,
     .     pcoo,stfrm,skymap,ans,nfitpsf,ifun)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Not currently working.  Attempts to fit a plane to sky.
      implicit none
      integer naxes(2),xmax,ymax,stmax,stn,stfrm,i,j,npt,ma,nfit,mp,np,
     .     ifun,nfitpsf
      parameter(xmax=600,ymax=600,stmax=20,nfit=5)
      integer stx(stmax,2),sty(stmax,2),k,niter
      real fitsdata(xmax,ymax),pcoo(stmax,2),datamin,datamax,
     .     skymap(xmax,ymax),x(xmax*ymax),y(xmax*ymax),z(xmax*ymax),
     .     a(nfit),u(xmax*ymax,nfit),v(nfit,nfit),w(nfit),
     .     sig(xmax*ymax),chisq,sky,mode,modecal,sdev,stdev,zsky,dis,
     .     ans(nfitpsf),gpsf

      

C     Collect data in STFRM
      npt=0
      do 10 i=1,naxes(1)
         do 11 j=1,naxes(2)
            if((i.ge.stx(stfrm,1)).and.(i.le.stx(stfrm,2)).and.
     .           (j.ge.sty(stfrm,1)).and.(j.le.sty(stfrm,2))) then
               if((fitsdata(i,j).ge.datamin).and.
     .              (fitsdata(i,j).le.datamax)) then
                  npt=npt+1
                  z(npt)=fitsdata(i,j)
               endif
            endif
 11      continue
 10   continue

C     Calculate mode and standard deviation.
      mode=modecal(npt,z,datamin,datamax)
      sdev=stdev(npt,z,mode)

c      write(6,*) "datarange", mode-sdev, mode+sdev

C     Put data into format for linear-fit routine
      npt=0
      do 30 i=1,naxes(1)
         do 31 j=1,naxes(2)
            if((i.ge.stx(stfrm,1)).and.(i.le.stx(stfrm,2)).and.
     .           (j.ge.sty(stfrm,1)).and.(j.le.sty(stfrm,2))) then
               dis=(pcoo(stfrm,1)-real(i))**2.0+
     .              (pcoo(stfrm,2)-real(j))**2.0
               dis=dis**0.5
               if(dis.gt.6.0) then
                  npt=npt+1
                  x(npt)=real(i)
                  y(npt)=real(j)
                  z(npt)=fitsdata(i,j)!-
C     .                 gpsf(x(npt),y(npt),ans,nfitpsf,ifun)
                  sig(npt)=sqrt(abs(z(npt)))+1.1
               endif
            endif
 31      continue
 30   continue
      
      ma=nfit
      mp=xmax*ymax
      np=nfit
c      write(6,*) "npt,ma",npt,ma

      niter=5
      zsky=0
      do 52 j=1,niter
         
         call svdfit(x,y,z,sig,npt,a,ma,u,v,w,mp,np,chisq)      
         
         do 50 i=1,npt
            sky=a(1)
            do 51 k=2,ma,2
               sky=sky+a(k)*real(x(i))**(k/2)+
     .              a(k+1)*real(y(i))**(k/2)
 51         continue
            if(j.eq.niter) zsky=zsky+sky
            sig(i)=abs(z(i)-sky)
            if(sig(i).lt.1.0)sig(i)=1.0
 50      continue
         
 52   continue
      zsky=zsky/npt


c      write(6,*) (a(i),i=1,ma)

      do 20 i=1,naxes(1)
         do 21 j=1,naxes(2)
            if((i.ge.stx(stfrm,1)).and.(i.le.stx(stfrm,2)).and.
     .           (j.ge.sty(stfrm,1)).and.(j.le.sty(stfrm,2))) then
               sky=a(1)
               do 22 k=2,ma,2 
                  sky=sky+a(k)*real(i)**(k/2)+
     .                 a(k+1)*real(j)**(k/2)
 22            continue
               skymap(i,j)=sky
               fitsdata(i,j)=fitsdata(i,j)-sky+zsky
            endif
 21      continue
 20   continue
         
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine apphot(mag,merr,stfrm,rad,skyrad,ans,fitsdata,naxes,
     .     datamin,datamax,ifun,nfitpsf,gain,rnoise,sky,stx,sty,xcoo,
     .     ycoo,rflux,rstx,rsty,itime,nfl,ftotal,fstar,nstack)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nsky,nx,ny,naxes(2),stmax,stfrm,ifun,xmax,ymax,nfitm,i,j,
     .     k,nfitpsf,nstack
      parameter(stmax=20,nfitm=200,xmax=600,ymax=600)
      integer stx(stmax,2),sty(stmax,2),rstx(stmax,2),rsty(stmax,2)
      real flux,sky,nfl,rnfl,rflux,prad,srad,xcoo,ycoo,ans(nfitm),
     .     crit,datamin,datamax,model,dis,gpsf,skypix(xmax,ymax),
     .     fitsdata(xmax,ymax),skyw(xmax,ymax),rad,skyrad,findsky,
     .     fluxw(xmax,ymax),skys(xmax*ymax),ddata(xmax,ymax),stdsky,
     .     gain,merr,rnoise,mag,z1,z2,itime,dgain,ftotal,fstar,
     .     pc(4),qp(4),qcoo(4,2),g,s,h,a,b,pdis,pcrit
     
      h=0.8 !hessian in quadpx algorithm

      flux=0.
      nsky=0
      nfl=0.0
      rnfl=0.0
      nx=0
      ny=0
      rflux=0.

      do 10 i=1,naxes(1)
         do 11 j=1,naxes(2)
            if((i.ge.stx(stfrm,1)).and.(i.le.stx(stfrm,2)).and.
     .           (j.ge.sty(stfrm,1)).and.(j.le.sty(stfrm,2))) then
               if((fitsdata(i,j).ge.datamin).and.
     .              (fitsdata(i,j).le.datamax)) then
                  prad=rad!ypsf(real(i),real(j),ans,nfitpsf,rad)
                  srad=skyrad!(prad/rad)*skyrad
                  dis= ((real(i-stx(stfrm,1)+1)-xcoo)**2.0+
     .                 (real(j-sty(stfrm,1)+1)-ycoo)**2.0)**0.5
                  crit=prad-dis
                  model=gpsf(real(i),real(j),ans,nfitpsf,ifun)
                  if(crit.gt.0.5) then
                     flux=flux+fitsdata(i,j)
                     nfl=nfl+1.0
                     fluxw(i,j)=1.0
                     skypix(i,j)=fitsdata(i,j)-model
                     skyw(i,j)=1.0/(0.9*(dis/skyrad)+0.1)
                  elseif((crit.ge.-0.5).and.(crit.le.0.5))then
                     fluxw(i,j)=1.0*(crit+0.5)
C                    old partial pixel scheme
c                     flux=flux+fitsdata(i,j)*(crit+0.5)
c                     nfl=nfl+1.0*(crit+0.5)
c                     fluxw(i,j)=1.0*(crit+0.5)
c                     skypix(i,j)=fitsdata(i,j)-model
c                     skyw(i,j)=1.0/(0.9*(dis/skyrad)+0.1)
                     if(j+1.le.sty(stfrm,2))then
                        pc(1)=fitsdata(i,j+1)
                     else
                        pc(1)=fitsdata(i,j)
                     endif
                     if(i+1.le.stx(stfrm,2))then
                        pc(2)=fitsdata(i+1,j)
                     else
                        pc(2)=fitsdata(i,j)
                     endif
                     if(j-1.ge.sty(stfrm,1))then
                        pc(3)=fitsdata(i,j-1)
                     else
                        pc(3)=fitsdata(i,j)
                     endif
                     if(i-1.ge.stx(stfrm,1))then
                        pc(4)=fitsdata(i-1,j)
                     else
                        pc(4)=fitsdata(i,j)
                     endif
                     g=min(max(h,0.0),1.0)
                     s=pc(1)+pc(2)+pc(3)+pc(4)
                     if(s.eq.0.0) g=0.0
                     a=(1.0-g)*(fitsdata(i,j)/4.0)
                     b=(g/s)*(fitsdata(i,j)/2.0)
                     qp(1)=(pc(1)+pc(4))*b+a
                     qcoo(1,1)=real(i-0.25)
                     qcoo(1,2)=real(j+0.25)
                     qp(2)=(pc(2)+pc(1))*b+a
                     qcoo(2,1)=real(i+0.25)
                     qcoo(2,2)=real(j+0.25)
                     qp(3)=(pc(4)+pc(3))*b+a
                     qcoo(3,1)=real(i-0.25)
                     qcoo(3,2)=real(j-0.25)
                     qp(4)=(pc(3)+pc(2))*b+a
                     qcoo(4,1)=real(i+0.25)
                     qcoo(4,2)=real(j-0.25)
                     do 12 k=1,4
                        pdis=((qcoo(k,1)-real(stx(stfrm,1))+1.0-xcoo)
     .                     **2.0+((qcoo(k,2)-real(sty(stfrm,1))+1.0-
     .                     ycoo)**2.0))**0.5
                        crit=prad-pdis
                        if((crit.ge.-0.5).or.(crit.le.0.5))then
                           flux=flux+qp(k)*(crit+0.5)
                           nfl=nfl+1.0*(crit+0.5)/4.0
                        elseif(crit.gt.0.5)then
                           flux=flux+qp(k)
                           nfl=nfl+1.0/4.0
                        endif
 12                  continue                     
                     skypix(i,j)=fitsdata(i,j)-model
                     skyw(i,j)=1.0/(0.9*(dis/skyrad)+0.1)
c                     write(6,*) flux,fitsdata(i,j),nfl,(crit+0.5)
                  elseif(dis.lt.srad)then
                     skypix(i,j)=fitsdata(i,j)-model
                     skyw(i,j)=1.0/(0.9*(dis/skyrad)+0.1)
                  endif
                  if(dis.gt.srad)then
                     skypix(i,j)=fitsdata(i,j)
                     skyw(i,j)=1.0
                     nsky=nsky+1
                     skys(nsky)=fitsdata(i,j)
                  endif
C     NEED TO MAKE THIS A PARAMETER
                  if(dis.lt.6.0) then
                     rflux=rflux+fitsdata(i,j)-model
                     rnfl=rnfl+1.0
                  endif
               endif      
            endif
 11      continue
 10   continue

C      sky=findsky(nsky,skys,stdsky)
C     This should be zero as gpsf does not tack on Ib

C     Xi Gem (110-130)
C     Procyon (0-5)
cffff FLATCODE
cf      if((sky.gt.1000.0).and.(sky.lt.1500.).and.(ans(5).lt.2.0).and.
cf     .     (ans(6).lt.2.0).and.(rad.eq.10.0)) then
cf         tt=tt+1
cf         write(6,*) tt,(stx(stfrm,i),sty(stfrm,i),i=1,2)
cf         call makeflat(stfrm,stmax,stx,sty,skypix)
cf      endif
cf      if((sky.gt.0.0).and.(stfrm.le.2)) then
cf         do 30 i=stx(stfrm,1),stx(stfrm,2)
cf            do 31 j=sty(stfrm,1),sty(stfrm,2)
cf               write(14,501) i-stx(stfrm,1)+rstx(stfrm,1),
cf     .              j-sty(stfrm,1)+rsty(stfrm,1),sky,fitsdata(i,j)!skypix(i,j)
cf 31         continue
cf 30      continue
cf      endif
 501  format(2(I4,1X),2(F8.2,1X))
cffff FLATCODE

C     Save total flux
      ftotal=flux
C     calculate flux from star
      flux=flux-sky*nfl
      fstar=flux
C     Rflux is for PSF corrections.
      rflux=rflux-sky*rnfl

      mag=25.0-2.5*log10(flux/itime*gain)
      merr=1.09/(flux*gain)*sqrt(flux*gain+
     .     nfl*((nstack*rnoise*gain)**2+sky*gain))
c      write(6,*) "ap:",flux,nfl,sky,mag,merr
c      read(5,*)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine apphot2(mag,merr,stfrm,rad,skyrad,ans,fitsdata,naxes,
     .     datamin,datamax,ifun,nfitpsf,gain,rnoise,sky,stx,sty,xcoo,
     .     ycoo,rflux,rstx,rsty,itime,nfl,ftotal,fstar,nstack)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nsky,nx,ny,naxes(2),stmax,stfrm,ifun,xmax,ymax,nfitm,i,j,
     .     k,nfitpsf,nstack,ii
      parameter(stmax=20,nfitm=200,xmax=600,ymax=600)
      integer stx(stmax,2),sty(stmax,2),rstx(stmax,2),rsty(stmax,2)
      real flux,sky,nfl,rnfl,rflux,prad,srad,xcoo,ycoo,ans(nfitm),
     .     crit,datamin,datamax,model,dis,gpsf,skypix(xmax,ymax),
     .     fitsdata(xmax,ymax),skyw(xmax,ymax),rad,skyrad,findsky,
     .     fluxw(xmax,ymax),skys(xmax*ymax),ddata(xmax,ymax),stdsky,
     .     gain,merr,rnoise,mag,z1,z2,itime,dgain,ftotal,fstar,fill,
     .     cirint,cirint2,ri,rj,coo(2,2),coot(2,2)

      flux=0.
      nsky=0
      nfl=0.0
      rnfl=0.0
      nx=0
      ny=0
      rflux=0.

      do 10 i=1,naxes(1)
         ri=real(i-stx(stfrm,1)+1)
         do 11 j=1,naxes(2)
            if((i.ge.stx(stfrm,1)).and.(i.le.stx(stfrm,2)).and.
     .           (j.ge.sty(stfrm,1)).and.(j.le.sty(stfrm,2))) then
               if((fitsdata(i,j).ge.datamin).and.
     .              (fitsdata(i,j).le.datamax)) then
                  prad=rad!ypsf(real(i),real(j),ans,nfitpsf,rad)
                  srad=skyrad!(prad/rad)*skyrad
                  dis= ((real(i-stx(stfrm,1)+1)-xcoo)**2.0+
     .                 (real(j-sty(stfrm,1)+1)-ycoo)**2.0)**0.5
                  crit=prad-dis
                  model=gpsf(real(i),real(j),ans,nfitpsf,ifun)
                  ii=0     !counts number of intersections for p.pixels
                  rj=real(j-sty(stfrm,1)+1) !turn counter into real
                  call getintc(ii,prad,ri,rj,xcoo,ycoo,coo)!find p.pixs
                  if((ii.eq.0).and.(crit.ge.0.5))then!inside ap
                     flux=flux+fitsdata(i,j)
                     nfl=nfl+1.0
                     fluxw(i,j)=1.0
                     skypix(i,j)=fitsdata(i,j)-model
                     skyw(i,j)=1.0/(0.9*(dis/skyrad)+0.1)
                  elseif((ii.eq.0).and.(crit.lt.-0.5))then !outside ap
                     fill=0.
                  elseif(ii.eq.2)then    !p. pixels
                     if((coo(1,2).ge.ycoo).and.(coo(2,2).ge.ycoo))then
                        fill=cirint(prad,xcoo,ycoo,coo,ri,rj)!top half
                     elseif((coo(1,2).lt.ycoo).and.
     .                    (coo(2,2).lt.ycoo))then
                        fill=cirint2(prad,xcoo,ycoo,coo,ri,rj)!bottom h
                     else
                        coot(1,1)=coo(1,2)  !in middle, rotate co-ords
                        coot(1,2)=coo(1,1)
                        coot(2,1)=coo(2,2)
                        coot(2,2)=coo(2,1)
                        if((coo(1,1).gt.xcoo).and.
     .                       (coo(2,1).gt.xcoo))then
                           fill=cirint(prad,ycoo,xcoo,coot,rj,ri)!top h
                        elseif((coo(1,1).lt.xcoo).and.
     .                       (coo(2,1).lt.xcoo))then
                           fill=cirint2(prad,ycoo,xcoo,coot,rj,ri)!b.h.
                        else
                           write(6,*) "err2:",ii
                           fill=crit+0.5  !mark errors here
                        endif
                     endif
                     flux=flux+fitsdata(i,j)*fill
                     nfl=nfl+1.0*fill
                     fluxw(i,j)=1.0*fill
                     skypix(i,j)=fitsdata(i,j)-model
                     skyw(i,j)=1.0/(0.9*(dis/skyrad)+0.1)               
                  else
                     write(6,*) "err:",ii
                     fill=crit+0.5  !mark errors here
                     if((fill.ge.0).and.(fill.le.1.0))then
                        flux=flux+fitsdata(i,j)*fill
                        nfl=nfl+1.0*fill
                        fluxw(i,j)=1.0*fill
                        skypix(i,j)=fitsdata(i,j)-model
                        skyw(i,j)=1.0/(0.9*(dis/skyrad)+0.1)
                     endif
                  endif
                  if(dis.gt.srad)then
                     skypix(i,j)=fitsdata(i,j)
                     skyw(i,j)=1.0
                     nsky=nsky+1
                     skys(nsky)=fitsdata(i,j)
                  endif
C     NEED TO MAKE THIS A PARAMETER
                  if(dis.lt.2.0) then
                     rflux=rflux+fitsdata(i,j)-model
                     rnfl=rnfl+1.0
                  endif
                  skypix(i,j)=fitsdata(i,j)
               endif      
            endif
 11      continue
 10   continue

C      sky=findsky(nsky,skys,stdsky)
C     This should be zero as gpsf does not tack on Ib

C     Xi Gem (110-130)
C     Procyon (0-5)
cffff FLATCODE
cf      if((sky.gt.1000.0).and.(sky.lt.1500.).and.(ans(5).lt.2.0).and.
cf     .     (ans(6).lt.2.0).and.(rad.eq.10.0)) then
cf         tt=tt+1
cf         write(6,*) tt,(stx(stfrm,i),sty(stfrm,i),i=1,2)
cf         call makeflat(stfrm,stmax,stx,sty,skypix)
cf      endif
cf      if((sky.gt.0.0).and.(stfrm.le.2)) then
cf         do 30 i=stx(stfrm,1),stx(stfrm,2)
cf            do 31 j=sty(stfrm,1),sty(stfrm,2)
cf               write(14,501) i-stx(stfrm,1)+rstx(stfrm,1),
cf     .              j-sty(stfrm,1)+rsty(stfrm,1),
cf     .              sky*9.0/real(nstack),fitsdata(i,j)*9.0/real(nstack)
cf 31         continue
cf 30      continue
cf      endif
 501  format(2(I4,1X),2(F8.2,1X))
cffff FLATCODE

C     Save total flux
      ftotal=flux
C     calculate flux from star
      flux=flux-sky*nfl
      fstar=flux
C     Rflux is for PSF corrections.
      rflux=rflux-sky*rnfl

      mag=25.0-2.5*log10(flux/itime*gain)
      merr=1.09/(flux*gain)*sqrt(flux*gain+
     .     nfl*((nstack*rnoise*gain)**2+sky*gain))
c      write(6,*) "ap:",flux,nfl,sky,mag,merr
c      read(5,*)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getintc(ii,ap,ri,rj,xc,yc,coo)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ii,k
      real ap,ri,rj,xc,yc,coo(2,2),discxp,discyp,discxm,discym,x(2),y(2)
     
      discxp=ap*ap-(ri+0.5-xc)**2
      if(discxp.ge.0)then 
         y(1)=sqrt(discxp)+yc
         y(2)=-sqrt(discxp)+yc
         do 13 k=1,2
            if((ii.lt.2).and.(abs(y(k)-rj).le.0.5))then
               ii=ii+1
               coo(ii,1)=ri+0.5
               coo(ii,2)=y(k)
            endif
 13      continue
      endif
      discyp=ap*ap-(rj+0.5-yc)**2
      if(discyp.ge.0)then
         x(1)=sqrt(discyp)+xc
         x(2)=-sqrt(discyp)+xc
         do 14 k=1,2
            if((ii.lt.2).and.(abs(x(k)-ri).le.0.5))then
               ii=ii+1
               coo(ii,1)=x(k)
               coo(ii,2)=rj+0.5
            endif
 14      continue
      endif
      discxm=ap*ap-(ri-0.5-xc)**2
      if(discxm.ge.0)then
         y(1)=sqrt(discxm)+yc
         y(2)=-sqrt(discxm)+yc
         do 15 k=1,2
            if((ii.lt.2).and.(abs(y(k)-rj).le.0.5))then
               ii=ii+1
               coo(ii,1)=ri-0.5
               coo(ii,2)=y(k)
            endif
 15      continue
      endif 
      discym=ap*ap-(rj-0.5-yc)**2
      if(discym.ge.0)then
         x(1)=sqrt(discym)+xc
         x(2)=-sqrt(discym)+xc
         do 16 k=1,2
            if((ii.lt.2).and.(abs(x(k)-ri).le.0.5))then
               ii=ii+1
               coo(ii,1)=x(k)
               coo(ii,2)=rj-0.5
            endif
 16      continue
      endif     
    
      return
      end  
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function cirint2(ap,xo,yo,coo,ri,rj)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real xo,yo,x1,x2,x,int1,int2,dx,ap,ri,rj,coo(2,2),y1,y2,corr
      
      x1=coo(1,1)
      x2=coo(2,1)
      y1=coo(1,2)
      y2=coo(2,2)
      
      x=max(x1,x2)
      int1=(Sqrt(ap**2 - (x - xo)**2)*(-x + xo) + 2*x*yo + 
     -    ap**2*ATan((-x + xo)/Sqrt(ap**2 - (x - xo)**2)))/2.
      x=min(x1,x2)
      int2=(Sqrt(ap**2 - (x - xo)**2)*(-x + xo) + 2*x*yo + 
     -    ap**2*ATan((-x + xo)/Sqrt(ap**2 - (x - xo)**2)))/2.
     
      dx=abs(x2-x1)
      cirint2=dx-(int1-int2-(rj-0.5)*dx)
      
      corr=0.
      if(dx.lt.1.0)then
         if(min(x1,x2).gt.xo)then
            corr=corr+min(x1,x2)-(ri-0.5)
         elseif(max(x1,x2).lt.xo)then
            corr=corr+(ri+0.5)-max(x1,x2)
         endif
      endif
      cirint2=cirint2+corr
      
      return
      end
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function cirint(ap,xo,yo,coo,ri,rj)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real xo,yo,x1,x2,x,int1,int2,dx,ap,ri,rj,y1,y2,coo(2,2),corr
      
      x1=coo(1,1)
      x2=coo(2,1)
      y1=coo(1,2)
      y2=coo(2,2)
      
      x=max(x1,x2)
      int1=(Sqrt(ap**2 - (x - xo)**2)*(x - xo))/2. + x*yo - 
     -  (ap**2*ATan((-x + xo)/Sqrt(ap**2 - (x - xo)**2)))/2.
c      write(6,*) ap**2 - (x - xo)**2
      x=min(x1,x2)
      int2=(Sqrt(ap**2 - (x - xo)**2)*(x - xo))/2. + x*yo - 
     -  (ap**2*ATan((-x + xo)/Sqrt(ap**2 - (x - xo)**2)))/2.

      dx=abs(x2-x1)
      cirint=int1-int2-(rj-0.5)*dx

      corr=0.
      if(dx.lt.1.0)then
         if(min(x1,x2).gt.xo)then
            corr=corr+min(x1,x2)-(ri-0.5)
         elseif(max(x1,x2).lt.xo)then
            corr=corr+(ri+0.5)-max(x1,x2)
         endif
      endif
      cirint=cirint+corr
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine makeflat(stfrm,stmax,stx,sty,skypix)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Uses star-subtracted images to estimate a flatfield.
      integer stfrm,stmax,stx(stmax,2),sty(stmax,2),xmax,ymax,i,j,fmax,
     .     ff,fm,tt,k,ttold
      parameter(xmax=600,ymax=600,fmax=20)
      real skypix(xmax,ymax),flatmedian(xmax,ymax),pts(fmax),median,
     .     flats(fmax,xmax,ymax),flat(xmax,ymax),flatm(fmax,xmax,ymax),
     .     findmed,stdev
      save flatmedian,flats,ff,fm,ttold

c      if(fm.gt.fmax) return

      if(tt.eq.1) then
         ttold=0
         do 5 i=1,xmax
            do 6 j=1,ymax
               flat(i,j)=0.0
 6          continue
 5       continue
         ff=1
         fm=1
      else
         if(tt.eq.ttold) return
         if(tt.gt.ttold) ttold=tt
         ff=ff+1
      endif
      
c      write(6,*) ff,fm

      if(fm.gt.fmax) then
         fm=fm-1
         open(unit=11,file="/home/rowe/MOST/flat.dat")
         do 40 i=stx(stfrm,1),stx(stfrm,2)
            do 41 j=sty(stfrm,1),sty(stfrm,2)
               flat(i,j)=0.0
               do 42 k=1,fm
                  flat(i,j)=flat(i,j)+flatm(k,i,j)
                  pts(k)=flatm(k,i,j)
 42            continue
c              if(i.eq.16) write(6,500)(int(pts(k)),k=1,fm)
               flat(i,j)=flat(i,j)/real(fm)
               write(11,*)i,j,flat(i,j),stdev(fm,pts,flat(i,j))
 41         continue
 40      continue
         close(11)
         stop
      endif

      if(ff.gt.fmax) then
         ff=ff-1
         do 30 i=stx(stfrm,1),stx(stfrm,2)
            do 31 j=sty(stfrm,1),sty(stfrm,2)
               do 32 k=1,ff
                  pts(k)=flats(k,i,j)
 32            continue
c               if(i.eq.16) write(6,500)(int(pts(k)),k=1,ff)
               median=findmed(ff,pts)
               flatmedian(i,j)=median
 31         continue
 30      continue
         do 20 i=stx(stfrm,1),stx(stfrm,2)
            do 21 j=sty(stfrm,1),sty(stfrm,2)
               flatm(fm,i,j)=flatmedian(i,j)
 21         continue
 20      continue
         fm=fm+1
         ff=1
      endif

      if(ff.le.fmax) then
         do 10 i=stx(stfrm,1),stx(stfrm,2)
            do 11 j=sty(stfrm,1),sty(stfrm,2)
C               if((i.eq.16).and.(j.eq.1)) then
C                  write(6,500) (int(skypix(15,k)),k=1,16)
C                  write(6,500) (int(skypix(i,k)),k=1,16)
Cc                  read(5,*)
C               endif
 500           format(20(I5,1X))
               flats(ff,i,j)=skypix(i,j)
 11         continue
 10      continue
      endif
            
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function exptime(header,nrecmax,nkeys,nstack)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrecmax,nkeys,i,nstack
      real temp,temp2
      character*80 header(nrecmax),record,headint,headint2

C     2nd generation data format
      headint="LEN_EXP"
C     3rd generation data format
      headint2="TOT_EXP"
      headint2="GS_T0000"


      temp=0.
      temp2=0.

      do 10 i=1,nkeys
         record=header(i)
         if(record(1:8).eq.headint(1:8))then
            read(record(10:80),*) temp
         endif
         if(record(1:8).eq.headint2(1:8))then
            read(record(10:80),*) temp2
         endif
 10   continue

      if(temp2.eq.0.0) then
         exptime=temp/1000.0
         nstack=1
      else
         exptime=temp2/1000.0
         nstack=temp2/temp
      endif
c      write(6,*) temp,temp2

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function magfield(header,nrecmax,nkeys)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrecmax,nkeys,i
      real temp
      character*80 header(nrecmax),record,headint

      headint="MAG_FLD"

      temp=20000.0

      do 10 i=1,nkeys
         record=header(i)
         if(record(1:8).eq.headint(1:8))then
            read(record(10:80),*) temp
         endif
 10   continue

      magfield=temp

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function boardtemp(header,nrecmax,nkeys)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrecmax,nkeys,i
      real temp
      character*80 header(nrecmax),record,headint

      headint="T_BOARD"

      temp=0

      do 10 i=1,nkeys
         record=header(i)
         if(record(1:8).eq.headint(1:8))then
            read(record(10:80),*) temp
         endif
 10   continue

      boardtemp=temp

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getgsmag(mag,merr,fstar,nfl,itime,header,nrecmax,nkeys,
     .  igns,gain,rnoise,sky,nstack)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrecmax,nkeys,gns,nstack,i,igns
      real flux,exptime,mag,merr,gain,sky,rnoise,fstar,itime,nfl
      character*80 header(nrecmax),headint,record,headint2
      
      gns=igns-1 !FITS headers start at entry 0
      
      if(gns.lt.10)then
        write(headint,500) "GS_I000",gns !header to look for
        write(headint2,500) "GS_T000",gns
 500    format(A7,I1)
      elseif(gns.lt.100)then
        write(headint,501) "GS_I00",gns
        write(headint2,501) "GS_T00",gns
 501    format(A6,I2)
      elseif(gns.lt.1000)then
        write(headint,502) "GS_I0",gns
        write(headint2,502) "GS_T0",gns
 502    format(A6,I3)
      else
        write(headint,503) "GS_I",gns
        write(headint2,503) "GS_T",gns
 503    format(A6,I4)
      endif 

      flux=0
      exptime=1
      do 10 i=1,nkeys
        record=header(i)
        if(record(1:8).eq.headint(1:8))then
            read(record(10:80),*) flux
        endif
        if(record(1:8).eq.headint2(1:8))then
            read(record(10:80),*) exptime
            exptime=exptime/1000.0 !convert from ms to sec
        endif
 10   continue
 
      nfl=9.0 !assume there are 9 pixels in aperture
      fstar=flux
      itime=exptime
 
      mag=25.0-2.5*log10(flux*gain/exptime)
      merr=1.09/(flux*gain)*sqrt(flux*gain+
     .     nfl*((nstack*rnoise*gain)**2+sky*gain))
     
c      write(0,*) gns,mag,merr
      

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function gsfields(header,nrecmax,nkeys)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrecmax,nkeys,gns,i,nf
      character*80 header(nrecmax),headint,record
      logical loop
      
      loop=.true.
      gns=0 !number of guidestars found
      do while(loop)
        nf=0 !flag to determine whether a header was found
                
        if(gns.lt.10)then
            write(headint,500) "GS_I000",gns !header to look for
 500        format(A7,I1)
        elseif(gns.lt.100)then
            write(headint,501) "GS_I00",gns
 501        format(A6,I2)
        elseif(gns.lt.1000)then
            write(headint,502) "GS_I0",gns
 502        format(A6,I3)
        else
            write(headint,503) "GS_I",gns
 503        format(A6,I4)
        endif

        do 10 i=1,nkeys
            record=header(i)
            if(record(1:8).eq.headint(1:8))then
                nf=1
            endif
 10     continue

        if(nf.eq.0)then  
            loop=.false.  !if we didn't find the header, then exit
        else
            gns=gns+1     !if we found the header, search for next
        endif
      enddo

      gsfields=gns
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function readdate(header,nrecmax,nkeys,datasty)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrecmax,nkeys,i,datasty,timecorr
      real date
      double precision ddate,hdcor
      character*80 header(nrecmax),record,headint

      if(datasty.eq.1) then
         headint="DATE-CAL"
      else
         timecorr=2
         headint="TIMECORR"
         do 20 i=1,nkeys
            record=header(i)
            if(record(1:8).eq.headint(1:8))then
               read(record(10:80),*) timecorr
            endif
 20      continue 
         if(timecorr.eq.0) then
            headint="JD-OBS"
         elseif(timecorr.eq.1) then
            headint="JD-CAL"
         elseif(timecorr.ge.2) then
            headint="JD-ACC"
         endif
      endif

      ddate=0.
      date=0.
      hdcor=0.

      do 10 i=1,nkeys
         record=header(i)
         if(record(1:8).eq.headint(1:8))then
            read(record(10:80),*) ddate
         endif
         if(record(1:8).eq."HJD-CORR")then
            read(record(10:80),*) hdcor
         endif
 10   continue
      write(17,*) ddate,hdcor,ddate+hdcor/8.64d4
      
      ddate=ddate+hdcor/8.64d4
c      if(hdcor.eq.0.) ddate=0.
      readdate=ddate
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine magfiles(ns,munit,basename)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ns,nsmax,i
      parameter(nsmax=100)
      integer munit(nsmax),uu,status,ln
      character*80 basename,filename
      character*2 cnum

c      ln=lnblnk(basename)
      ln=len_trim(basename)

      do 10 i=1,ns
         status=0
         call ftgiou(uu,status)
         if(i.lt.10) then
            write(cnum,500) "0",i
         elseif(i.lt.100) then
            write(cnum,501) i
         endif
         filename=basename(1:ln)//"_"//cnum//".dat"
         if(status.ne.0) pause "problem in magfiles"
         open(unit=uu,file=filename)
         munit(i)=uu
 10   continue

 500  format(A1,I1)
 501  format(I2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gsfiles(gns,gsunit,basename)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer gns,nsmax,i
      parameter(nsmax=100)
      integer gsunit(nsmax),uu,status,ln
      character*80 basename,filename
      character*2 cnum

c      ln=lnblnk(basename)
      ln=len_trim(basename)

      do 10 i=1,gns
         status=0
         call ftgiou(uu,status)
         if(i.lt.10) then
            write(cnum,500) "0",i
         elseif(i.lt.100) then
            write(cnum,501) i
         endif
         filename=basename(1:ln)//"_gs_"//cnum//".dat"
         if(status.ne.0) pause "problem in gsfiles"
         open(unit=uu,file=filename)
         gsunit(i)=uu
 10   continue

 500  format(A1,I1)
 501  format(I2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cmagfiles(ns,munit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ns,nsmax,i
      parameter(nsmax=100)
      integer munit(nsmax)

      do 10 i=1,ns
         close(munit(i))
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cgsfiles(gns,gsunit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer gns,nsmax,i
      parameter(nsmax=100)
      integer gsunit(nsmax)

      do 10 i=1,gns
         close(gsunit(i))
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine coooffset(stn,pcoo,pcool,peaks,bpeak,hasstar,coott)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer stn,stmax,nsmax,i,j,stc,coott
      parameter(stmax=20,nsmax=100)
      integer bpeak(stmax),hasstar(stmax)
      real peaks(stmax,2),pcoo(stmax,nsmax,2),ox(stmax),oy(stmax),
     .     pcool(stmax,nsmax,2),sigcut,meanx,meany,avg,std,stdev

      sigcut=1.0

      meanx=0.
      meany=0.
      stc=0
      do 10 i=1,stn
         if(hasstar(i).gt.0) then
c            if((i.ne.5).and.(i.ne.4))then
               stc=stc+1
               ox(stc)=peaks(i,1)-pcool(i,bpeak(i),1)
               oy(stc)=peaks(i,2)-pcool(i,bpeak(i),2)
c            endif
CNext two lines are proc
            if(coott.eq.1) then
               pcool(i,bpeak(i),1)=peaks(i,1)
               pcool(i,bpeak(i),2)=peaks(i,2)
            endif
            meanx=meanx+ox(stc)
            meany=meany+oy(stc)
C           write(6,*) "ox:",i,ox(stc),peaks(i,1),pcool(i,bpeak(i),1)
         endif
 10   continue
      meanx=meanx/real(stc)
      meany=meany/real(stc)

      if(stc.gt.2)then
         std=stdev(stc,ox,meanx)
         j=0
         avg=0.
         do 30 i=1,stc
            if(abs(ox(i)-meanx).lt.sigcut*std) then
               j=j+1
               avg=avg+ox(i)
            endif
 30      continue
         meanx=avg/real(j)

         std=stdev(stc,oy,meany)
         j=0
         avg=0.
         do 31 i=1,stc
            if(abs(oy(i)-meany).lt.sigcut*std) then
               j=j+1
               avg=avg+oy(i)
            endif
 31      continue
         meany=avg/real(j)

      endif

c      if((abs(meanx).lt.50.).and.(abs(meany).lt.50.))then
c         continue
c      else
c         meanx=0.
c         meany=0.
c      endif

C     Next two lines are proc
      if(coott.eq.1) then
         meanx=0.
         meany=0.
      endif

      do 20 i=1,stn
         do 21 j=1,hasstar(i)
            pcoo(i,j,1)=pcool(i,j,1)+meanx
            pcoo(i,j,2)=pcool(i,j,2)+meany
c            if(i.eq.2)write(6,*)"coo:",pcoo(i,j,1),pcoo(i,j,2)
 21      continue
 20   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fbpeak(stn,pmag,hasstar,bpeak,pcoo,pcool)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer stn,stmax,nsmax,i,j
      parameter(stmax=20,nsmax=100)
      integer hasstar(stmax),bpeak(stmax)
      real pmag(stmax,nsmax),bmag,pcoo(stmax,nsmax,2),
     .     pcool(stmax,nsmax,2)

      do 10 i=1,stn
         bmag=99.9
         bpeak(i)=1
         do 11 j=1,hasstar(i)
            pcool(i,j,1)=pcoo(i,j,1)
            pcool(i,j,2)=pcoo(i,j,2)
            if(pmag(i,j).lt.bmag) then
               bmag=pmag(i,j)
               bpeak(i)=j
            endif
 11      continue
C         write(6,*) "bpeak",i,bpeak(i)
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findpeaks(fitsdata,naxes,datamin,datamax,
     .     stn,stx,sty,pcoo,hasstar,pval)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     star detection and centroiding algorithim
C     Assumes there is ONE star on sub-raster
C     initial center is guessed from peak value
C     centroiding using "centre of mass" technique
C     since each star moves the same rate, positions are
C     further corrected for using offsets
      implicit none
      integer naxes(2),xmax,ymax,stmax,i,stn,j,k,ns,xc,yc
      parameter(xmax=600,ymax=600,stmax=20,ns=5)
      integer stx(stmax,2),sty(stmax,2),nx,ny,mx,my,hasstar(stmax),iter,
     .     maxiter
      real fitsdata(xmax,ymax),datamin,datamax,sumx,sumy,fsum,
     .     pcoo(stmax,2),mpeak,oldcoo(2),crit,pval(stmax),
     .     rdata(xmax,ymax),xcoo,ycoo

C     Avoid infinite loops
      maxiter=50 !maximum number of iterations for find centroid
      

      nx=naxes(1)  !number of x-pixels
      ny=naxes(2)  !number of y-pixels

      do 10 i=1,stn
         mpeak=-1.0e30 !some low number
         sumx=0.0 !summation in x-direction
         sumy=0.0 !summation in y-direction
         fsum=0.0 !sum over all pixels
         do 20 j=1,naxes(1) !loop over all pixels
            do 30 k=1,naxes(2)
C              Now we check to see that the pixels are in the subraster
               if((j.ge.stx(i,1)).and.(j.le.stx(i,2)).and.
     .              (k.ge.sty(i,1)).and.(k.le.sty(i,2))) then
C              Check that fits values are inside datamin/datamax bounds
                  if((fitsdata(j,k).gt.datamin).and.
     .                 (fitsdata(j,k).lt.datamax))then
                     if(fitsdata(j,k).gt.mpeak) then !find max-pixel
                        mpeak=fitsdata(j,k)  !store peak pixel value
                        mx=j
                        my=k
                     endif
C                    All of this is pretty much useless.
                     sumx=sumx+fitsdata(j,k)*real(j)
                     sumy=sumy+fitsdata(j,k)*real(k)
                     fsum=fsum+fitsdata(j,k)
                  endif
               endif
 30         continue
 20      continue
         pval(i)=mpeak
         pcoo(i,1)=real(mx)!sumx/fsum
         pcoo(i,2)=real(my)!sumy/fsum

         iter=0
 42      oldcoo(1)=pcoo(i,1)
         oldcoo(2)=pcoo(i,2)
         iter=iter+1

         sumx=0
         sumy=0
         fsum=0
         do 40 j=1,ns
            do 41 k=1,ns
               xc=int(pcoo(i,1))-3+j
               yc=int(pcoo(i,2))-3+k
               if((xc.ge.stx(i,1)).and.(xc.le.stx(i,2)).and.
     .              (yc.ge.sty(i,1)).and.(yc.le.sty(i,2))) then
                  if((fitsdata(xc,yc).gt.datamin).and.
     .                 (fitsdata(xc,yc).lt.datamax))then
                     sumx=sumx+fitsdata(xc,yc)*real(xc)
                     sumy=sumy+fitsdata(xc,yc)*real(yc)
                     fsum=fsum+fitsdata(xc,yc)
                  endif
               endif
 41         continue
 40      continue
         
         pcoo(i,1)=sumx/fsum
         pcoo(i,2)=sumy/fsum

         if(iter.gt.maxiter) goto 43
         crit=((pcoo(i,1)-oldcoo(1))**2.0 + 
     .      (pcoo(i,2)-oldcoo(2))**2.0)**0.5
         if(crit.gt.1.0) goto 42
 43      continue
         pcoo(i,1)=pcoo(i,1)-real(stx(i,1))+1.0
         pcoo(i,2)=pcoo(i,2)-real(sty(i,1))+1.0

c         write(6,*) i, "peak =",pcoo(i,1),pcoo(i,2),pval(i)

 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findstar(nx,ny,rdata,ns,xcoo,ycoo)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nx,ny,xmax,ymax,i,j,mx,my,ns,xc,yc,n1,n2,n3,k,nsmax
      parameter(xmax=600,ymax=600,n1=64,n3=1,nsmax=100)
      real rdata(xmax,ymax),xcoo,ycoo,sumx,sumy,fsum,
     .     oldx,oldy,crit,mpeak,data1(n1,n1),data2(n1,n1),fac,dat,
     .     pdata(xmax,ymax),z1,z2,pts(xmax*ymax),std,findmed,stdev,
     .     median,sigcut
      complex spec1(n1/2,n1),spec2(n1/2,n1),speq1(n1),speq2(n1),
     .     zpec1(n1*n1/2),zpec2(n1*n1/2),zpeq1(n1),zpeq2(n1)
      equivalence (data1,spec1,zpec1), (data2,spec2,zpec2),
     .     (speq1,zpeq1), (speq2,zpeq2)

      sigcut=3.0

      open(unit=12,file="/home/rowe/MOST/medpsf.dat",status='old')

      do 5 i=1,n1
         do 6 j=1,n1
            data1(i,j)=0.
            data2(i,j)=0.
 6       continue
 5    continue

 20   read(12,*,end=21) xc,yc,dat
      i=xc-8
      if(i.lt.1) i=i+nx
      j=yc-26
      if(j.lt.1) j=j+ny
      data2(i,j)=dat
      goto 20
 21   continue
      close(12)

      do 10 i=1,nx
         do 11 j=1,ny
            data1(i,j)=rdata(i,j)
            if(data2(i,j).eq.0.)data2(i,j)=116.0
 11      continue
 10   continue

      call rlft3(data1,speq1,N1,N1,N3,1)
      call rlft3(data2,speq2,N1,N1,N3,1)
      fac=2./(N1*N1*N3)

      do 12 j=1,n1*n1/2
         zpec1(j)=fac*zpec1(j)/zpec2(j)
 12   continue
      do 13 j=1,n1
         zpeq1(j)=fac*zpeq1(j)/zpeq2(j)
 13   continue
      call rlft3(data1,speq1,n1,n1,n3,-1)

      z1= 1.0e30
      z2=-1.0e30
      k=0
      do 25 i=1,n1
         do 26 j=1,n1
            k=k+1
            pdata(i,j)=data1(i,j)
            pts(k)=data1(i,j)
            if(pdata(i,j).lt.z1) z1=pdata(i,j)
            if(pdata(i,j).gt.z2) z2=pdata(i,j)
 26      continue
c         write(6,*) (pdata(i,j),j=1,n1)
 25   continue

      call sparewindow(n1,n1,pdata,z1,z2)
      read(5,*)

      if(k.eq.0) then
         write(6,*) "pts is zero!"
         pause
      endif
      median=findmed(k,pts)
      std=stdev(k,pts,median)

      do 27 i=1,n1
         do 28 j=1,n1
            if(data1(i,j).gt.median+std*sigcut) then
               write(6,*) "star ",i,j,data1(i,j)
               ns=ns+1
               if(ns.gt.nsmax) pause "nsmax is not big enough!"
               xcoo=real(i)
               ycoo=real(j)
            endif
 28      continue
 27   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine psfint(pmag,rad,nfit,ans,sky,rflux,gain,rnoise,ifun,
     .     hasstar,stfrm,itime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nfitm,ngrid,i,j,ifun,stmax,nsmax,k,id
      parameter(nfitm=30,ngrid=100,stmax=20,nsmax=100)
      integer hasstar(stmax),stfrm,ifun1
      real gmag,ans(nfitm),gain,rnoise,rad,x,y,xcen,ycen,step,flux,
     .     model,crit,dis,prad,nfl,gpsf,ypsf,sky,rflux,
     .     pmag(stmax,nsmax),ans1(nfitm),itime

      ans1(1)=ans(1)
      ans1(2)=ans(2)
      ans1(3)=ans(3)
      ans1(4)=ans(4)
      ifun1=ifun-4
      do 20 k=1,hasstar(stfrm)
         id=3*k
         xcen=ans(id+3)
         ycen=ans(id+4)
         ans1(5)=ans(id+2)
         ans1(6)=xcen
         ans1(7)=ycen
         step=2.*rad/real(ngrid)
         flux=0.
         do 10 i=1,ngrid
            do 11 j=1,ngrid
               x=xcen+step*real(i-ngrid/2)
               y=ycen+step*real(j-ngrid/2)
               prad=ypsf(x,y,ans,nfit,rad)
               dis=( (x-xcen)**2.0 + (y-ycen)**2.0 )**0.5
               crit=prad-dis
               model=gpsf(x,y,ans1,nfit,ifun1)
               if(crit.gt.0.5) then
                  flux=flux+model*step*step
                  nfl=nfl+1.0
               elseif((crit.ge.-0.5).and.(crit.le.0.5)) then
                  flux=flux+model*step*step*(crit+0.5)
                  nfl=nfl+1.0*(crit+0.5)
               endif
 11         continue
 10      continue

         gmag=25.0-2.5*log10(flux/itime*gain)
C         pmag(stfrm,k)=gmag
         pmag(stfrm,k)=flux
c         write(6,*) k,gmag
 20   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function ypsf(x,y,a,nfitpsf,rad)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Finding isophotal radii from Gaussian fit.. need to add other
C     PSF profiles, although Gaussian seems to work fine.
      implicit none
      integer nfitpsf,i
      real th,z,a(nfitpsf),temp(17),Sx,Sy,Sxy,xcoo,ycoo,fwhmx,x,y,
     .     Ib,Ic,rad,m,b,xd,yd,zt,ztold,pt,ptold,xt,yt,sign,p,pw,
     .     stepsize,expt1,expt2,expt3,expt
      logical loop

      Ib=a(1)
      Ic=a(5)
      xcoo=a(6)
      ycoo=a(7)
      Sx=a(2)
      Sy=a(3)
      Sxy=a(4)
      fwhmx=2.0*Sx*Sqrt(log(2.0))
     
      if((fwhmx.gt.10.0).or.(fwhmx.lt.0.)) then
         ypsf=rad
         return
      endif

      m=(y-ycoo)/(x-xcoo)
      b=y-x*m
      
C     Find flux at Fwhm_x
      z=Ib+Ic*exp(-(fwhmx*fwhmx)/(Sx*Sx))

      stepsize=0.1
      ztold=Ic+Ib
      ptold=0.
      xt=xcoo
      yt=ycoo
      sign=(x-xcoo)/abs(x-xcoo)
      
      i=0
      loop=.true.
      do while(loop)
         i=i+1
         xt=xt+sign*stepsize
         yt=m*xt+b
         pt=((xt-xcoo)**2.0+(yt-ycoo)**2.0)**0.5
         xd=xt-xcoo
         yd=yt-ycoo
         expt1=-xd*xd/(Sx*Sx)
         expt2=-yd*yd/(Sy*Sy)
         expt3=2*Sxy*xd*yd/(Sx*Sy)
         zt=Ib+Ic*exp(expt1+expt2+expt3)
         if(zt.lt.z) then
            loop=.false.
         else
            ztold=zt
            ptold=pt
         endif
         if(i.gt.200) then
            ypsf=rad
            return
            loop=.false.
         endif
      enddo

      p=pt/(zt-z)**2.0 + ptold/(ztold-z)**2.0
      pw=1.0/(zt-z)**2.0 + 1.0/(ztold-z)**2.0
      ypsf=p/pw*rad/fwhmx
      
c      write(6,*) i,p/pw,fwhmx
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function findsky(n,skys,stdsky)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,i,xmax,ymax,j,niter
      parameter(xmax=600,ymax=600)
      real skys(n),mean,cskys(xmax*ymax),stdev,findmed,clip,med,std,
     .     stdsky,osky


      if(n.eq.0) then
         write(6,*) "no pixels for sky determination!"
         findsky=0.
         goto 999
      endif

C     Sigma clip for rejecting extreme pixel values
      clip=2.5

      mean=0.
      do 11 i=1,n
         mean=mean+skys(i)
 11   continue
      mean=mean/real(n)

      med=findmed(n,skys)
      std=stdev(n,skys,mean)

      osky=mean

      niter=0
 13   j=0
      niter=niter+1
      do 10 i=1,n
         if(abs(skys(i)-med).lt.(clip*std)) then
            j=j+1
            cskys(j)=skys(i)
         endif
 10   continue

      med=findmed(j,cskys)
      std=stdev(j,cskys,mean)

      mean=0.
      do 12 i=1,j
         mean=mean+cskys(i)
 12   continue
      mean=mean/real(j)

      stdsky=std

c      if(med.lt.mean) then
c         findsky=3*med-2*mean
c      else
         findsky=mean
c      endif

      if(niter.gt.5) goto 14
      if(findsky.ne.osky) then
         osky=findsky
         goto 13
      endif

 14   continue

 999  return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine nstar(fitsdata,naxes,datamin,datamax,stn,stx,sty,
     .     hasstar,ns,a,ifun,nfit,zsky,stdsky,sat,fitsky,rad,
     .     skyrad,pcoo,rstx,rsty,stfrm,sflag,sint)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer xmax,ymax,naxes(2),stfrm,stn,npt,nfitm,ma,niter,ifun,
     .     stmax,nsmax,i,ns,nfit,fitsky,nsky,j,k,ii,sint
      parameter(xmax=600,ymax=600,stmax=20,nsmax=100,nfitm=200)
      integer stx(stmax,2),sty(stmax,2),hasstar(stmax),ia(nfitm),
     .     arej(nfitm),rstx(stmax,2),rsty(stmax,2),nst,sflag(stmax),
     .     nsat(stmax)
      real fitsdata(xmax,ymax),datamin,datamax,
     .     pcoo(stmax,nsmax,2),x(xmax*ymax),y(xmax*ymax),z(xmax*ymax),
     .     sig(xmax*ymax),a(nfitm),covar(nfitm,nfitm),
     .     alpha(nfitm,nfitm),chisq,alamda,chi2old,diff,ratio,zsky,
     .     sigcut,stdsky,gpsf,xc,yc,sat(stmax),pv,crit,
     .     skys(xmax*ymax),rad,skyrad,dis,xcoo,ycoo,xc2,yc2,findsky,
     .     mdis,xt,yt,itime
      logical loop

C     saturation flag is sflag
C     number of saturated pixels is contained in nsat

C     setting number of fitting parameters.
CS      nfit=4+3*(ns)
      nfit=4+3*hasstar(stfrm)
c      write(6,*) "nfit:",nfit,ns
      if (nfit.gt.nfitm) pause "Must increase nfitm to match nfit"

C     NPT is the number of pixels included in fit.
      npt=0
C     PV is the peak value
      pv=-1.0e30
C     NSKY is the number of skypixels
      nsky=0
      nsat(stfrm)=0
      do 10 i=1,naxes(1)
         do 20 j=1,naxes(2)
            k=stfrm
CS            do 30 k=1,stn
               if((i.ge.stx(k,1)).and.(i.le.stx(k,2)).and.
     .              (j.ge.sty(k,1)).and.(j.le.sty(k,2))) then
                  xc=i+rstx(k,1)-stx(k,1)
                  yc=j+rsty(k,1)-sty(k,1)
c                  if(k.eq.3) write(6,*) k,rstx(k,1),rsty(k,1),xc,yc

                  crit=-99.9e30
                  nst=hasstar(k)
                  do 11 ii=1,nst
                     dis= ((real(i-stx(k,1)+1)-pcoo(k,ii,1))**2.0+
     .                    (real(j-sty(k,1)+1)-pcoo(k,ii,2))**2.0)**0.5
                     if(rad-dis.ge.crit) then
                        crit=rad-dis
                        mdis=dis
                     endif
 11               continue
                  dis=mdis

                  if(crit.ge.0) then
                     if(fitsdata(i,j).lt.sat(stfrm))then
                        npt=npt+1
                        x(npt)=xc
                        y(npt)=yc
                        z(npt)=fitsdata(i,j)
c                     sig(npt)=100.0*(sqrt(abs(z(npt)))/z(npt))
                        sig(npt)=sqrt(abs(z(npt)))
                     else
                        nsat(stfrm)=nsat(stfrm)+1
                     endif
c                     write(6,*) "gg",k,x(npt),y(npt),z(npt),sig(npt)
                     if(fitsdata(i,j).gt.pv) pv=fitsdata(i,j)
                  endif
                  if(dis.gt.skyrad) then
                     nsky=nsky+1
                     skys(nsky)=fitsdata(i,j)
c                     write(6,*) "sky",i,j,fitsdata(i,j),nsky,skys(nsky)
                  endif
               endif
CS 30         continue
 20      continue
 10   continue

C     Find the skyvalue from sky pixels.
      zsky=findsky(nsky,skys,stdsky)

C     Saturation warning 
c      if(stfrm.eq.1) write(6,*) stfrm,pv,sat(stfrm)
      if((pv.ge.sat(stfrm)).or.(pv.le.-1.0e10)) then
         sflag(stfrm)=1
         if(nsat(stfrm).gt.1) sflag(stfrm)=2
         write(6,502) "Warning: saturation on frame ",stfrm," pv: ",
     .        pv,sat(1),zsky,nsat(stfrm),sflag(stfrm),npt
C     if we have more than one saturated pixel, forget about it.
      endif
 502  format(A29,I2,A5,3(F12.2,1X),I2,1X,I2,1X,I4)
c      if(npt.lt.nfit) then
c         write(6,*) "Warning: npt is less than nfit",npt,nfit
c      endif

C     when the image saturates, there may be no pixels for fitting
C     so let just get out of the routine before segmentation
C     problems arise.
      if(npt.eq.0) then
         sflag(stfrm)=2
         goto 101
      endif

C     initial guesses..
      
      ma=nfit
      
      k=0
CS      do 40 j=1,stn
      j=stfrm
         nst=hasstar(j)
         do 41 i=1,nst
            k=k+1
C     Convert to real co-ordinates
            xcoo=pcoo(j,i,1)+rstx(j,1)-1
            ycoo=pcoo(j,i,2)+rsty(j,1)-1
            a(3*(k-1)+5)=fitsdata(pcoo(j,i,1)+stx(j,1)-1,
     .           pcoo(j,i,2)+sty(j,1)-1)-zsky
            a(3*(k-1)+6)=xcoo
            a(3*(k-1)+7)=ycoo
c            write(6,*) "pp",a(3*(k-1)+5),a(3*(k-1)+6),
c     .           a(3*(k-1)+7),pcoo(j,i,1)+stx(j,1)-1,
c     .           pcoo(j,i,2)+sty(j,1)-1
 41      continue
CS 40   continue

C     tells minization routine which variables are held fixed.
C     seting ia(i)=0 fixes initial guesses at input values.
      do 13 i=1,ma
         ia(i)=1
 13   continue
      if(fitsky.eq.0) ia(1)=0

C     give a(1) value from findsky routine
      a(1)=zsky

C     Initial guess for shape parameters 
      a(2)=1.8
      a(3)=1.6
      a(4)=-0.08

C     if saturated.. skip minimization process
      if((sflag(stfrm).eq.2).or.((sflag(stfrm).eq.1).and.(sint.eq.0))) 
     .     goto 101

C     initialize fitting routine
C      loop=.true.
C      do while(loop)
      alamda=-1.0
      niter=0
      call mrqmin(x,y,z,sig,npt,a,ia,ma,covar,alpha,nfitm,chisq,
     .        alamda,ifun)
C      write(6,*) "mrq:",npt,ma,nfitm,chisq,alamda,ifun


C     Now we try to converge to local min. -> solution
C     This section needs improvement.  chisq is way too big.
 100  niter=niter+1

c      write(6,*) "niter:",niter,alamda,chisq
      chi2old=chisq
      call mrqmin(x,y,z,sig,npt,a,ia,ma,covar,alpha,nfitm,chisq,
     .     alamda,ifun)
      diff =abs(chi2old-chisq)
      ratio=diff/chisq
C     if(stfrm.eq.2) write(6,*) niter,alamda,diff,ratio,chisq
      
      if ((chi2old .gt. chisq .and. ratio .lt. 1.0E-05 .and.
     .     diff .lt. 1.0E-04) .or. (niter .ge. 400) .or.
     .     (alamda .lt. 1.0E-08) .or. (alamda .gt.1.0E08))
     .     then
         alamda=0
         call mrqmin(x,y,z,sig,npt,a,ia,ma,covar,alpha,nfitm,
     .        chisq,alamda,ifun)
      else
         goto 100
      end if

 101  continue

 500  format(I3,1X,F8.2,1X,F8.2,6(1X,F7.2))
 501  format(10(I7,1X))

C      nst=0
C      write(6,*) (k,":",a(k),k=1,4)
C      j=stfrm
CCS      do 51 j=1,stn
C         do 50 i=1,hasstar(j)
C            nst=nst+1
C            write(6,*) j,(3*(nst-1)+4+k,":",a(3*(nst-1)+4+k),k=1,3) 
C 50      continue
CCS 51   continue

      zsky=a(1)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine funcs2(x,y,a,z,dzda,na,ifun)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Contains the profile formulae and derivatives with respect to the
C     fitting parameters.
C     Functional form is controlled by ifun
C     1 - Gauss
C     2 - Moffat with fixed beta (current 3.5)
C     3 - Lorentz
C     4 - Penny (Gauss + Lorentz)
C     5 - N Gauss - fit more than one centre
C     6 - N Moffat - ditto
C     7 - heavy step function - used for Fabry photometry
      implicit none
      integer na,i,ifun,nstar
      real x,y,a(na),z,dzda(na),xd,yd,expt1,expt2,expt3,expt,Ib,Ic,
     .     xcoo,ycoo,sx,sy,sxy,tmp,beta,temp(2),Gz,Lz,gexpt,pc,Ip,Slxy,
     .     Sgxy,sech
      
      pc=-0.693

      beta=3.5
c      beta=1.5

      Ib=a(1)
      Ic=a(5)
      xcoo=a(6)
      ycoo=a(7)
      sx=a(2)
      sy=a(3)
      sxy=a(4)

      if(ifun.eq.4) then
         slxy=a(7)
         sgxy=a(8)
         Ip=a(9)
      endif


      if(ifun.gt.4)then
C     Double Gaussian
         if(na.eq.7) then
            nstar=1
         else
            nstar=(na-7)/3+1
         endif
         if(ifun.eq.5) call ngaussfit(nstar,x,y,na,a,z,dzda)
         if(ifun.eq.6) call nmoffatfit(nstar,x,y,na,a,z,dzda)
      endif

         
      xd=x-xcoo
      yd=y-ycoo

      if(ifun.eq.1) then
C     Gaussian Fit
         expt1=-((xd/Sx)**2.0)
         expt2=-((yd/Sy)**2.0)
         expt3=2.0*Sxy*(xd*yd)/(Sx*Sy)
         
         expt=exp(expt1+expt2+expt3)
         
         z= Ic*expt+Ib
         
         dzda(1)=1.0
         dzda(2)=expt
         tmp=2*xd/(Sx*Sx)-2*Sxy*yd/(Sx*Sy)
         dzda(3)=expt*Ic*tmp
         tmp=-2*Sxy*xd/(Sx*Sy)+2*yd/(Sy*Sy)
         dzda(4)=expt*Ic*tmp
         tmp=2*(xd*xd)/(sx*sx*sx) - 2*Sxy*(xd*yd)/(sx*sx*sy)
         dzda(5)=expt*Ic*tmp
         tmp=-2*Sxy*xd*yd/(sx*sy*sy) + 2*(yd*yd)/(sy*sy*sy)
         dzda(6)=expt*Ic*tmp
         tmp=Ic*xd*yd/(Sx*Sy)
         dzda(7)=2*expt*tmp

c      write(6,*) x,y,(dzda(i),i=1,na)
c      read(5,*)

      elseif(ifun.eq.2) then
C     Moffat Fit
         
         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Sxy
         expt=1.0+expt1+expt2+expt3
         
         z=Ib+Ic/(expt**beta)

         dzda(1)=1.0
         dzda(2)=1.0/(expt**beta)
         tmp=-2.0*xd/(Sx*Sx)-Sxy*yd
         dzda(3)=-Ic*beta*tmp*(expt**(-1.0-beta))
         tmp=-Sxy*xd-2.0*yd/(Sy*Sy)
         dzda(4)=-Ic*beta*tmp*(expt**(-1.0-beta))
         tmp=2*xd*xd/(Sx*Sx*Sx)
         dzda(5)=Ic*beta*tmp*(expt**(-1.0-beta))
         tmp=2*yd*yd/(Sy*Sy*Sy)
         dzda(6)=Ic*beta*tmp*(expt**(-1.0-beta))
         dzda(7)=-Ic*xd*yd*beta*(expt**(-1.0-beta))
C     for reference 8 is the derivative with respect to beta
C     incase someone wants to include beta as a fitting parameter.
C     Its easier to just measure it a few time with something like IRAF
C     and then pick a good choice.
c         dzda(8)=-Ic*(expt**(-beta))*Log(expt)

      elseif(ifun.eq.3) then
C     Lorentzian Fit

         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Sxy
         expt=1.0+expt1+expt2+expt3

         z=Ib+Ic/(expt)

         dzda(1)=1.0
         dzda(2)=1.0/expt
         tmp=-2*xd/(Sx*Sx)-Sxy*yd
         dzda(3)=-Ic*tmp/(expt**2.0)
         tmp=-Sxy*xd-2*yd/(Sy*Sy)
         dzda(4)=-Ic*tmp/(expt**2.0)
         tmp=2*xd*xd/(Sx*Sx*Sx)
         dzda(5)=Ic*tmp/(expt**2.0)
         tmp=2*yd*yd/(Sy*Sy*Sy)
         dzda(6)=Ic*tmp/(expt**2.0)
         dzda(7)=-Ic*xd*yd/(expt**2.0)


      elseif(ifun.eq.4) then
C     Penny "Lane" Function
C     The pretty nurse is selling poppies from a tray

         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Slxy
         expt=expt1+expt2+expt3
         
         Lz=(1.0-Ip)/expt

         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Sgxy
         gexpt=expt1+expt2+expt3

         Gz=Ip*Exp(pc*gexpt)

         z=Ib+Ic*(Lz+Gz)

         dzda(1)=1.0
         dzda(2)=Gz+Lz
         temp(1)=-2*xd/(Sx*Sx)-Sgxy*yd
         temp(2)=-2*xd/(Sx*Sx)-Slxy*yd
         dzda(3)=Ic*(temp(1)*pc*Gz-temp(2)*Lz/expt)
         temp(1)=-Sgxy*xd-2*yd/(Sy*Sy)
         temp(2)=-Slxy*xd-2*yd/(Sy*Sy)
         dzda(4)=Ic*(temp(1)*pc*Gz-temp(2)*Lz/expt)
         temp(1)=xd*xd/(Sx*Sx*Sx)
         temp(2)=2.0*temp(1)
         dzda(5)=Ic*(pc*pc*Gz*temp(1)+temp(2)*Lz/expt)
         temp(1)=yd*yd/(Sy*Sy*Sy)
         temp(2)=2.0*temp(2)
         dzda(6)=Ic*(pc*pc*Gz*temp(1)+temp(2)*Lz/expt)
         dzda(7)=Ic*xd*yd*Lz/expt
         dzda(8)=Ic*xd*yd*Gz*pc
         dzda(9)=Ic*(Exp(pc*gexpt)-1.0/expt)

      endif
      
C     Heavy Step function
C     f(x)=a(1)+a(2)*tanh((x-a(3))/a(4))
      if(ifun.eq.7) then      
         tmp=(x-a(3))/a(4)
         z=a(1)+a(2)*tanh(tmp)
         dzda(1)=1.0
         dzda(2)=tanh(tmp)
         dzda(3)=-a(2)*sech(tmp)*sech(tmp)/a(4)
         dzda(4)=dzda(3)*tmp
c         write(6,*) "dzda",(dzda(i),i=1,4)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function sech(z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real z
      
      sech=1.0/cosh(z)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine nmoffatfit(nstar,x,y,na,a,z,dzda)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nstar,na,nmax,i,offset,id
      parameter(nmax=2)
      real x,y,a(na),z,dzda(na),Ic,xcoo,ycoo,sx,sy,sxy,expt1,expt2,
     .     expt3,expt,Ib,tmp,xd,yd,beta

      beta=3.5

      z=0
      do 5 i=1,na
         dzda(i)=0
 5    continue

      Ib=a(1)
      sx=a(2)
      sy=a(3)
      sxy=a(4)

      dzda(1)=1.0

c     nst=0
c     write(6,*) (k,":",a(k),k=1,4)
c     do 51 j=1,stn
c        do 50 i=1,hasstar(j)
c           nst=nst+1
c           write(6,*) j,(3*(nst-1)+4+k,":",a(3*(nst-1)+4+k),k=1,3) 
c50      continue
c51   continue

      do 10 i=1,nstar
         id=3*i
         Ic=a(id+2)
         xcoo=a(id+3)
         ycoo=a(id+4)

         xd=x-xcoo
         yd=y-ycoo    

         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Sxy
         expt=1.0+expt1+expt2+expt3
         
         z=Ic/(expt**beta) + z

         dzda(id+2)=1.0/(expt**beta)
         tmp=-2.0*xd/(Sx*Sx)-Sxy*yd
         dzda(id+3)=-Ic*beta*tmp*(expt**(-1.0-beta))
         tmp=-Sxy*xd-2.0*yd/(Sy*Sy)
         dzda(id+4)=-Ic*beta*tmp*(expt**(-1.0-beta))
         tmp=2*xd*xd/(Sx*Sx*Sx)
         dzda(2)=Ic*beta*tmp*(expt**(-1.0-beta))+dzda(2)
         tmp=2*yd*yd/(Sy*Sy*Sy)
         dzda(3)=Ic*beta*tmp*(expt**(-1.0-beta))+dzda(3)
         dzda(4)=-Ic*xd*yd*beta*(expt**(-1.0-beta))+dzda(4)
 10   continue

      z=z+Ib

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ngaussfit(nstar,x,y,na,a,z,dzda)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nstar,na,nmax,i,offset,id
      parameter(nmax=2)
      real x,y,a(na),z,dzda(na),Ic,xcoo,ycoo,sx,sy,sxy,expt1,expt2,
     .     expt3,expt,Ib,tmp,xd,yd

      z=0
      do 5 i=1,na
         dzda(i)=0
 5    continue

      Ib=a(1)
      sx=a(2)
      sy=a(3)
      sxy=a(4)

      dzda(1)=1.0
      
      do 10 i=1,nstar
         id=3*i
         Ic=a(id+2)
         xcoo=a(id+3)
         ycoo=a(id+4)

         xd=x-xcoo
         yd=y-ycoo

         expt1=-((xd/Sx)**2.0)
         expt2=-((yd/Sy)**2.0)
         expt3=2.0*Sxy*(xd*yd)/(Sx*Sy)
         
         expt=exp(expt1+expt2+expt3)
         
         z= Ic*expt+z
         
         dzda(id+2)=expt
         tmp=2*xd/(Sx*Sx)-2*Sxy*yd/(Sx*Sy)
         dzda(id+3)=expt*Ic*tmp
         tmp=-2*Sxy*xd/(Sx*Sy)+2*yd/(Sy*Sy)
         dzda(id+4)=expt*Ic*tmp
         tmp=2*(xd*xd)/(sx*sx*sx) - 2*Sxy*(xd*yd)/(sx*sx*sy)
         dzda(2)=expt*Ic*tmp+dzda(2)
         tmp=-2*Sxy*xd*yd/(sx*sy*sy) + 2*(yd*yd)/(sy*sy*sy)
         dzda(3)=expt*Ic*tmp+dzda(3)
         tmp=Ic*xd*yd/(Sx*Sy)
         dzda(4)=2*expt*tmp+dzda(4)        

 10   continue
      
      z=z+Ib
c      write(6,*) x,y,z,(dzda(i),i=1,na)
c      read(5,*) 

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function gpsf(x,y,a,nfitpsf,ifun)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Returns value of PSF given offset from center (x,y)
      implicit none
      integer nfitpsf,ifun,i,nstar,offset
      real x,y,a(nfitpsf),Ib,Ic,xcoo,ycoo,Sx,Sy,Sxy,xd,yd,expt1,expt2,
     .     expt3,beta,temp(4),expt,Ip,Slxy,Sgxy,Gz,Lz,gexpt,pc

      
      pc=-0.693
      beta=3.5

      Ib=a(1)
      Ic=a(5)
      xcoo=a(6)
      ycoo=a(7)
      Sx=a(2)
      Sy=a(3)
      Sxy=a(4)   

      if(ifun.eq.4) then
         Slxy=a(7)
         Sgxy=a(8)
         Ip=a(9)
      endif

      if(ifun.gt.4) then
         nstar=(nfitpsf-7)/3+1
      endif

      xd=x-xcoo
      yd=y-ycoo

      if(ifun.eq.1) then
         
         expt1=-xd*xd/(Sx*Sx)
         expt2=-yd*yd/(Sy*Sy)
         expt3=2*Sxy*xd*yd/(Sx*Sy)
         gpsf=Ic*exp(expt1+expt2+expt3)
      elseif(ifun.eq.2) then
         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Sxy
         expt=1.0+expt1+expt2+expt3
         gpsf=Ic/(expt**beta)

c         temp(1)=beta-1
c         temp(2)=Sx*Sy
c         temp(3)=1+(xd/Sx)**2.0+(yd/Sy)**2.0+xd*yd*Sxy
c         temp(4)=temp(3)**beta
c         expt1=temp(1)/(temp(2)*temp(4))
c         gpsf=Ic*expt1+Ib
      elseif(ifun.eq.3) then
         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Sxy
         gpsf=Ic/(1.0+expt1+expt2+expt3)
      elseif(ifun.eq.4) then

         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Slxy
         expt=expt1+expt2+expt3
         
         Lz=(1.0-Ip)/expt

         expt1=xd*xd/(Sx*Sx)
         expt2=yd*yd/(Sy*Sy)
         expt3=xd*yd*Sgxy
         gexpt=expt1+expt2+expt3

         Gz=Ip*Exp(pc*gexpt)

         gpsf=Ic*(Lz+Gz)
      elseif(ifun.eq.5) then
         gpsf=0
         do 10 i=1,nstar            
            Ic=a(3*(i-1)+2)
            xcoo=a(3*(i-1)+3)
            ycoo=a(3*(i-1)+4)
            
            xd=x-xcoo
            yd=y-ycoo

            expt1=-((xd/Sx)**2.0)
            expt2=-((yd/Sy)**2.0)
            expt3=2.0*Sxy*(xd*yd)/(Sx*Sy)
            
            expt=exp(expt1+expt2+expt3)
            
            gpsf= Ic*expt+gpsf
 10      continue

      elseif(ifun.eq.6) then
         gpsf=0
         do 20 i=1,nstar
            Ic=a(3*(i-1)+2)
            xcoo=a(3*(i-1)+3)
            ycoo=a(3*(i-1)+4)
            
            xd=x-xcoo
            yd=y-ycoo    
            
            expt1=xd*xd/(Sx*Sx)
            expt2=yd*yd/(Sy*Sy)
            expt3=xd*yd*Sxy
            expt=1.0+expt1+expt2+expt3
            
            gpsf= Ic/(expt**beta) + gpsf
 20      continue

      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine starident(fitsdata,naxes,datamin,datamax,stn,stx,sty,
     .     pcool,hasstar,rstx,rsty,ns)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This routine allows the user to select the stars that are 
C     interesting for photometry - better than automated methods.
      implicit none
      integer naxes(2),xmax,ymax,stn,stmax,nsmax,i,j,k,nx,ny,ns
      parameter(xmax=600,ymax=600,stmax=20,nsmax=100)
      integer stx(stmax,2),sty(stmax,2),hasstar(stmax),rstx(stmax,2),
     .     rsty(stmax,2),nstot,xct,yct
      real fitsdata(xmax,ymax),datamin,datamax,pcool(stmax,nsmax,2),
     .     ddata(xmax,ymax),z1,z2,xc,yc,xc2,yc2,sumx,sumy,fsum
      character ch

c      pcool(1,1,1)=12.54
c      pcool(1,1,2)=11.79
c      pcool(2,1,1)=13.57
c      pcool(2,1,2)=14.27
c      pcool(3,1,1)=12.88
c      pcool(3,1,2)=12.22
c      hasstar(1)=1
c      hasstar(2)=1
c      hasstar(3)=1
c      goto 999
      


      nstot=0
      do 10 i=1,stn
         hasstar(i)=0
         nstot=nstot+ns
         ns=0
         do 11 j=stx(i,1),stx(i,2)
            do 12 k=sty(i,1),sty(i,2)
               nx=j-stx(i,1)+1
               ny=k-sty(i,1)+1
               ddata(nx,ny)=fitsdata(j,k)
c               write(6,*) nx,ny,ddata(nx,ny)
 12         continue
 11      continue

         z1=-1.
         z2=-1.
c         z1=0.
c         z2=450.
         nx=stx(i,2)-stx(i,1)+1
         ny=sty(i,2)-sty(i,1)+1

c         write(6,*) i,nx,ny,stx(i,1),stx(i,2),sty(i,1),sty(i,2)
         call sparewindow(nx,ny,ddata,z1,z2)        
C     now the user can pick stars from the frame.
         write(6,*) "Please begin to select stars, hit e to end."
C     The mouse buttons are: A,D,X
 13      call pgband(7,0,0.0,0.0,xc,yc,ch)
C     increase counter for total number of stars in frame.
           ns=ns+1
C     convert to fits co-ordinates
           xc2=xc+stx(i,1)-1
           yc2=yc+sty(i,1)-1
           if((ch.eq."E").or.(ch.eq."e")) goto 14
           if((xc2.lt.stx(i,1)).or.(xc2.gt.stx(i,2)).or.
     .          (yc2.lt.sty(i,1)).or.(yc2.gt.sty(i,2))) then
              write(6,*)
     .             "Error: Those co-ordinates are outside the raster"
              goto 13
           endif
C     convert to CCD co-ordinates
c         xc=xc+rstx(i,1)-1
c         yc=yc+rsty(i,1)-1
            sumx=0.
            sumy=0.
            fsum=0.
            do 40 j=1,5
               do 41 k=1,5
                  xct=int(xc2)-3+j
                  yct=int(yc2)-3+k
                  if((xct.ge.stx(i,1)).and.(xct.le.stx(i,2)).and.
     .             (yct.ge.sty(i,1)).and.(yct.le.sty(i,2))) then
                     sumx=sumx+fitsdata(xct,yct)*real(xct)
                     sumy=sumy+fitsdata(xct,yct)*real(yct)
                     fsum=fsum+fitsdata(xct,yct)
                  endif
 41            continue
 40         continue
            xc=sumx/fsum+real(1-stx(i,1))
            yc=sumy/fsum+real(1-sty(i,1))

           write(6,500) xc,yc   !,ch
 500       format(2(F6.2,1X),A1)

           pcool(i,ns,1)=xc
           pcool(i,ns,2)=yc
C     do psf subtraction stuff here.
         if((ch.ne."e").and.(ch.ne."E")) goto 13
 14      ns=ns-1
         hasstar(i)=ns

 10   continue
      ns=nstot

 999  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine displaypeaks(pcoo,stn,naxes,hasstar,stx,sty)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Marks location of stars on frames
      implicit none
      integer i,ny,nx,stn,stmax,naxes(2),nsmax,j
      parameter(stmax=20,nsmax=100)
      integer hasstar(stmax),stx(stmax,2),sty(stmax,2)
      real pcoo(stmax,nsmax,2),x1,x2,y1,y2

      call pgsvp(0.05,0.95,0.05,0.95)
      call pgwnad(0.0,1.0,0.0,1.0)

      nx=naxes(1)
      ny=naxes(2)

C     Get bounds of window.
      call pgqvp(0,x1,x2,y1,y2)
c      write(6,*) "window",x1,x2,y1,y2

C     Set up port size
c      call pgvport(x1,x1+(x2-x1)*real(nx)/real(ny),y1,y2)
      call pgvport(x1,x1+(x2-x1)*real(nx)/real(max(nx,ny)),y1,
     .   y1+(y2-y1)*real(ny)/real(max(nx,ny)))
C     Set scale
      call pgwindow(1.0,real(nx)+1.0,1.0,real(ny)+1.0)
C     mark boundaries and tick marks
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
c      CALL PGLABEL("X Pixel","Y Pixel"," ")
C     Draw the dots!
      do 10 i=1,stn
         do 20 j=1,hasstar(i)
            call pgsci(3)
            call pgpt1(pcoo(i,j,1)+stx(i,1)-0.5,
     .           pcoo(i,j,2)+sty(i,1)-0.5,17)
            call pgsci(1)
 20      continue
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sparewindow(nx,ny,ddata,z1,z2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Displays pixel values in a colour frame in device.
C     Handy for looking at processed data, like residuals.
      implicit none
      integer nx,ny,xmax,ymax,i,j,npt,ncol
      parameter(xmax=600,ymax=600)
      parameter(ncol=64)
      integer ia(xmax,ymax)
      real ddata(xmax,ymax),pts(xmax*ymax),z1,z2,dcut,x1,x2,y1,y2,R,G,B,
     .     z1n,z2n

C     dead parameter.  To be removed.
      dcut=0.98
      
      npt=0
      do 20 i=1,nx
         do 21 j=1,ny
            npt=npt+1
            pts(npt)=ddata(i,j)
 21         continue
 20   continue

C     Finds range for colour mapping of data.
      if(z1.eq.z2)call autoscale(npt,pts,z1,z2,dcut)
c      write(6,*) z1,z2


C     Map data onto colour map.
      do 10 i=1,nx
         do 11 j=1,ny
            IA(i,j)=int((ddata(i,j)-z1)/(z2-z1)*(NCOL-1)+16)
c            if(yn.eq.0) then
c               if(ddata(i,j).lt.z1) ia(i,j)=0
c            else
               if(ddata(i,j).lt.z1) ia(i,j)=16
c            endif
            if(ddata(i,j).gt.z2) ia(i,j)=ncol+15
 11      continue
C     write(6,*) (IA(i,k),k=1,ny)
 10   continue


C     set up pgplot window
      call pgscr(0,0.0,0.3,0.2)
      call pgsvp(0.65,0.95,0.07,0.37)
      call pgwnad(0.0,1.0,0.0,1.0)
      
C     Setup colour-map
      do 300 i=1,ncol
         R = REAL(I-1)/REAL(NCOL-1)*0.8 + 0.2
         G = MAX(0.0, 2.0*REAL(I-1-NCOL/2)/REAL(NCOL-1))
         B = 0.2 + 0.4*REAL(NCOL-I)/REAL(NCOL)
         CALL PGSCR(I+15, R, G, B)
 300  CONTINUE

C     Draw the pretty picture
      call pgpixl(ia,xmax,ymax,1,nx,1,ny,
     .     0.0,real(nx)/real(max(nx,ny)),
     .     0.0,real(ny)/real(max(nx,ny)))
      
C     get the bounds of the window so we can add labels.
      call pgqvp(0,x1,x2,y1,y2)
      write(6,*) "window",x1,x2,y1,y2
      
C     Define port size
c      call pgvport(x1,x1+(x2-x1)*real(nx)/real(ny),y1,y2)
      call pgvport(x1,x1+(x2-x1)*real(nx)/real(max(nx,ny)),y1,
     .   y1+(y2-y1)*real(ny)/real(max(nx,ny)))
C     define scale
      call pgwindow(1.0,real(nx)+1.0,1.0,real(ny)+1.0)
C     add tick marks, labels and borders
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)

c      write(6,*) "auto",z1,z2
c      read(5,*)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine displayfits(fitsdata,naxes,datamin,datamax,tix,tiy,
     .     stn,stx,sty,nstack)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     displays the fits data using pgplot libs
      implicit none
      integer stmax
      parameter(stmax=20)
      integer naxes(2),xmax,ymax,i,j,ncol,nx,ny,k,tix(2),tiy(2),
     .   stn,stx(stmax,2),sty(stmax,2),zs1(stmax),zs2(stmax),npt,nstack
      parameter(xmax=600,ymax=600,ncol=64)
      integer ia(xmax,ymax),yn,dumi
      real fitsdata(xmax,ymax),R,G,B,z1,z2,datamin,datamax,zt1,zt2,
     .     pts(xmax*ymax),dcut(stmax),dumr,satlev


      satlev=real(nstack)*2.0**14
      do 4 i=1,stmax
         dcut(i)=1.00
 4    continue

c      dcut(1)=1.00
c      dcut(2)=1.00
c      dcut(3)=1.00
c      dcut(4)=1.00
c      dcut(5)=1.00
c      dcut(6)=1.00

      do 5 i=1,stn
         npt=0
         do 6 j=1,naxes(1)
            do 7 k=1,naxes(2)
               if((j.ge.stx(i,1)).and.(j.le.stx(i,2)).and.
     .              (k.ge.sty(i,1)).and.(k.le.sty(i,2))) then
                  if((fitsdata(j,k).gt.-900.0).and.
     .                 (fitsdata(j,k).lt.satlev)) then
                     npt=npt+1
                     pts(npt)=fitsdata(j,k)
c                     write(6,*) i,j,k,fitsdata(j,k)
c                     write(6,*) pts(npt)
c                  else
c                     write(6,*) "oop?",j,k,fitsdata(j,k)
                  endif
               endif
 7          continue
 6       continue
         call autoscale(npt,pts,z1,z2,dcut(i))
         zs1(i)=z1
         zs2(i)=z2
c         write(6,*) "as:",i,z1,z2,npt,
c     .        stx(i,1),stx(i,2),sty(i,1),sty(i,2)
 5    continue
      
c      zs1(1)=0.0
c      zs2(1)=2000.0

      npt=0
      do 8 j=1,naxes(1)
         do 9 k=1,naxes(2)
            if((j.ge.tix(1)).and.(j.le.tix(2)).and.
     .           (k.ge.tiy(1)).and.(k.le.tiy(2))) then
               if((fitsdata(j,k).ge.datamin)
     .              .and.(fitsdata(j,k).le.datamax)) then
                  npt=npt+1
                  pts(npt)=fitsdata(j,k)
               endif
            endif
 9       continue
 8    continue
      call autoscale(npt,pts,z1,z2,1.00)

      zt1=z1
      zt2=z2
      nx=naxes(1)
      ny=naxes(2)

C      write(6,*) nx,ny

      do 10 i=1,nx
         do 11 j=1,ny
            yn=0
            if((i.ge.tix(1)).and.(i.le.tix(2)).and.(j.ge.tiy(1)).and.
     .           (j.le.tiy(2))) then
               z1=zt1
               z2=zt2
               yn=1
            else
               z1=0.0
               z1=65000.0
               do 12 k=1,stn
                  if((i.ge.stx(k,1)).and.(i.le.stx(k,2)).and.
     .                 (j.ge.sty(k,1)).and.(j.le.sty(k,2))) then
                     z1=zs1(k)
                     z2=zs2(k)
                     yn=1
                  endif
 12            continue
            endif
            IA(i,j)=int((fitsdata(i,j)-z1)/(z2-z1)*(NCOL-1))+16
            if(yn.eq.0) then
               if(fitsdata(i,j).lt.z1) ia(i,j)=0
            else
               if(fitsdata(i,j).lt.z1) ia(i,j)=16
            endif
            if(fitsdata(i,j).gt.z2) ia(i,j)=ncol+15
 11      continue
C     write(6,*) (IA(i,k),k=1,ny)
 10   continue

C     set up pgplot window
      call pgscr(0,0.0,0.3,0.2)
      call pgsvp(0.05,0.95,0.05,0.95)
      call pgwnad(0.0,1.0,0.0,1.0)

      open(unit=31,file="/iraf/iraf/unix/sun/heat.lut",
     .     status='old')

      read(31,*) dumi

C     Setup colour-map
      do 300 i=1,ncol
         read(31,*) r,g,b
         read(31,*) dumr,dumr,dumr
         read(31,*) dumr,dumr,dumr
         read(31,*) dumr,dumr,dumr
c         R = REAL(I-1)/REAL(NCOL-1)*0.8 + 0.2
c         G = MAX(0.0, 2.0*REAL(I-1-NCOL/2)/REAL(NCOL-1))
c         B = 0.2 + 0.4*REAL(NCOL-I)/REAL(NCOL)
         CALL PGSCR(ncol-I+16, R, G, B)
 300  CONTINUE

      close(31)

      
C     Setup colour-map
c      do 300 i=1,ncol
c         R = REAL(I-1)/REAL(NCOL-1)*0.8 + 0.2
c         G = MAX(0.0, 2.0*REAL(I-1-NCOL/2)/REAL(NCOL-1))
c         B = 0.2 + 0.4*REAL(NCOL-I)/REAL(NCOL)
c         CALL PGSCR(I+15, R, G, B)
c 300  CONTINUE

      call pgpixl(ia,xmax,ymax,1,nx,1,ny,
     .     0.0,real(nx)/real(max(nx,ny)),
     .     0.0,real(ny)/real(max(nx,ny)))

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine autoscale(npt,pts,z1,z2,dcut)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Selects display range of data
      implicit none
      integer npt,hin,lown,nc
      real pts(npt),z1,z2,dcut,ave,adev,sdev,var,skew,curt,modecal,mode,
     .     stdev

      if(npt.eq.0) then
         z1=0.
         z2=0.
         goto 999
      endif

c      write(6,*) "dcut",dcut
      call sort(npt,pts)

      nc=int(real(npt)*(1.0-dcut)/2.0)
c      write(6,*) "min,max",pts(1),pts(npt)
      lown=nc+1
      hin=npt-nc
c      write(6,*) "hin",hin,npt

      mode=modecal(npt,pts,pts(lown),pts(hin))
      sdev=stdev(npt,pts,mode)
c      write(6,*) "sdev",sdev,mode

      z1=mode-1.5*sdev
c      write(6,*) "z1",z1
      z2=pts(hin)
c      write(6,*) "z2",z2
c      write(6,*) "z1,z2",z1,z2

 999  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE sort(n,arr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n,M,NSTACK
      real arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      real a,temp
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then
         do 12 j=l+1,ir
            a=arr(j)
            do 11 i=j-1,l,-1
               if(arr(i).le.a)goto 2
               arr(i+1)=arr(i)
 11         enddo
            i=l-1
 2          arr(i+1)=a
 12      enddo
         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         temp=arr(k)
         arr(k)=arr(l+1)
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
         jstack=jstack+2
         if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function modecal(npt,pts,dmin,dmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbins,i,j
      parameter(nbins=20)
      integer bins(nbins),maxbin,maxi
      real pts(npt),dmin,dmax,bspace,mode,x(nbins),y(nbins),minx,maxx,
     .     miny,maxy

      if(npt.eq.0) then
         write(6,*) "modecal: no points!"
         pause
      endif

      do 10 i=1,nbins
         bins(i)=0
 10   continue
      maxbin=-99

      modecal=300.0

      do 20 j=1,npt
C     calculating the bspacing between each bin
         bspace=(dmax-dmin)/real(nbins)
C     calculating the bin number
         i=(pts(j)-dmin)/bspace
         if((i.gt.0).and.(i.le.nbins)) then
c            write(6,*) i
            bins(i)=bins(i)+1
            if(bins(i).gt.maxbin) then 
               maxbin=bins(i)
               maxi=i
               mode=real(i)*bspace+dmin
            endif
         endif
 20   continue

      modecal=mode

c      write(6,*) "mode,nbin",modecal,maxbin
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function findmed(n,pts)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calcuates median.  Sorts data and picks middle value.
      integer n,i
      real pts(n),med

      if(n.le.1) then
         write(6,*) "no points in findmed"
         findmed=0.
         goto 999
      endif

      call sort(n,pts)

      i=n
c 10   continue
c      if(i.eq.2) goto 20
c      if(pts(i).gt.160) then
c         i=i-1
c         goto 10
c      endif
 
 20   med=pts(i/2)

      findmed=med

 999  return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function stdev(npt,pts,mean)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i
      real pts(npt),mean,s,ep,adev,p,var,sdev,mean2

      if(npt.eq.0) then
         write(6,*) "no points in stdev"
         stdev=0.
         goto 999
      endif

      s=0.
c      do 11 i=1,npt
c         s=s+pts(i)
c 11   continue
c      mean=s/npt

      ep=0.
      var=0.
      do 10 i=1,npt
         s=pts(i)-mean
         ep=ep+s
         p=s*s
         var=var+p
 10   continue
      var=(var-ep**2/npt)/(npt-1)
      sdev=sqrt(var)

      stdev=sdev

 999  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readst(header,nrecmax,nkeys,stn,stx,sty,rstx,rsty)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in info about secondary image subraster positions
      implicit none
      integer stmax
      parameter(stmax=20)
      integer nrecmax,stn,stx(stmax,2),sty(stmax,2),i,j,dumi,nkeys,rst,
     .     rstx(stmax,2),rsty(stmax,2)
      character*80 header(nrecmax),record,headint1(stmax),
     .     headint2(stmax),headint3(stmax),headint4(stmax),record2,
     .     rhead1(stmax),rhead2(stmax),rhead3(stmax),rhead4(stmax)
      character dumc
      character*2 dumc2

      headint1(1)="NUM_SI  "
      headint1(2)="RST00_ID"
      
      do 5 i=1,stmax
         rstx(i,1)=0
         rstx(i,2)=0
         rsty(i,1)=0
         rsty(i,2)=0
 5    continue

      rst=0
      do 10 i=1,nkeys
         record=header(i)
C         write(6,*) "nkeys",i,nkeys,record
         if(headint1(1)(1:8).eq.record(1:8)) read(record(10:80),*) stn
         if(rst.eq.0) then
            record2="                 "
            if(headint1(2)(1:3).eq.record(1:3).and.
     .           (headint1(2)(6:8).eq.record(6:8))) then
               read(record(12:20),*,err=902,end=902) record2
 902           continue
               if(record2.eq."Secondary") then
                  read(record(4:5),*) rst
               endif
            endif
         endif
 10   continue

      do 20 i=1,stn
         dumi=i-1
         if(dumi.lt.10) then
            write(dumc,500) dumi
 500        format(i1)
            headint1(i)="SI000"//dumc//"_X"
            headint2(i)="SI000"//dumc//"_Y"
            headint3(i)="SI000"//dumc//"_W"
            headint4(i)="SI000"//dumc//"_H"
         elseif((dumi.ge.10).and.(dumi.le.99)) then
            write(dumc2,501) dumi
 501        format(i2)
            headint1(i)="SI00"//dumc2//"_X"
            headint2(i)="SI00"//dumc2//"_Y"
            headint3(i)="SI00"//dumc2//"_W"
            headint4(i)="SI00"//dumc2//"_H"
         endif
         dumi=rst+i-1
         if(dumi.lt.10) then
            write(dumc,500) rst+i-1
            rhead1(i)="RST0"//dumc//"_X1"
            rhead2(i)="RST0"//dumc//"_Y1"
            rhead3(i)="RST0"//dumc//"_X2"
            rhead4(i)="RST0"//dumc//"_Y2"
         elseif((dumi.ge.10).and.(dumi.le.99)) then
            write(dumc2,501) rst+i-1
            rhead1(i)="RST"//dumc2//"_X1"
            rhead2(i)="RST"//dumc2//"_Y1"
            rhead3(i)="RST"//dumc2//"_X2"
            rhead4(i)="RST"//dumc2//"_Y2"
         endif
 20   continue

      do 30 i=1,nkeys
         record=header(i)
         do 31 j=1,stn
            record2=headint1(j)
            if(record2(1:8).eq.record(1:8)) then
               read(record(10:80),*) stx(j,1)
            endif
            record2=headint2(j)
            if(record2(1:8).eq.record(1:8)) then
               read(record(10:80),*) sty(j,1)
            endif
            record2=headint3(j)
            if(record2(1:8).eq.record(1:8))then
               read(record(10:80),*) dumi
               stx(j,2)=stx(j,1)+dumi-1
            endif
            record2=headint4(j)
            if(record2(1:8).eq.record(1:8))then
               read(record(10:80),*) dumi
               sty(j,2)=sty(j,1)+dumi-1
            endif
            record2=rhead1(j)
            if(record2(1:8).eq.record(1:8)) then 
               read(record(10:80),*) rstx(j,1)
c               write(6,*) rstx(j,1)
            endif
            record2=rhead2(j)
            if(record2(1:8).eq.record(1:8)) then
               read(record(10:80),*) rsty(j,1)
c               write(6,*) rsty(j,1)
            endif
            record2=rhead3(j)
            if(record2(1:8).eq.record(1:8)) then 
               read(record(10:80),*) rstx(j,2)
c               write(6,*) rstx(j,2)
            endif
            record2=rhead4(j)
            if(record2(1:8).eq.record(1:8)) then
               read(record(10:80),*) rsty(j,2)
c               write(6,*) rsty(j,2)
            endif
 31      continue
 30   continue

      do 40 i=1,stn
         if(rstx(i,1).eq.rstx(i,2)) then
            rstx(i,1)=stx(i,1)
            rstx(i,2)=stx(i,2)
         endif
         if(rsty(i,1).eq.rsty(i,2)) then
            rsty(i,1)=sty(i,1)
            rsty(i,2)=sty(i,2)
         endif
 40   continue

c      do 32 j=1,stn
c         write(6,*) "rst:",rstx(j,1),rstx(j,2),rsty(j,1),rsty(j,2)
c 32   continue
c      read(5,*)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readti(header,nrecmax,nkeys,x,y)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in where the science image is.
      implicit none
      integer nrecmax,x(2),y(2),tmp(4),i,nkeys
      character*80 header(nrecmax),headint1,headint2,headint3,headint4,
     .     record

      headint1="TI_X"
      headint2="TI_Y"
      headint3="TI_W"
      headint4="TI_H"
      
      do 10 i=1,nkeys
         record=header(i)
         if(record(1:4).eq.headint1(1:4)) read(record(10:80),*) tmp(1)
         if(record(1:4).eq.headint2(1:4)) read(record(10:80),*) tmp(2)
         if(record(1:4).eq.headint3(1:4)) read(record(10:80),*) tmp(3)
         if(record(1:4).eq.headint4(1:4)) read(record(10:80),*) tmp(4)
 10   continue

      x(1)=tmp(1)
      x(2)=tmp(1)+tmp(3)-1
      y(1)=tmp(2)
      y(2)=tmp(2)+tmp(4)-1

c      write(6,*) x(1),y(1),x(2),y(2)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine flatmap(flat,dpix,flatfile)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,j,xmax,ymax
      parameter(xmax=600,ymax=600)
      real flat(xmax,ymax),dpix(xmax,ymax),dt,ft
      character*80 flatfile


      do 1 i=1,xmax
         do 2 j=1,ymax
            flat(i,j)=1.0
            dpix(i,j)=0.0
 2       continue
 1    continue

      open(unit=15,file=flatfile,status='old',err=901)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This has been changed to REAL X,Y pixel system for flatfield map
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     btw. the smaller xmax and ymax are, the faster the code is..
C     sucks because you have to initizalize some variables to zero
C     and it is really slow.
CB PUT THIS CODE BACK IN
 40   read(15,*,end=50) i,j,dt,ft
      if((i.lt.xmax).and.(j.lt.ymax)) then
         dpix(i,j)=dt
         flat(i,j)=ft
      endif
      goto 40 
 50   continue
      close(15)

C     procyon
c      dpix(4,28)=3.37
c      dpix(482,118)=-2000.


 502  format(2(I3,1X),2(F9.4,1X))

      goto 999
 901  write(6,*) "WARNING: No flat/dark pixel map: ",flatfile
      goto 999
 999  return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fabrymap(flat,dpix,flatfile)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,j,xmax,ymax
      parameter(xmax=600,ymax=600)
      real flat(xmax,ymax),dpix(xmax,ymax),dt,ft
      character*80 flatfile


      do 1 i=1,xmax
         do 2 j=1,ymax
            flat(i,j)=1.0
            dpix(i,j)=0.0
 2       continue
 1    continue

      open(unit=15,file=flatfile,status='old',err=901)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This has been changed to REAL X,Y pixel system for flatfield map
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     btw. the smaller xmax and ymax are, the faster the code is..
C     sucks because you have to initizalize some variables to zero
C     and it is really slow.
C     update: this is now only done once at the start of the program.
C      so the slowdown on performance is not an issue anymore.
CB PUT THIS CODE BACK IN
 40   read(15,502,end=50) i,j,dt,ft
      if((i.lt.xmax).and.(j.lt.ymax)) then
         dpix(i,j)=dt
         flat(i,j)=ft
c         if(flat(i,j).lt.0.01) flat(i,j)=flat(i,j)+1.0
c         if(dpix(i,j).gt.2300.0) dpix(i,j)=dpix(i,j)-2340.0
      endif
      goto 40 
 50   continue
      close(15)

C     procyon
c      dpix(4,28)=3.37
c      dpix(482,118)=-2000.


 502  format(2(I3,1X),2(F9.4,1X))

      goto 999
 901  write(6,*) "WARNING: No Fabry flat field information"
      goto 999
 999  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bdcor(header,nrecmax,nkeys,fitsdata,naxes,datamin,
     .     datamax,sat,stmax,stn,stx,sty,rstx,rsty,flat,dpix,sflag,
     .      nstack)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer xmax,ymax,naxes(2),i,j,nrecmax,nkeys,stmax,stn,k,xc,yc
      parameter(xmax=600,ymax=600)
      integer stx(stmax,2),sty(stmax,2),rstx(stmax,2),rsty(stmax,2),
     .     sflag(stmax),nstack
      real fitsdata(xmax,ymax),datamin,datamax,bias,dark,readbias,
     .     readdark,sat(stmax),flat(xmax,ymax),zcor,dpix(xmax,ymax),
     .     dt,ft
      character*80 header(nrecmax)

C     This is the correction info from the header.
      bias=readbias(header,nrecmax,nkeys)
      dark=readdark(header,nrecmax,nkeys,nstack)
C     M67
c      zcor=0.0
      zcor=dark

C     correct for global dark/bias.
      datamin= 1.0e30
      datamax=-1.0e30
      do 10 i=1,naxes(1)
         do 20 j=1,naxes(2)
            if(fitsdata(i,j).le.0.0) then
               fitsdata(i,j)=-999.9
            else
               fitsdata(i,j)=fitsdata(i,j)-zcor
c               fitsdata(i,j)=((fitsdata(i,j)-zcor)-dpix(i,j))/flat(i,j)
            endif
 20      continue
 10   continue

c      datamin=0.0-zcor
c      datamax=2.0**16-zcor
      do 30 i=1,stn
         sat(i)=5000000.0!1.0*(real(nstack)*2.0**14-zcor)
         sflag(i)=0
 30   continue
c      write(6,*) "zcor,datamin,datamax",zcor,datamin,datamax

      do 41 i=1,naxes(1)
         do 42 j=1,naxes(2)
            do 43 k=1,stn
               if((i.ge.stx(k,1)).and.(i.le.stx(k,2)).and.
     .              (j.ge.sty(k,1)).and.(j.le.sty(k,2))) then
                  xc=rstx(k,1)+i-stx(k,1)
                  yc=rsty(k,1)+j-sty(k,1)
c                  write(6,*) rstx(k,1),stx(k,1),rsty(k,1),sty(k,1)
                  if((xc.lt.xmax).and.(yc.lt.ymax)) then
c                     write(6,500) k,xc,yc,fitsdata(i,j),dpix(xc,yc),
c     .                    flat(xc,yc)
 500                 format(I2,1X,I4,1X,I4,3(1X,F10.4))
                     fitsdata(i,j)=fitsdata(i,j)*flat(xc,yc)+dpix(xc,yc)
c                     fitsdata(i,j)=fitsdata(i,j)+dpix(xc,yc)
c                     write(6,*) i,j,flat(xc,yc),dpix(xc,yc)
                  endif
               endif
 43         continue
            datamin=min(fitsdata(i,j),datamin)
            datamax=max(fitsdata(i,j),datamax)
c            read(5,*)
 42      continue
 41   continue


C     for 10 Aql
c      do 50 i=1,20
c         do 51 j=81,83
c            fitsdata(i,j)=-99.9e10 !bad data
c 51      continue
c 50   continue
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine closfits(unitfits,status)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     closes fits file and releases unit number
      implicit none
      integer unitfits,status

      call ftclos(unitfits,status)
      call ftfiou(unitfits,status)

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine openfits(readwrite,unitfits,filename,status)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     opens up a fits file and return unitnumber
C     readwrite sets read/write privilages 0=readonly
C     unitfits is the unit number assigned to the fits file
C     filename is the name of the fits file
C     status is required by cfitsio
      implicit none
      integer status,unitfits,readwrite,dumi
      character*80 filename

C     cfitsio wants this initialized
      status=0

C     gets an unused unit number to open fits file
      call ftgiou(unitfits,status)

C     setting to zero makes fits file readwrite
      readwrite=0

C     open this fits file
      call ftopen(unitfits,filename,readwrite,dumi,status)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function readbias(header,nrecmax,nkeys)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in bias info from headers and calculates a weighted mean
      implicit none
      integer nrecmax,nkeys,nbias,i,j,biasmax
      parameter(biasmax=9)
      integer bweight(biasmax)
      real bias(biasmax),weight
      character dumc
      character*80 header(nrecmax),record,headint1,headint2

      do 1 i=1,biasmax
         bias(i)=0.0
         bweight(i)=0.0
 1    continue


      headint1="NUM_BIAS"
      do 10 i=1,nkeys
         record=header(i)
         if(record(1:8).eq.headint1(1:8))then
C            write(6,*) record(10:80)
            read(record(10:80),*) nbias
         endif
 10   continue

      if(nbias.gt.biasmax) write(6,*) "nbias is greater than biasmax"

      do 20 i=1,nbias
         write(dumc,500) i-1
 500     format(I1)
         headint1="BIAS000"//dumc
         headint2="BSNP000"//dumc
C         write(6,*) headint1(1:8)," ",headint2(1:8)
         do 21 j=1,nkeys
            record=header(j)
            if(record(1:8).eq.headint1(1:8)) then
               read(record(10:80),*) bias(i)
c               write(6,*) record(10:80)
            endif
            if(record(1:8).eq.headint2(1:8)) then
               read(record(10:80),*) bweight(i)
c               write(6,*) record(10:80)
            endif
 21      continue
 20   continue

c      write(6,*) nbias,(bias(i),bweight(i),i=1,nbias)

      readbias=0.0
      weight=0.0
      do 30 i=1,nbias
         readbias=readbias+bias(i)*real(bweight(i))
         weight=weight+real(bweight(i))
 30   continue

      readbias=readbias/weight
      if(nbias.eq.0)readbias=0.
c      write(6,*) readbias

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function readdark(header,nrecmax,nkeys,nstack)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in dark current information and calcs weighted mean
      implicit none
      integer nrecmax,nkeys,ndark,i,j,darkmax
      parameter(darkmax=20)
      integer dweight(darkmax),nstack
      real dark(darkmax),weight
      character dumc
      character*80 header(nrecmax),record,headint1,headint2

      do 1 i=1,darkmax
         dark(i)=0.0
         dweight(i)=0.0
 1    continue


      headint1="NUM_DARK"
      do 10 i=1,nkeys
         record=header(i)
         if(record(1:8).eq.headint1(1:8))then
C            write(6,*) record(10:80)
            read(record(10:80),*) ndark
         endif
 10   continue

      if(ndark.gt.darkmax) write(6,*) "ndark is greater than darkmax"

      do 20 i=1,ndark
         write(dumc,500) i-1
 500     format(I1)
         headint1="DARK000"//dumc
         headint2="DKNP000"//dumc
C         write(6,*) headint1(1:8)," ",headint2(1:8)
         do 21 j=1,nkeys
            record=header(j)
            if(record(1:8).eq.headint1(1:8)) then
               read(record(10:80),*) dark(i)
c               write(6,*) record(10:80)
            endif
            if(record(1:8).eq.headint2(1:8)) then
               read(record(10:80),*) dweight(i)
c               write(6,*) record(10:80)
            endif
 21      continue
 20   continue

c      write(6,*) ndark,(dark(i),dweight(i),i=1,ndark)

      readdark=0.0
      weight=0.0
      do 30 i=1,ndark
C     The if statement is for M67 data only
c         if(i.eq.3) then
            readdark=readdark+dark(i)*real(dweight(i))
            weight=weight+real(dweight(i))
c         endif
 30   continue

      readdark=readdark/weight
      if(ndark.eq.0) readdark=real(nstack)*480.0  !guess for dark
c      write(6,*) readdark
C     is is for procyon data
c      readdark=dark(4)

      return
      end   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mrqmin(x,y,z,sig,ndata,a,ia,ma,covar,alpha,nca, 
     *     chisq,alamda,ifun) 
ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER ma,nca,ndata,ia(ma),MMAX,ifun 
      REAL alamda,chisq,funcs2,a(ma),alpha(nca,nca),covar(nca,nca), 
     *     sig(ndata),x(ndata),y(ndata),z(ndata) 
      PARAMETER (MMAX=200) 
      INTEGER j,k,l,mfit 
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX) 
      SAVE ochisq,atry,beta,da,mfit 
      if(alamda.lt.0.)then 
         mfit=0 
         do 11 j=1,ma 
            if (ia(j).ne.0) mfit=mfit+1 
 11      enddo  
         alamda=0.001 
C     mod line
         call mrqcof(x,y,z,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,ifun)
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
C     mod line
      call mrqcof(x,y,z,sig,ndata,atry,ia,ma,covar,da,nca,chisq,ifun)
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
      SUBROUTINE mrqcof(x,y,z,sig,ndata,a,ia,ma,alpha,beta,nalp, 
     *     chisq,ifun) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER ma,nalp,ndata,ia(ma),MMAX,ifun 
      REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata), 
     *     y(ndata),z(ndata)
      PARAMETER (MMAX=200) 
      INTEGER mfit,i,j,k,l,m 
      REAL dz,sig2i,wt,zmod,dyda(MMAX) 
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
         call funcs2(x(i),y(i),a,zmod,dyda,ma,ifun)!mod 
         sig2i=1./(sig(i)*sig(i)) 
         dz=z(i)-zmod 
         j=0 
         do 15 l=1,ma 
            if(ia(l).ne.0) then 
               j=j+1 
               wt=dyda(l)*sig2i 
c               write(6,*) "crap1",wt,dyda(l),sig2i
               k=0 
               do 14 m=1,l 
                  if(ia(m).ne.0) then 
                     k=k+1 
                     alpha(j,k)=alpha(j,k)+wt*dyda(m)
c                     write(6,*) "crap2",alpha(j,j),wt,dyda(m)
                  endif 
 14            enddo 
               beta(j)=beta(j)+dz*wt 
            endif 
 15      enddo 
         chisq=chisq+dz*dz*sig2i  
 16   enddo
      do 18 j=2,mfit 
         do 17 k=1,j-1 
            alpha(k,j)=alpha(j,k) 
 17      enddo
 18   enddo
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER isign,nn1,nn2,nn3 
      COMPLEX data(nn1/2,nn2,nn3),speq(nn2,nn3) 
C     USES fourn 
      INTEGER i1,i2,i3,j1,j2,j3,nn(3) 
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp 
      COMPLEX c1,c2,h1,h2,w 
      c1=cmplx(0.5,0.0) 
      c2=cmplx(0.0,-0.5*isign) 
      theta=6.28318530717959d0/dble(isign*nn1) 
      wpr=-2.0d0*sin(0.5d0*theta)**2 
      wpi=sin(theta) 
      nn(1)=nn1/2 
      nn(2)=nn2 
      nn(3)=nn3 
      if(isign.eq.1)then 
         call fourn(data,nn,3,isign) 
         do 12 i3=1,nn3 
            do 11 i2=1,nn2 
               speq(i2,i3)=data(1,i2,i3) 
 11         continue
 12      continue
      endif
      do 15 i3=1,nn3 
         j3=1 
         if (i3.ne.1) j3=nn3-i3+2 
         wr=1.0d0 
         wi=0.0d0 
         do 14 i1=1,nn1/4+1 
            j1=nn1/2-i1+2 
            do 13 i2=1,nn2 
               j2=1 
               if (i2.ne.1) j2=nn2-i2+2 
               if(i1.eq.1)then 
                  h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3))) 
                  h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
                  data(1,i2,i3)=h1+h2 
                  speq(j2,j3)=conjg(h1-h2) 
               else 
                  h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3))) 
                  h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3))) 
                  data(i1,i2,i3)=h1+w*h2 
                  data(j1,j2,j3)=conjg(h1-w*h2) 
               endif 
 13         continue
            wtemp=wr 
            wr=wr*wpr-wi*wpi+wr 
            wi=wi*wpr+wtemp*wpi+wi 
            w=cmplx(sngl(wr),sngl(wi)) 
 14      continue 
 15   continue
      if(isign.eq.-1)then 
         call fourn(data,nn,3,isign) 
      endif 
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE fourn(data,nn,ndim,isign) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER isign,ndim,nn(ndim) 
      REAL data(*) 
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2, 
     *     ip3,k1,k2,n,nprev,nrem,ntot 
      REAL tempi,tempr 
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp 
      ntot=1 
      do 11 idim=1,ndim 
         ntot=ntot*nn(idim) 
 11   continue
      nprev=1 
      do 18 idim=1,ndim 
         n=nn(idim) 
         nrem=ntot/(n*nprev) 
         ip1=2*nprev 
         ip2=ip1*n 
         ip3=ip2*nrem 
         i2rev=1 
         do 14 i2=1,ip2,ip1 
            if(i2.lt.i2rev)then 
               do 13 i1=i2,i2+ip1-2,2 
                  do 12 i3=i1,ip3,ip2 
                     i3rev=i2rev+i3-i2 
                     tempr=data(i3) 
                     tempi=data(i3+1) 
                     data(i3)=data(i3rev) 
                     data(i3+1)=data(i3rev+1) 
                     data(i3rev)=tempr 
                     data(i3rev+1)=tempi
 12               continue
 13            continue
            endif 
            ibit=ip2/2
 1          if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then 
               i2rev=i2rev-ibit 
               ibit=ibit/2 
               goto 1 
            endif 
            i2rev=i2rev+ibit 
 14      continue
         ifp1=ip1 
 2       if(ifp1.lt.ip2)then 
            ifp2=2*ifp1 
            theta=isign*6.28318530717959d0/(ifp2/ip1) 
            wpr=-2.d0*sin(0.5d0*theta)**2 
            wpi=sin(theta) 
            wr=1.d0 
            wi=0.d0 
            do 17 i3=1,ifp1,ip1 
               do 16 i1=i3,i3+ip1-2,2 
                  do 15 i2=i1,ip3,ifp2 
                     k1=i2 
                     k2=k1+ifp1 
                     tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1) 
                     tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2) 
                     data(k2)=data(k1)-tempr 
                     data(k2+1)=data(k1+1)-tempi
                     data(k1)=data(k1)+tempr 
                     data(k1+1)=data(k1+1)+tempi 
 15               continue
 16            continue
               wtemp=wr 
               wr=wr*wpr-wi*wpi+wr 
               wi=wi*wpr+wtemp*wpi+wi 
 17         continue
            ifp1=ifp2 
            goto 2 
         endif 
         nprev=n*nprev 
 18   continue
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE gaussj(a,n,np,b,m,mp) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
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
c                     pause 'singular matrix in gaussj' 
                     return
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
         if (a(icol,icol).eq.0.) then!pause 'singular matrix in gaussj' 
            return
         endif
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
      SUBROUTINE SVDFIT(X,Y,Z,SIG,NDATA,A,MA,U,V,W,MP,NP,CHISQ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
C     IN: REAL X(NDATA) X data values for fiting.
C     IN: REAL Y(NDATA) Y data values for fiting.
C     IN: REAL SIG(NDATA) standard deviations for data.
C     IN: INTEGER NDATA amount of data.
C     OUT: REAL A(MA) coefficients of fitting function.
C     IN: INTEGER MA order of fit
C     OUT: REAL U(MP,NP) work space for matrix solution
C     OUT: REAL V(NP,NP) work space for matrix solution, also used in
C                        SVDVAR to calculate covariance matrix
C     OUT: REAL W(NP)    work space for matrix solution, also used in
C                        SVDVAR to calculate covariance matrix
C     IN: INTEGER MP physical size of row index in U that must equal 
C                    physical size of U in routines that call this 
C                    routine MP must be greater that or equal to NDATA.
C     IN: INTEGER NP physical size for U,V,W that must equal physical 
C                    size of U,V,W in routines that call this routine.
C                    NP must be greater than or equal to MA.
C     OUT: REAL CHISQ value to determine whether we have a good fit.
C
C     Parameter NMAX should equal S from main and represent maximum 
C     amout of data, MMAX represent maximum order of fit and again 
C     should equal value from main.
C
      PARAMETER(NMAX=7000,MMAX=50,TOL=1.E-5)
      DIMENSION X(NDATA),Y(NDATA),Z(NDATA),SIG(NDATA),A(MA),V(NP,NP),
     *     U(MP,NP),W(NP),B(NMAX),AFUNC(MMAX)
      DO 12 I=1,NDATA
         CALL FUNCS(X(I),Y(I),AFUNC,MA)  !mod add y
         TMP=1./SIG(I)
         DO 11 J=1,MA
            U(I,J)=AFUNC(J)*TMP
 11      CONTINUE
         B(I)=Z(I)*TMP  !mod y->z
 12   CONTINUE
      CALL SVDCMP(U,NDATA,MA,MP,NP,W,V)
      WMAX=0.
      DO 13 J=1,MA
         IF(W(J).GT.WMAX)WMAX=W(J)
 13   CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,MA
         IF(W(J).LT.THRESH)W(J)=0.
 14   CONTINUE
      CALL SVBKSB(U,W,V,NDATA,MA,MP,NP,B,A)
      CHISQ=0.
      DO 16 I=1,NDATA
         CALL FUNCS(X(I),Y(I),AFUNC,MA)  !mod add y
         SUM=0.
         DO 15 J=1,MA
            SUM=SUM+A(J)*AFUNC(J)
 15      CONTINUE
         CHISQ=CHISQ+((Z(I)-SUM)/SIG(I))**2  !mod y-> z
 16   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SVDVAR(V,MA,NP,W,CVM,NCVM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
C     IN: REAL V(NP,NP) value returned from SVDFIT
C     IN: INTEGER MA order of fit
C     IN: INTEGER NP physical dimension of V and W as seen in calling
C                    routine
C     IN: REAL W(NP) value returned from SVDFIT
C     OUT: REAL CVM(NCVM,NCVM) covariance matrix generated from V and W
C     IN: INTEGER NCVM physical dimension of NCVM as seen in calling 
C         routine
C
      PARAMETER (MMAX=50)
      DIMENSION V(NP,NP),W(NP),CVM(NCVM,NCVM),WTI(MMAX)
      DO 11 I=1,MA
         WTI(I)=0.
         IF(W(I).NE.0.)  WTI(I)=1./(W(I)*W(I))
 11   CONTINUE
      DO 14 I=1,MA
         DO 13 J=1,I
            SUM=0.
            DO 12 K=1,MA
               SUM=SUM+V(I,K)*V(J,K)*WTI(K)
 12         CONTINUE
            CVM(I,J)=SUM
            CVM(J,I)=SUM
 13      CONTINUE
 14   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
      PARAMETER(NMAX=100)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      IF(M.LT.N) return
C     PAUSE 'You must augment A with extra zero rows.'
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
         L=I+1
         RV1(I)=SCALE*G
         G=0.0
         S=0.0
         SCALE=0.0
         IF (I.LE.M) THEN
            DO 11 K=I,M
               SCALE=SCALE+ABS(A(K,I))
 11         CONTINUE
            IF (SCALE.NE.0.0) THEN
               DO 12 K=I,M
                  A(K,I)=A(K,I)/SCALE
                  S=S+A(K,I)*A(K,I)
 12            CONTINUE
               F=A(I,I)
               G=-SIGN(SQRT(S),F)
               H=F*G-S
               A(I,I)=F-G
               IF (I.NE.N) THEN
                  DO 15 J=L,N
                     S=0.0
                     DO 13 K=I,M
                        S=S+A(K,I)*A(K,J)
 13                  CONTINUE
                     F=S/H
                     DO 14 K=I,M
                        A(K,J)=A(K,J)+F*A(K,I)
 14                  CONTINUE
 15               CONTINUE
               ENDIF
               DO 16 K=I,M
                  A(K,I)=SCALE*A(K,I)
 16            CONTINUE
            ENDIF
         ENDIF
         W(I)=SCALE*G
         G=0.0
         S=0.0
         SCALE=0.0
         IF ((I.LE.M).AND.(I.NE.N))THEN
            DO 17 K=L,N
               SCALE=SCALE+ABS(A(I,K))
 17         CONTINUE
            IF (SCALE.NE.0.0) THEN
               DO 18 K=L,N
                  A(I,K)=A(I,K)/SCALE
                  S=S+A(I,K)*A(I,K)
 18            CONTINUE
               F=A(I,L)
               G=-SIGN(SQRT(S),F)
               H=F*G-S
               A(I,L)=F-G
               DO 19 K=L,N
                  RV1(K)=A(I,K)/H
 19            CONTINUE
               IF (I.NE.M) THEN
                  DO 23 J=L,M
                     S=0.0
                     DO 21 K=L,N
                        S=S+A(J,K)*A(I,K)
 21                  CONTINUE
                     DO 22 K=L,N
                        A(J,K)=A(J,K)+S*RV1(K)
 22                  CONTINUE
 23               CONTINUE
               ENDIF
               DO 24 K=L,N
                  A(I,K)=SCALE*A(I,K)
 24            CONTINUE
            ENDIF
         ENDIF
         ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
 25   CONTINUE
      DO 32 I=N,1,-1
         IF(I.LT.N) THEN
            IF(G.NE.0.0) THEN
               DO 26 J=L,N
                  V(J,I)=(A(I,J)/A(I,L))/G
 26            CONTINUE
               DO 29 J=L,N
                  S=0.0
                  DO 27 K=L,N
                     S=S+A(I,K)*V(K,J)
 27               CONTINUE
                  DO 28 K=L,N
                     V(K,J)=V(K,J)+S*V(K,I)
 28               CONTINUE
 29            CONTINUE
            ENDIF
            DO 31 J=L,N
               V(I,J)=0.0
               V(J,I)=0.0
 31         CONTINUE
         ENDIF
         V(I,I)=1.0
         G=RV1(I)
         L=I
 32   CONTINUE
      DO 39 I=N,1,-1
         L=I+1
         G=W(I)
         IF (I.LT.N) THEN
            DO 33 J=L,N
               A(I,J)=0.0
 33         CONTINUE
         ENDIF
         IF (G.NE.0.0) THEN
            G=1.0/G
            IF (I.NE.N) THEN
               DO 36 J=L,N
                  S=0.0
                  DO 34 K=L,M
                     S=S+A(K,I)*A(K,J)
 34               CONTINUE
                  F=(S/A(I,I))*G
                  DO 35 K=I,M
                     A(K,J)=A(K,J)+F*A(K,I)
 35               CONTINUE
 36            CONTINUE
            ENDIF
            DO 37 J=I,M
               A(J,I)=A(J,I)*G
 37         CONTINUE
         ELSE
            DO 38 J=I,M
               A(J,I)=0.0
 38         CONTINUE
         ENDIF
         A(I,I)=A(I,I)+1.0
 39   CONTINUE
      DO 49 K=N,1,-1
         DO 48 ITS=1,30
            DO 41 L=K,1,-1
               NM=L-1
               IF ((ABS(RV1(L))+ANORM).EQ.ANORM) GO TO 2
               IF ((ABS(W(NM))+ANORM).EQ.ANORM) GO TO 1
 41         CONTINUE
 1          C=0.0
            S=1.0
            DO 43 I=L,K
               F=S*RV1(I)
               IF ((ABS(F)+ANORM).NE.ANORM) THEN
                  G=W(I)
                  H=SQRT(F*F+G*G)
                  W(I)=H
                  H=1.0/H
                  C=(G*H)
                  S=-(F*H)
                  DO 42 J=1,M
                     Y=A(J,NM)
                     Z=A(J,I)
                     A(J,NM)=(Y*C)+(Z*S)
                     A(J,I)=-(Y*S)+(Z*C)
 42               CONTINUE
               ENDIF
 43         CONTINUE
 2          Z=W(K)
            IF (L.EQ.K) THEN
               IF (Z.LT.0.0) THEN
                  W(K)=-Z
                  DO 44 J=1,N
                     V(J,K)=-V(J,K)
 44               CONTINUE
               ENDIF
               GO TO 3
            ENDIF
c            IF (ITS.EQ.30) PAUSE 'No covergence in 30 iterations '
            if(ITS.eq.30) return
            X=W(L)
            NM=K-1
            Y=W(NM)
            G=RV1(NM)
            H=RV1(K)
            F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
            G=SQRT(F*F+1.0)
            F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
            C=1.0
            S=1.0
            DO 47 J=L,NM
               I=J+1
               G=RV1(I)
               Y=W(I)
               H=S*G
               G=C*G
               Z=SQRT(F*F+H*H)
               RV1(J)=Z
               C=F/Z
               S=H/Z
               F= (X*C)+(G*S)
               G=-(X*S)+(G*C)
               H=Y*S
               Y=Y*C
               DO 45 NM=1,N
                  X=V(NM,J)
                  Z=V(NM,I)
                  V(NM,J)= (X*C)+(Z*S)
                  V(NM,I)=-(X*S)+(Z*C)
 45            CONTINUE
               Z=SQRT(F*F+H*H)
               W(J)=Z
               IF (Z.NE.0.0) THEN
                  Z=1.0/Z
                  C=F*Z
                  S=H*Z
               ENDIF
               F= (C*G)+(S*Y)
               X=-(S*G)+(C*Y)
               DO 46 NM=1,M
                  Y=A(NM,J)
                  Z=A(NM,I)
                  A(NM,J)= (Y*C)+(Z*S)
                  A(NM,I)=-(Y*S)+(Z*C)
 46            CONTINUE
 47         CONTINUE
            RV1(L)=0.0
            RV1(K)=F
            W(K)=X
 48      CONTINUE
 3       CONTINUE
 49   CONTINUE
      RETURN
      END    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
      PARAMETER(NMAX=100)
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      DO 12 J=1,N
         S=0.
         IF(W(J).NE.0.)THEN
            DO 11 I=1,M
               S=S+U(I,J)*B(I)
 11         CONTINUE
            S=S/W(J)
         ENDIF
         TMP(J)=S
 12   CONTINUE
      DO 14 J=1,N
         S=0.
         DO 13 JJ=1,N
            S=S+V(J,JJ)*TMP(JJ)
 13      CONTINUE
         X(J)=S
 14   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FUNCS(X,Y,P,NP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
      DIMENSION P(NP)
      P(1)=1.0
      DO 11 J=2,NP,2
         P(J)=X**(J/2)
         P(J+1)=Y**(J/2)
 11   CONTINUE
      RETURN
      END
