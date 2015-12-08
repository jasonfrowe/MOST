C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fitsreader(filename,stn,stx,sty,rstx,rsty,tix,tiy,
     .   nkeys,header,fitsdata,datamin,datamax,naxes)
      implicit none
C     this program tests out the new fits format for MOST data.  Will
C     incorporated into allMOST then become obs.
      integer unitfits,status,readwrite,dumi,i,clen,nrecmax,nkeys,j,
     .    xmax,ymax,naxes(2),k,naxesin(2),stmax,hud,stn,blocksize,
     .    hdutype,nhuds
      parameter(nrecmax=700,xmax=600,ymax=600,stmax=20)
      integer stx(stmax,2),sty(stmax,2),tix(2),tiy(2),rstx(stmax,2),
     .   rsty(stmax,2),tcheck
      real datamin,datamax,fitsdatain(xmax,ymax),fitsdata(xmax,ymax),
     .   dmin,dmax
      character*80 filename,header(nrecmax)
      character*8 imagetype,extname
      character ci1
      character*2 ci2

C     a check to see if the table information is in the fits files
      tcheck=0     
     
c      filename="1824_33240740_SDS2.fits"
C     if you have an error with obtaining the string length, try
C     using the g77 command lnblnk.  len_trim works with the intel
C     compiler.
      clen=len_trim(filename)

C     first we open the fits file
C     cfitsio wants this initialized
      status=0
      
C     setting to zero makes fits file readwrite
      readwrite=0

C     gets an unused unit number to open fits file
      call ftgiou(unitfits,status)
C     open the fits file
      call ftopen(unitfits,filename,readwrite,blocksize,status)
C     get the number of the extensions (huds) in the image
      call ftthdu(unitfits,nhuds,status)

C     the first HUD (i=1) is always the fabry lens.
      hud=1
C     stn counts the number of secondary images.
      stn=0
C     naxes keeps the size of the final image
      naxes(1)=0
      naxes(2)=0
C     see what the bounds of the data are, so initally set
C     to silly values
      datamin= 99.9e30
      datamax=-99.9e30 
C     open this fits file
      if(status.ne.0) then
         write(6,*) "whoops...",status
      endif
C     with fabry data, start this at hud=2
      do 10 hud=1,nhuds

         if(hud.gt.stmax) write(6,*) "overflow!! increase STMAX!!!"

C     if status is returned at 107 then we have reached the end
C     of the file

C        This command moves us to the next HUD
         call ftmahd(unitfits,hud,hdutype,status)

         if(status.eq.0) then
C     Read in all the header information
            call readheader(unitfits,status,header,nkeys)
            do 11 j=1,nkeys
C     The XTENSION keyword tells us if we have an image or table.
               if(header(j)(1:8).eq."XTENSION") then
                  read(header(j)(12:19),503) imagetype
 503              format(A8)
               endif
C     The extname is useful for identifying the correct table
               if(header(j)(1:8).eq."EXTNAME") then
                  read(header(j)(12:19),503) extname
c                  write(6,*) "extname:",extname
               endif
 11         continue
c            write(6,*) "imagetype:",imagetype,hud
CCCC        This line is for the new data set
            if(hud.eq.1)imagetype="IMAGE"           
 502        format(A80)
 
c            write(6,*) imagetype
            if(imagetype.eq."IMAGE") then
C     if we have an image extension, then read in the pixels and
C     add them to the master array (fitsdata)
               call readfits(unitfits,status,fitsdatain,naxesin,dmin,
     .         dmax)
C     this routine adds the read pixels (fitsdatain) to the master
C     array
               call asmfits(naxes,fitsdata,naxesin,fitsdatain,stn,stx,
     .            sty,datamin,datamax)
C     if we have a binary table, then check which one it is, then
C     call the routine that reads in the subraster locations from
C     the CCD
            elseif((imagetype.eq."BINTABLE").and.
     .            (extname.eq."SUBRSTRS")) then
               call rdfitbtab(unitfits,status,rstx,rsty)
               tcheck=1
            endif
         endif

C        check for any FITSIO error codes.
         if(status.ne.0) then
            write(6,*) "whoops check FITSIO error code:",status
            write(6,*) "filename:",filename(1:clen),"[",HUD,"]"
            pause
            status=0
         endif
   
 10   enddo

C     if the binary tables are missing, then we fill in the rstx/y 
C     raster information - or else the flatfielding maps will not work
C     and pieces of the subraster may get mapped to 0.
      tcheck=0
      if(tcheck.eq.0) then
         do 20 i=1,2
            do 21 j=1,stn
c              write(6,*) j,i,rstx(j,i)
c              write(6,*) j,i,rsty(j,i)
              rstx(j,i)=stx(j,i)
              rsty(j,i)=sty(j,i)
 21         continue
 20      continue
      endif      
c      read(5,*)

C     take on the primary image (not really necessary, but might as well
C     be complete
C     it is also useful
      status=0
C     the Fabry target is always supposed to be in HUD 1.
      hud=1
C     select the correct extension
      call ftmahd(unitfits,hud,hdutype,status)
C     read in header information
      call readheader(unitfits,status,header,nkeys)
C     read in the image data
      call readfits(unitfits,status,fitsdatain,naxesin,dmin,dmax)
C     close the fits file
      call ftclos(unitfits,status)
C     return unit number to unused list.
      call ftfiou(unitfits,status)
C     put the fabry lens info into main array      
      call priasm(naxes,fitsdata,naxesin,fitsdatain,tix,tiy)

C     I believe this format strings are obs.
 500  format(I1)
 501  format(I2)    

      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rdfitbtab(unitfits,status,rstx,rsty)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads binary tables in FITS file and get positional information
C     for subrasters.
      implicit none
      integer unitfits,status,nrows,ncols,i,datacode,repeat,width,dmax,
     .   stmax
      parameter(dmax=20,stmax=20)
      integer ival(dmax),j,sidx(dmax),stnc,rstx(stmax,2),rsty(stmax,2)
      real rval(dmax) 
      double precision dval(dmax)
      character*14 cval(dmax),ts
      logical anyf
    
C     get the number of rows and columns in the table.
      call FTGNRW(unitfits,nrows, status)
      call FTGNCL(unitfits,ncols, status)
      
c      write(6,*) "table:",nrows,ncols
            
C     loop over all the columns (have to do it this way with FITSIO
      do 10 i=1,ncols
C     see what kind of data string we are dealing with.
C     in the secondary target table, it should only be characters and
C     integers.
         call ftgtcl(unitfits,i,datacode,repeat,width,status)
c         write(6,*) "cc:",i,datacode,repeat,width,status
C     the first meanful column is in string format to label which
C     contain the subraster locations.
         if(datacode.eq.16) then
            do 11 j=1,nrows
C     read in the character string from the table.
               call ftgcvs(unitfits,i,j,1,1," ",ts,anyf,status)
               cval(j)=ts
C     label rows with secondary information.  If cval=1, then we care
C     about info in the next columns, or else we flag it 0 to skip.
               if(cval(j).eq."Secondary") then
                  sidx(j)=1
               else
                  sidx(j)=0
               endif
 11         continue
         elseif(datacode.eq.41) then
C     read in the integer values from the table.
            call ftgcvj(unitfits,i,1,1,nrows,0,ival,anyf,status)
c            write(6,*) (ival(j),j=1,nrows)
C     stnc counts the number of secondary images.  They should be
C     read out in the correct order and there is NO METHOD to check!
            stnc=0
            do 12 j=1,nrows
               if(sidx(j).eq.1)then
                  stnc=stnc+1
C     columns 3,4,5,6 contain the co-ordinates that map out the 
C     subraster on the CCD.
C     if you want the binning information, get it from column 7
                  if(i.eq.3) rstx(stnc,1)=ival(j)
                  if(i.eq.4) rsty(stnc,1)=ival(j)
                  if(i.eq.5) rstx(stnc,2)=ival(j)
                  if(i.eq.6) rsty(stnc,2)=ival(j)
               endif
 12         continue
         endif
           
 10   continue
 
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine priasm(naxes,fitsdata,naxesin,fitsdatain,tix,tiy)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this routine adds the Fabry pixels to the main array.
      implicit none
      integer naxes(2),naxesin(2),tix(2),tiy(2),xmax,ymax,i,j
      parameter(xmax=600,ymax=600)
      real fitsdata(xmax,ymax),fitsdatain(xmax,ymax)
      
C     set up these arrays to mark the relevent datasets in the array.
      tix(1)=naxes(1)+1
      tix(2)=naxes(1)+naxesin(1)
      tiy(1)=1
      tiy(2)=naxesin(2)

C     move the data into the correct positions in the array
      do 10 i=1,naxesin(1)
         do 11 j=1,naxesin(2)
            fitsdata(i+tix(1)-1,j+tiy(1)-1)=fitsdatain(i,j)
 11      continue
 10   continue
      
C     update the bounds of useful data.
      naxes(1)=max(naxes(1),tix(2))
      naxes(2)=max(naxes(2),tiy(2))
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine asmfits(naxes,fitsdata,naxesin,fitsdatain,stn,stx,sty,
     .   datamin,datamax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this routine maps the direct imaging data into the main data array
      implicit none
      integer naxes(2),naxesin(2),stn,stmax,xmax,ymax,i,j
      parameter(stmax=20,xmax=600,ymax=600)
      integer stx(stmax,2),sty(stmax,2)
      real fitsdata(xmax,ymax),fitsdatain(xmax,ymax),datamin,datamax
      
C     increase the number of secondary images
      stn=stn+1
C     update markers (basically recreating original SDS2 file)
      stx(stn,1)=1
      stx(stn,2)=naxesin(1)
      sty(stn,1)=naxes(2)+1
      sty(stn,2)=naxes(2)+naxesin(2)
     
C     loop over the data stream and move the data into the proper
C     positions.
      do 10 i=1,naxesin(1)
         do 11 j=1,naxesin(2)
            fitsdata(i+stx(stn,1)-1,j+sty(stn,1)-1)=fitsdatain(i,j)
            datamin=min(datamin,fitsdatain(i,j))
            datamax=max(datamax,fitsdatain(i,j))
 11      continue
 10   continue
c      write(0,*) "data:",datamin,datamax
      
C     update the arrays that map out useful data in the master array
      naxes(1)=max(naxes(1),stx(stn,2))
      naxes(2)=max(naxes(2),sty(stn,2))
      
      return
      end