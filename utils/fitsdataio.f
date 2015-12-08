CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readheader(unitfits,status,header,nkeys)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in the fits header
C     unitsfits is the unit number assigned to the fits file
C     status is from openfits and is used by cfitsio
      implicit none
      integer nrecmax
      parameter(nrecmax=700)
      integer status,unitfits,nkeys,nspace,i
      character*80 record,header(nrecmax)

      nkeys=0

C     get number of headers in image
      call ftghsp(unitfits,nkeys,nspace,status)

      if (nkeys.gt.nrecmax) then
         write(6,*) "WARNING: nkeys is more than nrecmax!!!!"
         write(6,*) "nkeys, nrecmax", nkeys, nrecmax
        pause
      endif

C     read in each header and move it into the master list.
      do 10 i=1,nkeys
         call ftgrec(unitfits,i,record,status)
         header(i)=record
 10   continue

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readfits(unitfits,status,fitsdata,naxes,datamin,
     .     datamax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in fits data
C     unitfies is the opened fits file
C     status is required by cfitsio for error reporting
C     fitsdata contains the image data
C     naxes is the size of each axes
C     datamin and datamax are the min and max of the data
      implicit none
      integer unitfits,status,naxes(2),nfound,xmax,ymax,npixels,group,
     .     firstpix,nullval,nbuf,nbufmax,nbuffer,i,j
      parameter(xmax=600,ymax=600,nbufmax=100)
      real datamin,datamax,buffer(nbufmax),fitsdata(xmax,ymax)
      logical anynull

C     determine size of fits image
      call ftgknj(unitfits,'NAXIS',1,2,naxes,nfound,status)
C      write(6,*) naxes,nfound,status,naxes(1)*naxes(2)

C     Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound.ne.2)then
         write(6,*) 'READIMAGE failed to read the NAXISn keywords.'
         pause
      end if

C     set up the total number of pixels to read      
      npixels=naxes(1)*naxes(2)
C     this is stuff required by FITSIO to work.
      group=1
      firstpix=1
      nullval=-999
C     for finding the maximum and minimum data values
      datamin= 1.0E30
      datamax=-1.0E30
      nbuf=naxes(1)
      j=0

C     keep looping until all pixels have been read.
      do while (npixels.gt.0)
C     read in 1 column at a time
         nbuffer=min(nbuf,npixels)
         
         call ftgpve(unitfits,group,firstpix,nbuffer,nullval,buffer,
     .        anynull,status)

         j=j+1
C     find max and min values
         do 10 i=1,nbuffer
            datamin=min(datamin,buffer(i))
            datamax=max(datamax,buffer(i))
            fitsdata(i,j)=buffer(i)
c            write(6,*) i,j,fitsdata(i,j)
 10      continue

C     update pointers and counters

         npixels=npixels-nbuffer
         firstpix=firstpix+nbuffer
      enddo

c      write(6,*) "datacheck:",datamin,datamax,j


      return
      end