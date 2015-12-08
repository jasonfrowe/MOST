      program hdjcalc
C     This program calculates the HJD for MOST data.
      implicit none
      integer nmax,npt,nunit,eunit
      parameter(nmax=500000)
      real aflux(nmax),ftotal(nmax),fx(nmax),fxy(nmax),fy(nmax),
     .   itime(nmax),mfield(nmax),nap(nmax),rmavg,skystd(nmax),
     .   tboard(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),
     .   yc(nmax),time(nmax)
      double precision dtime(nmax),flux(nmax),mjd(nmax),hmjd(nmax),
     .   hjd(nmax),oldhjd(nmax)
      character*80 filename
      
      rmavg=0.0
      
      nunit=10 !unit for file I/o
      eunit=11 !unit for file i/O
      
      write(6,*) "Enter filename containing photometry: "
      read(5,500) filename
 500  format(a80) !maximum length of filename is 80 characters
 
C     Now we open the filename, if we encounter an error we exit via
C     line 901
      open(unit=nunit,file=filename,status='old',err=901)

      write(6,*) "Enter filename for photometry output: "
      read(5,500) filename
      
      open(unit=eunit,file=filename)

      npt=nmax !pass along maximum storage to subroutine
      call readmostjdcor(nunit,npt,dtime,oldhjd)
c      call readdatasimp(nunit,npt,dtime,flux) !read in the data

c      call readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,
c     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)


      write(6,*) "Number of Data points: ",npt     
 
C     Get modified Julian date
      call modjuldate(npt,dtime,mjd)
C     Get modified heliocentric Julian date
      call gethjd(npt,mjd,hmjd)
C     Get Heliocentric Julian date in MOST format
      call hjdcor(npt,hmjd,hjd)
C     Export our results
      call exportdatamostjdcor(eunit,npt,hjd,dtime,oldhjd)
c      call exportdatasimp(eunit,npt,hjd,flux)
c      call exportdata(npt,hjd,mag,merr,filename,sky,xc,yc,
c     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,rmavg)
 
      close(nunit) !don't forget to close the files
      close(eunit)
      goto 999 !exit program error free
      
C     Error Section
 901  write(6,*) "Cannot open ",filename
      goto 999!exit program
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportdatamostjdcor(eunit,npt,hjd,dtime,oldhjd) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer eunit,npt,i
      double precision hjd(npt),dtime(npt),oldhjd(npt)
      
      do 10 i=1,npt
         write(eunit,510) oldhjd(i),dtime(i),hjd(i)
c         write(6,*) oldhjd(i),dtime(i),hjd(i)
 10   continue
      
 510  format(3(F13.8,1X))
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportdatasimp(eunit,npt,hjd,flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer eunit,npt,i
      double precision hjd(npt),flux(npt)
      
      do 10 i=1,npt
        write(eunit,*) hjd(i),flux(i)
 10   continue
 
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine hjdcor(npt,hmjd,hjd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      double precision hmjd(npt),hjd(npt)
      
      do 10 i=1,npt
         hjd(i)=hmjd(i)+2400000.5 !undo Modified HJD
         hjd(i)=hjd(i)-2451545.0 !MOST offset
 10   continue
      
      return
      end 
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gethjd(npt,mjd,hmjd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,ra_h,ra_m,dec_d,dec_m,nmax,i
      parameter(nmax=500000)
      double precision ausec,degtorad,ra,dec,ra_s,dec_s,correction_secs,
     .   cel,mjd(npt),hmjd(npt),source_ra,source_dec,earth_ra(nmax),
     .   earth_dec(nmax),earth_x,earth_y,earth_z,source_x,source_y,
     .   source_z
      
      ausec=499.01265d0 !Time for light to travel 1 AU
      degtorad=0.01745329251 !Conversion from degrees to radians
      
      write(6,*) "All co-ordinates should be J2000"
      write(6,*) " "
 10   write(6,*) "Enter RA (hh mm ss.s)"
      read(5,*,err=10) ra_h,ra_m,ra_s
      write(6,500) "RA: ",ra_h,":",ra_m,":",ra_s 
 500  format(A4,I3,A1,I2,A1,F5.2)
 11   write(6,*) "Enter DEC (dd mm ss,s)"
      read(5,*,err=11) dec_d,dec_m,dec_s
      write(6,500) "DEC:",dec_d,":",dec_m,":",dec_s
      
C     calculate RA and DEC in decimal degree notation
      source_ra=dble(ra_h)
      source_ra=source_ra+dble(ra_m)/60.0d0
      source_ra=source_ra+ra_s/3600.0d0      
      source_ra=source_ra*15.0d0
      
      source_dec=dble(dec_d)
      source_dec=source_dec+dble(dec_m)/60.0d0
      source_dec=source_dec+dec_s/3600.0d0
     
C     Get RA and DEC of Earth using Astronomical Calculator book
      call heliocentric_ra_dec(npt,mjd,earth_ra,earth_dec)
      
      do 12 i=1,npt
C     Calculate the heliocentric co-ordinates as X, Y and Z terms
         cel =cos(earth_dec(i)*(degtorad))
         earth_x=cos(earth_ra(i)*degtorad)*cel
         earth_y=sin(earth_ra(i)*degtorad)*cel
         earth_z=sin(earth_dec(i)*degtorad)

C     Calculate the X,Y,Z co-ordinate of your source
         cel =cos(source_ra*degtorad)
         source_x= cos(source_ra*degtorad)*cel
         source_y= sin(source_ra*degtorad)*cel
         source_z= sin(source_dec*degtorad)

C     Calculate the correction in seconds for the light travel time
C     between the source and the Earth vectors in a heliocentric 
C     reference frame
         correction_secs=ausec*
     .      (earth_x*source_x + earth_y*source_y + earth_z*source_z)
     
         hmjd(i)=mjd(i)+(correction_secs/(24.0d0*3600.0d0))
 12   continue      
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine heliocentric_ra_dec(npt,mjd,earth_ra,earth_dec)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,i
      parameter(nmax=500000)
      double precision mjd(npt),earth_ra(npt),earth_dec(npt),
     .   eccentricity,ecliptic_long,perigee_long,Pi,deg_to_rad,
     .   tropical_year,obliquity,mjd_1980,mean_anomoly(nmax),
     .   days_from_1980(nmax),solar_longitude(nmax),number_of_deg(nmax),
     .   equation_of_centres(nmax),x,y,beta,number_of_rotations
      Pi=acos(-1.0d0) !get Pi
     
C     Some Constants      
      eccentricity=0.016718 !Ecc of Earth's orbit
      ecliptic_long=278.833540 ! The longitude of the ecliptic at
                               ! 1 Jan 1980 0:00 UT
      perigee_long=282.596403 !The longitude of perigee at 
                              ! 1 Jan 1980 0:00 UT
      deg_to_rad = Pi/180.0d0 ! degrees to radians conversion
      tropical_year = 365.24219572 ! The length of the tropical year 
                                   ! in days
      obliquity = 23.441884 ! The obliquity of the orbit
      mjd_1980 = 44238.0 ! The MJD on 1 Jan 1980 0:00 UT
      
C     Calculate the number of days since 1 Jan 1980
      do 10 i=1,npt
         days_from_1980(i) = mjd(i) - mjd_1980   
 10   continue

C     Calculate the number of degrees around in the orbit travelled in
C     this time
      do 11 i=1,npt
         number_of_deg(i)=(360.0d0/tropical_year)*days_from_1980(i)
 11   continue
 
C     Adjust so the number of degrees in between 0 and 360.
      do 12 i=1,npt
         if((number_of_deg(i).lt.0.0d0).or.
     .      (number_of_deg(i).gt.360.0d0))then
            number_of_rotations = number_of_deg(i)/360.0
            number_of_rotations = dble(int(number_of_rotations))
            number_of_deg(i)=number_of_deg(i)-
     .         number_of_rotations*360.0d0
         endif 
 12   continue

C     Calculate the mean anomoly
      do 13 i=1,npt
         mean_anomoly(i)=number_of_deg(i)-perigee_long+ecliptic_long
 13   continue

C     Calculate equation of centres
      do 14 i=1,npt
         equation_of_centres(i)=(360.0d0/Pi)*eccentricity*
     .      sin(mean_anomoly(i)*deg_to_rad)      
 14   continue

C     Calculate the solar longitude
      do 15 i=1,npt
         solar_longitude(i)=number_of_deg(i) + equation_of_centres(i) +
     .      ecliptic_long
         if(solar_longitude(i).gt.360.0d0) 
     .      solar_longitude(i)=solar_longitude(i)-360.0d0      
 15   continue

C     The ecliptic latitude is zero for the Sun
      beta = 0.0

C     Calculate RA and DEC of the Sun
      do 16 i=1,npt      
         earth_dec(i)=asin((sin(beta*deg_to_rad)*
     .      cos(obliquity*deg_to_rad))+
     .      (cos(beta*deg_to_rad)*sin(obliquity*deg_to_rad)*
     .      sin(solar_longitude(i)*deg_to_rad)) )
         earth_dec(i)=earth_dec(i)/deg_to_rad
         x=cos(solar_longitude(i)*deg_to_rad)
         y=(sin(solar_longitude(i)*deg_to_rad) *
     .      cos(obliquity*deg_to_rad)) - 
     .      (tan(beta*deg_to_rad)*sin(obliquity*deg_to_rad))
         earth_ra(i)=atan(y/x)
         earth_ra(i)=earth_ra(i)/deg_to_rad
C     Convert from geocentric to helicentric co-ordinates for the Earth
         earth_dec(i)=-1.0d0*earth_dec(i)
         earth_ra(i)=earth_ra(i)-12.0d0
         if(earth_ra(i) .lt. 0.0d0) earth_ra(i)=earth_ra(i)+24.0
 16   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine modjuldate(npt,dtime,mjd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      double precision dtime(npt),mjd(npt)
      
      do 10 i=1,npt
         mjd(i)=dtime(i)+2451545.0d0 !correct for MOST offset
         mjd(i)=mjd(i)-2400000.5d0 !offset to a known Earth pos.
 10   continue
 
      return
      end
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readmostjdcor(nunit,npt,dtime,oldhjd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,i,npt
      double precision dtime(npt),oldhjd(npt),hjdcor,mostoffset,sec2day

      mostoffset=2451545.0d0
      sec2day=60.0*60.0*24.0

      i=1
 10   read(nunit,*,end=11,err=901) oldhjd(i),hjdcor
         oldhjd(i)=oldhjd(i)-mostoffset
         dtime(i)=oldhjd(i)-hjdcor/sec2day
         i=i+1
         goto 10
 11   continue
      npt=i-1

      goto 999 !exit error free
C     Error section
 901  write(6,*) "Error on line ",i, "in file"
      pause
      goto 10  !try and read next line

 999  return
      end    
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readdatasimp(nunit,npt,dtime,flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,npt,i
      double precision dtime(npt),flux(npt)
      
      i=1
C     read in file contents. when EOF break from loop.
C     read errors result in termination of program
 10   read(nunit,*,end=11,err=901) dtime(i),flux(i)
         i=i+1
         goto 10
 11   continue
      npt=i-1  !store number of read points
           
      goto 999 !exit error free
C     Error section
 901  write(6,*) "Error on line ",i, "in file"
      pause
      goto 10  !try and read next line
 999  return
      end
    
