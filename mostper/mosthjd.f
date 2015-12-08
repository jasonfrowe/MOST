      program MOSTHJD
      implicit none
      integer iday,imonth,iyear,ierr,ihour,iminute,irah,iram,idecd,
     .   idecm,oldhjd
      double precision jd,temp,sec,hjd,gethjd,ras,decs,mostoffset,
     .   hjdold,hjdcorr,sec2day
      
      sec2day=60.0d0*60.0d0*24.0d0
      jd=2454143.67587717d0
      irah=8
      iram=50
      ras=24.0d0
      idecd=11
      idecm=49
      decs=0.0d0
      mostoffset=2451545.0d0
      
 10   open(unit=10,file="hjd_dates.txt")
      
         read(10,*,end=11) hjdold,hjdcorr
         jd=hjdold-hjdcorr/sec2day
      
         call jd2day(jd,iday,imonth,iyear,ierr)
         temp=24.0d0*(jd-floor(jd))
         ihour=int(temp)+12
         temp=60.0d0*(temp-floor(temp))
         iminute=int(temp)
         sec=60.0d0*(temp-floor(temp))
         if(ihour.gt.24)then
            ihour=ihour-24
            iday=iday+1
         endif
      
c      write(6,*) "dd/mm/yyyy:",iday,imonth,iyear
c      write(6,*) "hh:mm:ss.s", ihour,iminute,sec
      
c      call getposition(iyear,imonth,iday,ihour,iminute,sec,3)
      
         hjd=gethjd(imonth,iday,iyear,ihour,iminute,sec,irah,iram,ras,
     .      idecd,idecm,decs)
         write(6,510) hjdold-mostoffset,hjd-mostoffset 
 510     format(3(F13.8,1X))
      goto 10
 11   continue
      close(10)
      
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function gethjd(imonth,iday,iyear,ihour,iminute,
     .   sec,irah,iram,ras,idecd,idecm,decs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer imonth,iday,iyear,ihour,iminute,iplanet,irah,iram,idecd,
     .   idecm
      double precision sec,ra,dec,rah,ram,decd,decm,ras,decs,ORAH,
     .   ORAM,ORAS,ODECD,ODECM,ODECS,ORA,ODEC,cel,earthx,earthy,earthz,
     .   pi,objectx,objecty,objectz,ausec,correction,D
      common /values/ ra,dec,D

      pi = 3.141592654d0
      gethjd=0.0d0

c  // Sun 
      iplanet=3
      call getposition(iyear,imonth,iday,ihour,iminute,sec,
     .   iplanet)
      rah=floor(ra)
      ram=floor((ra-floor(ra))*60.0d0)
      decd=floor(abs(dec))
      if (dec.lt.0) decd=-1.0d0*decd
      decm=floor((abs(dec)- floor(abs(dec)))*60.0d0)
c    ASDF=MM + "/" + DD + "/" + YY + "   " + HR + ":" + MN + ":" + SC;
c    ASDF=ASDF + "\nSun -- RA: "+rah+"h "+ram+"m   DEC: "+decd+"° "+decm+"'";
c  // Earth
      dec=-1.0d0*dec
      ra=ra+12.0d0
      if(ra.gt.24.0d0) ra=ra-24.0d0
      rah=floor(ra)
      ram=floor((ra - floor(ra))*60.0d0)
      decd=floor(abs(dec))
      if (dec.lt.0.0d0) decd=-1.0d0*decd
      decm=floor((abs(dec)- floor(abs(dec)))*60.0d0)
c    ASDF=ASDF + "\nEarth -- RA: "+rah+"h "+ram+"m   DEC: "+decd+"° "+decm+"'";

c //Object    
      ORAH=dble(irah)
      ORAM=dble(iram)
      ORAS=ras
      ORA = (ORAH + ORAM/60.0d0 + ORAS/3600.0d0)*15.0d0
      ODECD=dble(idecd)
      ODECM=dble(idecm)
      ODECS=decs
      ODEC = abs(ODECD) + ODECM/60.0d0 + ODECS/3600.0d0
      if (ODECD.lt.0) ODEC=-1.0d0*ODEC
c
c    //ORA=ra;
c    //ODEC=dec;
c
c    ASDF=ASDF + "\nObject -- RA: " + ORA + "°   DEC: " + ODEC + "° ";

c //Earth XYZ
      cel = cos(dec * pi/180.0d0)
      earthx = cos(ra * pi/12.0d0) * cel
      earthy = sin(ra * pi/12.0d0) * cel
      earthz = sin(dec * pi/180)
c    ASDF=ASDF + "\nEarth -- X,Y,Z: " + earthx + "," + earthy + "," + earthz;
c//Object XYZ
      cel = cos(ODEC * pi/180.0d0)
      objectx = cos(ORA * pi/180.0d0) * cel
      objecty = sin(ORA * pi/180.0d0) * cel
      objectz = sin(ODEC * pi/180.0d0)
c    ASDF=ASDF + "\nObject -- X,Y,Z: " + objectx + "," + objecty + "," + objectz;
     
c//Light Time (Minutes per AU)
      ausec=8.3168775d0
      correction = ausec*(earthx*objectx+earthy*objecty+earthz*objectz)
c      write(6,*) "Correction",correction
c    ASDF=ASDF + "\nLight Time: " + correction + " minutes";
      D=D+2451545
c      write(6,*) "JD:",D
c    ASDF=ASDF + "\nJulian Date: " + D;
      D=D+correction/(24.0d0*60.0d0)
c      write(6,*) "HJD:",D
      gethjd=D
c    ASDF=ASDF + "\nHeliocentric Julian Date: " + D;
c
c    form.coordinates.value=ASDF;
c}
c
    
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getposition(iyear,imonth,iday,ihour,iminute,sec,
     .   iplanet)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer iyear,imonth,iday,ihour,iminute,iplanet,p,q,xpix,ypix
      double precision sec,y,m,zDay,h,mins,secs,D,pi,Rads,el(54),ip,op,
     .   pp,ap,ep,lp,ie,oe,pe,ae,ee,le,Me1,B,e,delta,ve,re,xe,ye,ze,mp,
     .   vp,rp,xh,yh,zh,xg,yg,zg,ecl,xeq,yeq,zeq,ra,dec,rae,dece,raa
      common /values/ ra,dec,D

      pi = 3.141592654d0
      Rads = pi / 180.0d0

      y=dble(iyear)
      m=dble(imonth)
      zDay=dble(iday)
      h=dble(ihour)
      mins=dble(iminute)
      secs=sec
      h=h+mins/60.0d0+secs/3600.0d0
      p=iPlanet

      D=367.0d0*y-floor(7.0d0*(y+floor((m+9.0d0)/12.0d0))/4.0d0)+
     .   floor(275.0d0*m/9.0d0)+zDay-730531.5d0+h/24.0d0
c      write(6,*) "D:",D

c  // Mercury
      el(1) = (7.00487d0 - 0.000000178797d0 * D) * Rads
      el(2) = (48.33167d0 - 0.0000033942d0 * D) * Rads
      el(3) = (77.45645d0 + 0.00000436208d0 * D) * Rads
      el(4) = 0.38709893d0 + 1.80698d-11 * D
      el(5) = 0.20563069d0 + 0.000000000691855d0 * D
      el(6) = (Rads * (252.25084d0 + 4.092338796d0 * D))
c  // Venus
      el(7) = (3.39471d0 - 0.0000000217507d0 * D) * Rads
      el(8) = (76.68069d0 - 0.0000075815d0 * D) * Rads
      el(9) = (131.53298d0 - 0.000000827439d0 * D) * Rads
      el(10) = 0.72333199d0 + 2.51882d-11 * D
      el(11) = 0.00677323d0 - 0.00000000135195d0 * D
      el(12) = (Rads * (181.97973d0 + 1.602130474d0 * D))
c  // Earth
      el(13) = (0.00005d0 - 0.356985d-6 * D) * Rads
      el(14) = (-11.26064d0 - 0.00013863d0 * D) * Rads
      el(15) = (102.94719d0 + 0.00000911309d0 * D) * Rads
      el(16) = 1.00000011d0 - 1.36893d-12 * D
      el(17) = 0.01671022d0 - 0.104148d-8 * D
      el(18) = (Rads * (100.46435d0 + 0.985609101d0 * D))
c  // Mars
      el(19) = (1.85061d0 - 0.000000193703d0 * D) * Rads
      el(20) = (49.57854d0 - 0.0000077587d0 * D) * Rads
      el(21) = (336.04084d0 + 0.00001187d0 * D) * Rads
      el(22) = 1.52366231d0 - 0.000000001977d0 * D
      el(23) = 0.09341233d0 - 0.00000000325859d0 * D
      el(24) = (Rads * (355.45332d0 + 0.524033035d0 * D))
c  // Jupiter
      el(25) = (1.3053d0 - 0.0000000315613d0 * D) * Rads
      el(26) = (100.55615d0 + 0.00000925675d0 * D) * Rads
      el(27) = (14.75385d0 + 0.00000638779d0 * D) * Rads
      el(28) = 5.20336301d0 + 0.0000000166289d0 * D
      el(29) = 0.04839266d0 - 0.00000000352635d0 * D
      el(30) = (Rads * (34.40438d0 + 0.083086762d0 * D))
c  // Saturn
      el(31) = (2.48446d0 + 0.0000000464674d0 * D) * Rads
      el(32) = (113.71504d0 - 0.0000121d0 * D) * Rads
      el(33) = (92.43194d0 - 0.0000148216d0 * D) * Rads
      el(34) = 9.53707032d0 - 0.0000000825544d0 * D
      el(35) = 0.0541506d0 - 0.0000000100649d0 * D
      el(36) = (Rads * (49.94432d0 + 0.033470629d0 * D))
c  // Uranus
      el(37) = (0.76986d0 - 0.0000000158947d0 * D) * Rads
      el(38) = (74.22988d0 + 0.0000127873d0 * D) * Rads
      el(39) = (170.96424d0 + 0.0000099822d0 * D) * Rads
      el(40) = 19.19126393d0 + 0.0000000416222d0 * D
      el(41) = 0.04716771d0 - 0.00000000524298d0 * D
      el(42) = (Rads * (313.23218d0 + 0.011731294d0 * D))
c  // Neptune
      el(43) = (1.76917d0 - 0.0000000276827d0 * D) * Rads
      el(44) = (131.72169d0 - 0.0000011503d0 * D) * Rads
      el(45) = (44.97135d0 - 0.00000642201d0 * D) * Rads
      el(46) = 30.06896348d0 - 0.0000000342768d0 * D
      el(47) = 0.00858587d0 + 0.000000000688296d0 * D
      el(48) = (Rads * (304.88003d0 + 0.0059810572d0 * D))
c  // Pluto
      el(49) = (17.14175d0 + 0.0000000841889d0 * D) * Rads
      el(50) = (110.30347d0 - 0.0000002839d0 * D) * Rads
      el(51) = (224.06676d0 - 0.00000100578d0 * D) * Rads
      el(52) = 39.48168677d0 - 0.0000000210574d0 * D
      el(53) = 0.24880766d0 + 0.00000000177002d0 * D
      el(54) = (Rads * (238.92881d0 + 0.003931834d0 * D))      

      q = 6 * (p - 1)
      ip = el(q + 1)
      op = el(q + 2)
      pp = el(q + 3)
      ap = el(q + 4)
      ep = el(q + 5)
      lp = el(q + 6)
      ie = el(13)
      oe = el(14)
      pe = el(15)
      ae = el(16)
      ee = el(17)
      le = el(18)    

c  //Get Earths position using Keplers equation
      Me1 = (le - pe)
      B = Me1 / (2.0d0 * pi)
      Me1 = 2.0d0 * pi * (B - floor(abs(B)))
      if (B.lt.0.0d0) Me1 = 2.0d0 * pi * (B + floor(abs(B)))
      if (Me1.lt.0.0d0) Me1 = 2.0d0 * pi + Me1
      e = Me1
      delta = 0.05d0
      do while (abs(delta).ge.10.0d0**-12)
         delta = e - ee * sin(e) - Me1
         e = e - delta / (1.0d0 - ee * cos(e))
      enddo
      ve = 2.0d0 * atan(sqrt((1 + ee) / (1 - ee)) * tan(0.5d0 * e));
      if (ve.lt.0.0d0) ve = ve + 2.0d0 * pi
      re = ae * (1.0d0 - ee * ee) / (1.0d0 + ee * cos(ve))
      xe = re * cos(ve + pe)
      ye = re * sin(ve + pe)
      ze = 0.0d0

c  //Get planets position using Keplers equation
      mp = (lp - pp)
      B = mp / (2.0d0 * pi)
      mp = 2 * pi * (B - floor(abs(B)))
      if (B.lt.0.0d0) mp = 2.0d0 * pi * (B + floor(abs(B)))
      if (mp.lt.0) mp = 2.0d0 * pi + mp
      e = mp
      delta = 0.05
      do while (abs(delta).ge.10.0d0**-12)
         delta = e - ep * sin(e) - mp
         e = e - delta / (1.0d0 - ep * cos(e))
      enddo
      vp = 2 * atan(sqrt((1.0d0 + ep) / (1 - ep)) * tan(0.5d0 * e))
      if (vp.lt.0.0d0) vp = vp + 2.0d0 * pi
      rp = ap * (1.0d0 - ep * ep) / (1.0d0 + ep * cos(vp))
      xh = rp*(cos(op)*cos(vp+pp-op)-sin(op)*sin(vp+pp-op)*cos(ip))
      yh = rp*(sin(op)*cos(vp+pp-op)+cos(op)*sin(vp+pp-op)*cos(ip))
      zh = rp * (sin(vp + pp - op) * sin(ip))
      xg = xh - xe
      yg = yh - ye
      zg = zh
c  //compute RA and DEC
      ecl = 23.429292d0 * Rads
      xeq = xg
      yeq = yg * cos(ecl) - zg * sin(ecl)
      zeq = yg * sin(ecl) + zg * cos(ecl)
      ra = atan(yeq/ xeq)*12.0d0/pi
      if (xeq.lt.0.0d0) ra = ra + 12.0d0
      if (yeq.lt.0.0d0) then
         if (xeq.gt.0) ra = ra + 24.0d0
      endif
      dec = 180.0d0*atan(zeq / sqrt(xeq * xeq + yeq * yeq))/pi
c  // Sun Coodinates
      xeq = xe
      yeq = ye * cos(ecl) - ze * sin(ecl)
      zeq = ye * sin(ecl) + ze * cos(ecl)
      rae = 12.0d0 + atan(yeq/ xeq)*12.0d0/pi
      if (xe.lt.0.0) rae = rae + 12.0d0
      if (ye.lt.0) then
         if (xe.gt.0) rae = rae + 24
      endif
      dece = -180.0d0*atan(zeq / sqrt(xeq * xeq + yeq * yeq))/pi
      if (p.eq.3)then
         ra=rae
         dec=dece
      endif
      if (ra.lt.12.0d0) then
         raa=12.0d0-ra
      else 
         raa=36.0d0-ra
      endif
      xpix=32+floor(raa*965.0d0/24.0d0)
      ypix=201+floor(dec*-3.0d0)
      
      return
      end


c   //the Planets
c      nm[1] = "Mercury"
c      nm[2] = "Venus"
c      nm[3] = "Sun"
c      nm[4] = "Mars"
c      nm[5] = "Jupiter"
c      nm[6] = "Saturn"
c      nm[7] = "Uranus"
c      nm[8] = "Neptune"
c      nm[9] = "Pluto"

     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jd2day(jdate,iday,imonth,iyear,ierr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer iday,imonth,iyear,ierr,igreg,ia,ja,jb,jc,jd,je,julian
      parameter(igreg=2299161)
      double precision jdate,xc
      
      julian=int(jdate)

      if(julian.lt.0) then
        ierr = -1
        return
      else
        ierr = 0
      endif
      
      if (julian.ge.igreg) then
        ia = (real(julian-1867216)-0.25d0)/36524.25d0
        ja = julian + 1+ia-int(0.25*ia)
      else
        ja = julian
      end if      

      jb = ja + 1524
      xc = (dble(jb-2439870)-122.1d0)/365.25d0
      jc = 6680.0d0 + xc
      jd = 365*jc + int(0.25d0*real(jc))
      je = int(real(jb-jd)/30.6001d0)

      iday = jb - jd - int(30.6001d0*real(je))

      imonth = je - 1
      if (imonth.gt.12) imonth = imonth - 12

      iyear = jc - 4715
      if (imonth.gt.2) iyear = iyear - 1
      if (iyear.le.0) iyear = iyear - 1
      
      return
      end
      