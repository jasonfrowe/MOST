      program widthfix
      implicit none
      integer npt,nmax,nunit,i
      parameter(nmax=400000)
      real rmavg,time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),
     .   yc(nmax),fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),skystd(nmax)
      double precision dtime(nmax)
      character*80 filename
      
      
      write(6,*) "Enter MOST photometry filename"
      read(5,*) filename
      
      nunit=10 !set unit number for file read
      
      open(unit=nunit,file=filename,status='old',err=901)
      
C     read in photometry data  
      call readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
C     we have photometry, can close the file now.
      close(nunit)
      
      call widthfit(npt,time,dtime,mag,merr,sky,xc,yc,
     .         fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
      
      filename="widthcor.dat"
      rmavg=0.
      call exportdata(npt,dtime,mag,merr,filename,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,rmavg)
      
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine widthfit(npt,time,dtime,mag,merr,sky,xc,yc,
     .         fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,nfit,ma,np,mp,i,k,npt2
      parameter(nmax=400000,nfit=5)
      real time(npt),mag(npt),merr(npt),sky(npt),xc(npt),yc(npt),
     .   fx(npt),fy(npt),fxy(npt),tboard(npt),mfield(npt),ftotal(npt),
     .   aflux(npt),nap(npt),itime(npt),skystd(npt),a(nfit),
     .   u(nmax,nfit),v(nfit,nfit),w(nfit),chisq,magcor,x(nmax),y(nmax),
     .   z(nmax),fxoff,fyoff,magoff,sig(nmax)
      double precision dtime(npt)
      
      
C     offset to bring all data to co-ords 0,0,0 for fitting
      call avgrm(npt,fx,fxoff)
      call avgrm(npt,fy,fyoff)
      call avgrm(npt,mag,magoff)
      write(6,*) "offsets: ",fxoff,fyoff,magoff
      
      npt2=0
      do 10 i=1,npt
         if((fx(i).gt.-0.4).and.(fx(i).lt.0.4).and.(fy(i).gt.-0.4).and.
     .      (fy(i).lt.0.4).and.(mag(i).gt.-0.4).and.
     .      (mag(i).lt.0.4))then
            npt2=npt2+1
            x(npt2)=fx(i)
            y(npt2)=fy(i)
            z(npt2)=mag(i)
            sig(npt2)=merr(i)
         endif
 10   continue
      
      ma=nfit
      np=nfit
      mp=nmax
      call svdfit(x,y,z,sig,npt2,a,ma,u,v,w,mp,np,chisq)
      
      write(6,*) (a(i),i=1,nfit),chisq/real(npt2-1),npt2
      
      do 21 i=1,npt
         magcor=0.
         do 22 k=2,ma,2 
            magcor=magcor+a(k)*fx(i)**(k/2)+
     .         a(k+1)*fy(i)**(k/2)
 22      continue
         mag(i)=mag(i)-magcor!+magoff
         fx(i)=fx(i)+fxoff
         fy(i)=fy(i)+fyoff
c         write(6,*) magcor
 21   continue
         
      return
      end
      
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
      PARAMETER(NMAX=400000,MMAX=50,TOL=1.E-5)
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