§      PROGRAM GALROT
C
      DOUBLE PRECISION R0,DIST,DDIST,PI,AL,R2,D,R
      DOUBLE PRECISION VR,B,CL,SL,CB,CB2,THETA0,THETA
      DOUBLE PRECISION A1,A2,A3
      CHARACTER*30 FILEOUT
C
      write(6,*)'Program to calculate LSR radial velocity assuming' 
      write(6,*)'circular rotation and the Galactic rotation curve'
      write(6,*)'of Brand & Blitz (1993, A&A, 275, 67) with R_0=8.5 kpc'
      write(6,*)'and Theta_0=220 km/s. The derived values of V_LSR as a'
      write(6,*)'function of distance are written to galrot.dat'
      write(6,*)
      write(6,*)'Enter name output file'
      read(*,'(A30)') FILEOUT
C
      open(unit=10,file=FILEOUT,status='unknown')
C
      r0=8.5D0
      theta0=220.D0
C
      write(*,*)'enter galactic longitude and latitude: '
      read(*,*)al,b
      write(*,*)'enter max distance (kpc), dist. increment (kpc): '
      read(*,*)dist,ddist
C
      pi=3.14159265D0
      al=al*pi/180.D0
      b=b*pi/180.D0
      cl=dcos(al)
      sl=dsin(al)
      cb=dcos(b)
      cb2=cb*cb
C
      a1=1.00767D0
      a2=0.0394D0
      a3=0.00712D0
C
      do d=0,dist,ddist
         r2=r0**2+d**2*cb2-2*r0*d*cl*cb
         r=sqrt(r2)
         theta=(a1*(r/r0)**a2+a3)*theta0
         vr=(theta*r0/r - theta0)*sl*cb
         write(10,'(F6.2,F10.3)')d,vr
      enddo
      end
