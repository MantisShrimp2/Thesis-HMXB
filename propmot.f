C Program propmot.f
C Computes proper motion in the local rest frame (i.e. corrected
C for differential galactic rotation and solar peculiar motion) at
C 3 distances (d/1.4, d, and 1.4d). The assumption is made that
C the target is in the galactic plane (b=0). The input parameters
C are read from a list containing 5 columns: 
C (1) Target name; 
C (2) galactic longitude l (degrees); 
C (3) observed proper motion mu_l; 
C (4) mu_b; 
C (5) distance d (in kpc). 
C Based on formulae in Comeron et al. 1998, A&A 330, 975.
       DIMENSION D(3),CORMUL(3),CORMUB(3),VTL(3),VTB(3),VT(3)
       REAL*4 LONG,LATT,MUL(3),MUB(3)
       CHARACTER*20 NAME,FILEIN,FILEOUT
       PI = 3.1415926535898
C
C Read input from file
C
       WRITE(*,*) 'Enter name input file'
       READ(*,*) FILEIN
       WRITE(*,*) 'Enter name output file'
       READ(*,*) FILEOUT
       OPEN(10,FILE=FILEIN,FORM='FORMATTED',
     +      STATUS='OLD')
       OPEN(11,FILE=FILEOUT,FORM='FORMATTED',
     +      STATUS='UNKNOWN')
       DO K = 1,200
          READ(10,500,END=999) NAME,LONG,OBSMUL,OBSMUB,DIST
C          READ(10,500,END=999) NAME,LONG,LATT,DIST,Z,OBSMUL,OBSMUB
 500      FORMAT(A15,4(F10.2))
C 500      FORMAT(A19,F12.8,F13.8,F6.2,I6,F8.2,F7.2)
C
C Compute MUL and MUB due to galactic rotation and solar motion
C          
          RADLONG = PI*LONG/180.0
          D(1) = DIST/1.4
          D(2) = DIST
          D(3) = 1.4*DIST
          DO L = 1,3
             TERM1 = 1.56*SIN(RADLONG) - 3.22*COS(RADLONG)
             TERM2 = 2.74*COS(2.0*RADLONG) - 2.74
             MUL(L) = TERM1/D(L) + TERM2
             MUB(L) = 8.0/(4.74*D(L))
C
C     Calculate proper motion in LSR and convert
C     proper motions ("/yr) in velocities (km/s) at
C     given distance: 1 pc/yr := 977767.6376 km/s
C
             CORMUL(L) = OBSMUL - MUL(L)
             CORMUB(L) = OBSMUB - MUB(L)
             PMLDEG = CORMUL(L)/3600.0
             PMLRAD = PI*PMLDEG/180.0
             PMBDEG = CORMUB(L)/3600.0
             PMBRAD = PI*PMBDEG/180.0
             VTL(L) = D(L)*TAN(PMLRAD)*977767.6376             
             VTB(L) = D(L)*TAN(PMBRAD)*977767.6376
             VT(L) = SQRT(VTL(L)**2 + VTB(L)**2)
          ENDDO
          WRITE(11,1000) NAME,'Distance (kpc):',D(1),D(2),D(3)
          WRITE(11,1000) NAME,'MUL ("/yr):',MUL(1),MUL(2),MUL(3)
          WRITE(11,1000) NAME,'PM long ("/yr):',CORMUL(1),
     +                   CORMUL(2),CORMUL(3)
          WRITE(11,1000) NAME,'PM latt ("/yr):',CORMUB(1),
     +                   CORMUB(2),CORMUB(3)
          WRITE(11,1000) NAME,'Vtan long (km/s):',VTL(1),
     +                   VTL(2),VTL(3)          
          WRITE(11,1000) NAME,'Vtan latt (km/s):',VTB(1),
     +                   VTB(2),VTB(3)
          WRITE(11,1000) NAME,'Vtan      (km/s):',VT(1),
     +                   VT(2),VT(3)
          WRITE(11,*) ' '
 1000     FORMAT(A10,A20,3(F12.3))
       ENDDO
 999   CLOSE(10)
       CLOSE(11)
       END
