      PROGRAM PENDULA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C     PENDULA IS A PROGRAM THAT NUMERICALLY SOLVES THE EQUATIONS OF   C
C     MOTION FOR A COUPLED SET OF PENDULA WITH PHYSICAL PARAMETERS    C
C     STORED IN COMMON ARRAY PHYSPARAM AND USER'S CHOICE OF IC'S.     C
C*********************************************************************C
      PARAMETER(NMAX=50)  
C
      COMMON/PHYSPARAM/GRAV,P1M,P2M,P1L,P2L
C
      DIMENSION YINP(NMAX),YOUT(NMAX),DYDT(NMAX)
C
      DATA PI/3.141592653589D0/
C
C     COEFFICIENT FOR GRAVITY
      GRAV = 9.8D+00
C
C     MASSES OF PENDULA (KG)
      P1M = 1.0D+00
      P2M = 2.0D+00
C
C     LENGTHS OF RODS (M)
      P1L = 1.0D+00
      P2L = 2.0D+00
C
C     INITIAL ANGLES (RADIANS OF ELEVATION FROM DIRECTLY DOWN, ANTICLOCKWISE)
      YINP(1) = PI/3.0D0
      YINP(3) =-PI/3.0D0
C
C     INITIAL ANGULAR VELOCITIES (RADIANS PER SECOND)
      YINP(2) = 9.0D+00
      YINP(4) =15.0D+00
C
C     INITAL AND FINAL TIME VALUES (SECONDS), NUMBER OF STEPS BETWEEN
      T1    = 0.0D+00
      T2    = 1.0D+01
C
C     GNUPLOT COMMAND ASKS FOR 50 FRAMES PER SECOND
      NSTEP = INT(50*(T2-T1))
C     NSTEP = 2000
C
C     CALL THE ODE SOLVER
      CALL DERIVS(T1,YINP,DYDT)
      CALL RKDUMB(T1,T2,YINP,YOUT,NSTEP,DERIVS)
C
C     ENERGY DOMAIN BOUNDARIES
      EDOM1 = 2.0D0*P2L*P2M*GRAV
      EDOM2 = 2.0D0*P1L*P1M*GRAV + EDOM1
C
C     TELL THE USER WHICH ENERGY DOMAIN WE ARE IN
      IF(EINIT.LT.EDOM1) THEN
        WRITE(*,*) 'Low energy domain (neither pendulum can flip)'
      ELSEIF(EINIT.GE.EDOM1.AND.EINIT.LT.EDOM2) THEN
        WRITE(*,*) 'Medium energy domain (one pendulum can flip)'
      ELSEIF(EINIT.GE.EDOM2) THEN
        WRITE(*,*) 'High energy domain (both pendula can flip)' 
      ENDIF
C
C     OUTPUT RESULTS
      WRITE(*,*) 'TIME    ',T1,T2
      WRITE(*,*) 'THETA(1)',YINP(1),YOUT(1)
      WRITE(*,*) 'OMEGA(1)',YINP(2),YOUT(2)
      WRITE(*,*) 'THETA(2)',YINP(3),YOUT(3)
      WRITE(*,*) 'OMEGA(2)',YINP(4),YOUT(4)
C
20    FORMAT(F6.3)
21    FORMAT(I5)
C
C     PRINT PHYSICAL PARAMETERS TO A FILE FOR PLOTTING PURPOSES
      OPEN(UNIT=17,FILE='plots/physparam',STATUS='UNKNOWN')
      REWIND(UNIT=17)
        WRITE(17,20) GRAV
        WRITE(17,20) P1M
        WRITE(17,20) P2M
        WRITE(17,20) P1L
        WRITE(17,20) P2L
        WRITE(17,20) EINIT
        WRITE(17,21) NSTEP
      CLOSE(UNIT=17)   
C     
      STOP
      END
C
      SUBROUTINE RKDUMB(T1,T2,YINP,YOUT,NSTEP,DERIVS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C     RKDUMB STARTS FROM INITIAL VALUES YINP(1:N) KNOWN AT T1, AND    C
C     USES THE FOURTH-ORDER RUNGE-KUTTA METHOD TO ADVANCE IN NSTEP    C
C     EQUAL INCREMENTS TO T2. RESULTS ARE STORED IN A COMMON BLOCK.   C
C ------------------------------------------------------------------- C
C  INPUT:                                                             C
C       N - NUMBER OF COUPLED ORDINARY DIFFERENTIAL EQUATIONS.        C
C    YINP - INITIAL VALUES AT LOCATION T1.                            C
C      T1 - STARTING POINT OF SOLVER.                                 C
C      T2 - TERMINATING POINT OF SOLVER.                              C
C   NSTEP - NUMBER OF EQUALLY-SPACED INCREMENTS.                      C
C  DERIVS - USER-SUPPLIED SUBROUTINE FOR CALCULATING RHS DERIVATIVE.  C
C  OUTPUT:                                                            C
C      TT - STORAGE OF INTERMEDIATE STEP LOCATIONS.                   C
C    YOUT - STORAGE OF INTERMEDIATE FUNCTION VALUES.                  C
C*********************************************************************C
      PARAMETER(NMAX=50)
C
      COMMON/PHYSPARAM/GRAV,P1M,P2M,P1L,P2L
C
      DIMENSION YINP(NMAX),YOUT(NMAX),Y(NMAX),DV(NMAX),V(NMAX)
C
      DATA PI/3.141592653589D0/
C
30    FORMAT(11F9.4)
C
C     REWIND THE ODE DATA FILE AND START RECORDING RESULTS
      OPEN(UNIT=15,FILE='plots/odesolver',STATUS='UNKNOWN')
      REWIND(UNIT=15)   
C
C       LOAD STARTING VALUES
        DO I=1,4
          V(I) = YINP(I)
          Y(I) = V(I)
        ENDDO
C
C       FOR PENDULUM: EVALUATE GENERALISED MOMENTA AND JACOBI INTEGRAL
        PTM = P1M+P2M
        P1 = PTM*P1L*P1L*Y(2) + P2M*P1L*P2L*Y(4)*DCOS(Y(3)-Y(1))
        P2 = P2M*P2L*P2L*Y(4) + P2M*P1L*P2L*Y(2)*DCOS(Y(3)-Y(1))
        E1 = P1L*PTM*(0.5D0*P1L*Y(2)*Y(2)-GRAV*(DCOS(Y(1))-1.0D0))
        E2 = P2L*P2M*(0.5D0*P2L*Y(4)*Y(4)-GRAV*(DCOS(Y(3))-1.0D0))
        EB = P2M*P1L*P2L*Y(2)*Y(4)*DCOS(Y(3)-Y(1))
        ET = E1 + E2 + EB
C
        T = T1
        H = (T2-T1)/NSTEP
C
C       WRITE STARTING VALUES
        WRITE(15,30) T,(Y(I),I=1,4),P1,P2,E1,E2,EB,ET
C
C       TAKE NSTEP STEPS
        DO K=1,NSTEP
C
          CALL DERIVS(X,V,DV)
          CALL RK4(T,V,DV,H,DERIVS,V)
C
C         VARIABLES ARE RADIANS SO KEEP ARGUMENTS WITHIN WINDOW [-PI,+PI]
C          DO I=1,4
C           IF(V(1).GT.+PI) V(1) = V(1) - 2.0D0*PI
C           IF(V(1).LT.-PI) V(1) = V(1) + 2.0D0*PI
C           IF(V(3).GT.+PI) V(3) = V(3) - 2.0D0*PI
C           IF(V(3).LT.-PI) V(3) = V(3) + 2.0D0*PI
CC          ENDDO
C         
C         STEP SIZE TOO SMALL
          IF(T+H.EQ.T) THEN
            WRITE(*,*) 'In RKDUMB: step size not significant.'
            STOP
          ENDIF
C
C         FOR PENDULUM: EVALUATE GENERALISED MOMENTA AND JACOBI INTEGRAL
          PTM = P1M+P2M
          P1 = PTM*P1L*P1L*Y(2) + P2M*P1L*P2L*Y(4)*DCOS(Y(3)-Y(1))
          P2 = P2M*P2L*P2L*Y(4) + P2M*P1L*P2L*Y(2)*DCOS(Y(3)-Y(1))
          E1 = P1L*PTM*(0.5D0*P1L*Y(2)*Y(2)-GRAV*(DCOS(Y(1))-1.0D0))
          E2 = P2L*P2M*(0.5D0*P2L*Y(4)*Y(4)-GRAV*(DCOS(Y(3))-1.0D0))
          EB = P2M*P1L*P2L*Y(2)*Y(4)*DCOS(Y(3)-Y(1))
          ET = E1 + E2 + EB
C
C         STEP UP T VALUE AND STORE INTERMEDIATE VALUES
          T = T + H
          DO I=1,4
            Y(I) = V(I)
          ENDDO
C
          WRITE(15,30) T,(Y(I),I=1,4),P1,P2,E1,E2,EB,ET
C          
        ENDDO
C
C     CLOSE THE DATA FILE
      CLOSE(UNIT=15)
C
      DO I=1,4
        YOUT(I) = Y(I)
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE RK4(T,YIN,DYDT,H,DERIVS,YOUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C     RK4 USES THE FOURTH-ORDER RUNGE-KUTTA METHOD TO INCREMENT A     C
C     SYSTEM OF N COUPLED ORDINARY DIFFERENTIAL EQUATIONS BY STEP H,  C
C     GIVEN VALUES Y(1:N) AND THEIR DERIVATIVES DYDT(1:N) AT T.       C 
C ------------------------------------------------------------------- C
C  INPUT:                                                             C
C    4 - NUMBER OF COUPLED ORDINARY DIFFERENTIAL EQUATIONS (PARAM).C
C       T - STARTING POINT OF SOLVER (T + H IS THE ENDING POINT).     C
C     YIN - ARRAY OF DIMENSION 4 CONTAINING STARTING VALUES.       C
C    DYDT - ARRAY OF DIMENSION 4 CONTAINING STARTING DERIVATIVES.  C
C  DERIVS - USER-SUPPLIED SUBROUTINE FOR CALCULATING RHS DERIVATIVE.  C
C       H - STEP SIZE INTERVAL.                                       C
C  OUTPUT:                                                            C
C    YOUT - ARRAY OF DIMENSION 4 CONTAINING INCREMENTED VARIABLES. C
C*********************************************************************C
      PARAMETER(NMAX=50)
C
      COMMON/PHYSPARAM/GRAV,P1M,P2M,P1L,P2L
C
      DIMENSION YIN(NMAX),DYDT(NMAX),YOUT(NMAX),DYM(NMAX),DYT(NMAX),
     &          YT(NMAX)
C
      HH   = H*0.5D0
      H6   = H/6.0D0
      TH   = T + HH
C
C     FIRST STEP
      DO I=1,4
        YT(I) = YIN(I) + HH*DYDT(I)
      ENDDO
C
C     SECOND STEP
      CALL DERIVS(TH,YT,DYT)
      DO I=1,4
        YT(I) = YIN(I) + HH*DYT(I)
      ENDDO
C
C     THIRD STEP
      CALL DERIVS(TH,YT,DYM)
      DO I=1,4
        YT(I)  = YIN(I) + H*DYM(I)
        DYM(I) = DYT(I) +   DYM(I)
      ENDDO
C
C     FOURTH STEP
      CALL DERIVS(T+H,YT,DYT)
      DO I=1,4
        YOUT(I) = YIN(I) + H6*(DYDT(I) + DYT(I) + 2.0D0*DYM(I))
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE DERIVS(T,Y,DYDT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C     DERIVS CONTAINS THE INFORMATION FOR DERIVATIVES OF FUNCTIONS.   C
C     ODD ENTRIES ARE ANGULAR VELOCITIES, EVEN ENTRIES ARE ANGLES.    C
C*********************************************************************C
      PARAMETER(NMAX=50)
C
      DIMENSION Y(NMAX),DYDT(NMAX)
C 
      COMMON/PHYSPARAM/GRAV,P1M,P2M,P1L,P2L
C
      PTM = P1M + P2M
      TDF = Y(1) - Y(3)
      TUN = Y(1) - 2.0D0*Y(3)

      XA = -GRAV*(2.0D0*P1M+P2M)*DSIN(Y(1))
      XB = -P2M*GRAV*DSIN(TUN)
      XC = -2.0D0*DSIN(TDF)*P2M*(Y(4)*Y(4)*P2L+Y(2)*Y(2)*P1L*DCOS(TDF))
      XD =  P1L*(2.0D0*P1M + P2M - P2M*DCOS(2.0D0*TDF))

      XE = Y(2)*Y(2)*P1L*PTM + GRAV*PTM*DCOS(Y(1)) 
      XF = Y(4)*Y(4)*P2L*P2M*DCOS(TDF)
      XG = 2.0D0*DSIN(TDF)*(XE+XF)
      XH = P2L*(2.0D0*P1M + P2M - P2M*DCOS(2.0D0*TDF))
C
C     DEFINE DERIVATIVE VALUES AT T
      DYDT(1) = Y(2)
      DYDT(2) =(XA+XB+XC)/XD
      DYDT(3) = Y(4)
      DYDT(4) = XG/XH
C
      RETURN
      END
