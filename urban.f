      PROGRAM URBAN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             UU    UU RRRRRRR  BBBBBBB     AA    NN    NN             C
C             UU    UU RR    RR BB    BB   AAAA   NNN   NN             C
C             UU    UU RR    RR BB    BB  AA  AA  NNNN  NN             C
C             UU    UU RR    RR BBBBBBB  AA    AA NN NN NN             C
C             UU    UU RRRRRRR  BB    BB AAAAAAAA NN  NNNN             C
C             UU    UU RR    RR BB    BB AA    AA NN   NNN             C
C              UUUUUU  RR    RR BBBBBBB  AA    AA NN    NN             C
C -------------------------------------------------------------------- C
C            PROGRAM FOR ANALYSING URBAN POPULATION DENSITY.           C
C**********************************************************************C
      PARAMETER(NDIM=50,NMAX=40)
C
      CHARACTER*80 TITLE
C
      DIMENSION NMAP(NDIM,NDIM),NTAL(0:NMAX)
C
C     READ MATRIX RESULTS (INTEGERS)
      OPEN(UNIT=10,FILE='4.txt',STATUS='UNKNOWN')
      REWIND(UNIT=10)
      DO I=1,NDIM
        READ(10,*) (NMAP(I,J),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=10)
C
      TITLE = 'citypop'
      CALL IGNUMTRX(NMAP,TITLE,NDIM)
C
C     TOTAL POPULATION OF CITY
      NTOT = NPOP(NMAP,NDIM)
      WRITE(*,*) 'Total population of city = ',NTOT
C
C     HISTOGRAM OF POPULATION DENSITIES
      CALL NHIST(NMAP,NDIM,NTAL,NMAX)
      WRITE(*,*) 'Histogram of cell populations:'
      DO N=0,NMAX
        WRITE(*,*) N,NTAL(N)
      ENDDO
C
C     HISTOGRAM INDICATING REGION TYPES
      CALL REGIONS(NMAP,NDIM,NTAL,NMAX)
      WRITE(*,*) 'Histogram of region types:'
      DO N=0,NMAX
        WRITE(*,*) NTAL(N)
      ENDDO
C
C     BLUR AND PLOT THE MATRIX
      CALL BLUR(NMAP,NDIM)
C     
      STOP
      END
C
C
      FUNCTION NPOP(MTRX,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C     NPOP COUNTS THE TOTAL POPULATION OF THE CITY.                    C
C**********************************************************************C
C
      DIMENSION MTRX(NDIM,NDIM)
C
      NCOUNT = 0
C
      DO M=1,NDIM
        DO N=1,NDIM
          NCOUNT = NCOUNT + MTRX(M,N)
        ENDDO
      ENDDO
C
      NPOP = NCOUNT
C
      RETURN
      END
C
C
      SUBROUTINE NHIST(MTRX,NDIM,NTAL,NMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C     NPOP COUNTS THE TOTAL POPULATION OF THE CITY.                    C
C**********************************************************************C
C
      DIMENSION MTRX(NDIM,NDIM),NTAL(0:NMAX)
C
C     INITIALISE GRID
      DO N=0,NMAX
        NTAL(N) = 0
      ENDDO
C
C     SORT CELL POPULATIONS
      DO I=1,NDIM
        DO J=1,NDIM
          NCELL = MTRX(I,J)
          NTAL(NCELL) = NTAL(NCELL)+1
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE REGIONS(MTRX,NDIM,NTAL,NMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C     REGIONS PROVIDES A SLIDING 'CITY'/'SUBURBS' LABEL TO CELLS.      C
C**********************************************************************C
C
      CHARACTER*80 TITLE
C
      DIMENSION CITY(NDIM,NDIM)
      DIMENSION MTRX(NDIM,NDIM),NTAL(0:NMAX)
C
C     INITIALISE MAXIMUM REGION COUNTER
      RMAX = 0.0D0
C
C     LOOP OVER ALL CELL ADDRESSES
      DO I=1,NDIM
        DO J=1,NDIM
C
C         INITIALISE REGION COUNTER
          RCOUNT = 0.0D0
C
C         LOOP OVER ALL NEIGHBOURING CELLS
          DO K=1,NDIM
            DO L=1,NDIM
C
C             EXCLUDING (I,J), ADD PROXIMITY CONTRIBUTION TO COUNTER
              IF(K.NE.I.AND.L.NE.J) THEN
C
C               SHORTEST DISTANCE FROM I TO K AND J TO L
                IK = MIN(IABS(K-I),NDIM-IABS(K-I))
                JL = MIN(IABS(L-J),NDIM-IABS(L-J))
C
                DIST = DFLOAT(IK*IK + JL*JL)
                POP  = DFLOAT(MTRX(K,L))
C
                RCOUNT = RCOUNT + POP/DIST
C
              ENDIF
C
C           END LOOP OVER NEIGHBOURING CELLS
            ENDDO
          ENDDO
C
C         UPDATE MAXIMUM REGION COUNTER
          IF(RCOUNT.GT.RMAX) RMAX = RCOUNT
C
C         ADD REGION COUNTER TO ARRAY
          CITY(I,J) = RCOUNT
C
C       END LOOP OVER CELL ADDRESSES
        ENDDO
      ENDDO
C
C     NORMALISE CITY INDEX
c      DO I=1,NDIM
c        DO J=1,NDIM
c          CITY(I,J) = CITY(I,J)/RMAX
c        ENDDO
c      ENDDO
C
C     INITIALISE HISTOGRAM TALLY
      DO N=0,NMAX
        NTAL(N) = 0
      ENDDO

C     SORT RESULTS INTO THE TALLY NTAL
      DO I=1,NDIM
        DO J=1,NDIM
          NIND = INT(CITY(I,J)*NMAX/RMAX)
          NTAL(NIND) = NTAL(NIND)+1
        ENDDO
      ENDDO
      
      TITLE = "regions"
      CALL GNUMTRX(CITY,TITLE,NDIM)
C
      RETURN
      END
C
C
      SUBROUTINE BLUR(MTRX,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C     BLUR AVERAGES OVER 3x3 CELL CITES AND PLOTS THE RESULT.          C
C**********************************************************************C
C
      CHARACTER*80 TITLE
      DIMENSION MTRX(NDIM,NDIM),MBLR(NDIM,NDIM)
C
C     INITIALISE BLURRED GRID
      DO I=1,NDIM
        DO J=1,NDIM
          MBLR(I,J) = 0
        ENDDO
      ENDDO
C
C     ADD CONTRIBUTIONS
      DO I=1,NDIM
        DO J=1,NDIM
          IF(I.EQ.1) THEN
            IL = NDIM
          ELSE
            IL = I-1
          ENDIF
          IF(I.EQ.NDIM) THEN
            IR = 1
          ELSE
            IR = I+1
          ENDIF
C
          IF(J.EQ.1) THEN
            JL = NDIM
          ELSE
            JL = J-1
          ENDIF
          IF(J.EQ.NDIM) THEN
            JR = 1
          ELSE
            JR = J+1
          ENDIF
C          
          MBLR(I,J) = MTRX(IL,JL) + MTRX(IL,J ) + MTRX(IL,JR)
     &              + MTRX(I ,JL) + MTRX(I ,J ) + MTRX(I ,JR)
     &              + MTRX(IR,JL) + MTRX(IR,J ) + MTRX(IR,JR)
C
        ENDDO
      ENDDO
C
C     DIVIDE BY 4 AND PUT BACK INTO MTRX
      DO I=1,NDIM
        DO J=1,NDIM
          MTRX(I,J) = MBLR(I,J)/9
        ENDDO
      ENDDO
C
C     PLOT THE RESULT
      TITLE = 'blur'
      CALL IGNUMTRX(MTRX,TITLE,NDIM)
C
      RETURN
      END
C
C
      SUBROUTINE IGNUMTRX(IARRAY,TITLE,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   GGGGGG  NN    NN UU    UU MM       MM TTTTTTTT RRRRRRR  XX     XX  C
C  GG    GG NNN   NN UU    UU MMM     MMM    TT    RR    RR  XX   XX   C
C  GG       NNNN  NN UU    UU MMMM   MMMM    TT    RR    RR   XX XX    C
C  GG       NN NN NN UU    UU MM MM MM MM    TT    RR    RR    XXX     C
C  GG   GGG NN  NNNN UU    UU MM  MMM  MM    TT    RRRRRRR    XX XX    C
C  GG    GG NN   NNN UU    UU MM   M   MM    TT    RR    RR  XX   XX   C
C   GGGGGG  NN    NN  UUUUUU  MM       MM    TT    RR    RR XX     XX  C
C                                                                      C
C -------------------------------------------------------------------- C
C  IGNUMTRX EXPORTS AN IARRAY TO AN EXTERNAL DATA FILE AND PLOTS IT.   C
C**********************************************************************C
C
      CHARACTER*80 TITLE
C
      DIMENSION IARRAY(NDIM,NDIM)
C
C     PRINT TO EXTERNAL DATA FILE
      OPEN(UNIT=8,FILE="plots/"//TRIM(TITLE)//".dat",STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO I=1,NDIM
        WRITE(8, *) (IARRAY(I,J),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=8)
C
      XEND = DFLOAT(NDIM)-0.5D0
      YEND = DFLOAT(NDIM)-0.5D0
C
C     WRITE GNUPLOT MAKE FILE
      OPEN(UNIT=9,FILE='plots/'//TRIM(TITLE)//'.gnuplot',
     &                                                 STATUS='REPLACE')
      WRITE(9,'(A)') '#'//TRIM(TITLE)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#  Usage:'
      WRITE(9,'(A)') '#  gnuplot < '//TRIM(TITLE)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Terminal output specs'
      WRITE(9,'(A)') 'set terminal pdf size 4,4'
      WRITE(9,'(A)') 'set output "plots/'//TRIM(TITLE)//'.pdf"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') 'load "plots/pals/jet.pal"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Axes and title'
      WRITE(9,'(A)') 'set title sprintf("'//TRIM(TITLE)//'")'
      WRITE(9,'(A,F6.1,A)') 'set xrange [-0.5:',XEND,']'
      WRITE(9,'(A,F6.1,A)') 'set yrange [',YEND,':-0.5] reverse'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plot data to file'
      WRITE(9,'(A,I2,A,I2,A)') 'plot "plots/'//TRIM(TITLE)//'.dat"'
     &                                    //' matrix with image notitle'
      CLOSE(UNIT=9)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(TITLE)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(TITLE)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE GNUMTRX(ARRAY,TITLE,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   GGGGGG  NN    NN UU    UU MM       MM TTTTTTTT RRRRRRR  XX     XX  C
C  GG    GG NNN   NN UU    UU MMM     MMM    TT    RR    RR  XX   XX   C
C  GG       NNNN  NN UU    UU MMMM   MMMM    TT    RR    RR   XX XX    C
C  GG       NN NN NN UU    UU MM MM MM MM    TT    RR    RR    XXX     C
C  GG   GGG NN  NNNN UU    UU MM  MMM  MM    TT    RRRRRRR    XX XX    C
C  GG    GG NN   NNN UU    UU MM   M   MM    TT    RR    RR  XX   XX   C
C   GGGGGG  NN    NN  UUUUUU  MM       MM    TT    RR    RR XX     XX  C
C                                                                      C
C -------------------------------------------------------------------- C
C  GNUMTRX EXPORTS AN IARRAY TO AN EXTERNAL DATA FILE AND PLOTS IT.    C
C**********************************************************************C
C
      CHARACTER*80 TITLE
C
      DIMENSION ARRAY(NDIM,NDIM)
C
C     PRINT TO EXTERNAL DATA FILE
      OPEN(UNIT=8,FILE="plots/"//TRIM(TITLE)//".dat",STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO I=1,NDIM
        WRITE(8, *) (ARRAY(I,J),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=8)
C
      XEND = DFLOAT(NDIM)-0.5D0
      YEND = DFLOAT(NDIM)-0.5D0
C
C     WRITE GNUPLOT MAKE FILE
      OPEN(UNIT=9,FILE='plots/'//TRIM(TITLE)//'.gnuplot',
     &                                                 STATUS='REPLACE')
      WRITE(9,'(A)') '#'//TRIM(TITLE)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#  Usage:'
      WRITE(9,'(A)') '#  gnuplot < '//TRIM(TITLE)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Terminal output specs'
      WRITE(9,'(A)') 'set terminal pdf size 4,4'
      WRITE(9,'(A)') 'set output "plots/'//TRIM(TITLE)//'.pdf"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') 'load "plots/pals/jet.pal"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Axes and title'
      WRITE(9,'(A)') 'set title sprintf("'//TRIM(TITLE)//'")'
      WRITE(9,'(A,F6.1,A)') 'set xrange [-0.5:',XEND,']'
      WRITE(9,'(A,F6.1,A)') 'set yrange [',YEND,':-0.5] reverse'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plot data to file'
      WRITE(9,'(A,I2,A,I2,A)') 'plot "plots/'//TRIM(TITLE)//'.dat"'
     &                                    //' matrix with image notitle'
      CLOSE(UNIT=9)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(TITLE)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(TITLE)//'.pdf')
C
      RETURN
      END

