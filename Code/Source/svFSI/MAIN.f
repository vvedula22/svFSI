!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------------------------------------------------------------------
!
!     Main routine that contains the calls to all major routines and
!     general structure of the code.
!
!--------------------------------------------------------------------

      PROGRAM MAIN
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL l1, l2, l3, kflag, lnsrch
      INTEGER(KIND=IKIND) i, iM, iBc, ierr, iEqOld, stopTS
      REAL(KIND=RKIND) timeP(3)

      INTEGER(KIND=IKIND), ALLOCATABLE :: incL(:)
      REAL(KIND=RKIND), ALLOCATABLE :: Ag(:,:), Yg(:,:), Dg(:,:), res(:)

      IF (IKIND.NE.LSIP .OR. RKIND.NE.LSRP) THEN
         STOP "Incompatible datatype precision between solver and FSILS"
      END IF

      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.

      savedOnce = .FALSE.
      CALL MPI_INIT(i)
      CALL cm%new(MPI_COMM_WORLD)

!     Initiating the exception tracing
      CALL EXCEPTIONS

      resetSim  = .FALSE.
      rmsh%cntr = 0

!     Reading the user-defined parameters from foo.inp
 101  CALL READFILES

!     Doing the partitioning and distributing the data to the all
!     Processors
      CALL DISTRIBUTE

!     Initializing the solution vectors and constructing LHS matrix
!     format
      CALL INITIALIZE(timeP)
      stopTS = nTS

      dbg = 'Allocating intermediate variables'
      ALLOCATE(Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo),
     2   res(nFacesLS), incL(nFacesLS))

!--------------------------------------------------------------------
!     Outer loop for marching in time. When entring this loop, all old
!     variables are completely set and satisfy BCs.
      IF (cTS .LE. nITS) dt = dt/10._RKIND
      DO
!     Adjusting the time step size once initialization stage is over
         IF (cTS .EQ. nITS) THEN
            dt = dt*10._RKIND
            std = " New time step size: "//dt
         END IF
!     Incrementing time step, hence cTS will be associated with new
!     variables, i.e. An, Yn, and Dn
         cTS    = cTS + 1
         time   = time + dt
         cEq    = 1
         eq%itr = 0
         eq%ok  = .FALSE.

!     Compute mesh properties to check if remeshing is required
         IF (mvMsh .AND. rmsh%isReqd) THEN
            CALL CALCMESHPROPS(nMsh, msh)
            IF (resetSim) EXIT
         END IF

!     Predictor step
!     Predicts new quantities (An, Yn, Dn) from old quantities (Ao, Yo, Do) using
!     gen-alpha relations
         CALL PICP

!     Apply Dirichlet BCs strongly
         CALL SETBCDIR(An, Yn, Dn)

!     Inner loop for Newton iteration (with line search)
         DO
            iEqOld = cEq

!           Increment Newton iteration counter
            eq(cEq)%itr = eq(cEq)%itr + 1


!           BEGIN CALCKR.
            kflag = .TRUE.
            CALL CALCKR(kflag, eq(cEq), incL, res)
!           END CALCKR

!           Solve linear system K*DA = R. The solution is stored in R. 
!           Note that the correct linear system for Newtons method is
!           K*DA = -R, so after this function, we treat -R as the Newton increment.
            dbg = "Solving equation <"//eq(cEq)%sym//">"
            CALL LSSOLVE(eq(cEq), incL, res)

!           ADD LINE SEARCH LOOP HERE
            lnsrch = .TRUE.
            IF (lnsrch) THEN
               CALL LINESEARCH
            ELSE
!              Solution is obtained, now updating (Corrector)
!              Note the corrector step inside the NR loop, which is
!              slightly different from https://www.scorec.rpi.edu/~kjansen/genalf.pdf
!   
!              REMOVE convergence checks and equation increment from PICC
               CALL PICC
            END IF

!        Increment equation (taken from PICC, since in LINESEARCH, we may call 
!        PICC multiple times.)
         IF (eq(cEq)%coupled) THEN
            cEq = cEq + 1
            IF (ALL(.NOT.eq%coupled .OR. eq%ok)) THEN
               DO WHILE (cEq .LE. nEq)
                  IF (.NOT.eq(cEq)%coupled) EXIT
                  cEq = cEq + 1
               END DO
            ELSE
               IF (cEq .GT. nEq) cEq = 1
               DO WHILE (.NOT.eq(cEq)%coupled)
                  cEq = cEq + 1
                  IF (cEq .GT. nEq) cEq = 1
               END DO
            END IF
         ELSE
            IF (eq(cEq)%ok) cEq = cEq + 1
         END IF

!        Checking for exceptions
            CALL EXCEPTIONS

!        Writing out the time passed, residual, and etc.
            IF (ALL(eq%ok)) EXIT
            CALL OUTRESULT(timeP, 2, iEqOld)
         END DO
!     End of inner loop

!     IB treatment: interpolate flow data on IB mesh from background
!     fluid mesh for explicit coupling, update old solution for implicit
!     coupling
         IF (ibFlag) THEN
            CALL IB_INTERPYU(Yn, Dn)
            IF (ib%cpld .EQ. ibCpld_I) THEN
               ib%Auo = ib%Aun
               ib%Ubo = ib%Ubn
            END IF
         END IF

!     Saving the TXT files containing average and fluxes
         CALL TXT(.FALSE.)

         IF (rmsh%isReqd) THEN
            l1 = MOD(cTS,rmsh%cpVar) .EQ. 0
            IF (l1) THEN
               rmsh%rTS = cTS-1
               rmsh%time = time-dt
               rmsh%iNorm(:) = eq(:)%iNorm
               rmsh%A0(:,:) = Ao(:,:)
               rmsh%Y0(:,:) = Yo(:,:)
               rmsh%D0(:,:) = Do(:,:)
            END IF
         END IF

         IF (cm%mas()) THEN
            INQUIRE(FILE=stopTrigName, EXIST=l1)
            IF (l1) THEN
               OPEN(664,FILE=TRIM(stopTrigName))
               READ(664,'(I10)',ADVANCE='NO',IOSTAT=ierr) stopTS
               CLOSE(664)
               IF (ierr .EQ. -1) stopTS = cTS
            ELSE
               stopTS = nTS
            END IF
         END IF
         CALL cm%bcast(stopTS)
         l1 = cTS .GE. stopTS
         l2 = MOD(cTS,stFileIncr) .EQ. 0
!     Saving the result to restart bin file
         IF (l1 .OR. l2) CALL WRITERESTART(timeP)

!     Writing results into the disk with VTU format
         IF (saveVTK) THEN
            l2 = MOD(cTS,saveIncr) .EQ. 0
            l3 = cTS .GE. saveATS
            IF (l2 .AND. l3) THEN
               CALL OUTRESULT(timeP, 3, iEqOld)
               CALL WRITEVTUS(An, Yn, Dn, .FALSE.)
               IF (ibFlag) CALL IB_WRITEVTUS(ib%Yb, ib%Ubo)
            ELSE
               CALL OUTRESULT(timeP, 2, iEqOld)
            END IF
         ELSE
            CALL OUTRESULT(timeP, 2, iEqOld)
         END IF
         IF (pstEq) CALL OUTDNORM()

         IF (ibFlag) CALL IB_OUTCPUT()

!     Exiting outer loop if l1
         IF (l1) EXIT

!     Solution is stored here before replacing it at next time step
         Ao = An
         Yo = Yn
         IF (dFlag) Do = Dn
         cplBC%xo = cplBC%xn
      END DO
!     End of outer loop

      IF (resetSim) THEN
         CALL REMESHRESTART(timeP)
         DEALLOCATE(Ag, Yg, Dg, incL, res)
         IF (ALLOCATED(tls)) THEN
            DEALLOCATE(tls%ltg, tls%W, tls%R)
            DEALLOCATE(tls)
         END IF
         GOTO 101
      END IF

      IF (l1 .AND. saveAve) CALL CALCAVE

      DEALLOCATE(Ag, Yg, Dg, incL, res)
      CALL FINALIZE()
      CALL MPI_FINALIZE(ierr)

      END PROGRAM MAIN
!####################################################################
      SUBROUTINE STOPSIM()

      CALL FINALIZE
      STOP "MPI is forced to stop by a fatal error"

      END SUBROUTINE STOPSIM
!####################################################################
