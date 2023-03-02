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
!     Subroutines related to initializing linear solver arrays and
!     function calls to svFSILS and Trilinos solver library
!
!--------------------------------------------------------------------

      SUBROUTINE LSALLOC(lEq)
      USE COMMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq

      IF (ALLOCATED(R)) DEALLOCATE(R)
      ALLOCATE (R(dof,tnNo))
      R(:,:) = 0._RKIND

      IF (.NOT.lEq%assmTLS) THEN
         IF (ALLOCATED(Val)) DEALLOCATE(Val)
         ALLOCATE (Val(dof*dof,lhs%nnz))
         Val(:,:) = 0._RKIND
      END IF

#ifdef WITH_TRILINOS
      IF (lEq%useTLS) THEN
         IF (ALLOCATED(tls%W)) THEN
            DEALLOCATE(tls%W, tls%R)
            CALL TRILINOS_LHS_FREE()
         END IF
         ALLOCATE(tls%W(dof,tnNo), tls%R(dof,tnNo))
         CALL TRILINOS_LHS_CREATE(gtnNo, lhs%mynNo, tnNo, lhs%nnz,
     2      tls%ltg, ltg, rowPtr, colPtr, dof)
      END IF
#endif

      RETURN
      END SUBROUTINE LSALLOC
!####################################################################
!     The solution of the linear system is return in the global variable R, which
!     contains the residual before the linear solve
      SUBROUTINE LSSOLVE(lEq, incL, res)
      USE COMMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER(KIND=IKIND), INTENT(IN) :: incL(nFacesLS)
      REAL(KIND=RKIND), INTENT(IN) :: res(nFacesLS) ! Neumann cplBC resistance

#ifdef WITH_TRILINOS
      INTEGER(KIND=IKIND) a

      IF (lEq%useTLS) CALL INIT_DIR_AND_COUPNEU_BC(incL, res)

      IF (lEq%assmTLS) THEN
         lEq%FSILS%RI%suc = .FALSE.
         CALL TRILINOS_SOLVE(tls%R, tls%W, lEq%FSILS%RI%fNorm,
     2      lEq%FSILS%RI%iNorm, lEq%FSILS%RI%itr, lEq%FSILS%RI%callD,
     3      lEq%FSILS%RI%dB, lEq%FSILS%RI%suc, lEq%ls%LS_type,
     4      lEq%FSILS%RI%reltol, lEq%FSILS%RI%mItr, lEq%FSILS%RI%sD,
     5      lEq%ls%PREC_Type, lEq%assmTLS)

      ELSE IF(lEq%useTLS) THEN
         CALL TRILINOS_GLOBAL_SOLVE(Val, R, tls%R, tls%W,
     2      lEq%FSILS%RI%fNorm, lEq%FSILS%RI%iNorm, lEq%FSILS%RI%itr,
     3      lEq%FSILS%RI%callD, lEq%FSILS%RI%dB, lEq%FSILS%RI%suc,
     4      lEq%ls%LS_type, lEq%FSILS%RI%reltol, lEq%FSILS%RI%mItr,
     5      lEq%FSILS%RI%sD, lEq%ls%PREC_Type)

      ELSE
#endif
         CALL FSILS_SOLVE(lhs, lEq%FSILS, dof, R, Val,
     2      lEq%ls%PREC_Type, incL=incL, res=res)
#ifdef WITH_TRILINOS
      END IF

      IF (lEq%useTLS) THEN
         DO a=1, tnNo
            R(:,a) = tls%R(:,lhs%map(a))
         END DO
      END IF
#endif

      RETURN
      END SUBROUTINE LSSOLVE
!--------------------------------------------------------------------
      SUBROUTINE INIT_DIR_AND_COUPNEU_BC(incL, res)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: incL(lhs%nFaces)
      REAL(KIND=RKIND), INTENT(IN) :: res(lhs%nFaces)

      REAL(KIND=RKIND), ALLOCATABLE :: v(:,:)

      INTEGER(KIND=IKIND) a, i, Ac, faIn, faDof
      LOGICAL flag, isCoupledBC

      IF (lhs%nFaces .NE. 0) THEN
         lhs%face%incFlag = .TRUE.
         DO faIn=1, lhs%nFaces
            IF (incL(faIn) .EQ. 0)  lhs%face(faIn)%incFlag = .FALSE.
         END DO
         DO faIn=1, lhs%nFaces
            lhs%face(faIn)%coupledFlag = .FALSE.
            IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
            flag = lhs%face(faIn)%bGrp .EQ. BC_TYPE_Neu
            IF (flag .AND. res(faIn).NE.0._RKIND) THEN
               lhs%face(faIn)%res = res(faIn)
               lhs%face(faIn)%coupledFlag = .TRUE.
            END IF
         END DO
      END IF

      tls%W = 1._RKIND
      DO faIn=1, lhs%nFaces
         IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
         faDof = MIN(lhs%face(faIn)%dof,dof)
         IF (lhs%face(faIn)%bGrp .EQ. BC_TYPE_Dir) THEN
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               DO i=1, faDof
                  tls%W(i,Ac) = tls%W(i,Ac) * lhs%face(faIn)%val(i,a)
               END DO
            END DO
         END IF
      END DO

      ALLOCATE(v(dof,tnNo))
      v = 0._RKIND
      isCoupledBC = .FALSE.
      DO faIn=1, lhs%nFaces
         IF (lhs%face(faIn)%coupledFlag) THEN
            isCoupledBC = .TRUE.
            faDof = MIN(lhs%face(faIn)%dof,dof)
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               DO i=1, faDof
                  v(i,Ac) = v(i,Ac) +
     2               SQRT(ABS(res(faIn)))*lhs%face(faIn)%val(i,a)
               END DO
            END DO
         END IF
      END DO

#ifdef WITH_TRILINOS
      CALL TRILINOS_BC_CREATE(v, isCoupledBC)
#endif
      DEALLOCATE(v)

      END SUBROUTINE INIT_DIR_AND_COUPNEU_BC

!####################################################################
!     Compute the residual with nodal DOFs An
! 
!     How to compute the residual? Copy most of the functions from the inner loop
!     in MAIN(). But don't want to compute tangent or solve linear system.
!
!     Compute residual and take dot product at element level, then sum to compute g.
!     See GLOBALEQASSEM

      SUBROUTINE CALCRES(Phin)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Phin
      REAL(KIND=LSRP) g, FSILS_DOTV

!     SKETCH. Functions taken from MAIN(). But want to eliminate functions that
!     compute the tangent and do other unnecessary things.

      An_old = An
      Yn_old = Yn
      Dn_old = Dn

      An = Phin


      CALL PICC(Ag, Yg, Dg)

      IF (cplBC%coupled .AND. cEq.EQ.1) THEN
         CALL SETBCCPL
         CALL SETBCDIR(An, Yn, Dn)
      END IF

!        Initiator step (quantities at n+am, n+af)
      CALL PICI(Ag, Yg, Dg)
      IF (ALLOCATED(Rd)) THEN
         Rd = 0._RKIND
         Kd = 0._RKIND
      END IF

      dbg = 'Allocating the RHS and LHS'
      CALL LSALLOC(eq(cEq))

!        Compute body forces. If phys is shells or CMM (init), apply
!        contribution from body forces (pressure) to residue
      CALL SETBF(Dg)

      dbg = "Assembling equation <"//eq(cEq)%sym//">"
      DO iM=1, nMsh
         CALL GLOBALEQASSEM(msh(iM), Ag, Yg, Dg)
         dbg = "Mesh "//iM//" is assembled"
      END DO

!        Treatment of boundary conditions on faces
!        Apply Neumman or Traction boundary conditions
      CALL SETBCNEU(Yg, Dg)

!        Apply CMM BC conditions
      IF (.NOT.cmmInit) CALL SETBCCMM(Ag, Dg)

!        Apply weakly applied Dirichlet BCs
      CALL SETBCDIRW(Yg, Dg)

!        Apply contact model and add its contribution to residue
      IF (iCntct) CALL CONTACTFORCES(Dg)

!        Synchronize R across processes. Note: that it is important
!        to synchronize residue, R before treating immersed bodies as
!        ib%R is already communicated across processes
      IF (.NOT.eq(cEq)%assmTLS) CALL COMMU(R)

!        Update residue in displacement equation for USTRUCT phys.
!        Note that this step is done only first iteration. Residue
!        will be 0 for subsequent iterations
      IF (sstEq) CALL USTRUCTR(Yg)

      IF ((eq(cEq)%phys .EQ. phys_stokes) .OR.
2          (eq(cEq)%phys .EQ. phys_fluid)  .OR.
3          (eq(cEq)%phys .EQ. phys_ustruct).OR.
4          (eq(cEq)%phys .EQ. phys_fsi)) THEN
         CALL THOOD_ValRC()
      END IF

      CALL SETBCUNDEFNEU()

!        IB treatment: for explicit coupling, simply construct residue.
      IF (ibFlag) THEN
         IF (ib%cpld .EQ. ibCpld_I) THEN
            CALL IB_IMPLICIT(Ag, Yg, Dg)
         END IF
         CALL IB_CONSTRUCT()
      END IF

      incL = 0
      IF (eq(cEq)%phys .EQ. phys_mesh) incL(nFacesLS) = 1
      IF (cmmInit) incL(nFacesLS) = 1
      DO iBc=1, eq(cEq)%nBc
         i = eq(cEq)%bc(iBc)%lsPtr
         IF (i .NE. 0) THEN
!                 scaled resistance value for Neumann surface, to be "added" to stiffness matrix in LSSOLVE
            res(i) = eq(cEq)%gam*dt*eq(cEq)%bc(iBc)%r 
!                 For DEBUGGING
!                  IF (cm%mas()) THEN
!                     PRINT*, "iBc: ", iBc, 'i: ', i, 'res(i): ', res(i)
!                  END IF
            incL(i) = 1
         END IF
      END DO

      An = An_old 
      Yn = Yn_old
      Dn = Dn_old
      END SUBROUTINE CALCRES
!####################################################################

!     Compute the scalar function g(alpha) = (Delta A)^T * G
!     where Delta A is the Newton increment in dof and G is the residual, as 
!     part of the Newton method linear system K * (Delta A) = - G
!     Used in LINESEARCH below
      SUBROUTINE CALCG(alpha)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: alpha
      REAL(KIND=LSRP) g, FSILS_DOTV

!     The full Newton increment in nodal DOFs
      DA = R(a:b)

!     The residual vector at this alpha. G(A + alpha * DA). 
      G = CALCRES(A)

!     Compute dot product between vectors, then MPI ALL REDUCE
      g = FSILS_DOTV(dof,lhs%mynNo, lhs%commu, DA, G)

!     Reset R to the full Newton increment, since the rest of the code thinks
!     R contains the Newton increment, not the residual
      R = DA
      
      END SUBROUTINE CALCG
!####################################################################
!     Perform a line search to determine the Newton damping parameter 
!     alpha = [0,1] that reduces |g(alpha)| <= 0.8 |g(0)|, where
!     g(alpha) = (Delta A)^T * G
!     where Delta A is the Newton increment in dof and G is the residual, as 
!     part of the Newton method linear system K * (Delta A) = - G
!     In this code, the variable R contains the residual G before the linear solve
!     but contains Delta A after the linear solve.
!     This subroutine uses the secant method to determine alpha.
!     See FEM book by Wriggers pg. 158
      SUBROUTINE LINESEARCH
      USE COMMOD
      IMPLICIT NONE

!     Compute g(alpha = 0)
      alpha_km1 = 0
      g_km1 = CALCG(alpha_km1)
      g_0 = g_km1

!     Compute g(alpha = 1)
      alpha_k = 0
      g_k = CALCG(alpha_k)

!     If there is a zero in 0 <= alpha <= 1, do secant method. Otherwise, use alpha = 1
!     Check this by checking g(alpha = 0) * g(alpha = 1) < 0
      IF (g_km1 * g_k > 0) THEN 
         alpha = 1._RKIND
      ELSE
         DO
!        Get new alpha from secant method
         alpha_kp1 = alpha_k - g_k * (alpha_k - alpha_km1) / (g_k - g_km1)

!        Compute g(alpha_k+1)
         g_kp1 = CALCG(alpha_kp1)

!        Check if |g(alpha_k+1)| < 0.8 |g(alpha_0)|
         IF (ABS(g_kp1) < 0.8 * ABS(g_0)) THEN
            alpha = alpha_kp1
            EXIT
         END IF
         END DO
      END IF

!     Multiply Newton increment with damping parameter alpha
      R = R * alpha

      END SUBROUTINE LINESEARCH

