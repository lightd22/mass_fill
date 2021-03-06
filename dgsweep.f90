! --------------------------------------------------------------------
! Nodal DG sweep method for DG/FV hybrid 1D Advection
! By: Devin Light
! --------------------------------------------------------------------

SUBROUTINE DGSWEEP(rhoq,num_elem,elemdx,nnodes,nodes,wghts,u,N,DG_D,DG_C,DG_CINV,dodghybrid,doposlimit,oddstep,dt)

    USE DGmod ! Contains useful functions, parameters, and subroutines
    
    IMPLICIT NONE

    INTEGER, PARAMETER :: DOUBLE = KIND(1D0) ! Specification of DOUBLE

    ! ---
    ! Inputs
    ! ---
    INTEGER, INTENT(IN) :: num_elem,N,nnodes
    REAL(KIND=DOUBLE), INTENT(IN) :: elemdx,dt
    REAL(KIND=DOUBLE), DIMENSION(1:N), INTENT(IN) :: u
    REAL(KIND=DOUBLE), DIMENSION(1:N), INTENT(INOUT) :: rhoq
    REAL(KIND=DOUBLE), DIMENSION(0:nnodes,0:nnodes), INTENT(IN) :: DG_D,DG_C,DG_CINV
    REAL(KIND=DOUBLE), DIMENSION(0:nnodes), INTENT(IN) :: nodes, wghts
	LOGICAL, INTENT(IN) :: dodghybrid,doposlimit,oddstep
    
    ! ---
    ! Local variables
    ! ---
    INTEGER :: i,j,k
    REAL(KIND=DOUBLE), DIMENSION(0:nnodes,0:num_elem+1) :: A,A1,A2, utild
    REAL(KIND=DOUBLE), DIMENSION(0:nnodes,0:num_elem+1) :: B
    REAL(KIND=DOUBLE), DIMENSION(0:nnodes,1:num_elem) :: rqBAR
    REAL(KIND=DOUBLE) :: PI
	REAL(KIND=DOUBLE), DIMENSION(0:nnodes) :: HOLDER1,HOLDER2


    PI = DACOS(-1D0)
    HOLDER1 = 0D0
    HOLDER2 = 0D0
    
    ! ########################################
    ! A(k,j) gives a_k(t) in the jth element for rhoq 
    ! B(k,j) gives a_k(t) in the jth element for rhop
    ! ########################################

    ! Reform incoming nodal values to be more convienent
    DO j=1,num_elem
        rqBAR(:,j) = rhoq(1+(nnodes+1)*(j-1) : (nnodes+1)*j)
        utild(:,j) = u(1+(nnodes+1)*(j-1) : (nnodes+1)*j)
    END DO

    ! For ppmdg hybrid, values incoming assumed to be cell averages, invert with CmatINV, giving coefficents a_k(t) 
	! in series expansion which when averaged over the evenly spaced subcells in each element, give the incoming values.
    ! Use these values to initialize the A and B matricies
	IF(dodghybrid) THEN
	    DO j=1,num_elem
    		    DO i=0,nnodes
    		        DO k = 0,nnodes
    		            HOLDER1(k) = DG_CINV(i,k)*rqBAR(k,j)
    		        END DO
    		        A(i,j) = SUM(HOLDER1)
    		    END DO
    		END DO
	ELSE
	! Otherwise, just use reshaped values
	A(:,1:num_elem) = rqBAR(:,1:num_elem)
	END IF

	! -- Enforce periodicity
	A(:,0) = A(:,num_elem)
	A(:,num_elem+1) = A(:,1)
	utild(:,0) = utild(:,num_elem)
	utild(:,num_elem+1) = utild(:,1)

    ! #######################
    ! Time step using SSPRK3
    ! #######################

    DO j = 1, num_elem
       DO k = 0, nnodes
           A1(k,j) = A(k,j) + dt*rk3rhs(A,utild,k,j,nnodes,num_elem,DG_D,elemdx,wghts)  
       END DO
    END DO
	! -- Enforce periodicity
    A1(:,0) = A1(:,num_elem)
	A1(:,num_elem+1) = A1(:,1)

    DO j = 1, num_elem
      DO k = 0, nnodes
           A2(k,j) = (3D0/4D0)*A(k,j) + (1D0/4D0)*(A1(k,j) + dt*rk3rhs(A1,utild,k,j,nnodes,num_elem,DG_D,elemdx,wghts))
      END DO
    END DO
	! -- Enforce periodicity
	A2(:,0) = A2(:,num_elem)
	A2(:,num_elem+1) = A2(:,1)

   DO j = 1,num_elem
       DO k = 0, nnodes
           A(k,j) = (1D0/3D0)*A(k,j) + (2D0/3D0)*(A2(k,j) + dt*rk3rhs(A2,utild,k,j,nnodes,num_elem,DG_D,elemdx,wghts))
       END DO
   END DO

	IF(dodghybrid) THEN
    ! After time stepping is complete, use DG_C to re-average the series expansion to get back to cell-averaged values
    ! on the evenly spaced grid to send back to PPM
	    HOLDER1 = 0D0
	    HOLDER2 = 0D0
	    DO j = 1,num_elem
	        DO i=0,nnodes    
	            DO k=0,nnodes
	                HOLDER1(k) = DG_C(i,k)*A(k,j)
	            END DO
	            rqBAR(i,j) = SUM(HOLDER1)
	       END DO
	    END DO
	ELSE
	! Otherwise, just send back DG nodal values (read: coefficents) 
		DO j=1,num_elem
			rqBAR(:,j) = A(:,j)
		END DO
	END IF

	! Do mass-filling based positivity preservation on each element
	IF(doposlimit) THEN
		CALL MFILL(rqBAR,nnodes,num_elem,oddstep)
	END IF
	
	write(*,*)
	DO j=1,num_elem
		write(*,*) SUM(rqBAR(:,j))
   	END DO

    ! Reform original rp and rq vectors to send back
	! -- Note: There are 2 ghost cells that we dont update
    DO j=1,num_elem
        rhoq(1+(nnodes+1)*(j-1) : (nnodes+1)*j) = rqBAR(:,j)
    END DO
    CONTAINS

	REAL(KIND=DOUBLE) FUNCTION Gflux(A,jleft,jright,u,nnodes,nelems)
		IMPLICIT NONE
		! -- Inputs
		INTEGER, INTENT(IN) :: nnodes,nelems,jleft,jright
		REAL(KIND=DOUBLE), DIMENSION(0:nnodes,0:nelems+1), INTENT(IN) :: A,u

		! -- Note that it is assumed that u(nnodes,jleft) = u(0,jright)
		IF(u(nnodes,jleft)>0D0) THEN
			Gflux = u(nnodes,jleft)*A(nnodes,jleft)
		ELSE
			Gflux = u(0,jright)*A(0,jright)
		END IF
	END FUNCTION Gflux

    REAL(KIND = DOUBLE) FUNCTION rk3rhs(Ain,utild,k,j,nnodes,num_elem,DG_D,elemdx,wghts)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: k,j,num_elem, nnodes
        REAL(KIND = DOUBLE), DIMENSION(0:nnodes,0:num_elem+1), INTENT(IN) :: Ain,utild
        REAL(KIND = DOUBLE), DIMENSION(0:nnodes,0:nnodes), INTENT(IN) :: DG_D
        REAL(KIND = DOUBLE), DIMENSION(0:nnodes), INTENT(IN) ::wghts
        REAL(KIND = DOUBLE), INTENT(IN) :: elemdx

		! -- Local variables            
        INTEGER :: n
		REAL(KIND=DOUBLE), DIMENSION(0:nnodes) :: HOLDER

        HOLDER = 0D0
        DO n = 0, nnodes
			HOLDER(n) = utild(n,j)*Ain(n,j)*DG_D(k,n)*wghts(n)
        END DO
		rk3rhs = SUM(HOLDER)

		! -- Fluxes
		IF(k .eq. nnodes) THEN
			rk3rhs = rk3rhs - Gflux(Ain,j,j+1,utild,nnodes,num_elem) 
		END IF
		IF(k .eq. 0) THEN
			rk3rhs = rk3rhs + Gflux(Ain,j-1,j,utild,nnodes,num_elem) 
		END IF
           
        rk3rhs = rk3rhs*(2D0/(wghts(k)*elemdx))

        END FUNCTION rk3rhs

	SUBROUTINE MFILL(rhoq,nnodes,nelems,oddstep)
		IMPLICIT NONE
		! Inputs
		! --
		INTEGER, INTENT(IN) :: nnodes,nelems
		LOGICAL, INTENT(IN) :: oddstep
		REAL(KIND=DOUBLE), DIMENSION(0:nnodes,1:nelems), INTENT(INOUT) :: rhoq

		! Local Variables
		! --
		INTEGER :: max_i,min_i,i,j,k,step
		REAL(KIND=DOUBLE) :: mxfer

		IF(oddstep) THEN
			max_i = nnodes
			min_i = 0
			step = 1
		ELSE
			max_i = 0
			min_i = nnodes
			step = -1
		END IF

		DO j=1,nelems
			DO i=min_i,max_i,step
				IF(rhoq(i,j) .lt. 0D0) THEN
!						write(*,*) '(i,j)',i,j
					DO k=i+step,max_i,step
						IF(rhoq(i,j) .ge. 0D0) EXIT
						mxfer = DMIN1(DMAX1(0D0,rhoq(k,j)),DABS(rhoq(i,j)))
!						write(*,*) 'mxfer=',mxfer
						rhoq(i,j) = rhoq(i,j) + mxfer
						rhoq(k,j) = rhoq(k,j) - mxfer
					END DO
					DO k=i-step,min_i,-1*step
						IF(rhoq(i,j) .ge. 0D0) EXIT
						mxfer = DMIN1(DMAX1(0D0,rhoq(k,j)),DABS(rhoq(i,j)))
!						write(*,*) 'mxfer=',mxfer
						rhoq(i,j) = rhoq(i,j) + mxfer
						rhoq(k,j) = rhoq(k,j) - mxfer
					END DO
				END IF
			END DO
		END DO
	END SUBROUTINE MFILL

END SUBROUTINE DGSWEEP
