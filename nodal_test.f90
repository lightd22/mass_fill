! 1D Nodal Advection Test 2
! Checking constant advection with variable nodes in 1d 
! -------------------------------------------------------

PROGRAM EXECUTE
    USE DGmod
	USE netcdf

    IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	INTEGER :: nelems,nnodes
	REAL(KIND=DOUBLE) :: mu,tfinal
	LOGICAL :: dodghybrid, doposlimit

	nnodes = 5
	nelems = 16
	mu = 0.3D0
	tfinal = 1.0D0
	dodghybrid = .true.
	doposlimit = .true.

	WRITE(*,*)
	WRITE(*,*) '1d Square Wave Advection w/ Mass Filling'
	CALL DRIVER(nnodes,nelems,mu,tfinal,dodghybrid,doposlimit)

	WRITE(*,*)
	WRITE(*,*) '1d Square Wave Advection w/out Mass Filling'
	doposlimit = .false.
	CALL DRIVER(nnodes,nelems,mu,tfinal,dodghybrid,doposlimit)
	
	WRITE(*,*)
	WRITE(*,*) 'SUBROUTINE COMPLETE!'

    CONTAINS

    SUBROUTINE DRIVER(nnodes,nelems,mu,tfinal,dodghybrid,doposlimit)
        IMPLICIT NONE
		! -- Inputs
		INTEGER, INTENT(IN) :: nnodes,nelems
		REAL(KIND=DOUBLE), INTENT(IN) :: mu,tfinal
		LOGICAL, INTENT(IN) :: dodghybrid,doposlimit

		! -- Local variables
		INTEGER :: nsteps, nxp, nxppm
		REAL(KIND=DOUBLE), DIMENSION(1:nelems) :: ecent
		REAL(KIND=DOUBLE) :: dxel, PI, dt,dxp,t,dxmin, dxppm
		REAL(KIND=DOUBLE), DIMENSION(0:nnodes,1:nelems) :: xdg
		REAL(KIND=DOUBLE), DIMENSION(0:nnodes,0:nelems+1) :: A,u
		REAL(KIND=DOUBLE), DIMENSION(0:nnodes,0:nnodes) :: D,C,Cinv
		REAL(KIND=DOUBLE), DIMENSION(0:nnodes) :: nodes,wghts
		REAL(KIND=DOUBLE), ALLOCATABLE, DIMENSION(:) :: tvals,xppm,inplt,soln,rhoq,utild
		REAL(KIND=DOUBLE), DIMENSION(0:nnodes) :: HOLDER
		REAL(KIND=DOUBLE), DIMENSION(1:nnodes) :: nodespace
		LOGICAL :: oddstep
		REAL(KIND=4), DIMENSION(2) :: tstart,tend,thold
		REAL(KIND=4) :: t0,tf
		REAL(KIND=DOUBLE) :: qmin

		INTEGER :: i,j,k,n ! Looping variables

		! ---- Stuff For NETCDF 
        CHARACTER(LEN=13) :: cdf_out
		INTEGER :: cdfid ! ID for netCDF file
		INTEGER, PARAMETER :: NDIMS = 2, NDIMS2 = 3
		INTEGER :: NX, NE,NT, ierr, DGN,DGE
	    INTEGER :: idinit,idrhoq,idt,idx, dimids(NDIMS),dimids2(NDIMS2),idA,idnodes,idCPUt
	    INTEGER :: x_dimid, t_dimid, dgn_dimid,dge_dimid,cpu_dimid
		INTEGER, DIMENSION(1:NDIMS) :: start, count
		INTEGER, DIMENSION(1:NDIMS2) :: start2,count2
		
        PI = DACOS(-1D0)
		
		IF(doposlimit) THEN
			cdf_out = 'dg1d_mfill.nc'
		ELSE
			cdf_out = 'dg1d_nofil.nc'
		END IF

		! -- Start real time calculation
		t0 = etime(tstart)

		! --  Compute GLL nodes, weights, and the D-matrix
		CALL gllNewton(nnodes,nodes)
		CALL weights(nnodes,nodes,wghts)
		CALL Dmat(nnodes,nodes, D)

		! -- Set up x-grid via elements
		dxel = 1D0/nelems
		ecent(1) = dxel/2D0
		DO j=2,nelems
			ecent(j) = ecent(j-1)+dxel
		END DO
		DO j=1,nelems
			xdg(:,j) = ecent(j) + 0.5D0*dxel*nodes(:)
		END DO

		DO k=1,nnodes
			nodespace(k) = xdg(k,1)-xdg(k-1,1)
		END DO
!		write(*,*)
!		write(*,*) nodespace

		dxmin = MINVAL(nodespace)

		! -- Set up ppm grid
		nxppm = (nnodes+1)*nelems
		dxppm = 1D0/nxppm

		ALLOCATE(xppm(1:nxppm),rhoq(1:nxppm),utild(1:nxppm))

		xppm(1) = dxppm/2D0
		DO j=2,nxppm
			xppm(j) = xppm(j-1)+dxppm
		END DO

		! -- Initialize rhoq on ppm grid for dgsweep
!		DO j=1,nxppm
!			rhoq(j) = DSIN(6*PI*xppm(j)) + DSIN(8*PI*xppm(j))	
!		END DO
		rhoq = 0D0
		where (xppm .lt. 0.75D0 .and. xppm .gt. 0.25D0)
			rhoq = 1D0
		end where

		! -- Initialize u (For now, just constant velocity)
		u = 1D0
		! --  Reshape for dgsweep
		DO j=1,nelems
			utild(1+(j-1)*(nnodes+1):(j)*(nnodes+1)) = u(:,j)
		END DO
		
		! -- Set up dt and time array
		dt = (mu*dxmin)/MAXVAL(DABS(u))
		nsteps = INT(tfinal/dt)
		ALLOCATE(tvals(1:nsteps))

		! -- Compute inital conditions on plotting grid
		ALLOCATE(inplt(1:nxppm))
		inplt = 0D0
!		DO j=1,nxppm
!			inplt(j) = DSIN(6*PI*xppm(j)) + DSIN(8*PI*xppm(j))
!		END DO
		where (xppm .lt. 0.75D0 .and. xppm .gt. 0.25D0)
			inplt = 1D0
		end where

	
		! Initialize C and Cinv
		! --
		CALL Cmat(nnodes,nodes,wghts,dxppm,dxel,C)
		CALL GEPP_INV(C,nnodes+1,Cinv)

		! -- Create netCDF file for output
		! --
		NX = (nnodes+1)*nelems
		NT = nsteps
		DGN = (nnodes+1)
		DGE = nelems

		ierr = nf90_create(cdf_out,NF90_CLOBBER,cdfid)
		write(*,*) 'Creating dimensions..'
		ierr = nf90_def_dim(cdfid, "nx", NX, x_dimid)
		ierr = nf90_def_dim(cdfid, "nt", NT, t_dimid)
		ierr = nf90_def_dim(cdfid, "dgn",DGN,dgn_dimid)
		ierr = nf90_def_dim(cdfid, "dge",DGE,dge_dimid)
		ierr = nf90_def_dim(cdfid,"cputime",2,cpu_dimid)

		dimids = (/ t_dimid, x_dimid /)
		dimids2 = (/ t_dimid,dge_dimid,dgn_dimid /)
		write(*,*) 'Creating variables..'
		ierr = NF90_DEF_VAR(cdfid, "rhoq", NF90_FLOAT, dimids, idrhoq)
		ierr = NF90_DEF_VAR(cdfid, "x",NF90_FLOAT,x_dimid,idx)
		ierr = NF90_DEF_VAR(cdfid, "time", NF90_FLOAT, t_dimid,idt)
		ierr = NF90_DEF_VAR(cdfid, "nodes",NF90_FLOAT, dgn_dimid, idnodes)
		ierr = NF90_DEF_VAR(cdfid, "init",NF90_FLOAT,x_dimid,idinit)
		ierr = NF90_DEF_VAR(cdfid, "ctime",NF90_FLOAT,cpu_dimid,idcput)
		write(*,*) 'Exiting define mode..'
		ierr = nf90_enddef(cdfid)
		! --

		! -- Output plotting grid and initial conditions
		ierr = nf90_put_var(cdfid, idinit,inplt)
		ierr = nf90_put_var(cdfid, idx, xppm)
		ierr = nf90_put_var(cdfid, idnodes,nodes)

		qmin = 0D0
		! -- Time step using SSPRK3
		write(*,*) 'Time stepping..'
		ALLOCATE(soln(1:nxppm))
		t = 0D0	
		DO n=1,nsteps
          
		  IF(mod(n,2).eq.1) THEN
             oddstep = .true.
          ELSE
             oddstep = .false.
          END IF

		  CALL DGSWEEP(rhoq,nelems,dxel,nnodes,nodes,wghts,utild,nxppm,D,C,Cinv,dodghybrid,doposlimit,oddstep,dt)
		  t = t + dt
			! -- Fill in tvals
			tvals(n) = t
			soln = rhoq
			! -- Write solution at current time to output
			start = (/ n , 1 /)
			count = (/ 1 , NX /)
!			start2 = (/ n,1,1 /)
!			count2 = (/ 1,DGE,DGN/)
			ierr = NF90_PUT_VAR(cdfid,idrhoq,soln,start,count)
			
			! -- Keep track of minimum value of q
			qmin = DMIN1(MINVAL(rhoq),qmin)
			
		END DO
		! -- Write time values to output
		ierr = NF90_PUT_VAR(cdfid,idt,tvals)
		tf = etime(tend)-t0
		write(*,*) 'CPU time = ', tf
		thold(1) = t0
		thold(2) = tf

		write(*,*) 'Overall minimum value = ', qmin
		ierr = NF90_PUT_VAR(cdfid,idcput,thold)
		! -- Close netCDF file
		write(*,*) 'Closing netCDF file..'
		ierr = nf90_close(cdfid)
       	DEALLOCATE(tvals,xppm,inplt,soln)

    END SUBROUTINE DRIVER
END PROGRAM EXECUTE
