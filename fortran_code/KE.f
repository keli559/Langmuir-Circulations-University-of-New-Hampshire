c ------------------------------------------------------------------- c
c                      	     LC-IW code                               c
c ------------------------------------------------------------------- c
c                                                                     c
c  solve LC equations in u, psi, omega & T form               	    c
c                                                                     c
c  domain (non-dimensional):                                          c
c	y: [0,L] (fourier) x-wind direction                             c
c	z: [0,H] vertical direction                                     c
c	Z: [1,-1] (chebyshev) transformed vertical direction            c
c                                                                     c
c  time discretization:                                               c
c     s - present time                                                c
c     s+1 - new time                                                  c
c     s-1 - old time                                                  c
c                                                                     c
c  matrix structure:                                                  c
c     row  -->  varying y, Z = const                                  c
c     col  -->  varying Z, y = const                                  c
c	e.g. A(i,j) = A(Z,y): i [0,N], j [1,N]                          c
c                                                                     c
c  time structure:                                                    c
c	u(y,Z,t): t = 1, 2 or 3                                         c
c	  1: new time step (s+1)                                        c
c	  2: present time step (s)                                      c
c	  3: old time step (s-1)                                        c
c                                                                     c
c ------------------------------------------------------------------- c
       
	program KE

c	DECLARE AND DEFINE PARAMETERS

c	declare parameters variables
	integer*8 Ny, Nz, Mz, tstep, datagap, filtgap, energygap
	real*8 pi, H, L, LJ, PJ, QJ, Lat, eps, S,  
     ^		SD, xi, Res, ER, dt, A, B, C, D, flrad

c       declare plan variables
        integer*8  plany_f, plany_b, planz_f, planz_b, planzc_f,
     ^       planzc_b
       
c	parameters:
c	N  - number of modes (grid size)
c	M  - twice N
c	H  - nondimensional physical vertical domain, tilde H
C       L  - nondimensional physical horizontal domain, tilde L
c	LJ - Langmuir # in JB         	ER - Earth Rotation # (1 or 0)
c	PJ - Peclet # in JB		SD - Stokes-Drift damping constant
c	QJ - Gr*La^2 in JB		S  - Burger # 
c       Res- Renolds' #                 Lat- turbulence langmuir #
c       eps- Rossby #
c	dt - time step			xi - domain transformation factor
c	tstep - # of time steps		datagap - tsteps b/w data output 
C       A,B,C,D   -  B.C. value switches
	parameter (pi = 2.0d0*dacos(0.0d0))
	parameter (Ny = 1300, Nz = 128, Mz = Nz*2)
	parameter (H = 1.0d0, L = 2.0d0*pi)
	parameter (LJ = 0.0001d0, PJ = 4000.0d0, QJ = 0.15d0)
	parameter (SD = 6.0d0, xi = 2.0d0/H, S = 0.5d0/0.15d0, 
     ^		Res = 1.0d0, ER = 0.0d0, eps = 2.0d0/QJ, Lat = 1.0d0)
 	parameter (dt = 0.0005d0, tstep = 5, datagap = 1)
        parameter (A = 0.0d0, B = 0.0d0, C = 0.0d0, D = 0.0d0)
	parameter (flrad = 200.0d0, filtgap = 10)
        parameter (energygap = 1)


c	DECLARE VARIABLES FOR SUBROUTINES:

c	declare variables for dmn_setup
	real*8 k, y(1:Ny), Z(0:Nz), zz(0:Nz)

c	declare variables for base_st
	real*8 Us(0:Nz,1:Ny), dUs(0:Nz,1:Ny), Usy(0:Nz,1:Ny), 
     ^		Vs(0:Nz,1:Ny), 
     ^		dVs(0:Nz,1:Ny), Vsy(0:Nz,1:Ny), U_b(0:Nz,1:Ny),
     ^		T_b(0:Nz,1:Ny), psi_b(0:Nz,1:Ny), omega_b(0:Nz,1:Ny),
     ^		U_b_ZZ(0:Nz,1:Ny), T_b_ZZ(0:Nz,1:Ny)

c	declare variables for lc_IC
	real*8 u(0:Nz,1:Ny,1:3), T(0:Nz,1:Ny,1:3), psi(0:Nz,1:Ny,1:3),
     ^		omega(0:Nz,1:Ny,1:3)

c	declare variables for lc_BC
	complex*16 BC_u_hat(1:Ny,1:2), BC_psi_hat(1:Ny,1:2), 
     ^		BC_omega_hat(1:Ny,1:2), BC_T_hat(1:Ny,1:2)

c	declare variables for slv_setup
	real*8 lam(1:Ny), lam_T(1:Ny), gam(1:Ny), gam_T(1:Ny), 
     ^		eta(1:Ny)

c	declare variables for time loop
	integer*8 tt, BC_u_type(1:2), BC_o_type(1:2), BC_T_type(1:2),
     ^		BC_psi_type(1:2)

c	declare variables for write subroutines
	character*4 test*30, run_nm*2, fseq*3
	integer*8 r, ver

c	declare timing variables
	integer n, time
	real*8 before, after, cputime
	character*24 ctime, string

c	temp write
	real*8 temp(0:tstep, 1:Ny)

c	temp filter
	real*8 uftmp(0:Nz, 1:Ny), oftmp(0:Nz, 1:Ny), 
     ^		pftmp(0:Nz, 1:Ny), tftmp(0:Nz, 1:Ny)
	integer*8 ii, jj

c	declare norm variables
	real*8 nrm, E_u(1:tstep), E_psi(1:tstep), E_omega(1:tstep),
     ^		E_T(1:tstep)

c       declare variables for energy tracking
        integer enerstep, flag
        parameter (enerstep = tstep/energygap)
        real*8 v(0:Nz,1:Ny), w(0:Nz,1:Ny), uke, pke, poe, ps, pu,
     ^       bb, dh, dp, dpe2, dpe3, pke2,
     ^       uke_vec(1:enerstep), pke_vec(1:enerstep),
     ^       pke2_vec(1:enerstep),
     ^       poe_vec(1:enerstep), ps_vec(1:enerstep),
     ^       pu_vec(1:enerstep), bb_vec(1:enerstep), 
     ^       dh_vec(1:enerstep), dp_vec(1:enerstep),
     ^       dpe2_vec(1:enerstep), dpe3_vec(1:enerstep)


c	run name - !!! -- adjust run_nm length -- !!!
c	counter   1234567890123456789
	run_nm = "KE"

c	start time measurement
! 	before = second()
!       n = time()
!       string = ctime(n)
!       write(*,*) 'Start Time:'
!       write(*,*) string

c	BEGIN COMPUTATION:

c	setup domain
	call dmn_setup(Ny, Nz, L, k, xi, y, Z, zz)
	
c       initialize all the plans
        call plan_init(Ny, Nz, Mz, plany_f, plany_b,
     ^       planz_f, planz_b, planzc_f, planzc_b)

c	define basic state variables
	call base_st(Ny, Nz, Mz,  plany_f, plany_b, planz_f, planz_b, 
     ^       H, y, Z, zz, k, SD, xi, Us, dUs, Usy,Vs, dVs, Vsy, U_b,
     ^       T_b, psi_b, omega_b, U_b_ZZ, T_b_ZZ)


c ------------------------------------------------------------------- c
c --------------------------- DEFINE IC'S --------------------------- c
c ------------------------------------------------------------------- c
c 	use either lc_IC to specify IC's or use ld_ICs to load IC's 
c	file and comment out unused subroutine

cc	define IC's
c	call lc_IC(Ny, Nz, k, xi, U_b, T_b, psi_b, omega_b, u, T, psi, 
c     ^		omega)
c	load IC's
	call ld_ICs(Ny, Nz, u, omega, T, psi)
c ------------------------------------------------------------------- c
c	define BC's
	call lc_BC(Ny, Nz, plany_f, xi, A, B, C, D, Res, S, eps,   
     ^		BC_u_hat, BC_psi_hat, BC_omega_hat, BC_T_hat)
        BC_u_type(1) = 1
        BC_u_type(2) = 1
        BC_o_type(1) = 0
        BC_o_type(2) = 0
        BC_T_type(1) = 1
        BC_T_type(2) = 1
        BC_psi_type(1) = 0
        BC_psi_type(2) = 0

c	write parameters being used to screen
	call params(Ny, Nz, H, L, xi, LJ, PJ, QJ, S, ER, SD, Res,
     ^		eps, Lat, dt, tstep, Us, Vs, A, B, C, D,
     ^		datagap, energygap, run_nm, BC_u_hat, BC_psi_hat, BC_omega_hat, 
     ^		BC_T_hat, BC_u_type, BC_psi_type, BC_o_type, BC_T_type)

c	verify correct parameters
	call wait()

c	setup constants for solving 
	call slv_setup(Ny, k, xi, dt, LJ, PJ, lam, lam_T, gam, gam_T,
     ^		eta)

c	begin time loop
	r = 1
        flag = 1
	do tt = 1,tstep

c       define BC type: 0 refers to Direchlet BC, 1 refers to Norman BC 
        BC_u_type(1) = 1
        BC_u_type(2) = 1
        BC_o_type(1) = 0
        BC_o_type(2) = 0
        BC_T_type(1) = 1
        BC_T_type(2) = 1
        BC_psi_type(1) = 0
        BC_psi_type(2) = 0

c	  compute x velocity
c	  BC_u_type(1) = 1
c	  BC_u_type(2) = 1
	  call slv_u(Ny, Nz, Mz, plany_f, plany_b, planz_f, planz_b,
     ^         planzc_f, planzc_b, LJ, eps, ER, Lat, Vs, k, xi,  
     ^		lam, gam, BC_u_hat,BC_u_type, psi, u, U_b_ZZ)
c	  compute x vorticity
c	  BC_o_type(1) = 0
c	  BC_o_type(2) = 0

	  call slv_omega(Ny, Nz, Mz, plany_f, plany_b, planz_f, planz_b,
     ^         planzc_f, planzc_b, LJ, ER, Lat, eps, QJ, dUs,  
     ^		Usy, Vs, dVs, Vsy, k, xi, lam, gam, 
     ^		BC_omega_hat, BC_o_type, u, psi, T, omega)

c	  compute temperature
c 	  BC_T_type(1) = 0
c 	  BC_T_type(2) = 0
	  call slv_T(Ny, Nz, Mz, plany_f, plany_b, planz_f, planz_b,
     ^         planzc_f, planzc_b, PJ, eps, Lat, S, Vs, k, xi,  
     ^		lam_T, gam_T, BC_T_hat, BC_T_type, psi, T, T_b_ZZ)
c	  compute streamfunction
c	  BC_psi_type(1) = 0
c	  BC_psi_type(2) = 0
	
	  call slv_psi(Ny, Nz, plany_f, plany_b, planzc_f, planzc_b, 
     ^         xi, eta, BC_psi_hat, BC_psi_type, omega, psi)

c______________________________________________________________________
c	  if tt (time) == filtgap --> filter solutions once 
c	  if (mod(tt, filtgap) .eq. 0) then
c	  if (1 .eq. 0) then
c	     call flmanip(u(0:Nz, 1:Ny, 2), omega(0:Nz, 1:Ny, 2),    
c     ^		psi(0:Nz, 1:Ny, 2), T(0:Nz, 1:Ny, 2), Ny, Nz, L, 0.0d0,    
c     ^		0.0d0, flrad, uftmp, oftmp, pftmp, tftmp)
c	     do ii = 0, Nz
c		do jj = 1, Ny
c		   u(ii, jj, 2) = uftmp(ii, jj)
c		   omega(ii, jj, 2) = oftmp(ii, jj)
c		   psi(ii, jj, 2) = pftmp(ii, jj)
c		   T(ii, jj, 2) = tftmp(ii, jj)
c		end do
c	     end do
c	     call flmanip(u(0:Nz, 1:Ny, 3), omega(0:Nz, 1:Ny, 3),    
c     ^		psi(0:Nz, 1:Ny, 3), T(0:Nz, 1:Ny, 3), Ny, Nz, L, 0.0d0,    
c     ^		0.0d0, flrad, uftmp, oftmp, pftmp, tftmp)
c	     do ii = 0, Nz
c		do jj = 1, Ny
c		   u(ii, jj, 3) = uftmp(ii, jj)
c		   omega(ii, jj, 3) = oftmp(ii, jj)
c		   psi(ii, jj, 3) = pftmp(ii, jj)
c		   T(ii, jj, 3) = tftmp(ii, jj)
c		end do
c	     end do
c	  end if
c______________________________________________________________________
	
c	  if tt (time) == datagap --> write solutions to file 
	  if (mod(tt, datagap) .eq. 0) then
	    call num2str(r, fseq)
	    ver = 1
	    test = fseq// ".dat"
	    call writesoln(Ny, Nz, u, psi, omega, T, ver, test)
	    write(*,*) 'time step:'
	    write(*,*) tt
	    r = r+1
	  end if

c         if tt (time) == energygap --> save data to vector
          if (mod(tt, energygap) .eq. 0) then
c           calculate the energy
c           compute the cross-secional velocity
            ver=1
            call velocity(Ny, Nz, Mz, ver, plany_f, plany_b, 
     ^           planz_f, planz_b, xi, k, psi, v, w)
c           compute the y-z-averaged energy by Parseval's identity
            call comp_energy(Ny, Nz, Mz, ver, plany_f, plany_b, 
     ^           planz_f, planz_b, zz, k, xi, PJ, LJ, QJ, L, H, dUs, 
     ^           u, v, w, T, psi, omega, uke, pke, pke2, poe, ps, pu,
     ^           bb, dh, dp, dpe2, dpe3)

c           store to another variable for future output
               uke_vec(flag) = uke
               pke_vec(flag) = pke
               poe_vec(flag) = poe
               ps_vec(flag) = ps
               pu_vec(flag) = pu
               bb_vec(flag) = bb
               dp_vec(flag) = dp
               dh_vec(flag) = dh
               dpe2_vec(flag) = dpe2
               dpe3_vec(flag) = dpe3
               pke2_vec(flag) = pke2

            flag = flag + 1

          end if

c	  store T(y,Z = 0)
	  do i = 1,Ny
	    temp(tt, i) = T(Nz/2,i,2)
	  end do

c	  advance time step
	  do i = 0,Nz
	    do j = 1,Ny
c		x velocity
		u(i,j,3) = u(i,j,2)
		u(i,j,2) = u(i,j,1)
c		x vorticity
		omega(i,j,3) = omega(i,j,2)
		omega(i,j,2) = omega(i,j,1)
c		temperature
		T(i,j,3) = T(i,j,2)
		T(i,j,2) = T(i,j,1)
c		streamfunction
		psi(i,j,3) = psi(i,j,2)
		psi(i,j,2) = psi(i,j,1)
 	    end do
	  end do

! c	  determine norm of each variable
! 	  E_u(tt) = nrm(Ny, Nz, u)
! 	  E_psi(tt) = nrm(Ny, Nz, psi)
! 	  E_omega(tt) = nrm(Ny, Nz, omega)
! 	  E_T(tt) = nrm(Ny, Nz, T)

	end do

! c	stop time measurement
! 	after    = second()          
! c	compute elapsed time (secs)
!       cputime  = after - before  
! c	print time results to screen
! 	write(*,*) 'FORTRAN SOLUTION TIMING'
! 	write(*,*) 'Elapsed time:'
! 	write(*,'(1x,f15.2)') cputime
! 	write(*,*)
! 	write(*,*) 'End time:'
!       n = time()
!       string = ctime(n)
!       write(*,*) string

c	write old time step
	ver = 3
	test = run_nm // ".dat"
	call writesoln(Ny, Nz, u, psi, omega, T, ver, test)

! c	write E_u, E_psi, E_omega, E_T to file
! 	test = "E_u_" // run_nm // ".dat"
! 	call writevec(1, tstep, E_u, test)
! 	test = "E_psi_" // run_nm // ".dat"
! 	call writevec(1, tstep, E_psi, test)
! 	test = "E_omega_" // run_nm // ".dat"
! 	call writevec(1, tstep, E_omega, test)
! 	test = "E_T_" // run_nm // ".dat"
! 	call writevec(1, tstep, E_T, test)

c	write temp to file
c 	test = "T_Z0_" // run_nm // ".dat"
c 	call writemtrx(Ny, tstep, temp, test)

c	write y and Z
	call writevec(1, Ny, y, "y.dat")
	call writevec(0, Nz, zz, "z.dat")

c       write integrated energy varied in t
	call writevec(1, enerstep, uke_vec, "ukee.dat")
        call writevec(1, enerstep, pke_vec, "pke1.dat")
        call writevec(1, enerstep, poe_vec, "poee.dat")
        call writevec(1, enerstep, ps_vec, "psps.dat")
        call writevec(1, enerstep, pu_vec, "pupu.dat")
        call writevec(1, enerstep, pv_vec, "pvpv.dat")
        call writevec(1, enerstep, bb_vec, "bbbb.dat")
        call writevec(1, enerstep, dh_vec, "dhdh.dat")
        call writevec(1, enerstep, dp_vec, "dpdp.dat")
        call writevec(1, enerstep, dpe2_vec, "dpe2.dat")
        call writevec(1, enerstep, dpe3_vec, "dpe3.dat")
        call writevec(1, enerstep, pke2_vec, "pke2.dat")


! c	write physical vertical domain, z
! 	test = "z.dat"
! 	call writevec(0, Nz, zz, test)

	end

c *********************************************************************
c **************************** SUB ROUTINES ***************************
c *********************************************************************

c ******************************* SETUP *******************************
	subroutine dmn_setup(Ny, Nz, L, k, xi, y, Z, zz)
c	setup domain
c	input: 
c	  Ny, Nz, L
c	output:
c	  k, y, Z, zz

c	declare inputs/output variables
	integer*8 Ny, Nz
	real*8 L, k, xi, y(1:Ny), Z(0:Nz), zz(0:Nz)

c	declare subroutine specific variables
	real*8 pi, h

c	define pi
	pi = 2.0d0*dacos(0.0d0)

c	domain periodicity transformation constant: 2*pi -->  L
	k = (2.0d0*pi)/L

c	grid spacing
	h = (2.0d0*pi/dble(Ny))/k

c	define y - fourier grid
	do i = 1,Ny
	  y(i) = (dble(i)-1.0d0)*h
	end do

c	define z - cheb grid
	do i = 0,Nz
	  Z(i) = dcos(pi*dble(i)/dble(Nz))
	  zz(i) = (Z(i)-1.0d0)/xi
	end do

	end

	
c **************************** BASIC STATE ****************************
	subroutine base_st(Ny, Nz, Mz, plany_f, plany_b, planz_f, planz_b, 
     ^     H, y, Z, zz, k, SD, xi, Us, dUs, Usy,Vs, dVs, Vsy, U_b,
     ^     T_b, psi_b, omega_b, U_b_ZZ, T_b_ZZ)
c	setup basic state (BS)
c	input: 
c	  Ny, Nz, Mz, y, Z, k, SD, xi
c	output:
c	  Us, Vs, U_b, T_b, psi_b

c	declare inputs/output variables
	integer*8 Ny, Nz, Mz, ii, jj, plany_f, plany_b, planz_f, planz_b
	real*8 H, y(1:Ny), Z(0:Nz), zz(0:Nz), k, SD, xi, Us(0:Nz,1:Ny),
     ^		dUs(0:Nz,1:Ny), Usy(0:Nz,1:Ny), Vs(0:Nz,1:Ny),
     ^		dVs(0:Nz,1:Ny), Vsy(0:Nz,1:Ny), U_b(0:Nz,1:Ny),  
     ^		T_b(0:Nz,1:Ny), psi_b(0:Nz,1:Ny), omega_b(0:Nz,1:Ny),
     ^		U_b_ZZ(0:Nz,1:Ny), T_b_ZZ(0:Nz,1:Ny)
c       M --- perturbation amplitude
	real*8 M, upr(1:Ny), utemp(1:Ny)

c	declare subroutine specific variables
	real*8 pi, psi_y(0:Nz,1:Ny), psi_yy(0:Nz,1:Ny), 
     ^		psi_Z(0:Nz,1:Ny), psi_ZZ(0:Nz,1:Ny), 
     ^		temp1(0:Nz,1:Ny), A, c, T0, U0, B, sig, y1, y2, Z1

c	define pi
	pi = 2.0d0*dacos(0.0d0)
c	define perturbation amplitude M
	M = 0.005d0

c	define Stokes-Drift
	do i = 0,Nz
	  do j = 1,Ny
c	    Us(i,j) = dexp(SD*zz(i))
c	    dUs(i,j) = (SD/xi)*Us(i,j)
            Usy(i,j) = 0.0d0
            Vs(i,j) = 0.0d0
	    dVs(i,j) = 0.0d0
            Vsy(i,j) = 0.0d0
	    Us(i,j) = (Z(i) + 1.0d0)/xi
	    dUs(i,j) = 1.0d0/xi
c	    Us(i,j) = 0.0d0
c	    dUs(i,j) = 0.0d0
	  end do
	end do
c	difine perturbation upr(y)
	do jj = 1,Ny
	   upr(jj) = 0.0d0
	   utemp(jj) = 0.0d0
	end do

	ii = 15
	   do jj = 1,Ny
	      utemp(jj) = M*(dcos(dble(ii)*k*y(jj)) 
     ^		+dsin(dble(ii)*k*y(jj)))
	   end do
	   do jj = 1, Ny
	      upr(jj)=utemp(jj)+upr(jj)
	   end do



c	define U_b and T_b parameters
	A = 0.0d0
	c = 5.5d0
	T0 = 0.0d0
	U0 = 0.0d0
	offset = 1.0d0

c	define x-velocity, temperature and streamfunction BS
	do i = 0,Nz
	  do j = 1,Ny
	    U_b(i,j) = 0.0d0
c((1.0d0/xi)*(Z(i) + 1.0d0))
c     ^		*A*(dtanh(c*pi*(Z(i)-U0)) + offset)
	    T_b(i,j) = 0.0d0 
c	 	A*(dtanh(c*pi*(Z(i)-U0)) + offset)
c		A*(dtanh(c*pi*(Z(i)-T0)) + offset)
	    psi_b(i,j) = upr(j)*sin(pi*zz(i))
c            M*(((Z(i) - 1.0d0)/2.0d0)**4.0d0 +  
c     ^		2.0d0*((Z(i) - 1.0d0)/2.0d0)**3.0d0 -  
c     ^		((Z(i) - 1.0d0)/2.0d0)) + upr(j)
c                0.0d0
c		dsin(5.0d0*k*y(j))*dcos(pi/2.0d0*Z(i)) 
c               0.1d0*dsin(2.0d0*k*y(j))
c     ^		*dsin(2.5d0*pi*((Z(i)-1.0d0)/2.0d0))
c     ^		*A*(dtanh(c*pi*(Z(i)-U0)) + offset)
	  end do
	end do

c	compute U_b_ZZ and T_b_ZZ
	call diffZ(Ny, Nz, Mz, planz_f, planz_b, U_b, temp1)
	call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp1, U_b_ZZ)

	call diffZ(Ny, Nz, Mz, planz_f, planz_b, T_b, temp1)
	call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp1, T_b_ZZ)

c	compute  psi_yy and psi_ZZ
	call diffy(Ny, Nz, k, plany_f, plany_b, psi_b, psi_y)
	call diffy(Ny, Nz, k, plany_f, plany_b, psi_y, psi_yy)

	call diffZ(Ny, Nz, Mz, planz_f, planz_b, psi_b, psi_Z)
	call diffZ(Ny, Nz, Mz, planz_f, planz_b, psi_Z, psi_ZZ)
     
c	compute omega_b
	do i = 0,Nz
	  do j = 1,Ny
	    omega_b(i,j) = -(psi_yy(i,j) + xi*xi*psi_ZZ(i,j))
	  end do
	end do

	end


c ******************************* IC's ********************************
	subroutine lc_IC(Ny, Nz, k, xi, U_b, T_b, psi_b, omega_b, u, T,  
     ^		psi, omega)
c	setup IC's
c	input: 
c	  Ny, Nz, k, xi, U_b, T_b, psi_b
c	output:
c	  u, T, psi, omega

c	declare inputs/output variables
	integer*8 Ny, Nz
	real*8 k, xi, U_b(0:Nz,1:Ny), T_b(0:Nz,1:Ny), psi_b(0:Nz,1:Ny), 
     ^		omega_b(0:Nz,1:Ny), u(0:Nz,1:Ny,1:3),  
     ^		T(0:Nz,1:Ny,1:3), psi(0:Nz,1:Ny,1:3), 
     ^		omega(0:Nz,1:Ny,1:3)

c	declare subroutine specific variables
	integer*8 tt

c	define x-velocity, temperature and streamfunction IC
	do tt = 2,3
	  do i = 0,Nz
	    do j = 1,Ny
	      u(i,j,tt) = U_b(i,j)
	      T(i,j,tt) = T_b(i,j)
	      psi(i,j,tt) = psi_b(i,j)
	      omega(i,j,tt) = omega_b(i,j)
	    end do
	  end do
	end do

	do i = 0,Nz
	  do j = 1,Ny
	    u(i,j,1) = 0.0d0
	    T(i,j,1) = 0.0d0
	    psi(i,j,1) = 0.0d0
	    omega(i,j,1) = 0.0d0
	  end do
	end do
	end


c ************************* LOAD BASIC STATE **************************
	subroutine ld_ICs(Ny, Nz, u, omega, T, psi)
c	setup basic state (BS)
c	input: 
c	  Ny, Nz, Mz, y, Z, k, SD, xi
c	output:
c	  dUs, U_b, T_b, psi_b

c	declare inputs/output variables
	integer*8 Ny, Nz
	integer c
	real*8 u(0:Nz,1:Ny,1:3), T(0:Nz,1:Ny,1:3), psi(0:Nz,1:Ny,1:3), 
     ^		omega(0:Nz,1:Ny,1:3)
	
c	declare subroutine specific variables
	real*8 temp(0:Nz,1:Ny)
	character flname13*13, flname15*15, flname17*17, flname19*19, 
     ^		flname21*21

c	load u
	flname13 = "../data/u.dat"
	call readfle(Ny, Nz, temp, flname13)
	do i = 0,Nz
	  do j = 1,Ny
	    u(i,j,2) = temp(i,j)
	  end do 
	end do

	flname17 = "../data/u_old.dat"
	call readfle(Ny, Nz, temp, flname17)
	do i = 0,Nz
	  do j = 1,Ny
	    u(i,j,3) = temp(i,j)
	  end do 
	end do

c	load omega
	flname17 = "../data/omega.dat"
	call readfle(Ny, Nz, temp, flname17)
	do i = 0,Nz
	  do j = 1,Ny
	    omega(i,j,2) = temp(i,j)
	  end do 
	end do	

	flname21 = "../data/omega_old.dat"
	call readfle(Ny, Nz, temp, flname21)
	do i = 0,Nz
	  do j = 1,Ny
	    omega(i,j,3) = temp(i,j)
	  end do 
	end do

c	load T
	flname13 = "../data/T.dat"
	call readfle(Ny, Nz, temp, flname13)
	do i = 0,Nz
	  do j = 1,Ny
	    T(i,j,2) = temp(i,j)
	  end do 
	end do

	flname17 = "../data/T_old.dat"
	call readfle(Ny, Nz, temp, flname17)
	do i = 0,Nz
	  do j = 1,Ny
	    T(i,j,3) = temp(i,j)
	  end do 
	end do

c	load psi
	flname15 = "../data/psi.dat"
	call readfle(Ny, Nz, temp, flname15)
	do i = 0,Nz
	  do j = 1,Ny
	    psi(i,j,2) = temp(i,j)
	  end do 
	end do
	
	flname19 = "../data/psi_old.dat"
	call readfle(Ny, Nz, temp, flname19)
	do i = 0,Nz
	  do j = 1,Ny
	    psi(i,j,3) = temp(i,j)
	  end do 
	end do

c	set s+1 t.s. to zero
	do i = 0,Nz
	  do j = 1,Ny
	    u(i,j,1) = 0.0d0
	    T(i,j,1) = 0.0d0
	    psi(i,j,1) = 0.0d0
	    omega(i,j,1) = 0.0d0
	  end do
	end do
	end


c ******************************* BC's ********************************
	subroutine lc_BC(Ny, Nz, plany_f, xi, A, B, C, D, Res, S, eps,   
     ^		BC_u_hat, BC_psi_hat, BC_omega_hat, BC_T_hat)
c	setup BC's
c	input: 
c	  Ny, Nz, xi, Res, S, eps, A, B, C, D
c	output:
c	  BC_u_hat, BC_psi_hat, BC_omega_hat, BC_T_hat
c	structure
c	  BC_u_hat(:,i): i = 1 or 2
c	    1: Z = 1 BC
c	    2: Z = -1 BC

c	declare inputs/output variables
	integer*8 Ny, Nz, plany_f 
	real*8 xi, Res, S, eps
	complex*16 BC_u_hat(1:Ny,1:2), BC_psi_hat(1:Ny,1:2), 
     ^		BC_omega_hat(1:Ny,1:2), BC_T_hat(1:Ny,1:2)

c	declare subroutine specific variables
	real*8 BC_u(1:Ny), BC_psi(1:Ny), BC_omega(1:Ny), BC_T(1:Ny),  
     ^		A, B, C, D
	complex*16 temp1(1:Ny), temp2(1:Ny), temp3(1:Ny), temp4(1:Ny)

c *******************
c ***** Z = 1 ******
c********************
c	define BC's
	do i = 1,Ny
	  BC_u(i) = A*Res/xi
	  BC_psi(i) = 0.0d0
	  BC_omega(i) = 0.0d0
	  BC_T(i) = B/xi  
c     ^		-S/(2.0d0*eps)
	end do

c	define y transformed BC's
	call fft_f(Ny, plany_f, BC_u, temp1)
	call fft_f(Ny, plany_f, BC_psi, temp2)
	call fft_f(Ny, plany_f, BC_omega, temp3)
	call fft_f(Ny, plany_f, BC_T, temp4)

c	correct normalization factor
	do i = 1,Ny
	  BC_u_hat(i,1) = temp1(i)/dcmplx(dble(Ny), 0.0d0)
	  BC_psi_hat(i,1) = temp2(i)/dcmplx(dble(Ny), 0.0d0)
	  BC_omega_hat(i,1) = temp3(i)/dcmplx(dble(Ny), 0.0d0)
	  BC_T_hat(i,1) = temp4(i)/dcmplx(dble(Ny), 0.0d0)
	end do
c *******************
c ***** Z = -1 ******
c********************

c	define BC's
	do i = 1,Ny
	  BC_u(i) = A*Res/xi
	  BC_psi(i) = D*(0.5d0*Res*(1.0d0-A)/xi)
	  BC_omega(i) = 0.0d0
	  BC_T(i) = C/xi  
c     ^		-0.5d0*S/eps
	end do

c	define y transformed BC's
	call fft_f(Ny, plany_f, BC_u, temp1)
	call fft_f(Ny, plany_f, BC_psi, temp2)
	call fft_f(Ny, plany_f, BC_omega, temp3)
	call fft_f(Ny, plany_f, BC_T, temp4)

c	correct normalization factor
	do i = 1,Ny
	  BC_u_hat(i,2) = temp1(i)/dcmplx(dble(Ny), 0.0d0)
	  BC_psi_hat(i,2) = temp2(i)/dcmplx(dble(Ny), 0.0d0)
	  BC_omega_hat(i,2) = temp3(i)/dcmplx(dble(Ny), 0.0d0)
	  BC_T_hat(i,2) = temp4(i)/dcmplx(dble(Ny), 0.0d0)
	end do


	end


c *************************** SOLVER SETUP ****************************
	subroutine slv_setup(Ny, k, xi, dt, LJ, PJ, lam, lam_T, gam, 
     ^		gam_T, eta)
c	solver setup
c	input: 
c	  Ny, k, xi, dt, LJ, PJ
c	output:
c	  lam, lam_T, gam, gam_t, eta

c	declare inputs/output variables
	integer*8 Ny
	real*8 k, xi, dt, LJ, PJ, lam(1:Ny), lam_T(1:Ny), gam(1:Ny), 
     ^		gam_T(1:Ny), eta(1:Ny)

c	declare subroutine specific variables
	real*8 n_w(1:Ny)

c	define wavenumber vector - order (0:N/2-1 N/2 -N/2+1:-1)
	do i = 1,Ny/2+1
	  n_w(i) = dble(i) - 1.0d0
	end do
	j = 1
	do i = Ny/2+2,Ny
	  n_w(i) = dble(-Ny)/2.0d0 + dble(j)
	  j = j+1
	end do

	do i = 1,Ny
c	  helmholtz coefficient (new time step s+1)
	  lam(i) = (((n_w(i)*k)**2.0d0)/(xi*xi))  
     ^		+ (2.0d0/(xi*xi*dt*LJ))

c	  helmholtz coefficient (new time step s+1) - temperature
	  lam_T(i) = (((n_w(i)*k)**2.0d0)/(xi*xi)) 
     ^		+ ((2.0d0*PJ)/(xi*xi*dt))

c	  pseudo helmholz coefficient (present time step s)
	  gam(i) = (((n_w(i)*k)**2.0d0)/(xi*xi)) 
     ^		- (2.0d0/(xi*xi*dt*LJ))

c	  pseudo helmholz coefficient (present time step s) - temperature
	  gam_T(i) = (((n_w(i)*k)**2.0d0)/(xi*xi)) 
     ^		- ((2.0d0*PJ)/(xi*xi*dt))

c	  psi helmholz coefficient (present time step s+1)
	  eta(i) = ((n_w(i)*k)**2.0d0)/(xi*xi)
	end do

	end


c ***************************** u SOLVER ******************************
	subroutine slv_u(Ny, Nz, Mz, plany_f, plany_b, planz_f, planz_b,
     ^     planzc_f, planzc_b, LJ, eps, ER, Lat, Vs, k, xi,  
     ^		lam, gam, BC_hat,BC_type, psi, u, U_b_ZZ)
c	u solver 
c	input: 
c	  Ny, Nz, Mz, LJ, eps, ER, Lat, Vs, k, xi, lam, 
c         gam, BC_hat, BC_type, psi, u
c	  
c	output:
c	  u

c	declare inputs/output variables
	integer*8 Ny, Nz, Mz, BC_type(1:2),  plany_f, plany_b,
     ^       planz_f, planz_b, planzc_f, planzc_b
	real*8 LJ, eps, ER, Lat, Vs(0:Nz, 1:Ny), k, xi,   
     ^		u(0:Nz,1:Ny,1:3), 
     ^		psi(0:Nz,1:Ny,1:3), U_b_ZZ(0:Nz,1:Ny), lam(1:Ny), 
     ^		gam(1:Ny) 
	complex*16 BC_hat(1:Ny,1:2)
 
c	declare subroutine specific variables
	integer*8 tt
	real*8 u_y(0:Nz,1:Ny,2:3), u_Z(0:Nz,1:Ny,2:3), u_ZZ(0:Nz,1:Ny), 
     ^		psi_y(0:Nz,1:Ny,2:3), psi_Z(0:Nz,1:Ny,2:3), 
     ^		temp1(0:Nz,1:Ny), temp2(0:Nz,1:Ny), temp3(0:Nz,1:Ny)

	complex*16 u_hat(0:Nz,1:Ny,1:2), u_rhs(0:Nz,1:Ny),
     ^		slv_temp(0:Nz), rhs_temp(0:Nz), bc_temp(1:2),
     ^		temp1c(0:Nz,1:Ny), F(0:Nz,1:Ny,2:3)

c 	compute derivatives for present and old time step
	do tt = 2,3
	  do i = 0,Nz
	    do j = 1,Ny
	      temp1(i,j) = u(i,j,tt)
	      temp2(i,j) = psi(i,j,tt)
	    end do
	  end do

c	  compute u_y
	  call diffy(Ny, Nz, k, plany_f, plany_b, temp1, temp3)
	  do i = 0,Nz
	    do j = 1,Ny
	      u_y(i,j,tt) = temp3(i,j)
	    end do
	  end do

c	  compute u_z
	  call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp1, temp3)
	  do i = 0,Nz
	    do j = 1,Ny
	      u_Z(i,j,tt) = temp3(i,j)
	    end do
	  end do

c	  compute psi_y
	  call diffy(Ny, Nz, k, plany_f, plany_b, temp2, temp3)
	  do i = 0,Nz
	    do j = 1,Ny
	      psi_y(i,j,tt) = temp3(i,j)
	    end do
	  end do

c	  compute psi_Z
	  call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp2, temp3)
	  do i = 0,Nz
	    do j = 1,Ny
	      psi_Z(i,j,tt) = temp3(i,j)
	    end do
	  end do

c	  compute non-linear terms (F)
	  do i = 0,Nz
	    do j = 1,Ny
	      temp1(i,j) = xi*(psi_y(i,j,tt)*u_Z(i,j,tt)  
     ^		- psi_Z(i,j,tt)*u_y(i,j,tt)) + psi_y(i,j,tt)
     ^          + xi * psi_z(i,j,tt) / eps
c     ^		- ((xi*xi)*LJ)*U_b_ZZ(i,j)
	    end do
	  end do
C         compute earth rotating terms (ER)
c	  do i = 0,Nz
c	    do j = 1,Ny
c	      temp1(i,j) = temp1(i,j)  
c     ^		+ER*(xi*psi_Z(i,j,tt)+(1.0d0/(Lat*Lat))*Vs(i,j))  
c     ^		-(eps/(Lat*Lat))*Vs(i,j)*u_y(i,j,tt)
c	    end do
c	  end do

	  call trans_y(Ny, Nz, plany_f, temp1, temp1c)
	  do i = 0,Nz
	    do j = 1,Ny
	      F(i,j,tt) = temp1c(i,j)
	    end do 
	  end do
	end do


c	compute u_hat at present time step (s)
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = u(i,j,2)
	  end do
	end do 
	call trans_y(Ny, Nz, plany_f, temp1, temp1c)
	do i = 0,Nz
	  do j = 1,Ny
	    u_hat(i,j,2) = temp1c(i,j)
	  end do
	end do 
	
c	compute first term in u_rhs
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = u_Z(i,j,2)
	  end do
	end do
	call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp1, u_ZZ)
	call trans_y(Ny, Nz, plany_f, u_ZZ, temp1c)

c	compute sum of four terms in u_rhs
	do i = 0,Nz
	  do j = 1,Ny
	    u_rhs(i,j) = - temp1c(i,j)  
     ^		+ dcmplx(gam(j), 0.0d0)*u_hat(i,j,2) 
     ^		- dcmplx(1.0d0/(xi*xi*LJ), 0.0d0)
     ^		*(dcmplx(3.0d0, 0.0d0)*F(i,j,2) - F(i,j,3))
	  end do
	end do

c     solve helmholtz equation at each y location
	do j = 1,Ny
	  do i = 0,Nz
	    rhs_temp(i) = u_rhs(i,j)
	  end do

	  do i = 1,2
	   bc_temp(i) = BC_hat(j,i)
	  end do

	  call slv_helm2(Nz, planzc_f, planzc_b, lam(j), rhs_temp,
     ^         bc_temp, BC_type, slv_temp)

	  do i = 0,Nz
	  u_hat(i,j,1) = slv_temp(i)
	  end do
	end do

c	transform to get new u value
	do i = 0,Nz
	  do j = 1,Ny
	    temp1c(i,j) = u_hat(i,j,1)
	  end do
	end do
	call itrans_y(Ny, Nz, plany_b, temp1c, temp1)
	do i = 0,Nz
	  do j = 1,Ny
	    u(i,j,1) = temp1(i,j)
	  end do
	end do

	end


c *************************** omega SOLVER ****************************
	subroutine slv_omega(Ny, Nz, Mz, plany_f, plany_b, planz_f, 
     ^     planz_b, planzc_f, planzc_b, LJ, ER, Lat, eps, QJ, dUs,  
     ^		Usy, Vs, dVs, Vsy, k, xi, lam, gam, 
     ^		BC_hat, BC_type, u, psi, T, omega)
c	omega solver 
c	input: 
c	  Ny, Nz, Mz, LJ, ER, eps, QJ, k, xi, lam, gam, BC_hat,   
C         BC_type, u, psi, T, omega
c	  
c	  
c	output:
c	  omega

c	declare inputs/output variables
	integer*8 Ny, Nz, Mz, BC_type(1:2), plany_f, 
     ^     plany_b, planz_f, planz_b, planzc_f, planzc_b
	real*8 LJ, ER, Lat, eps, QJ, k, xi,
     ^		u(0:Nz,1:Ny,1:3), psi(0:Nz,1:Ny,1:3),
     ^		omega(0:Nz,1:Ny,1:3), T(0:Nz,1:Ny,1:3), 
     ^		Us(0:Nz,1:Ny), dUs(0:Nz,1:Ny), Vs(0:Nz,1:Ny),  
     ^		dVs(0:Nz,1:Ny), Usy(0:Nz,1:Ny), Vsy(0:Nz,1:Ny),  
     ^		lam(1:Ny), gam(1:Ny)
	complex*16 BC_hat(1:Ny,1:2)
 
c	declare subroutine specific variables
	integer*8 tt
	real*8 omega_y(0:Nz,1:Ny,2:3), omega_Z(0:Nz,1:Ny,2:3), 
     ^		omega_ZZ(0:Nz,1:Ny), psi_y(0:Nz,1:Ny,2:3), 
     ^		psi_Z(0:Nz,1:Ny,2:3), u_y(0:Nz,1:Ny,2:3), 
     ^		u_Z(0:Nz, 1:Ny, 2:3), 
     ^		T_y(0:Nz,1:Ny,2:3), temp1(0:Nz,1:Ny), 
     ^		temp2(0:Nz,1:Ny),temp3(0:Nz,1:Ny), temp4(0:Nz,1:Ny),
     ^		temp5(0:Nz,1:Ny)
	complex*16 omega_hat(0:Nz,1:Ny,1:2), G(0:Nz,1:Ny,2:3), 
     ^		omega_rhs(0:Nz,1:Ny), slv_temp(0:Nz), rhs_temp(0:Nz), 
     ^		bc_temp(1:2), temp1c(0:Nz,1:Ny)

c 	compute derivatives for present and old time step
	do tt = 2,3
	  do i = 0,Nz
	    do j = 1,Ny
	      temp1(i,j) = omega(i,j,tt)
	      temp2(i,j) = psi(i,j,tt)
	      temp3(i,j) = u(i,j,tt)
	      temp4(i,j) = T(i,j,tt)
	    end do
	  end do

c	  compute omega_y
	  call diffy(Ny, Nz, k, plany_f, plany_b, temp1, temp5)
	  do i = 0,Nz
	    do j = 1,Ny
	      omega_y(i,j,tt) = temp5(i,j)
	    end do
	  end do

c	  compute omega_z
	  call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp1, temp5)
	  do i = 0,Nz
	    do j = 1,Ny
	      omega_Z(i,j,tt) = temp5(i,j)
	    end do
	  end do

c	  compute psi_y
	  call diffy(Ny, Nz, k, plany_f, plany_b, temp2, temp5)
	  do i = 0,Nz
	    do j = 1,Ny
	      psi_y(i,j,tt) = temp5(i,j)
	    end do
	  end do

c	  compute psi_Z
	  call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp2, temp5)
	  do i = 0,Nz
	    do j = 1,Ny
	      psi_Z(i,j,tt) = temp5(i,j)
	    end do
	  end do

c	  compute u_y
	  call diffy(Ny, Nz, k, plany_f, plany_b, temp3, temp5)
	  do i = 0,Nz
	    do j = 1,Ny
	      u_y(i,j,tt) = temp5(i,j)
	    end do
	  end do

c	  compute u_z
	  call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp3, temp5)
	  do i = 0,Nz
	    do j = 1,Ny
	      u_Z(i,j,tt) = temp5(i,j)
	    end do
	  end do

c	  compute T_y
	  call diffy(Ny, Nz, k, plany_f, plany_b, temp4, temp5)
	  do i = 0,Nz
	    do j = 1,Ny
	      T_y(i,j,tt) = temp5(i,j)
	    end do
	  end do

c	  compute non-linear terms (G)
	  do i = 0,Nz
	    do j = 1,Ny
	      temp1(i,j) =  
     ^		xi*(psi_y(i,j,tt)*omega_Z(i,j,tt) 
     ^		- psi_Z(i,j,tt)*omega_y(i,j,tt)) 
     ^		+ QJ*T_y(i,j,tt) 
c     ^		- ER*xi*(u_Z(i,j,tt)+dUs(i,j)/(Lat*Lat)) 
c     ^		+ (eps/(Lat*Lat))*( 
     ^		- xi*dUs(i,j)*u_y(i,j,tt)
     ^          + (xi*u_Z(i,j,tt))/eps
c     ^		-xi*Usy(i,j)*u_Z(i,j,tt)
c     ^		-Vsy(i,j)*omega(i,j,tt) 
c     ^		-Vs(i,j)*omega_y(i,j,tt) 
c     ^		)
	    end do
	  end do

	  call trans_y(Ny, Nz, plany_f, temp1, temp1c)

	  do i = 0,Nz
	    do j = 1,Ny
	      G(i,j,tt) = temp1c(i,j)
	    end do 
	  end do
	end do

c	compute omega_hat at present time step (s)
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = omega(i,j,2)
	  end do
	end do 
	call trans_y(Ny, Nz, plany_f, temp1, temp1c)
	do i = 0,Nz
	  do j = 1,Ny
	    omega_hat(i,j,2) = temp1c(i,j)
	  end do
	end do 
	
c	compute first term in omega_rhs
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = omega_Z(i,j,2)
	  end do
	end do
	call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp1, omega_ZZ)
	call trans_y(Ny, Nz, plany_f, omega_ZZ, temp1c)

c	compute sum of four terms in omega_rhs
	do i = 0,Nz
	  do j = 1,Ny
	    omega_rhs(i,j) = - temp1c(i,j)      
     ^		+ dcmplx(gam(j), 0.0d0)*omega_hat(i,j,2) 
     ^		- dcmplx(1.0d0/(xi*xi*LJ), 0.0d0)
     ^		*(dcmplx(3.0d0, 0.0d0)*G(i,j,2)- G(i,j,3))
	  end do
	end do

c	solve helmholtz equation at each y location
	do j = 1,Ny
	  do i = 0,Nz
	    rhs_temp(i) = omega_rhs(i,j)
	  end do

	  do i = 1,2
	   bc_temp(i) = BC_hat(j,i)
	  end do

	  call slv_helm2(Nz, planzc_f, planzc_b, lam(j), rhs_temp,
     ^         bc_temp, BC_type, slv_temp)

	  do i = 0,Nz
	  omega_hat(i,j,1) = slv_temp(i)
	  end do
	end do

c	transform to get new omega value
	do i = 0,Nz
	  do j = 1,Ny
	    temp1c(i,j) = omega_hat(i,j,1)
	  end do
	end do
	call itrans_y(Ny, Nz, plany_b, temp1c, temp1)
	do i = 0,Nz
	  do j = 1,Ny
	    omega(i,j,1) = temp1(i,j)
	  end do
	end do

	end


c ***************************** T SOLVER ******************************
	subroutine slv_T(Ny, Nz, Mz,  plany_f, plany_b, planz_f, 
     ^     planz_b, planzc_f, planzc_b, PJ, eps, Lat, S, Vs, k, xi,  
     ^     lam, gam, BC_hat, BC_type, psi, T, T_b_ZZ)
c	T solver 
c	input: 
c	  Ny, Nz, PJ, Lat, S, Vs k, xi, lam, gam, BC_hat, BC_type, psi, T
c	  
c	output:
c	  T

c	declare inputs/output variables
	integer*8 Ny, Nz, Mz, BC_type(1:2), plany_f, 
     ^     plany_b, planz_f, planz_b, planzc_f, planzc_b
	real*8 PJ, Lat, S, k, xi, eps, Vs(0:Nz,1:Ny),  
     ^		psi(0:Nz,1:Ny,1:3),  
     ^		T(0:Nz,1:Ny,1:3), T_b_ZZ(0:Nz,1:Ny),   
     ^		lam(1:Ny), gam(1:Ny)
	complex*16 BC_hat(1:Ny,1:2) 
 
c	declare subroutine specific variables
	integer*8 tt
	real*8 T_y(0:Nz,1:Ny,2:3), T_Z(0:Nz,1:Ny,2:3), T_ZZ(0:Nz,1:Ny), 
     ^		psi_y(0:Nz,1:Ny,2:3), psi_Z(0:Nz,1:Ny,2:3), 
     ^		temp1(0:Nz,1:Ny), temp2(0:Nz,1:Ny), temp3(0:Nz,1:Ny)
	complex*16 H(0:Nz,1:Ny,2:3), T_hat(0:Nz,1:Ny,1:2), 
     ^		T_rhs(0:Nz,1:Ny), slv_temp(0:Nz), rhs_temp(0:Nz),  
     ^		bc_temp(1:2), temp1c(0:Nz,1:Ny)
c	compute derivatives for present and old time step
	do tt = 2,3
	  do i = 0,Nz
	    do j = 1,Ny
	      temp1(i,j) = T(i,j,tt)
	      temp2(i,j) = psi(i,j,tt)
	    end do
	  end do
c	  compute T_y
	  call diffy(Ny, Nz, k, plany_f, plany_b, temp1, temp3)
	  do i = 0,Nz
	    do j = 1,Ny
	      T_y(i,j,tt) = temp3(i,j)
	    end do
	  end do

c	  compute T_z
	  call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp1, temp3)
	  do i = 0,Nz
	    do j = 1,Ny
	      T_Z(i,j,tt) = temp3(i,j)
	    end do
	  end do

c	  compute psi_y
	  call diffy(Ny, Nz, k, plany_f, plany_b, temp2, temp3)
	  do i = 0,Nz
	    do j = 1,Ny
	      psi_y(i,j,tt) = temp3(i,j)
	    end do
	  end do
	  
c	  compute psi_Z
	  call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp2, temp3)
	  do i = 0,Nz
	    do j = 1,Ny
	      psi_Z(i,j,tt) = temp3(i,j)
	    end do
	  end do
c	  compute non-linear terms (H)
	  do i = 0,Nz
	    do j = 1,Ny
 	      temp1(i,j) =
     ^		xi*(psi_y(i,j,tt)*T_Z(i,j,tt)  
     ^		- psi_Z(i,j,tt)*T_y(i,j,tt)) 
c     ^		- (1.0d0/(Lat*Lat))*Vs(i,j)*T_y(i,j,tt) 
     ^		+ S * psi_y(i,j,tt) 
c     ^		- ((xi*xi)/PJ)*T_b_ZZ(i,j) 
     ^		+ xi * psi_Z(i,j,tt) 
            end do
	  end do
	  call trans_y(Ny, Nz, plany_f, temp1, temp1c)
	  
	  do i = 0,Nz
	    do j = 1,Ny
	      H(i,j,tt) = temp1c(i,j)
	    end do 
	  end do	  
	end do

c	compute T_hat at present time step (s)
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = T(i,j,2)
	  end do
	end do 
	call trans_y(Ny, Nz, plany_f, temp1, temp1c)
	do i = 0,Nz
	  do j = 1,Ny
	    T_hat(i,j,2) = temp1c(i,j)
	  end do
	end do 

c	compute first term in T_rhs
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = T_Z(i,j,2)
	  end do
	end do
	call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp1, T_ZZ)
	call trans_y(Ny, Nz, plany_f, T_ZZ, temp1c)


c	compute sum of four terms in T_rhs
	do i = 0,Nz
	  do j = 1,Ny
	    T_rhs(i,j) = - temp1c(i,j) 
     ^		+ dcmplx(gam(j), 0.0d0)*T_hat(i,j,2) 
     ^		- dcmplx((PJ)/(xi*xi), 0.0d0)
     ^		*(dcmplx(3.0d0, 0.0d0)*H(i,j,2) - H(i,j,3))
	  end do
	end do

c	solve helmholtz equation at each y location
	do j = 1,Ny
	  do i = 0,Nz
	    rhs_temp(i) = T_rhs(i,j)
	  end do
	  do i = 1,2
	   bc_temp(i) = BC_hat(j,i)
	  end do

	  call slv_helm2(Nz, planzc_f, planzc_b, lam(j), rhs_temp,
     ^         bc_temp, BC_type, slv_temp)

	  do i = 0,Nz
	    T_hat(i,j,1) = slv_temp(i)
	  end do
	end do
       
c	transform to get new T value
	do i = 0,Nz
	  do j = 1,Ny
	    temp1c(i,j) =T_hat(i,j,1)
	  end do
	end do

	call itrans_y(Ny, Nz, plany_b, temp1c, temp1)

	do i = 0,Nz
	  do j = 1,Ny
	    T(i,j,1) = temp1(i,j)
	  end do
	end do

	end


c **************************** psi SOLVER *****************************
	subroutine slv_psi(Ny, Nz, plany_f, plany_b, planzc_f, planzc_b, 
     ^     xi, eta, BC_hat, BC_type, omega, psi)
c	psi solver 
c	input: 
c	  Ny, Nz, xi, eta, BC_hat, BC_type, omega, psi
c	  
c	output:
c	  psi

c	declare inputs/output variables
	integer*8 Ny, Nz, BC_type(1:2), plany_f, plany_b,
     ^       planzc_f, planzc_b
	real*8 xi, omega(0:Nz,1:Ny,1:3), psi(0:Nz,1:Ny,1:3), eta(1:Ny)
	complex*16 BC_hat(1:Ny,1:2)
 
c	declare subroutine specific variables
	integer*8 tt
	real*8 temp1(0:Nz,1:Ny)
	complex*16 psi_rhs(0:Nz,1:Ny), slv_temp(0:Nz), rhs_temp(0:Nz), 
     ^		bc_temp(1:2), psi_hat(0:Nz,1:Ny), temp1c(0:Nz,1:Ny)

c 	compute transform of omega at new time step (s+1)
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = omega(i,j,1)
	  end do
	end do
	call trans_y(Ny, Nz, plany_f, temp1, temp1c)
	
c	compute psi_rhs
	do i = 0,Nz
	  do j = 1,Ny
	    psi_rhs(i,j) = -temp1c(i,j)/dcmplx(xi*xi, 0.0d0)
	  end do
	end do

c	solve helmholtz equation at each y location
	do j = 1,Ny
	  do i = 0,Nz
	    rhs_temp(i) = psi_rhs(i,j)
	  end do

	  do i = 1,2
	   bc_temp(i) = BC_hat(j,i)
	  end do

	  call slv_helm2(Nz, planzc_f, planzc_b, eta(j), rhs_temp,
     ^         bc_temp, BC_type, slv_temp)

	  do i = 0,Nz
	    psi_hat(i,j) = slv_temp(i)
	  end do
	end do

c	transform to get new psi value
	do i = 0,Nz
	  do j = 1,Ny
	    temp1c(i,j) = psi_hat(i,j)
	  end do
	end do
	call itrans_y(Ny, Nz, plany_b, temp1c, temp1)
	do i = 0,Nz
	  do j = 1,Ny
	    psi(i,j,1) = temp1(i,j)
	  end do
	end do

	end


c **********************CROSS-SECTIONAL VELOCITY V & W ************************
        subroutine velocity(Ny, Nz, Mz, ver, plany_f, plany_b, 
     ^     planz_f, planz_b, xi, k, psi, v, w)
c       calculate v & w by the streamfunction
        integer*8 Ny, Nz, Mz, ver, plany_f, plany_b, planz_f, 
     ^       planz_b
        real*8 xi, k, psi(0:Nz,1:Ny,1:3), v(0:Nz,1:Ny),
     ^         w(0:Nz,1:Ny), temp_psi(0:Nz,1:Ny),
     ^         temp_v(0:Nz,1:Ny), temp_w(0:Nz,1:Ny)

        do i = 0,Nz
           do j = 1,Ny
              temp_psi(i,j) = psi(i,j,ver)
           end do
        end do

        call diffZ(Ny, Nz, Mz, planz_f, planz_b, temp_psi, temp_v)    
        call diffy(Ny, Nz, k, plany_f, plany_b, temp_psi, temp_w)

        do i = 0,Nz
          do j = 1,Ny
             v(i,j) = xi*temp_v(i,j)
             w(i,j) = -temp_w(i,j)
          end do
        end do

        end

c *************************ENERGY CALCULATION *************************
       subroutine comp_energy(Ny, Nz, Mz, ver, plany_f, plany_b, 
     ^   planz_f, planz_b, zz, k, xi, Pe, LJ, QJ, L, H, dUs, u, v, w, T,
     ^  psi, omega, uke, pke, pke2, poe, ps, pu, bb, dh, dp, dpe2, dpe3)
c       calculate the energy at each x-location

c       define inputs/outputs variables
        integer*8 Ny, Nz, Mz, ver, planz_f, planz_b, plany_f, plany_b
        real*8 xi, L, H, k, Pe, LJ, QJ, zz(0:Nz), u(0:Nz,1:Ny,1:3),
     ^       T(0:Nz,1:Ny,1:3), v(0:Nz,1:Ny), w(0:Nz,1:Ny), 
     ^       psi(0:Nz,1:Ny,1:3), omega(0:Nz,1:Ny,1:3),
     ^       dUs(0:Nz,1:Ny), uke, pke, pke2, poe, dh, dp, bb, ps, pu,
     ^       dpe2, dpe3
        
c       define the subroutine variables
        real*8 uprime(0:Nz,1:Ny), uke_temp(0:Nz,1:Ny), 
     ^       pke_temp(0:Nz,1:Ny), poe_temp(0:Nz,1:Ny), 
     ^       pu_temp(0:Nz,1:Ny), ps_temp(0:Nz,1:Ny), 
     ^       dpe2_temp(0:Nz,1:Ny), dpe3_temp(0:Nz,1:Ny),
     ^       bb_temp(0:Nz,1:Ny), dh_temp(0:Nz,1:Ny), dp_temp(0:Nz,1:Ny),
     ^       u_z(0:Nz,1:Ny), u_y(0:Nz,1:Ny), v_z(0:Nz,1:Ny),
     ^       v_y(0:Nz,1:Ny), w_z(0:Nz,1:Ny), w_y(0:Nz,1:Ny),
     ^       vavg_temp(1:Ny), Tavg_temp(1:Ny), vavg(0:Nz), 
     ^       Tavg(0:Nz,1:Ny), Tavg_z(0:Nz,1:Ny),
     ^       pke2_temp(0:Nz,1:Ny)


c       subtract the base state of velocity to get the perturbation
        do i = 0,Nz
           do j = 1,Ny
              uprime(i,j) = u(i,j,ver)
c - (zz(i)+1.0d0)
           end do
        end do

c       store the y-average of v and T
        do i = 0,Nz
           do j = 1,Ny
              vavg_temp(j) = v(i,j)
              Tavg_temp(j) = T(i,j,ver)
           end do
           vavg(i) = sum(vavg_temp)/dble(Ny)
           do j = 1,Ny
              Tavg(i,j) = sum(Tavg_temp)/dble(Ny)
           end do
        end do


c       compute the derivative of u,v and w
        call diffZ(Ny, Nz, Mz, planz_f, planz_b, uprime, u_z)    
        call diffy(Ny, Nz, k, plany_f, plany_b, uprime, u_y)
        call diffZ(Ny, Nz, Mz, planz_f, planz_b, v, v_z)    
        call diffy(Ny, Nz, k, plany_f, plany_b, v, v_y)
        call diffZ(Ny, Nz, Mz, planz_f, planz_b, w, w_z)    
        call diffy(Ny, Nz, k, plany_f, plany_b, w, w_y)

c       compute the derivative of Tavg
        call diffZ(Ny, Nz, Mz, planz_f, planz_b, Tavg, Tavg_z)  

c       store the integrand matrix
        do i = 0,Nz
          do j = 1,Ny
               uke_temp(i,j) = 0.5d0*uprime(i,j)*uprime(i,j)
               pke_temp(i,j) = 0.5d0*(v(i,j)*v(i,j)+w(i,j)*w(i,j))
               poe_temp(i,j) = -QJ*T(i,j,ver)*zz(i) 
               pu_temp(i,j) = -uprime(i,j)*w(i,j)
               ps_temp(i,j) = -uprime(i,j)*w(i,j)*xi*dUs(i,j)
               bb_temp(i,j) = QJ*w(i,j)*T(i,j,ver)
               dh_temp(i,j) = LJ*(u_y(i,j)*u_y(i,j)
     ^                        +xi*xi*u_z(i,j)*u_z(i,j))
               dp_temp(i,j) = LJ*(v_y(i,j)*v_y(i,j)
     ^                        +xi*xi*v_z(i,j)*v_z(i,j)
     ^                        +w_y(i,j)*w_y(i,j)
     ^                        +xi*xi*w_z(i,j)*w_z(i,j))
               dpe2_temp(i,j) = QJ*zz(i)*vavg(i)
               dpe3_temp(i,j) = QJ*xi*Tavg_z(i,j)/Pe
               pke2_temp(i,j) = -0.5d0*psi(i,j,ver)*omega(i,j,ver)
          end do
        end do

c       Do the integration in y and Z
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, uke_temp, uke)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, pke_temp, pke)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, poe_temp, poe)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, pu_temp, pu)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, ps_temp, ps)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, bb_temp, bb)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, dh_temp, dh)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, dp_temp, dp)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, dpe2_temp, dpe2)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, dpe3_temp, dpe3)
        call inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, pke2_temp, pke2)


c       save the y-z-averaged energy which depends on the frequency in x.
          uke = uke/(xi*L*H)
          pke = pke/(xi*L*H)
          poe = poe/(xi*L*H)
          pu = pu/(xi*L*H)
          ps = ps/(xi*L*H)
          bb = bb/(xi*L*H)
          dh = dh/(xi*L*H)
          dp = dp/(xi*L*H)
          dpe2 = dpe2/(xi*L*H)
          dpe3 = dpe3/(xi*L*H)
          pke2 = pke2/(xi*L*H)
          
        end

c ************************* WRITE SOLN TO FILE ************************
	subroutine writesoln(Ny, Nz, u, psi, omega, T, ver, run)
c	write solution to file
	
	integer*8 Ny, Nz, ver
	real*8 temp1(0:Nz,1:Ny), u(0:Nz,1:Ny,1:3), psi(0:Nz,1:Ny,1:3),  
     ^		omega(0:Nz,1:Ny,1:3), T(0:Nz,1:Ny,1:3)
	character run*30, flname*40

c	write u to file	  
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = u(i,j,ver)
	  end do
	end do
	if (ver .eq. 1) then
	  flname = "u_" // run
	  call writemtrx(Ny, Nz, temp1, flname)
	else if (ver .eq. 3) then
	  flname = "u_old_" // run
	  call writemtrx(Ny, Nz, temp1, flname)
	end if

c	write psi to file	  
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = psi(i,j,ver)
	  end do
	end do
	if (ver .eq. 1) then
	  flname = "psi_" // run
	  call writemtrx(Ny, Nz, temp1, flname)
	else if (ver .eq. 3) then
	  flname = "psi_old_" // run
	  call writemtrx(Ny, Nz, temp1, flname)
	end if

c	write omega to file	  
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = omega(i,j,ver)
	  end do
	end do
	if (ver .eq. 1) then
	  flname = "omega_" // run
	  call writemtrx(Ny, Nz, temp1, flname)
	else if (ver .eq. 3) then
	  flname = "omega_old_" // run
	  call writemtrx(Ny, Nz, temp1, flname)
	end if

c	write T to file	  
	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = T(i,j,ver)
	  end do
	end do
	if (ver .eq. 1) then
	  flname = "T_" // run
	  call writemtrx(Ny, Nz, temp1, flname)
	else if (ver .eq. 3) then
	  flname = "T_old_" // run
	  call writemtrx(Ny, Nz, temp1, flname)
	end if

	end


c ********************** WRITE CONSTANT TO FILE ***********************
	subroutine writecnt(in, filename)
c	write in to file

	real*8 in
	character filename*25, location*50

	location = "../data/" // filename

c	write in to file
	open(unit=2,file=location,status='unknown')
	write(2,'(f25.5)') in
	close(2)

	end


c ************************ WRITE VECTOR TO FILE ***********************
	subroutine writevec(start, stp, in, filename)
c	write in to file

	integer start, stp
	real*8 in(start:stp)
	character filename*8, location*70

	location = "../data/" // filename

c	write in to file
	open(unit=2,file=location,status='unknown')
	write(2,'(1X,f25.15)') (in(j),j=start,stp)
	close(2)

	end


c ************************ WRITE MATRIX TO FILE ***********************
	subroutine writemtrx(Ny, Nz, in, filename)
c	write in to file

	integer*8 Ny, Nz
	real*8 in(0:Nz,1:Ny)
	character filename*40, location*70

	location = "../data/" // filename

c	write in to file
	open(unit=2,file=location,status='unknown')
	do i = 0,Nz
	  write(2,'(10000f35.25)') (in(i,j),j=1,Ny)
	end do
	close(2)

	end


c ************************ READ FILE TO MATRIX ************************
	subroutine readfle(Ny, Nz, in, filename)
c	write in to file

	integer*8 Ny, Nz
	integer e, c
	real*8 in(0:Nz,1:Ny)
	character filename*(*)

c	write in to file
	open(unit=2,file=filename, err=8, iostat=e, status='old')
	do i = 0,Nz
	  read(2,'(10000f35.25)', err=9, iostat=e) (in(i,j),j=1,Ny)
	end do
	close(2)
	return

8	write(*,*) 'open I/O error # '
	write(*,*) e
	write(*,*) filename
	stop

9	write(*,*) 'read I/O error # '
	write(*,*) e
	write(*,*) filename
	stop

	end


c ******************************* NORM ********************************
	function nrm(Ny, Nz, in)
c	determine norm of input --> max(abs(new-present))

	integer*8 Ny, Nz
	real*8 nrm, temp1(0:Nz,1:Ny), temp2(0:Nz), mx, in(0:Nz,1:Ny,1:3)

	do i = 0,Nz
	  do j = 1,Ny
	    temp1(i,j) = dabs(in(i,j,2) - in(i,j,3))
	  end do
	end do

	do i = 0,Nz
	  mx = temp1(i,1)
	  do j = 2,Ny
	    mx = dmax1(mx, temp1(i,j))
	  end do
	  temp2(i) = mx
	end do

	nrm = temp2(0)
	do i = 1,Nz
	  nrm = dmax1(nrm, temp2(i))
	end do
	
	end


c ******************************* MAXdD *******************************
	function maxdd(Ny, Nz, in)
c	determine maximum value of a matrix
	integer*8 Ny, Nz
	real*8 maxdd, in(0:Nz,1:Ny)

	real*8 mx, temp1(0:Nz)

	do i = 0,Nz
	  mx = in(i,2)
	  do j = 3,Ny
	    mx = dmax1(mx, in(i,j))
	  end do
	  temp1(i) = mx
	  
	end do

	maxdd = temp1(0)
	
	do i = 1,Nz
	  maxdd = dmax1(maxdd, temp1(i))
	end do

	end
c********************************NUM2STR*******************************
	subroutine num2str(num, str)
c	convert a 3 digit integer into a 3 digit string
	character*4 str*3
	integer*8 num, unit1, unit2, dig1, dig2, dig3
c	construct unit1 and unit2 representing hundred and teen digit
	unit1 = 100
	unit2 = 10
	dig1 = num/unit1
	dig2 = (num - dig1*unit1)/unit2
	dig3 = num - dig1 * unit1 - dig2 * unit2
c	num+48 is the asicii codes for numbers: 48 is '0' and 57 is '9'
	str = char(dig1+48)//char(dig2+48)//char(dig3+48)
	end


c ****************************** PARAMS *******************************
	subroutine wait()
c	wait for input to continue
	character in*1

	write(*,*) 'Continue (y or n):'
	read(*,*) in

	if (in .eq. 'y') then
	  write(*,*) 'Here we go ...'
	elseif (in .eq. 'n') then
	  stop
	endif

	end 


c ****************************** PARAMS *******************************
	subroutine params(Ny, Nz, H, L, xi, LJ, PJ, QJ, S, ER, SD, Res,
     ^		eps, Lat, dt, tstep, Us, Vs, A, B, C, D,
     ^		datagap, energygap, run_nm, BC_u_hat, BC_psi_hat, BC_omega_hat, 
     ^		BC_T_hat, BC_u_type, BC_psi_type, BC_o_type, BC_T_type)
c	write all parameters used to terminal

	integer*8 Ny, Nz, tstep, datagap, energygap, BC_u_type(1:2),
     ^       BC_o_type(1:2), BC_psi_type(1:2), BC_T_type(1:2)
	real*8 H, L, xi, LJ, PJ, QJ, S, ER, SD, Res, dt, eps, Lat, 
     ^		 Us(0:Nz,1:Ny), Vs(0:Nz, 1:Ny), A, B, C, D
	complex*16 BC_u_hat(1:Ny,1:2), BC_psi_hat(1:Ny,1:2), 
     ^		BC_omega_hat(1:Ny,1:2), BC_T_hat(1:Ny,1:2)
	character*4 run_nm*(*), location*40

	location = "../data/" // run_nm // ".txt"

	write(*,*) ''
	write(*,'(15x,a,a)') 'RUN NAME: ', run_nm 
	write(*,*) ''
        write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,1x,a,i3,2x,a,1x,a,i3,2x,a,1x,a,f4.3,2x,a,1x,a,
     ^f7.1,1x,a)') '|', 'Ny = ', Ny, '|', 'Nz = ', Nz, '|', 'LJ = ', 
     ^LJ,'|', 'PJ = ', PJ, '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,1x,a,f4.2,1x,a,2x,a,f4.1,2x,a,1x,a,
     ^f5.3,2x,a,2x,a,f6.3,1x,a)') '|', 'QJ = ', QJ, '|', 'S = ', S, 
     ^ '|', 'H = ', H, '|', 'L = ', L, '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,1x,a,f2.0,1x,a,1x,a,f4.1,2x,a,1x,a,
     ^f4.1,1x,a,2x,a,f4.2,2x,a)') '|', 'ER = ', ER, '|', 'Res = ', Res, 
     ^ '|', 'eps = ', eps, '|', 'Lat = ', Lat, '|'
        write(*,*) '----------------------------------------------------
     ^-'
        write(*,'(1x,a,1x,a,f4.1,1x,a,1x,a,f4.1,1x,a,1x,a,
     ^f4.1,3x,a,1x,a,i2,1x,a)') '|', ' A = ', A, '|', ' B = ', B, 
     ^ '|', ' C = ', C, '|', ' e.gap = ', energygap, '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,1x,a,1x,i7,1x,a,1x,a,1x,f5.4,1x,a,2x,a,
     ^1x,i5,1x,a)') '|','t.s = ', tstep, '|','dt = ', dt,'|',
     ^'t.s. gap = ', datagap, '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,24x,a,23x,a)') '|', 'U_s:', '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,9x,a,9x,a,8x,a,9x,a)') '|', '@ Z = 1', '|', 
     ^'@ Z = -1', '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,9x,f8.6,8x,a,8x,f8.6,9x,a)') '|', Us(0,1), '|', 
     ^Us(Nz,1), '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,24x,a,23x,a)') '|', 'V_s:', '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,9x,a,9x,a,8x,a,9x,a)') '|', '@ Z = 1', '|', 
     ^'@ Z = -1', '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,9x,f8.6,8x,a,8x,f8.6,9x,a)') '|', Vs(0,1), '|', 
     ^Vs(Nz,1), '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,23x,a,23x,a)') '|', 'BC''s:', '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,9x,a,9x,a,8x,a,9x,a)') '|', '@ Z = 1', '|', 
     ^'@ Z = -1', '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,4x,a,i1,a,f3.1,6x,a,4x,a,i1,a,f3.1,5x,a)') '|', 
     ^ 'u(y,1)(',BC_u_type(1),') = ', dreal(BC_u_hat(1,1)),
     ^ '|', 'u(y,-1)(',BC_u_type(2),') = ',dreal(BC_u_hat(1,2)),'|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,4x,a,i1,a,f3.1,4x,a,4x,a,i1,a,f3.1,3x,a)') '|', 
     ^'psi(y,1)(',BC_psi_type(1),') = ', dreal(BC_psi_hat(1,1)), '|',  
     ^'psi(y,-1)(',BC_psi_type(2),') = ', dreal(BC_psi_hat(1,2)),'|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,3x,a,i1,a,f3.1,3x,a,3x,a,i1,a,f3.1,2x,a)') '|', 
     ^'omega(y,1)(',BC_o_type(1),') = ', dreal(BC_omega_hat(1,1)), '|', 
     ^ 'omega(y,-1)(',BC_o_type(2),') = ', dreal(BC_omega_hat(1,2)), '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,'(1x,a,4x,a,i1,a,f3.1,6x,a,4x,a,i1,a,f3.1,5x,a)') '|', 
     ^'T(y,1)(',BC_T_type(1),') = ', dreal(BC_T_hat(1,1)), '|',  
     ^'T(y,-1)(',BC_T_type(2),') = ', dreal(BC_T_hat(1,2)), '|'
	write(*,*) '----------------------------------------------------
     ^-'
	write(*,*) ''

c	write in to file
	open(unit=2,file=location,status='unknown')

	write(2,*) ''
	write(2,'(15x,a,a)') 'RUN NAME: ', run_nm 
	write(2,*) ''
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,1x,a,i3,2x,a,1x,a,i3,2x,a,1x,a,f4.3,2x,a,1x,a,
     ^f7.1,1x,a)') '|', 'Ny = ', Ny, '|', 'Nz = ', Nz, '|', 'LJ = ', 
     ^LJ,'|', 'PJ = ', PJ, '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,1x,a,f4.2,1x,a,2x,a,f4.1,2x,a,1x,a,
     ^f5.3,2x,a,2x,a,f6.3,1x,a)') '|', 'QJ = ', QJ, '|', 'S = ', S, 
     ^ '|', 'H = ', H, '|', 'L = ', L, '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,1x,a,f2.0,1x,a,1x,a,f4.1,2x,a,1x,a,
     ^f4.1,1x,a,2x,a,f4.2,2x,a)') '|', 'ER = ', ER, '|', 'Res = ', Res, 
     ^ '|', 'eps = ', eps, '|', 'Lat = ', Lat, '|'
        write(2,*) '----------------------------------------------------
     ^-'
        write(2,'(1x,a,1x,a,f4.1,1x,a,1x,a,f4.1,1x,a,1x,a,
     ^f4.1,3x,a,1x,a,i2,1x,a)') '|', ' A = ', A, '|', ' B = ', B, 
     ^ '|', ' C = ', C, '|', ' e.gap = ', energygap, '|'
	write(2,*) '----------------------------------------------------
     ^-'
        write(2,'(1x,a,1x,a,1x,i7,1x,a,1x,a,1x,f5.4,1x,a,2x,a,
     ^1x,i5,1x,a)') '|','t.s = ', tstep, '|','dt = ', dt,'|',
     ^'t.s. gap = ', datagap, '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,24x,a,23x,a)') '|', 'U_s:', '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,9x,a,9x,a,8x,a,9x,a)') '|', '@ Z = 1', '|', 
     ^'@ Z = -1', '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,9x,f8.6,8x,a,8x,f8.6,9x,a)') '|', Us(0,1), '|', 
     ^Us(Nz,1), '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,24x,a,23x,a)') '|', 'V_s:', '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,9x,a,9x,a,8x,a,9x,a)') '|', '@ Z = 1', '|', 
     ^'@ Z = -1', '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,9x,f8.6,8x,a,8x,f8.6,9x,a)') '|', Vs(0,1), '|', 
     ^Vs(Nz,1), '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,23x,a,23x,a)') '|', 'BC''s:', '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,9x,a,9x,a,8x,a,9x,a)') '|', '@ Z = 1', '|', 
     ^'@ Z = -1', '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,4x,a,i1,a,f3.1,6x,a,4x,a,i1,a,f3.1,5x,a)') '|', 
     ^ 'u(y,1)(',BC_u_type(1),') = ', dreal(BC_u_hat(1,1)),
     ^ '|', 'u(y,-1)(',BC_u_type(2),') = ',dreal(BC_u_hat(1,2)),'|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,4x,a,i1,a,f3.1,4x,a,4x,a,i1,a,f3.1,3x,a)') '|', 
     ^'psi(y,1)(',BC_psi_type(1),') = ', dreal(BC_psi_hat(1,1)), '|',  
     ^'psi(y,-1)(',BC_psi_type(2),') = ', dreal(BC_psi_hat(1,2)),'|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,3x,a,i1,a,f3.1,3x,a,3x,a,i1,a,f3.1,2x,a)') '|', 
     ^'omega(y,1)(',BC_o_type(1),') = ', dreal(BC_omega_hat(1,1)), '|', 
     ^ 'omega(y,-1)(',BC_o_type(2),') = ', dreal(BC_omega_hat(1,2)), '|'
	write(2,*) '----------------------------------------------------
     ^-'
	write(2,'(1x,a,4x,a,i1,a,f3.1,6x,a,4x,a,i1,a,f3.1,5x,a)') '|', 
     ^'T(y,1)(',BC_T_type(1),') = ', dreal(BC_T_hat(1,1)), '|',  
     ^'T(y,-1)(',BC_T_type(2),') = ', dreal(BC_T_hat(1,2)), '|'
	write(2,*) '----------------------------------------------------
     ^-'	
        write(2,*) ''

	close(2)

	end




 
