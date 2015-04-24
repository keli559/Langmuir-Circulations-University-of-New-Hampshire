c ------------------------------------------------------------------- c
c                      Helmholtz Equation Solver                      c
c	       using Blanchard Gaussian Elimination Method            c
c ------------------------------------------------------------------- c
c                                                                     c
c  function to solve the helmholtz equation:                          c
c    u'' - lam*u = g                                                  c
c                                                                     c
c  inputs:                                                            c
c    	lam = coeffiecent                                             c
c    	g = rhs                                                       c
c    	bc(1:2): boundard condition values                            c
c	  bc(1) = BC @ Z = 1                                          c
c       bc(2) = BC @ Z = -1                                           c 
c    	b_t(1:2) = boundary condition type (int)                      c
c	  b_t(1) = BC type @ Z = 1                                    c
c	  b_t(2) = BC type @ Z = -1                                   c
c       options:                                                      c
c         0      -->     u(-1) = u(+1) = bc(1) or bc(2)               c
c         1      -->     u'(-1) = u'(+1) = bc(1) or bc(2)             c
c                                                                     c
c  internal subroutines used:		   			      c
c    	cheb_f(N, M, in, out)            			      c
c   	makeA(N, in, out)					      c
c    	makeb(N, in, out)                                             c
c	solveAx(N, A, b, x, info)				      c
c	cheb_b(N, M, in, out)                                         c
c								      c
c ------------------------------------------------------------------- c

	subroutine slv_helm2(Nz, planzc_f, planzc_b, lam, g, bc, b_t, u)

c	DECLARE CONSTANTS
c	matrix bounds
	integer*8 Nz, Mz, b_t(1:2), planzc_f, planzc_b

c	lam
	real*8 lam

c	problem variables
	complex*16 g(0:Nz), g_hat(0:Nz), bcl(1:2,0:Nz+1),  
     ^		dgl(1:Nz-1,1:4),G_vec(0:Nz), u_hat(0:Nz), 
     ^		bc(1:2), u(0:Nz)

c	DEFINITIONS
c	define M
	Mz = 2*Nz

c	SOLVE
c	cheb forward transform 
	call cheb_f(Nz, Mz, planzc_f, g, g_hat)

c	make A matrix
	call makeA(Nz, lam, b_t, bcl, dgl)

c	make G vector
	call makeb(Nz, g_hat, bc, G_vec)

c	solve [A]*u_hat = G_vec
	call solveAx(Nz, bcl, dgl, G_vec, u_hat)

c	cheb backward transform 
	call cheb_b(Nz, Mz, planzc_b, u_hat, u)

	end

c *********************************************************************
c **************************** SUB ROUTINES ***************************
c *********************************************************************

c ************************** CHEB FORWARD *****************************
	subroutine cheb_f(N, M, plan, in, out)
c	compute chebyshev transform

c	define variables
	integer*8 N, M, plan
	complex*16 in(0:N), in_t(M), in_t_hat(M), out(0:N)

c	extend in to in_t_2N-j = in_j, j = 1,2,...,N-1
	do i = 1,N
	  in_t(i) = in(i-1)
	  in_t(M-i+1) = in(i)
	end do

c	compute fft and remove plan
c	call dfftw_plan_dft_1d(plan,M,in_t,in_t_hat,forward,estimate)
	call dfftw_execute_dft(plan, in_t, in_t_hat)
c	call dfftw_destroy_plan(plan)

c	keep N+1 terms 
	do i = 0,N
	  out(i) = in_t_hat(i+1)
	end do

c 	correct for normalization factor
	out(0) = out(0)/dcmplx(2.0d0,0.0d0)
	out(N) = out(N)/dcmplx(2.0d0,0.0d0)
	do i = 0,N
	  out(i) = out(i)/dcmplx(dble(N),0.0d0)
	end do

	end


c ************************** A-MATRIX *********************************
	subroutine makeA(N, in, b_t, out1, out2)
c	create sparse A matrix
	
c	define variables
	integer*8 N, b_t(1:2)
	real*8 in, n_c(2:N), beta(2:N+2), c(2:N), d(2:N), e(2:N)
	complex*16 out1(1:2,0:N+1), out2(1:N-1,1:4) 

c	counter and zero-function
	do i = 2,N
	  n_c(i) = dble(i)
	  beta(i) = 1.0d0
	end do
	beta(N+2) = 0.d0
	beta(N+1) = 0.d0
	beta(N) = 0.0d0
	beta(N-1) = 0.0d0

c	coefficients 
	do i = 2,N
	  c(i) = (-in)/(4.0d0*n_c(i)*(n_c(i)-1.0d0))
	  d(i) = 1.0d0 + (in/(2.0d0*(n_c(i)**2.0d0
     ^         - 1.0d0)))*beta(i)
	  e(i) = ((-in)/(4.0d0*n_c(i)*(n_c(i)
     ^         + 1.0d0)))*beta(i+2)
	end do
	c(2) = c(2)*2.0d0

	
c	create out matrix
	do i = 0,N+1
	  out1(1,i) = dcmplx(0.0d0, 0.0d0)
	  out1(2,i) = dcmplx(0.0d0, 0.0d0)
	end do
	
	do i = 1,N-1
	  out2(i,1) = dcmplx(0.0d0, 0.0d0)
	  out2(i,2) = dcmplx(0.0d0, 0.0d0)
	  out2(i,3) = dcmplx(0.0d0, 0.0d0)
	  out2(i,4) = dcmplx(0.0d0, 0.0d0)
	end do
	
c	create boundary condition rows
	if (b_t(1) .eq. 0) then
	  do j = 0,N
	    out1(1,j) = dcmplx(1.0d0, 0.0d0)
	  end do
	else if (b_t(1) .eq. 1) then
 	  do j = 0,N
	    out1(1,j) = dcmplx(dble(j)**2.0d0, 0.0d0)
	  end do
	else
	  write(*,*) 'ERROR -- boundary condition type1'
	  return
	end if

	if (b_t(2) .eq. 0) then
	  do j = 0,N
	    out1(2,j) = dcmplx((-1.0d0)**dble(j), 0.0d0)
	  end do
	else if (b_t(2) .eq. 1) then
 	  do j = 0,N
	    out1(2,j) = dcmplx((-1.0d0)**dble(j+1)*dble(j)**2.0d0, 0.0d0)
	  end do
	else
	  write(*,*) 'ERROR -- boundary condition type2'
	  return
	end if

	do i = 2,N
	  out2(i-1,1) = dcmplx(c(i), 0.0d0)
	  out2(i-1,2) = dcmplx(d(i), 0.0d0)
	  out2(i-1,3) = dcmplx(e(i), 0.0d0)
	end do

	end


c ************************** b-VECTOR *********************************
	subroutine makeb(N, in, bc, out)
c	create b vector according to dirichlet BC's

c	define variables
	integer*8 N
	real*8 n_c_b(2:N), Cg(2:N), Dg(2:N-2), Eg(2:N-4)
	complex*16 in(0:N), out(0:N), bc(1:2), temp(0:N-2)

c	counter and zero-function
	do i = 2,N
	  n_c_b(i) = dble(i)
	end do

c	coefficients 
	do i = 2,N-4
	  Cg(i) = 1.0d0/(4.0d0*n_c_b(i)*(n_c_b(i)-1.0d0))
	  Dg(i) = (-1.0d0/(2.0d0*((n_c_b(i)**2.0d0) - 1.0d0)))
	  Eg(i) = (1.0d0/(4.0d0*n_c_b(i)*(n_c_b(i)+1.0d0)))
	end do
	do i = N-3,N
	  Cg(i) = 1.0d0/(4.0d0*n_c_b(i)*(n_c_b(i)-1.0d0))
	end do
	Dg(N-3) = (-1.0d0/(2.0d0*((n_c_b(N-3)**2.0d0) - 1.0d0)))
	Dg(N-2) = (-1.0d0/(2.0d0*((n_c_b(N-2)**2.0d0) - 1.0d0)))
	Cg(2) = Cg(2)*2.0d0

c	store in to temp, last 2 elements removed for BC's
	do i = 0,N-2
	  temp(i) = in(i)
	end do

c	create empty out vector
	do i = 0,N
	  out(i) = dcmplx(0.0d0, 0.0d0)
	end do

c 	add boundary conditions
	out(0) = bc(1)
	out(1) = bc(2)

c	compute out vector 
	do i = 2,N-4
	  out(i) = dcmplx(Cg(i), 0.0d0)*temp(i-2) 
     ^           + dcmplx(Dg(i), 0.0d0)*temp(i) 
     ^           + dcmplx(Eg(i), 0.0d0)*temp(i+2)
	end do
	out(N-3) = dcmplx(Cg(N-3), 0.0d0)*temp((N-3)-2) 
     ^           + dcmplx(Dg(N-3), 0.0d0)*temp(N-3)
	out(N-2) = dcmplx(Cg(N-2), 0.0d0)*temp((N-2)-2)
     ^	     + dcmplx(Dg(N-2), 0.0d0)*temp(N-2)
	out(N-1) = dcmplx(Cg(N-1), 0.0d0)*temp((N-1)-2)
	out(N) = dcmplx(Cg(N), 0.0d0)*temp(N-2)

	end


c ************************** SOLVE Ax=b *******************************
	subroutine solveAx(N, bcl, dgl, b, x)
c	solve Ax=b 

c	define variables
	integer*8 N
	complex*16 bcl(1:2,0:N+1), dgl(1:N-1,1:4), b(0:N), x(0:N) 

c	define routine variables
	complex*16 y1, y2, A(1:N+1,1:3)
	
c	combine rhs with bcl and dgl 
	do i = 1,2
	  bcl(i,N+1) = b(i-1)
	end do

	do i = 1,N-1
	  dgl(i,4) = b(i+1)
	end do

c	eliminate upper diagonal
	do i = N-3,3,-1
	  y1 = dgl(i-2,3)/dgl(i,2)
	  dgl(i-2,2) = dgl(i-2,2) - y1*dgl(i,1)
	  dgl(i-2,3) = dgl(i-2,3) - y1*dgl(i,2)
	  dgl(i-2,4) = dgl(i-2,4) - y1*dgl(i,4)
	end do

c	eliminate BC rows
	do i = N+1,3,-1
	  y1 = bcl(1,i-1)/dgl(i-2,2)
	  y2 = bcl(2,i-1)/dgl(i-2,2)
	  bcl(1,i-1) = bcl(1,i-1) - y1*dgl(i-2,2)
	  bcl(1,i-3) = bcl(1,i-3) - y1*dgl(i-2,1)
	  bcl(1,N+1) = bcl(1,N+1) - y1*dgl(i-2,4)
	  
	  bcl(2,i-1) = bcl(2,i-1) - y2*dgl(i-2,2)
	  bcl(2,i-3) = bcl(2,i-3) - y2*dgl(i-2,1)
	  bcl(2,N+1) = bcl(2,N+1) - y2*dgl(i-2,4)	
	end do
	y1 = bcl(1,1)/bcl(2,1)
	bcl(1,0) = bcl(1,0) - y1*bcl(2,0)
	bcl(1,1) = bcl(1,1) - y1*bcl(2,1)
	bcl(1,N+1) = bcl(1,N+1) - y1*bcl(2,N+1)

c	define temp A matrix for forward substitution	
	do i = 1,N+1
	  do j = 1,3
	    A(i,j) = dcmplx(0.0d0, 0.0d0)
	  end do 
	end do

c	store values to A
	do i = 1,2
	  do j = 0,1
	    A(i,j+1) = bcl(i,j)
	  end do
	  A(i,3) = bcl(i,N+1)
	end do
	do i = 3,N+1
	  do j = 1,2
	    A(i,j) = dgl(i-2,j)
	  end do
	  A(i,3) = dgl(i-2,4)
	end do

c	forward substitute 
	if (A(1,1) .eq. 0) then
	  x(0) = 0
	else
	  x(0) = A(1,3)/A(1,1)
	end if
	x(1) = (A(2,3) - A(2,1)*x(0))/A(2,2)
	do i = 3,N+1
	  x(i-1) = (A(i,3) - A(i,1)*x(i-3))/A(i,2)
	end do
	
	end


c ************************** CHEB BACKWARD ****************************
	subroutine cheb_b(N, M, plan, in, out)
c	compute inverse chebyshev transform

c	define variables
	integer*8 N, M, plan
	complex*16 in(0:N), out(0:N), in_t(M), out_t(M)

c	factor of 2 correction
	do i = 1,N-1
	  in(i) = in(i)/dcmplx(2.0d0,0.0d0)
	end do

c	extend data
	do i = 1,N
	  in_t(i) = in(i-1)
	  in_t(M-i+1) = in(i)
	end do

c	compute fft
c	call dfftw_plan_dft_1d(plan,M,in_t,out_t,backward,estimate)
	call dfftw_execute_dft(plan, in_t, out_t)
c	call dfftw_destroy_plan(plan)

c	only keep N+1 terms
	do i = 0,N
	  out(i) = out_t(i+1)
	end do

	end









