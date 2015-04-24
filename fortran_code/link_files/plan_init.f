c ------------------------------------------------------------------- c
c                 Initialization of all the plans                     c
c ------------------------------------------------------------------- c
c                                                                     c
c  inputs:                                                            c
c    	Ny, Nz, Mz, Nyf, Nzf, Mzf                                     c
c  outputs:                                                           c
c       plany_f, plany_b, planz_f, planz_b,                           c
c       planyf_f, planyf_b, planzf_f, planzf_b,                       c
c       planzc_f, planzc_b                                            c
c                                                                     c
c                                                                     c
c  internal subroutines used:		   		              c
c    	fft_f_init(N, in, out, plan)           	                      c
c    	fft_b_init(N, in, out, plan)     	                      c
c	chebyshev_f_init(N, M, in, out, plan)                         c
c	chebyshev_b_init(N, M, in, out, plan)                         c
c	chebyshev_c_f_init(N, M, in, out, plan)                       c
c	chebyshev_c_b_init(N, M, in, out, plan)                       c
c								      c
c ------------------------------------------------------------------- c

      subroutine plan_init(Ny, Nz, Mz, plany_f, plany_b,
     ^  planz_f, planz_b, planzc_f, planzc_b)

c     define variables
      integer*8 Ny, Nz, Mz, plany_f, plany_b, planz_f, planz_b, 
     ^     planzc_f, planzc_b

      real*8 pi, u(0:Nz,1:Ny), uy(1:Ny), uz(0:Nz), y(1:Ny), z(0:Nz)

      complex*16 uyout(1:Ny), uzout(0:Nz)

      parameter (pi = 2.0d0*dacos(0.0d0))

c     define normal grids
      do j=1,Ny
         y(j)=dble(j-1)/dble(Ny)
      end do
      do i = 0,Nz
         z(i) = dcos(pi*dble(i)/dble(Nz))
      end do

c     define function on normal grids
      do i=0,Nz
         do j=1,Ny
            u(i,j)=dcos(2.0d0*pi*y(j))*(z(i)+1)
         end do
      end do

c     initialize plans on normal grids
      do j=1,Ny
         uy(j)=u(10,j)
      end do

      do i=0,Nz
         uz(i)=u(i,10)
      end do

      call fft_f_init(Ny, uy, uyout, plany_f)
      call chebyshev_f_init(Nz, Mz, uz, uzout, planz_f)
      call chebyshev_c_f_init(Nz, Mz, uzout, uzout, planzc_f)

      call fft_b_init(Ny, uyout, uy, plany_b)
      call chebyshev_b_init(Nz, Mz, uzout, uz, planz_b)
      call chebyshev_c_b_init(Nz, Mz, uzout, uzout, planzc_b)

      end





c *********************************************************************
c **************************** SUB ROUTINES ***************************
c *********************************************************************

c *************************** fft_f_init*******************************
	subroutine fft_f_init(N, in, out, plan)
c	compute fft of in

c	define variables
c	rows --> z
c	cols --> y
	integer*8 N
	real*8 in(1:N)
	complex*16 out(1:N)

	integer*8 plan, estimate
	parameter (estimate = 64)
	complex*16 temp_out(1:N/2+1)

c	compute fft and remove plan
	call dfftw_plan_dft_r2c_1d(plan,N,in,temp_out,estimate)
	call dfftw_execute_dft_r2c(plan, in, temp_out)
c	call dfftw_destroy_plan(plan)

c	store temp_out to out
	do i = 1,N/2+1
	  out(i) = temp_out(i)
	end do

c	add redundent complex conjugate data
	do i = N/2+2,N
	  out(i) = dconjg(temp_out(N-i+2))  
	end do

	end


c *************************** fft_b_init*******************************
	subroutine fft_b_init(N, in, out, plan)
c	compute inverse fft of in

c	define variables
c	rows --> z
c	cols --> y
	integer*8 N
	real*8 out(1:N)
	complex*16 in(1:N)

	integer*8 plan, estimate
	parameter (estimate = 64)
	complex*16 temp_in(1:N/2+1)

	do i = 1,N/2+1
	  temp_in(i) = in(i)
	end do

	call dfftw_plan_dft_c2r_1d(plan,N,temp_in,out,estimate)
	call dfftw_execute_dft_c2r(plan, temp_in, out)
c	call dfftw_destroy_plan(plan)	

	end


c ************************** CHEB_r2c_FORWARD *****************************
	subroutine chebyshev_f_init(N, M, in, out, plan)
c	compute chebyshev transform

c	define variables
	integer*8 N, M, plan
	real*8 in(0:N), in_t(1:M)
        complex*16 in_t_hat(1:M), out(0:N)

c	extend in to in_t_2N-j = in_j, j = 1,2,...,N-1
	do i = 1,N
	  in_t(i) = in(i-1)
	  in_t(M-i+1) = in(i)
	end do

c	compute fft
        call fft_f_init(M, in_t, in_t_hat, plan)

c	keep N+1 terms 
	do i = 0,N
	  out(i) = in_t_hat(i+1)
	end do

c 	correct for normalization factor
	out(0) = out(0)*dcmplx(2.0d0,0.0d0)
	out(N) = out(N)*dcmplx(2.0d0,0.0d0)
	do i = 0,N
	  out(i) = out(i)/dcmplx(dble(M),0.0d0)
	end do

	end


c ************************** CHEB_c2r_BACKWARD *****************************
	subroutine chebyshev_b_init(N, M, in, out, plan)
c	compute inverse chebyshev transform

c	define variables
	integer*8 N, M, plan
	complex*16 in(0:N), in_t(1:M)
        real*8 out(0:N), out_t(1:M)

c	factor of 2 correction
	in(0) = in(0)/dcmplx(2.0d0,0.0d0)
	in(N) = in(N)/dcmplx(2.0d0,0.0d0)

c	extend data
	do i = 1,N
	  in_t(i) = in(i-1)
	  in_t(M-i+1) = in(i)
	end do

c	compute fft
        call fft_b_init(M, in_t, out_t, plan)

c	only keep N+1 terms
c	do i = 0,N
c	  out(i) = out_t(i+1)
c	end do

	end

c ************************** CHEB_c2c_FORWARD *****************************
	subroutine chebyshev_c_f_init(N, M, in, out, plan)
c	compute chebyshev transform

c	define variables
	integer*8 N, M, plan, forward, estimate
	complex*16 in(0:N), in_t(M), in_t_hat(M), out(0:N)

	parameter (forward = -1, estimate=64)

c	extend in to in_t_2N-j = in_j, j = 1,2,...,N-1
	do i = 1,N
	  in_t(i) = in(i-1)
	  in_t(M-i+1) = in(i)
	end do

c	compute fft and plan
	call dfftw_plan_dft_1d(plan,M,in_t,in_t_hat,forward,estimate)
	call dfftw_execute_dft(plan, in_t, in_t_hat)
c	call dfftw_destroy_plan(plan)

c	keep N+1 terms 
	do i = 0,N
	  out(i) = in_t_hat(i+1)
	end do

c 	correct for normalization factor
	out(0) = out(0)*dcmplx(2.0d0,0.0d0)
	out(N) = out(N)*dcmplx(2.0d0,0.0d0)
	do i = 0,N
	  out(i) = out(i)/dcmplx(dble(M),0.0d0)
	end do

	end

c ************************** CHEB_c2c_BACKWARD *****************************
	subroutine chebyshev_c_b_init(N, M, in, out, plan)
c	compute chebyshev transform

c	define variables
	integer*8 N, M, plan, backward, estimate
	parameter (backward = +1, estimate = 64)
	complex*16 in(0:N), out(0:N), in_t(M), out_t(M)

c	extend in to in_t_2N-j = in_j, j = 1,2,...,N-1
	do i = 1,N
	  in_t(i) = in(i-1)
	  in_t(M-i+1) = in(i)
	end do

c	compute fft and plan
	call dfftw_plan_dft_1d(plan,M,in_t,out_t,backward,estimate)
	call dfftw_execute_dft(plan, in_t, out_t)
c	call dfftw_destroy_plan(plan)

c	keep N+1 terms 
c	do i = 0,N
c	  out(i) = out_t(i+1)
c	end do

	end
