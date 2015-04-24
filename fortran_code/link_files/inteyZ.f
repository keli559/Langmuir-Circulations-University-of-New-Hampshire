c ------------------------------------------------------------------- c
c       Integration of Z from -1 to 1 using Chebyshev integral        c
c       and integration of y in a period using trapezoidal rule       c
c ------------------------------------------------------------------- c
c                                                                     c
c  FF(y) =  \int_-1^1 f(y,Z) dZ at each y location                    c
c  FFF = \int FF(y) dy                                                c
c                                                                     c
c  inputs:                                                            c
c    	f - the function on which we do the integral                  c
c                                                                     c
c  variables we are using:                                            c
c       a, chebyshev coefficient of integrand f                       c
c       b, chebyshev coefficient of integral FF                       c
c                                                                     c
c  internal subroutines used:		   		              c
c    	cheby_f(N, M, in, out)            	                      c
c  	make_b(N, a, b)           			              c
c       intey
c	cheby_b(N, M, in, out)                                         c
c								      c
c ------------------------------------------------------------------- c


        subroutine inteyZ(Ny, Nz, Mz, L, planz_f, planz_b, f, FFF)

c	declare inputs/output variables
	integer*8 Ny, Nz, Mz, planz_f, planz_b
        real*8 L, f(0:Nz,1:Ny), FFF

c	declare subroutine specific variables
        real*8 pi, temp_f(0:Nz), temp_FF(0:Nz), FF(1:Ny)
        complex*16 a(0:Nz), b(0:Nz)

        parameter (pi = 2.0d0*dacos(0.0d0))

        do j = 1,Ny

           do i = 0,Nz
              temp_f(i) = f(i,j)
           end do

c          cheb forward transform 
           call cheby_f(Nz, Mz, planz_f, temp_f, a)

c          make chebyshev coefficient for FF
           call make_b(Nz, a, b)

c	   cheb backward transform 
           call cheby_b(Nz, Mz, planz_b, b, temp_FF)
           
           FF(j) = temp_FF(0) - temp_FF(Nz)
           
        end do
        
        FFF = sum(FF)*L/dble(Ny)
        
        end
      

c *********************************************************************
c **************************** SUB ROUTINES ***************************
c *********************************************************************

c ************************** CHEB FORWARD *****************************
	subroutine cheby_f(N, M, plan, in, out)
c	compute chebyshev transform

c	define variables
	integer*8 N, M, planz_f
        real*8 in(0:N), in_t(M)
        complex*16 in_t_hat(M), out(0:N)

c	extend in to in_t_2N-j = in_j, j = 1,2,...,N-1
	do i = 1,N
	  in_t(i) = in(i-1)
	  in_t(M-i+1) = in(i)
	end do

c	compute fft and remove plan
	call dfftw_execute_dft_r2c(plan, in_t, in_t_hat)

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


c ****************************** makeb ********************************
        subroutine make_b(N, a, b)
c       create coefficient b

c	define variables
        integer*8 N
        complex*16 a(0:N), b(0:N)
        
c       create b using the relationship with a
        b(0) = dcmplx(0.0d0, 0.0d0)
        b(1) = a(0) - dcmplx(0.5d0,0.0d0)*a(2)
        b(N) = dcmplx(1.0d0/(2.0d0*dble(N)), 0.0d0)*a(N-1)

        do i = 2,N-1
          b(i) = dcmplx(1.0d0/(2.0d0*dble(i)), 0.0d0)*(a(i-1)-a(i+1))
        end do

        end

c ************************** CHEB BACKWARD ****************************
	subroutine cheby_b(N, M, plan, in, out)
c	compute inverse chebyshev transform

c	define variables
	integer*8 N, M, plan
	complex*16 in(0:N), in_t(M) 
        real*8 out(0:N), out_t(M)

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
	call dfftw_execute_dft_c2r(plan, in_t, out_t)

c	only keep N+1 terms
	do i = 0,N
	  out(i) = out_t(i+1)
	end do

	end
