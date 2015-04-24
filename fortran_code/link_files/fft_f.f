	
	subroutine fft_f(N, plan, in, out)
c	compute fft of in

c	define variables
c	rows --> z
c	cols --> y
	integer*8 N
	real*8 in(1:N)
	complex*16 out(1:N)

	integer*8 plan
	complex*16 temp_out(1:N/2+1)

c	compute fft and remove plan
c	call dfftw_plan_dft_r2c_1d(plan,N,in,temp_out,estimate)
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
