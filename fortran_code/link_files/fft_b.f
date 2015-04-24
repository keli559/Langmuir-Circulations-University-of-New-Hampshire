
	subroutine fft_b(N, plan, in, out)
c	compute inverse fft of in

c	define variables
c	rows --> z
c	cols --> y
	integer*8 N
	real*8 out(1:N)
	complex*16 in(1:N)

	integer*8 plan
	complex*16 temp_in(1:N/2+1)

	do i = 1,N/2+1
	  temp_in(i) = in(i)
	end do

c	call dfftw_plan_dft_c2r_1d(plan,N,temp_in,out,estimate)
	call dfftw_execute_dft_c2r(plan, temp_in, out)
c	call dfftw_destroy_plan(plan)	

	end
