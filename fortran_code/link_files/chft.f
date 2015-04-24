c ************************ Chebyshev transform ************************
	subroutine chft(N, in, out)

c	Chebyshev transform: real physical -> real spectrum
c	Input:
c		in		- colume vector with N+1 elements
c		N		- resolution of "in"
c	Output:
c		out		- vector after transform, 
c				  with N+1 elements
c  	external subroutines used:                                     
c   		fft_f(N, in, out) 
        integer*8 N
	real*8 in(0:N,1), out(0:N,1), inflip(0:N,1), temp(1:2*N,1) 
     ^		, ck(0:N,1)
	complex*16 temp1(1:2*N,1)
	integer*8 ii

c	flip "in" over, store it to "inflip"
	do ii = 0, N
		inflip(ii,1) = in(N-ii,1)
	end do

c	combine "in" and mirror image "inflip" into "temp"
	do ii = 0, N 
		temp(ii+1,1) = in(ii,1)
	end do
	do ii = 2, N
		temp(N+ii,1) = inflip(ii-1,1)
	end do

c	fator "ck"
	do ii = 0, N
		ck(ii,1) = 1.0d0
	end do
	ck(0,1) = ck(0,1)/2.0d0
	ck(N,1) = ck(N,1)/2.0d0

c	FFT on "temp"
	call fft_f(2*N, temp, temp1)

c	Get real part out of "temp1", store it into "temp"
	do ii = 1, 2*N
		temp(ii,1) = dreal(temp1(ii,1))/dble(N)
	end do

c	calculate final result "out"
	do ii = 0, N
		out(ii,1) = temp(ii+1,1)*ck(ii,1) 
	end do	
	end
