c ************************Inverse Chebyshev transform *****************
	subroutine chifft(N, in, out)

c	Chebyshev transform: real physical -> real spectrum
c	Input:
c		in		- colume vector with N+1 elements
c		N		- resolution of "in"
c	Output:
c		out		- vector after inverse transform, 
c				  with N+1 elements
c  	external subroutines used:                                     
c   		fft_b(N, in, out) 
        integer*8 N
	real*8 in(0:N,1), out(0:N,1), inflip(0:N,1), temp(1:2*N,1)
    	complex*16 temp1(1:2*N,1)
	integer*8 ii

	in(0,1) = 2.0d0 * in(0,1)
	in(N,1) = 2.0d0 * in(N,1)
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
c	transform "temp" into complex number "temp2"
	do ii = 1, 2*N
		temp1(ii,1) = dcmplx(temp(ii,1),0.0d0)
	end do

c	inverse fft of "temp1", store into "temp"
	call fft_b(2*N, temp1, temp)
c	calculate final result "out"
	do ii = 0, N
		out(ii,1) = temp(ii+1,1)/2.0d0 
	end do	
	end
