c ------------------------------------------------------------------- c
c  			    Z differentiation subroutine                    c
c  inputs:                                                            c
c    N - grid size                                                    c
c    k - periodicity                                                  c
c    u - input matrix                                                 c
c                                                                     c
c  outputs:                                                           c
c    w - y deriv matrix                                               c
c                                                                     c
c  matrix structure:                                                  c
c    row  -->  varying y, Z = const                                   c
c    coln  -->  varying Z, y = const                                  c
c    e.g.: A(i,j) = A(Z,y): i [0,N], j [1,N]                          c
c                                                                     c
c  method and notation adopted from:                                  c
c    Trefethen, L. N. Spectral Methods in Matlab. 2000                c
c                                                                     c
c  external subroutines used:                                         c
c    fft_f(N, in, out)                                                c
c    fft_b(N, in, out)                                                c
c                                                                     c
c ------------------------------------------------------------------- c

	subroutine diffZ(Ny, Nz, Mz, planz_f, planz_b, u, w)

c	declare inputs/output variables
	integer*8 Ny, Nz, Mz, planz_f, planz_b
	real*8 u(0:Nz,1:Ny), w(0:Nz,1:Ny)

c	declare subroutine specific variables
	real*8 pi, x(0:Nz), V(1:Mz), WW(1:Mz)
	complex*16 n_w(1:Mz), V_hat(1:Mz), W_hat(1:Mz)

c	define pi
	pi = 2.0d0*dacos(0.0d0)

c	define x - cheb points
	do i = 0,Nz
	  x(i) = dcos(pi*dble(i)/dble(Nz))
	end do

c	wavenumber vector - order (0:N/2 -N/2+1:-1)
	do i = 1,Mz/2+1
	  n_w(i) = dcmplx(dble(i)-1.0d0, 0.0d0)
	end do
	l = 1
 	do i = Mz/2+2,Mz
 	  n_w(i) = dcmplx(dble(-Mz)/2.0d0 + dble(l), 0.0d0)
 	  l = l+1
 	end do

	do i = 1,Mz
	  WW(i) = (0.0d0, 0.0d0)
	end do

c	compute w at each y location
	do j = 1,Ny
c	  extend v to V_2N-j = v_j, j = 1,2,...,N-1
	  do i = 1,Nz
	    V(i) = u(i-1,j)
	    V(Mz-i+1) = u(i,j)
	  end do
	  
c	  fft V to V_hat
	  call fft_f(Mz, planz_f, V, V_hat)
	  
c 	  compute W_hat
	  do i = 1,Mz
 	    W_hat(i) = (sqrt((-1.0d0,0.0d0))*n_w(i))*V_hat(i)
 	  end do
	  
c	  set W_hat(k=N/2) = 0
	  W_hat(Mz/2+1) = (0.0d0, 0.0d0)
	  
c	  ifft W_hat to W
	  call fft_b(Mz, planz_b, W_hat, WW)
	  
c	  correct fft factor
	  do i = 1,Mz
	    WW(i) = WW(i)/dble(Mz)
	  end do
	  
c	  transfrom from theta --> x
	  do i = 1,Nz-1
	    w(i,j) = -WW(i+1)/sqrt(1.0d0 - x(i)*x(i))
	  end do
	  
c	  correct end points
	  w(0,j) = 0.0d0
	  do i = 0,Nz-1
	    w(0,j) = w(0,j) + dble(i)*dble(i)*dreal(V_hat(i+1))
	  end do
	  w(0,j) = w(0,j) + 0.5d0*dble(Nz)*dble(Nz)*dreal(V_hat(Nz+1))
	  w(0,j) = w(0,j)/dble(Nz)
	   
	  w(Nz,j) = 0.0d0
	  do i = 0,Nz-1
	    w(Nz,j) = w(Nz,j) + ((-1.0d0)**(dble(i)+1.0d0))*dble(i)
     ^		*dble(i)*dreal(V_hat(i+1))
	  end do
	  w(Nz,j) = w(Nz,j) + dcmplx(0.5d0*((-1.0d0)**(dble(Nz)+1.0d0))
     ^		*dble(Nz)*dble(Nz), 0.0d0)*dreal(V_hat(Nz+1))
	  w(Nz,j) = w(Nz,j)/dble(Nz)
	end do

	end


