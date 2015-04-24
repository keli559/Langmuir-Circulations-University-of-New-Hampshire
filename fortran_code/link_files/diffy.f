c ------------------------------------------------------------------- c
c  			    y differentiation subroutine                    c
c  inputs:                                                            c
c    N - grid size                                                    c
c    k - periodicity                                                  c
c    v - input matrix                                                 c
c                                                                     c
c  outputs:                                                           c
c    w - y deriv matrix                                               c
c                                                                     c
c  matrix structure                                                   c
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

	subroutine diffy(Ny, Nz, k, plany_f, plany_b, v, w)

c	declare inputs/output variables
	integer*8 Ny, Nz, plany_f, plany_b
	real*8 k, v(0:Nz,1:Ny), w(0:Nz,1:Ny)

c	declare subroutine specific variables
	real*8 temp_v(1:Ny), temp_w(1:Ny)
	complex*16 n_w(1:Ny), k_c, v_hat(1:Ny), w_hat(1:Ny)

c	define wavenumber vector - order (0:N/2-1 N/2 -N/2+1:-1)
	do i = 1,Ny/2+1
	  n_w(i) = dcmplx(dble(i)-1.0d0, 0.0d0)
	end do
	j = 1;
	do i = Ny/2+2,Ny
	  n_w(i) = dcmplx(dble(-Ny)/2.0d0+dble(j), 0.0d0)
	  j = j+1
	end do

c	convert k --> k_c
	k_c = dcmplx(k, 0.0d0)

c	compute w each z location
	do i = 0,Nz
c	  store v vec at each z location
	  do j = 1,Ny
	    temp_v(j) = v(i,j)
	  end do

c	  fft v to v_hat
	  call fft_f(Ny, plany_f, temp_v, v_hat)

c	  determine w_hat = (i*n_w*k)*v_hat
	  do j = 1,Ny
	    w_hat(j) = (sqrt((-1.0d0, 0.0d0))*n_w(j)*k_c)*v_hat(j)
	  end do

c 	  correct N/2 mode 
	  w_hat(Ny/2+1) = (0.0d0, 0.0d0)

c	  ifft w_hat to temp_w 
	  call fft_b(Ny, plany_b, w_hat, temp_w)

c	  combine each temp_w vec to get w matrix
	  do j = 1,Ny
	    w(i,j) = temp_w(j)/dble(Ny)
	  end do
	end do 

	end

