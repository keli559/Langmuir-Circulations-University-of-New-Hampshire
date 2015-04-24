c ------------------------------------------------------------------- c
c  			  y fourier transform subroutine                    c
c                                                                     c
c  inputs:                                                            c
c    N - grid size                                                    c
c    in - input matrix                                                c
c    d - input dimension: 1 - 1d, 2 - 2d                              c
c                                                                     c
c  outputs:                                                           c
c    out - y transform matrix                                         c
c                                                                     c
c  matrix structure                                                   c
c    row  -->  varying y, Z = const                                   c
c    coln  -->  varying Z, y = const                                  c
c    e.g.: A(i,j) = A(Z,y): i [0,N], j [1,N]                          c
c                                                                     c
c  external subroutines used:                                         c
c    fft_f(N, in, out)                                                c
c    fft_b(N, in, out)                                                c
c                                                                     c
c ------------------------------------------------------------------- c

	subroutine trans_y(Ny, Nz, plany_f, in, out)

c	declare inputs/output variables
	integer*8 Ny, Nz, plany_f
	real*8 in(0:Nz,1:Ny)
	complex*16 out(0:Nz,1:Ny)

c	declare subroutine specific variables
	real*8 temp_in(1:Ny)
	complex*16 temp_out(1:Ny)

c	compute w each z location
	do i = 0,Nz
c	  store v vec at each z location
	  do j = 1,Ny
	    temp_in(j) = in(i,j)
	  end do

c	  fft temp_in to temp_out
	  call fft_f(Ny, plany_f, temp_in, temp_out)

c	  combine each temp_out to get out
	  do j = 1,Ny
	    out(i,j) = temp_out(j)/dcmplx(dble(Ny), 0.0d0)
	  end do
	end do

	end
