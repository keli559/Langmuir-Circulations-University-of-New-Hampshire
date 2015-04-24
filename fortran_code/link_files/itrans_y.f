c ------------------------------------------------------------------- c
c              y fourier inverse transform subroutine                 c
c                                                                     c
c  inputs:                                                            c
c    N - grid size                                                    c
c    in - input matrix                                                c
c    d - input dimension: 1 - 1d, 2 - 2d                              c
c                                                                     c
c  outputs:                                                           c
c    out - y inverse transformed matrix                               c
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

	subroutine itrans_y(Ny, Nz, plany_b, in, out)

c	declare inputs/output variables
	integer*8 Ny, Nz, plany_b
	real*8 out(0:Nz,1:Ny)
	complex*16 in(0:Nz,1:Ny)
	
c	declare subroutine specific variables
	real*8 temp_out(1:Ny)
	complex*16 temp_in(1:Ny)

c	compute w each z location
	do i = 0,Nz
c	  store v vec at each z location
	  do j = 1,Ny
	    temp_in(j) = in(i,j)
	  end do

c	  fft temp_in to temp_out
	  call fft_b(Ny, plany_b, temp_in, temp_out)

c	  combine temp_out to get out, correct for fft normalization
	  do j = 1,Ny
	    out(i,j) = temp_out(j)
	  end do
	end do

	end
