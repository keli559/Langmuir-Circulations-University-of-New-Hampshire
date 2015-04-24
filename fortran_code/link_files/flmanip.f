c*************************FILTER MANIPULATION FUNCTION****************C
	subroutine flmanip(ur, or, pr, tr, Ny, Nz, L, m0, n0, r0   
     ^		, u, o, p, t)
c	Input:
c		ur, or, pr, tr		- oiriginal physical data for filtering
c		L			- horizontal physical domain length
c		r0			- boundary on m (for y spectrum)
c		
c	Output;
c		u, o ,p, t 		- physical data after filtering
c*********************************************************************c

c	0. Declaration
	integer*8 Nz, Ny
	real*8 ur(0:Nz,1:Ny), or(0:Nz,1:Ny), pr(0:Nz,1:Ny),    
     ^		tr(0:Nz,1:Ny), u(0:Nz,1:Ny), o(0:Nz,1:Ny),    
     ^		p(0:Nz,1:Ny), t(0:Nz,1:Ny), filt(0:Nz,1:Ny)
	real*8 m0, n0, r0
	complex*16 ucf(0:Nz,1:Ny), ocf(0:Nz,1:Ny), pcf(0:Nz,1:Ny),    
     ^		tcf(0:Nz,1:Ny), up(0:Nz,1:Ny), op(0:Nz,1:Ny),    
     ^		pp(0:Nz,1:Ny), tp(0:Nz,1:Ny)
	integer*8 ii, jj, kk
c	1. Physical -> Spcetrum
	call spec_cf(ur, or, pr, tr, Ny, Nz, ucf, ocf, pcf, tcf)

c	2. flip over
c		from 	[0,...Ny/2, -Ny/2+1,...-1]
c		to	[-Ny/2+1,..., -1, 0, ...Ny/2]
	call spflip(ucf, ocf, pcf, tcf, Ny, Nz, up, op, pp, tp)

c	3. filter function construction
	call filtfun(Ny, Nz, m0, n0, r0, filt)

c	4. filtering
	do ii = 0, Nz
		do jj = 1, Ny
			up(ii, jj) = dcmplx(dreal(up(ii, jj))   
     ^		*filt(ii, jj), dimag(up(ii, jj))*filt(ii, jj))
			op(ii, jj) = dcmplx(dreal(op(ii, jj))   
     ^		*filt(ii, jj), dimag(op(ii, jj))*filt(ii, jj))
			pp(ii, jj) = dcmplx(dreal(pp(ii, jj))  
     ^		*filt(ii, jj), dimag(pp(ii, jj))*filt(ii, jj))
			tp(ii, jj) = dcmplx(dreal(tp(ii, jj))   
     ^		*filt(ii, jj), dimag(tp(ii, jj))*filt(ii, jj))
		end do
	end do

c	5. flip back to Matlab conventional sequence on m (for y spectrum)
	call ispflip(up, op, pp, tp, Ny, Nz, ucf, ocf, pcf, tcf)

c	6. Spectral -> Physical
	call spec_infc(ucf, ocf, pcf, tcf, Ny, Nz, u, o,  
     ^		 p, t)
	end
c *********************************************************************
c **************************** SUB ROUTINES ***************************
c *********************************************************************

c*************************Chebyshev-Fourier transform*****************C
	subroutine spec_cf(u, o, p, t, Ny, Nz, ucf, ocf, pcf, tcf)
c	Input:
c		u, o, p, t		- original physical data for filtering
c		Ny, Nz			- resolutions at y and Z
c		
c	Output;
c		ucf, ocf, pcf, tcf	- spectral variables
c*********************************************************************c

c	0. Declaration
	integer*8 Ny, Nz
	real*8 u(0:Nz,1:Ny), o(0:Nz,1:Ny), p(0:Nz,1:Ny), t(0:Nz,1:Ny)
	complex*16 uf(0:Nz,1:Ny), of(0:Nz,1:Ny),
     ^		   pf(0:Nz,1:Ny), tf(0:Nz,1:Ny), 
     ^		   ucf(0:Nz,1:Ny), ocf(0:Nz,1:Ny), 
     ^		   pcf(0:Nz,1:Ny), tcf(0:Nz,1:Ny)
	integer*8 ii, jj
c	1. FFT on y
	do ii = 0, Nz
		call fft_f(Ny, u(ii, 1:Ny), uf(ii, 1:Ny))
		call fft_f(Ny, o(ii, 1:Ny), of(ii, 1:Ny))
		call fft_f(Ny, p(ii, 1:Ny), pf(ii, 1:Ny))
		call fft_f(Ny, t(ii, 1:Ny), tf(ii, 1:Ny))
	end do

c	2. Chebyshev on Z
	do jj = 1, Ny
		call cheb(Nz, uf(0:Nz, jj), ucf(0:Nz, jj))
		call cheb(Nz, of(0:Nz, jj), ocf(0:Nz, jj))
		call cheb(Nz, pf(0:Nz, jj), pcf(0:Nz, jj))
		call cheb(Nz, tf(0:Nz, jj), tcf(0:Nz, jj))
	end do
	end

c*************************Spectrum Flipping on y**********************C
	subroutine spflip(u, o, p, t, Ny, Nz, up, op, pp, tp)
c	Input:
c		u, o, p, t		- oiriginal SPECTRAL data
c		Ny, Nz			- resolutions at y and Z
c		
c	Output;
c		up, op ,pp, tp 		- flipped spectral data
c*********************************************************************c
	integer*8 Ny, Nz
	complex*16 u(0:Nz,1:Ny), o(0:Nz,1:Ny),
     ^		   p(0:Nz,1:Ny), t(0:Nz,1:Ny),
     ^		   up(0:Nz,1:Ny), op(0:Nz,1:Ny),
     ^		   pp(0:Nz,1:Ny), tp(0:Nz,1:Ny),
     ^		   temp1(0:Nz, 1:Ny/2+1), temp2(0:Nz, 1:Ny/2-1)
	integer*8 ii, jj
c	******** u *********
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			temp1(ii, jj) = u(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			temp2(ii, jj) = u(ii, jj+Ny/2+1)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			up(ii, jj) = temp2(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			up(ii, jj+Ny/2-1) = temp1(ii, jj)
		end do
	end do
c	******** omega *********

	do ii = 0, Nz
		do jj = 1, Ny/2+1
			temp1(ii, jj) = o(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			temp2(ii, jj) = o(ii, jj+Ny/2+1)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			op(ii, jj) = temp2(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			op(ii, jj+Ny/2-1) = temp1(ii, jj)
		end do
	end do
c	******** psi *********

	do ii = 0, Nz
		do jj = 1, Ny/2+1
			temp1(ii, jj) = p(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			temp2(ii, jj) = p(ii, jj+Ny/2+1)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			pp(ii, jj) = temp2(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			pp(ii, jj+Ny/2-1) = temp1(ii, jj)
		end do
	end do
c	******** t *********

	do ii = 0, Nz
		do jj = 1, Ny/2+1
			temp1(ii, jj) = t(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			temp2(ii, jj) = t(ii, jj+Ny/2+1)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			tp(ii, jj) = temp2(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			tp(ii, jj+Ny/2-1) = temp1(ii, jj)
		end do
	end do
	end

c********************************FITLER FUNCTION**********************C
	subroutine filtfun(Ny, Nz, m0, n0, r0, filt)
c	Input:
c		Ny, Nz			- resolutions on y and Z
c		(m0, n0)		- center of spectrum with 
c					  respect of (y, Z)
c		r0 			- spectrum radius with 
c				  	  respect of y
c	Output;
c		filt			- filter function
c*********************************************************************c
	integer*8 Ny, Nz
	real*8 m0, n0, r0, filt(0:Nz, 1:Ny), m(1, 1:Ny), n(0:Nz,1)
	real*8 A, B
	integer*8 ii, jj
	call spset(Ny, Nz, m, n)
	
	A = 100.0d0
	B = 1.0d0

	do ii = 0, Nz
		do jj = 1, Ny
			filt(ii, jj) = 
     ^		   -0.5d0* derf(
     ^		   A/dble(Ny)*(
     ^		   dsqrt(
     ^		   (m(1, jj)-m0)**2.0d0
     ^		 + B*dble(Ny)**2.0d0/((dble(Nz)+1.0d0)**2.0d0)
     ^		 * (n(ii,1)-n0)**2.0d0
     ^		        )
     ^		   	    -r0
     ^		              )
     ^		               )+0.5d0
		end do
	end do
	end

c*****************Spectrum flipping back to Matlab Sequence***********C
	subroutine ispflip(u, o, p, t, Ny, Nz, up, op, pp, tp)
c	Input:
c		u, o, p, t		- spectral space variables 
c					  not physical
c		Ny, Nz			- resolutions on m (for y) 
c					  and n (for z) in spectral space
c		
c	Output;
c		up, op ,pp, tp 		- spectral data after flipping
c*********************************************************************c
	integer*8 Ny, Nz
	complex*16 u(0:Nz, 1:Ny), o(0:Nz, 1:Ny), p(0:Nz, 1:Ny), 
     ^		  t(0:Nz, 1:Ny), up(0:Nz, 1:Ny), op(0:Nz, 1:Ny), 
     ^		  pp(0:Nz, 1:Ny), tp(0:Nz, 1:Ny)
	complex*16 temp1(0:Nz, 1:Ny/2-1), temp2(0:Nz, 1:Ny/2+1)
	integer*8 ii, jj
c********** u ************
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			temp1(ii, jj) = u(ii, jj)
		end do
	end do 
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			temp2(ii, jj) = u(ii, jj+Ny/2-1)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			up(ii, jj) = temp2(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			up(ii, jj+Ny/2+1) = temp1(ii, jj)
		end do
	end do 
c********** omega ************
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			temp1(ii, jj) = o(ii, jj)
		end do
	end do 
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			temp2(ii, jj) = o(ii, jj+Ny/2-1)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			op(ii, jj) = temp2(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			op(ii, jj+Ny/2+1) = temp1(ii, jj)
		end do
	end do 
c********** psi ************
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			temp1(ii, jj) = p(ii, jj)
		end do
	end do 
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			temp2(ii, jj) = p(ii, jj+Ny/2-1)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			pp(ii, jj) = temp2(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			pp(ii, jj+Ny/2+1) = temp1(ii, jj)
		end do
	end do 
c********** t ************
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			temp1(ii, jj) = t(ii, jj)
		end do
	end do 
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			temp2(ii, jj) = t(ii, jj+Ny/2-1)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2+1
			tp(ii, jj) = temp2(ii, jj)
		end do
	end do
	do ii = 0, Nz
		do jj = 1, Ny/2-1
			tp(ii, jj+Ny/2+1) = temp1(ii, jj)
		end do
	end do 
	end

c******************Chebyshev-Fourier inverse transform****************C
	subroutine spec_infc(ufc, ofc, pfc, tfc, Ny, Nz, u, o,  
     ^		 p, t)
c	Input:
c		ufc, ofc, pfc, tfc	- spectral
c		Ny, Nz			- resolutions
c		
c	Output;
c		u, o ,p, t 		- physical data
c*********************************************************************c
	integer*8 Ny, Nz
	complex*16 ufc(0:Nz, 1:Ny), ofc(0:Nz, 1:Ny),  
     ^		 pfc(0:Nz, 1:Ny), tfc(0:Nz, 1:Ny), 
     ^		 uf(0:Nz, 1:Ny), of(0:Nz, 1:Ny),  
     ^		 pf(0:Nz, 1:Ny), tf(0:Nz, 1:Ny)
	real*8 u(0:Nz, 1:Ny), o(0:Nz, 1:Ny),  
     ^		 p(0:Nz, 1:Ny), t(0:Nz, 1:Ny)
	integer*8 ii, jj
	do jj = 1, Ny
		call icheb(Nz, ufc(0:Nz, jj), uf(0:Nz, jj))
		call icheb(Nz, ofc(0:Nz, jj), of(0:Nz, jj))
		call icheb(Nz, pfc(0:Nz, jj), pf(0:Nz, jj))
		call icheb(Nz, tfc(0:Nz, jj), tf(0:Nz, jj))
	end do

	do ii = 0, Nz
		call fft_b(Ny, uf(ii, 1:Ny), u(ii, 1:Ny))
		call fft_b(Ny, of(ii, 1:Ny), o(ii, 1:Ny))
		call fft_b(Ny, pf(ii, 1:Ny), p(ii, 1:Ny))
		call fft_b(Ny, tf(ii, 1:Ny), t(ii, 1:Ny))
	end do

	do ii = 0, Nz
		do jj = 1, Ny
			u(ii, jj) = u(ii, jj)/dble(Ny)
			o(ii, jj) = o(ii, jj)/dble(Ny)
			p(ii, jj) = p(ii, jj)/dble(Ny)
			t(ii, jj) = t(ii, jj)/dble(Ny)
		end do
	end do
	end

c********************Chebyshev Transform in complex regime************C
	subroutine cheb(N, v, vp)

c 	Input:
c		v			- complex colume vector 
c						in physical space
c		N			- resolution of v
c		
c	Output;
c		vp			- complex colume vector
c						in spectral space
c*********************************************************************c
	integer*8 N
	real*8 vre(0:N,1), vim(0:N,1), vpre(0:N,1), vpim(0:N,1)
	complex*16 v(0:N,1), vp(0:N,1)
	integer*8 ii
	
	do ii = 0, N
		vre(ii,1) = dble(v(ii,1))
		vim(ii,1) = dimag(v(ii,1))
	end do
	call chft(N, vre, vpre)
	call chft(N, vim, vpim)
	do ii = 0, N
		vp(ii,1) = dcmplx(vpre(ii,1), vpim(ii,1))
	end do
	end

c************Chebyshev inverse Transform in complex regime************C
	subroutine icheb(N, vp, v)

c 	Input:
c		vp			- complex colume vector 
c						in spectral space
c		N			- resolution of v
c		
c	Output;
c		v			- complex colume vector
c						in physical space
c*********************************************************************c
	integer*8 N
	real*8 vre(0:N,1), vim(0:N,1), vpre(0:N,1), vpim(0:N,1)
	complex*16 v(0:N,1), vp(0:N,1)
	integer*8 ii
	
	do ii = 0, N
		vpre(ii,1) = dble(vp(ii,1))
		vpim(ii,1) = dimag(vp(ii,1))
	end do
	call chifft(N, vpre, vre)
	call chifft(N, vpim, vim)
	do ii = 0, N
		v(ii,1) = dcmplx(vre(ii,1), vim(ii,1))
	end do
	end
c***************** SPECTRAL DOMAIN SETUP *****************************c
	subroutine spset(Ny, Nz, m, n)
	integer*8 Ny, Nz
	real*8 m(1, 1:Ny), n(0:Nz,1)
	integer*8 ii, jj
	do ii = 0, Nz
		n(ii,1) = dble(ii)
	end do
	do jj = 1, Ny
		m(1,jj) = dble(jj)-dble(Ny)/2.0d0
	end do
	end


