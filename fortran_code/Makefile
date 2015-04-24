# lc make file
#
# to compile:
#    $make run
# to execute:
#    $./run
# to clean up files
#    $make clean

# name
NAME = run

# compiler
FC = gfortran

# all files in program
OBJECTS = KE.o slv_helm2.o diffy.o diffZ.o fft_f.o fft_b.o itrans_y.o trans_y.o plan_init.o inteyZ.o

# library location
LIB = -L/home/iam/zrg2/lib/

# linking header files
#UMFLIB = /usr/local/UMFPACK/Lib/libumfpack.a
#AMDLIB = /usr/local/AMD/Lib/libamd.a
FFTLIB = -lfftw3
#LAPACKLIB = -llapack -lrefblas

# run program
$(NAME): $(OBJECTS)
#	$(FC) -o $(NAME) $(OBJECTS) $(LIB) $(FFTLIB) $(UMFLIB) $(AMDLIB) $(LAPACKLIB) -lm
	$(FC) -o $(NAME) $(OBJECTS) $(LIB) $(FFTLIB) -lm
	$ clear

KE.o: KE.f
	$(FC) -c KE.f

slv_helm2.o: link_files/slv_helm2.f
	$(FC) -c link_files/slv_helm2.f

diffy.o: link_files/diffy.f
	$(FC) -c link_files/diffy.f

diffZ.o: link_files/diffZ.f
	$(FC) -c link_files/diffZ.f

fft_f.o: link_files/fft_f.f
	$(FC) -c link_files/fft_f.f

fft_b.o: link_files/fft_b.f
	$(FC) -c link_files/fft_b.f

trans_y.o: link_files/trans_y.f
	$(FC) -c link_files/trans_y.f

itrans_y.o: link_files/itrans_y.f
	$(FC) -c link_files/itrans_y.f

plan_init.o: link_files/plan_init.f
	$(FC) -c link_files/plan_init.f

inteyZ.o: link_files/inteyZ.f
	$(FC) -c link_files/inteyZ.f

rerun:
	rm -f $(NAME)

data:
	rm -f ../data/*.*

clean:
	rm -f *.o *~* $(NAME)  
