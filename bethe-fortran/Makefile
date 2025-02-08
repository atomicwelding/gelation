SRC = main.f90 mod_utils.f90
LDFLAG = -framework Accelerate

.PHONY : all

all: $(SRC)
	gfortran -c mod_utils.f90
	gfortran -c main.f90
	gfortran $(SRC) -o main $(LDFLAG)


clean:
	rm *.o *.mod
	rm main

clean_data:
	rm *.txt
