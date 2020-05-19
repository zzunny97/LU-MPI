all:
	mpic++ -O3 -g -Wall -o LU_par LU.cpp -lm

seq:
	g++ LU_seq.cpp -o LU_seq

clean:
	rm -rf LU_seq LU_par 
