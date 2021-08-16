all : Rfuncs

Rfuncs:
	R CMD SHLIB glasso.c
	R CMD SHLIB SFGen.c -o SFGen.so

clean:
	rm -f *.o *.so

