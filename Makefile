CC=mpic++
CFLAGS=-O9 -std=c++11 -fopenmp

all: cliqueParallel cliqueShared

cliqueParallel : kclique_distributed.cpp
	$(CC) $(CFLAGS) kclique_distributed.cpp -o cliqueParallel

cliqueShared : kclique_shared.cpp
	$(CC) $(CFLAGS) kclique_shared.cpp -o cliqueShared 

clean:
	rm cliqueParallel
	rm cliqueShared
