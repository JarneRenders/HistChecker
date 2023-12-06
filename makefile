compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces -O3
profileflags=-std=gnu11 -march=native -Wall -fsanitize=address -g -pg

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.
64bit: histChecker.c readGraph/readGraph6.c bitset.h 
	$(compiler) -DUSE_64_BIT -o histChecker histChecker.c readGraph/readGraph6.c $(flags)

128bit: histChecker.c readGraph/readGraph6.c bitset.h 
	$(compiler) -DUSE_128_BIT -o histChecker-128 histChecker.c readGraph/readGraph6.c $(flags)

128bitarray: histChecker.c readGraph/readGraph6.c bitset.h 
	$(compiler) -DUSE_128_BIT_ARRAY -o histChecker-128a histChecker.c readGraph/readGraph6.c $(flags)

profile: histChecker.c readGraph/readGraph6.c bitset.h 
	$(compiler) -DUSE_64_BIT -o histChecker-pr histChecker.c readGraph/readGraph6.c $(profileflags)

all: 64bit 128bit 128bitarray

.PHONY: clean
clean:
	rm -f histChecker histChecker-pr histChecker-128 histChecker-128a gmon.out

