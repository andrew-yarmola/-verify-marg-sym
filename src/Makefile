CC := g++ -std=c++1y -Wall -mfpmath=sse -mieee-fp -mno-fma -mno-fma4
CXX := g++ -std=c++1y  -Wall -mfpmath=sse -mieee-fp -mno-fma -mno-fma4
COMMON := arithmetic/AJCC.o arithmetic/SL2.o box.o elimination.o arithmetic/Complex.o arithmetic/roundoff.o

verify: verify.o $(COMMON)
	$(CC) verify.o $(COMMON) -lm -o ../bin/verify

rootcat: rootcat.o
	$(CC) rootcat.o -o ../bin/rootcat

# never optimize here, since the timing of arithmetic is important
tests: test_underflow.o test_overflow.o test_powers.o roundoff.o
	$(CC) test_underflow.o roundoff.o -o ../bin/test_underflow
	$(CC) test_overflow.o roundoff.o -o ../bin/test_overflow
	$(CC) test_powers.o roundoff.o -o ../bin/test_powers

