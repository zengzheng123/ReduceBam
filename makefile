OBJS = ReduceBam.o VariantFile.o
CC = g++
CFLAGS = -Wall -c
LFLAGS = -Wall -o3
GZFLAGS = -I./gzstream/include/ -L./gzstream/lib/ -lgzstream -lz
BAMFLAGS = -I./bamtools-master/include/ -L./bamtools-master/lib/ -lbamtools

ReduceBam: $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -o ReduceBam $(GZFLAGS) $(BAMFLAGS)

ReduceBam.o: ReduceBam.cpp VariantFile.h
	$(CC) $(CFLAGS) $(GZFLAGS) $(BAMFLAGS) ReduceBam.cpp

VariantFile.o: VariantFile.cpp VariantFile.h
	$(CC) $(CFLAGS) $(GZFLAGS) VariantFile.cpp

clean:
	rm -f *.o ReduceBam

