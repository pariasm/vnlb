CC=c99
CFLAGS=-O3 -I../iio/
OMPFLAGS=-fopenmp

tvl1flow: main.c tvl1flow_lib.c bicubic_interpolation.c mask.c zoom.c iio.o
	$(CC) $(CFLAGS) $(OMPFLAGS) -o tvl1flow main.c iio.o -ltiff -lpng -ljpeg -lm

vqm: main_mesure.c tvl1flow_lib.c bicubic_interpolation.c mask.c zoom.c iio.o
	$(CC) $(CFLAGS) $(OMPFLAGS) -o vqm main_mesure.c iio.o -ltiff -lpng -ljpeg -lm

iio.o: ../iio/iio.c
	$(CC) $(CFLAGS) -DNDEBUG -D_GNU_SOURCE -c ../iio/iio.c

clean:
	rm -f iio.o main.o tvl1flow vqm
