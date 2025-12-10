

#use clang	mingw32-make x64
#use gcc	mingw32-make x64 USEGCC=1


ifeq ($(USEGCC), 1)
CC = gcc
else
CC = clang
endif

CFLAGS = -Wall -O3 -s

ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif


x64:	l1c-sse41  l1c-avx2  l1c-avx512
#l1c-sse41-linux  l1c-avx2-linux
x86:	l1c-ssse3-x86  l1c-sse41-x86


l1c-avx2:	l1c-avx2.c
	$(CC) $(CFLAGS) -mavx2 -mlzcnt -Dcodec_l1_avx2=main l1c-avx2.c -o l1c-avx2 $(LIBS)

l1c-sse41:	l1c-sse41.c
	$(CC) $(CFLAGS) -msse4.1 -Dcodec_l1_sse41=main l1c-sse41.c -o l1c-sse41 $(LIBS)

l1c-sse41-x86:	l1c-sse41.c
	$(CC) $(CFLAGS) -m32 -msse4.1 -Dcodec_l1_sse41=main l1c-sse41.c -o l1c-sse41-x86 $(LIBS)


l1c-sse41-linux:	l1c-sse41.c
	$(CC) $(CFLAGS) --target=x86_64-pc-linux-gnu -msse4.1 -Dcodec_l1_sse41=main l1c-sse41.c -o l1c-sse41 $(LIBS)

l1c-avx2-linux:		l1c-avx2.c
	$(CC) $(CFLAGS) --target=x86_64-pc-linux-gnu -mavx2 -mlzcnt -Dcodec_l1_avx2=main l1c-avx2.c -o l1c-avx2-linux $(LIBS)


#incompatible with the rest	TODO
l1c-ssse3-x86:	l1c-ssse3.c
	$(CC) $(CFLAGS) -m32 -mssse3 -Dcodec_l1_ssse3=main l1c-ssse3.c -o l1c-ssse3-x86 $(LIBS)


#incompatible with the rest	TODO
l1c-avx512:	l1c-avx512.c
	$(CC) $(CFLAGS) -mavx512f -mavx512bw -mavx512cd -mavx512dq -mavx512vl -mavx512vbmi2 -mbmi -mlzcnt -Dcodec_l1_avx512=main l1c-avx512.c -o l1c-avx512 $(LIBS)

