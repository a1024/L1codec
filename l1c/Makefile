

#use clang	mingw32-make x64
#use gcc	mingw32-make x64 USEGCC=1


ifeq ($(USEGCC), 1)
CC = gcc
else
CC = clang
endif

ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif


l1c:		l1c-avx2.c
	$(CC) -Wall -O3 -s -mavx2 -mlzcnt -Dcodec_l1_avx2=main l1c-avx2.c -o l1c $(LIBS)

l1cx86:		l1c-x86.c
	$(CC) -Wall -O3 -s -m32 -msse4.1 -Dcodec_l1_x86=main l1c-x86.c -o l1cx86 $(LIBS)

l1cavx512:	l1c-avx512.c
	$(CC) -Wall -O3 -s -mavx512f -mavx512bw -mavx512cd -mavx512dq -mavx512vl -mavx512vbmi2 -mbmi -mlzcnt -Dcodec_l1_avx512=main l1c-avx512.c -o l1cavx512 $(LIBS)

x64:	l1c  l1cavx512
x86:	l1cx86

