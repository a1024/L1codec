# release:	make l1c
# profiler:	make l1c PROF=1
# debug:	make l1c DEBUG=1


ifeq ($(DEBUG), 1)
FLAGS = -Wall -O3 -g -D_DEBUG
else
FLAGS = -Wall -O3 -s
endif

ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif


l1csse41:
	clang $(FLAGS) -msse4.1  l1c.c main.c util.c  -o l1csse41 $(LIBS)

l1c32:
	gcc $(FLAGS) -m32 -msse4.1  l1c.c main.c util.c  -o l1c32 $(LIBS)

l1c-gcc:
	gcc $(FLAGS) -msse4.1  l1c.c main.c util.c  -o l1csse41 $(LIBS)
