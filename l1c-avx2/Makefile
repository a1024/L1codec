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


l1c-clang:
	clang $(FLAGS) -mavx2 -mbmi -mlzcnt  l1c.c main.c util.c  -o l1c $(LIBS)

l1c-gcc:
	gcc $(FLAGS) -mavx2 -mbmi -mlzcnt  l1c.c main.c util.c  -o l1c $(LIBS)
