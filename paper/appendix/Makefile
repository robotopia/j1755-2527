CC = gcc
CFLAGS = -fPIC

TARGETS = erfcxinv.so

erfcxinv.so: erfcxinv.o
	$(LINK.c) -shared $^ -o $@

clean:
	$(RM) *.so *.o
