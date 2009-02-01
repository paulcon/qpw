CC = gcc
FLAGS = -Wall -O2 -g
LIB = -lm
INC = -I/home/paulcon/CLAPACK/F2CLIBS -I/home/paulcon/CLAPACK
SRC = ${wildcard *.c}
OBJS = ${SRC:%.c=%.o}
NAME = qpw

all: ${NAME}

${NAME}: ${OBJS}
	${CC} ${FLAGS} ${OBJS} -o ${NAME} ${LIB}

%.o: %.c
	${CC} ${FLAGS} -c $< -o $@

clean:
	rm -rf *.o ${NAME}