CC = gcc
FLAGS = -Wall -O2 -g
LIB = -lm -lcblas -llapack
SRC = main.c quad_functions.c spectral_functions.c
OBJS = ${SRC:%.c=%.o}
TEST_SRC = tester.c quad_functions.c spectral_functions.c unit_tests.c
TEST_OBJS = ${TEST_SRC:%.c=%.o} 
NAME = qpw
TEST = qpw_tester

all: ${NAME}

test: ${TEST}

${NAME}: ${OBJS}
	${CC} ${FLAGS} ${OBJS} -o ${NAME} ${LIB}

${TEST}: ${TEST_OBJS}
	${CC} ${FLAGS} ${TEST_OBJS} -o ${TEST} ${LIB}

%.o: %.c
	${CC} ${FLAGS} -c $< -o $@

clean:
	rm -rf *.o ${NAME} ${TEST}