VPATH=alglib-3.15/src
OBJ:= linalg.o alglibmisc.o alglibinternal.o ap.o solvers.o statistics.o specialfunctions.o arPLS.o main.o

arPLS: $(OBJ)
	g++ -g -Wall -o $@ $^

generateTestData: generateTestData.cpp ap.o
	g++ -g -Wall -Ialglib-3.15/src -o $@ $^

$(OBJ): %.o : %.cpp
	g++ -g -Wall -Ialglib-3.15/src -c -o $@ $<



