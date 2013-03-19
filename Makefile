export CPP=g++

INC = -I/usr/local/include/boost
C_FLAGS = -g -pthread -c -O3 $(INC)

OBJS = driver.o FitnessTranslation.o GAIndividual.o CGA.o EvaluationThread.o

garun: $(OBJS)
	$(CPP) -g -pthread -o ./bin/garun $(OBJS)

driver.o: driver.cpp
	$(CPP) $(C_FLAGS) driver.cpp

FitnessTranslation.o: FitnessTranslation.cpp
	$(CPP) $(C_FLAGS) FitnessTranslation.cpp

EvaluationThread.o: EvaluationThread.cpp
	$(CPP) $(C_FLAGS) EvaluationThread.cpp

GAIndividual.o: GAIndividual.cpp
	$(CPP) $(C_FLAGS) GAIndividual.cpp

CGA.o: CGA.cpp
	$(CPP) $(C_FLAGS) CGA.cpp

clean:
	rm -f *.o *~ ./bin/garun



