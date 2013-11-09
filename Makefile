a.out: main.o InputDeck.o SourceIteration.o Utilities.o
	g++ -g main.o InputDeck.o SourceIteration.o Utilities.o -o a.out 

main.o: main.cpp InputDeck.h SourceIteration.h
	g++ -g -c main.cpp -o main.o 

InputDeck.o: InputDeck.cpp InputDeck.h Utilities.h
	g++ -g -c InputDeck.cpp -o InputDeck.o 

SourceIteration.o: SourceIteration.cpp SourceIteration.h InputDeck.h Utilities.h
	g++ -g -c SourceIteration.cpp -o SourceIteration.o 

Utilities.o: Utilities.cpp Utilities.h
	g++ -g -c Utilities.cpp -o Utilities.o 

clean:
	rm -rf main.o InputDeck.o SourceIteration.o Utilities.o InputDeck.h.gch SourceIteration.h.gch Utilities.h.gch

cleanall:
	rm -rf main.o InputDeck.o SourceIteration.o Utilities.o InputDeck.h.gch SourceIteration.h.gch Utilities.h.gch a.out
