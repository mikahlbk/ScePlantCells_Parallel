CC=g++ -std=c++11

CFLAGS=-c -Wall

all: program

program: folder main.o coord.o node.o cell.o cell_calibration.o cell_div.o tissue.o rand.o
		$(CC) main.o coord.o node.o cell.o cell_calibration.o cell_div.o tissue.o rand.o -o program

folder: 
		mkdir -p ./DataOutput ./Animation

main.o: main.cpp
		$(CC) $(CFLAGS) main.cpp

coord.o: coord.cpp
		$(CC) $(CFLAGS) coord.cpp

node.o: node.cpp
		$(CC) $(CFLAGS) node.cpp

cell.o: cell.cpp
		$(CC) $(CFLAGS) cell.cpp

cell_div.o: cell_div.cpp
		$(CC) $(CFLAGS) cell_div.cpp

cell_calibration.o: cell_calibration.cpp
		$(CC) $(CFLAGS) cell_calibration.cpp

tissue.o: tissue.cpp
		$(CC) $(CFLAGS) tissue.cpp

rand.o: rand.cpp
		$(CC) $(CFLAGS) rand.cpp

clean: wipe
		rm -rf strain_vec.txt stress_vec.txt *o program

wipe:
		rm -rf ./Animation ./DataOutput
