# **************************************************
# Variables to control Makefile operation
CXX= -g++ -fopenmp -static-libstdc++
CXXFLAGS= -Wall -g 
LIBMAT1 := -L/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018b/bin/matlab/bin/glnxa64/
LIBMAT2 := -L/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018b/extern/bin/glnxa64/
LIBMAT3 := -I/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018b/extern/include/
LIBSFINAL := -lMatlabEngine -lMatlabDataArray
#*****************************************************
# Targets needed to bring the executable up to date
all: program 

program:folder main.o coord.o node.o cell.o cell_div.o tissue.o rand.o
	$(CXX) $(CXXFLAGS) $(LIBMAT1) $(LIBMAT2) -o program  main.o coord.o node.o cell.o cell_div.o tissue.o rand.o $(LIBSFINAL)
	#$(CXX) $(CXXFLAGS) -o program main.o coord.o node.o cell.o cell_div.o tissue.o rand.o
folder: 
		mkdir -p ./DataOutput ./Animation
main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(LIBMAT1) $(LIBMAT2) $(LIBMAT3) -c main.cpp $(LIBSFINAL)
	
	#$(CXX) $(CXXFLAGS) -c main.cpp

coord.o: coord.cpp
	$(CXX) $(CXXFLAGS) -c coord.cpp

node.o: node.cpp
	$(CXX) $(CXXFLAGS) -c node.cpp

cell.o: cell.cpp
	$(CXX) $(CXXFLAGS) -c cell.cpp
	
cell_div.o: cell_div.cpp
	$(CXX) $(CXXFLAGS) -c cell_div.cpp

tissue.o: tissue.cpp
	$(CXX) $(CXXFLAGS) -c tissue.cpp

rand.o: rand.cpp
	$(CXX) $(CXXFLAGS) -c rand.cpp

clean: wipe
		rm -rf strain_vec.txt stress_vec.txt *o program

wipe:
