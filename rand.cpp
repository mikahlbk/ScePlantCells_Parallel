#include <iostream>
#include <ctime>

using namespace std;

double unifRand() {
	return rand() / double (RAND_MAX);
}

double unifRand(double a, double b) {
	return (b-a)*unifRand() + a;
}

void seed() {
	srand(time(0));
}
