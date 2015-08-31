/*
 * Utils.cpp
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#include "Utils.h"
#include <sys/time.h>
#include <stdarg.h>

void Utils::log(const char* args, ...) {
	va_list va_ptr;
	//fprintf(stderr, "[P%d]", getpid());
	va_start(va_ptr, args);
	vfprintf(stderr, args, va_ptr);
	va_end(va_ptr);
}
void Utils::exit(const char* args, ...) {
	va_list va_ptr;
	//fprintf(stderr, "[P%d]", getpid());
	//print out the message
	va_start(va_ptr, args);
	vfprintf(stderr, args, va_ptr);
	va_end(va_ptr);

	//exit the program
	::exit(-1);
}
void Utils::cmpExit(bool v, const char* args, ...) {
	va_list va_ptr;
	//check the truth of v
	if (v) {
		//print out the message
		va_start(va_ptr, args);
		vfprintf(stderr, args, va_ptr);
		va_end(va_ptr);

		//exit the program
		::exit(-1);
	}
}
double Utils::getSysTime() {
	double dtime;
	struct timeval tv;

	gettimeofday(&tv, NULL);

	dtime = (double) tv.tv_sec;
	dtime += (double) (tv.tv_usec) / 1000000.0;

	return dtime;

}
uint64_t Utils::getTimeStamp() {
	uint64_t stamp;
	struct timeval tv;

	gettimeofday(&tv, NULL);

	stamp = tv.tv_sec * 1000;
	stamp += tv.tv_usec;

	return stamp;
}
bool Utils::exists(const char* fileName) {
	if (!fileName) {
		return false;
	}
	FILE* file = fopen(fileName, "rb");
	if (file == NULL) {
		return false;
	}
	fclose(file);
	return true;
}
string Utils::concaten(string& str, int index) {
	char buffer[1024];
	sprintf(buffer, "%d", index);

	return str + buffer;
}
