/*
 * Utils.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef UTILS_H_
#define UTILS_H_
#include "Macros.h"

class Utils {
public:
	static void log(const char* args, ...);
	static void exit(const char* args, ...);
	static void cmpExit(bool v, const char* args, ...);
	static double getSysTime();
	static uint64_t getTimeStamp();
	static bool exists(const char* fileName);
	static string concaten(string& str, int index);
};

#endif /* UTILS_H_ */
