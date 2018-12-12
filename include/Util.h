#ifndef _UTIL_h
#define _UTIL_h

#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <vector>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>

#include <stdlib.h>
#include <algorithm>
#include <numeric>


typedef struct {
    double  x,y,z;    
} Point_t;

typedef struct {
    vector<double> ax, ay, az;
} Curve_t;

void Fatal(const char *format, ...); 
vector<string> split(const string& str, const string& delim); 
void WashString(string &str);
int  GetValID(const vector<Var_t> &vals, const string &name);
int GetColIDFromTable(const Table_t &table, const string &name);
vector<double>  GenerateSequence(double from, double to, double meshSize);
void FoldBox(real8 boundMin[3], real8 boundMax[3], real8 *x, real8 *y, real8 *z);
real8 LinearInterpolation(const Curve_t &curve, real8 x, real8 min = -1, real8 max = -1);

#endif
