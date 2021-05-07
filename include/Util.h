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

#include "MD.h"

#define DotProduct(vec1,vec2)  ((vec1[0])*(vec2[0])+(vec1[1])*(vec2[1]) +(vec1[2])*(vec2[2]))

#define  NUMBER_DATA    0x0001
#define  CHAR_DATA      0x0002

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
void FoldBox(int pbc, real8 boundMin[3], real8 boundMax[3], real8 *x, real8 *y, real8 *z);
void ZImage(int pbc, real8 boundMin[3], real8 boundMax[3], real8 *x, real8 *y, real8 *z);
real8 LinearInterpolation(const Curve_t &curve, real8 x, real8 min = -1, real8 max = -1);
void SwapTable(Table_t &table);
void SwapLineList(LineList_t &list);
void CleanDump(Dump_t &dum);
void NormalizeVec(real8 vec[3]);
void StitchTecplotData(vector<Table_t> &tables, Table_t &table, int eigenID = 0);
int DataType(const string &str);
void SpecifyEquations_PLTDATA(InArgs_t *inArgs);
void HandleTecplotData(InArgs_t *inArgs);
void SpecifyEquations(Table_t &table);
bool Analysis(int pointID, double sf, Table_t &table, real8 &sigma, real8 &hard, real8 &thard ,real8 &twindef, real8 &crss);
bool Analysis(Table_t &table, real8 &crss);
void FormatVector(real8 vec[3], const char *msg);

void InitList(LineList_t &list);
void InitTable(Table_t &table);
void SortTable(Table_t &table, int sortColID);

void AnimateAuxData(InArgs_t *inArgs);
void AnimateCurve(InArgs_t *inArgs);
#endif
