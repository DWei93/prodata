#ifndef _MD_h
#define _MD_h

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

typedef struct{
    int     id, type;
    double  x, y, z;

    vector<double> vars;
}Atom_t;

typedef struct {
    int     timestep;
    real8   box[6];

    vector<string> bounds;
    vector<string> variables;
    vector<Atom_t> atom;
}MgData_t;

void HandleExtendedDislocation_MD(InArgs_t *inArgs);

#endif
