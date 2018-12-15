#ifndef _PRODATAIO_h
#define _PRODATAID_h

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

#include <iomanip>


int ReadTecplotNormalData(string &file, Table_t &table, string &secLine);
void    WriteTecplotNormalData(const LineList_t &list, const string &file,
                               double precision = 6);

#endif
