#ifndef _MATH_h
#define _MATH_h

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

#include <stdlib.h>
#include <algorithm>
#include <numeric>

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define SQUARE(a) ((a)*(a))
#define VECTOR_ADD(a,b)  {(a)[0] += (b)[0]; (a)[1] += (b)[1]; (a)[2] += (b)[2];}
#define VECTOR_COPY(a,b) {(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2];}
#define VECTOR_ZERO(a)   {(a)[0] = 0; (a)[1]  = 0; (a)[2] = 0;}
#define VECTOR_MULTIPLY(a,x)	{(a)[0] *= (x);    (a)[1] *= (x);    (a)[2] *= (x);}

#define NO_INTERSECTION     0
#define POINT_INTERSECTION  1
#define LINE_INTERSECTION   2
#define PLANE_INTERSECTION  3
#define EPS1 1.0E-3
#define EPS2 1.0E-1

int PlanePlaneIntersection(double *n1, double *p1, double *n2, double*p2,
                           double *vec1, double *vec2);
int LinePlaneIntersection(double *d, double *p1, double *n, double *p2, 
                          double *vec1, double *vec2);
int LineLineIntersection(double *d1, double *p1, double *d2, double *p2,
                         double *vec1, double *vec2, double *dis);
int PointPlaneIntersection(double *p1, double *n, double *p2, double *vec, double *dis);
int PointLineIntersection(double *p1, double *d, double *p2, double *vec, double *dis);
int PointPointIntersection(double *p1, double *p2, double *vec, double *dis);
void cross(real8 a[3], real8 b[3], real8 c[3]);
real8 Normal(real8 a[3]);

void AverageLines(InArgs_t *inArgs);

#endif
