#ifndef _HOME_h
#define _HOME_h

#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <cstring>

#include <string>
#include <vector>

#define real8 double
using namespace std;

/*
 *      Define a structure containing all items corresponding to
 *      all command line options that have associated values.  This
 *      gives us an easy way to pass the command line arg values around
 *      to various functions.
 */
typedef struct {
        double  burgMag;
        int     seed;
        int     type;

        std::vector<string>      inpFiles;
        std::vector<string>      outFiles;
}InArgs_t;

/*
 *      Define a structure to hold a nodal configuration type, name, and
 *      a pointer to the function to invoke to create that type of nodal
 *      configuration.
 */
typedef struct {
        int   funcType;
        const char  *funcName;
/*
        void  (* func)();
*/
} FuncData_t;

/*
 *      Define an integer identifier to be associated with each
 *      posible command line argument.  To be used to index the
 *      option-specific data in the optList array below.
 */
typedef	enum{
    OPT_HELP,
    OPT_INPFILE,
    OPT_OUTFILE,
    OPT_SEED,
    OPT_TYPE,
	OPT_MAX
}OPT_t;

typedef	enum{
	FTYPE_AVERAGE_LINES = 0,
    FTYPE_PROC_EXTEND_DIS,
	FTYPE_MAX
}FTYPE_t;

/*
 *      Define a structure to hold a command line option's id (type),
 *      name, the shortest possible unique abbreviation of the option
 *      name, and the number of parameters paired with the opt.
 */
typedef struct {
        int         optType;
        const char  *optName;
        int         optMinAbbrev;
        int         optPaired;
} Option_t;

typedef struct {
    double  x,y,z;    
} Point_t;

typedef struct {
    double  x,y,z;  
} Line_t;

typedef struct {
    vector<string>          variables;
    vector<vector<double> > data;
} Table_t;





#endif 
