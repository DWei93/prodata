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
#include <list>
#include <map>

#ifdef _OPENMP

#include <omp.h>

#define INIT_LOCK(a)    omp_init_lock((a))
#define LOCK(a)         omp_set_lock((a))
#define UNLOCK(a)       omp_unset_lock((a))
#define DESTROY_LOCK(a) omp_destroy_lock((a))

#else /* _OPENMP not defined */

#define INIT_LOCK(a)
#define LOCK(a)
#define UNLOCK(a)
#define DESTROY_LOCK(a)

#endif /* end ifdef _OPENMP */


#define real8 double
using namespace std;

typedef struct {
    string          name;
    vector<string>  vals;
}Var_t;

/*
 *      Define a structure containing all items corresponding to
 *      all command line options that have associated values.  This
 *      gives us an easy way to pass the command line arg values around
 *      to various functions.
 */
typedef struct {
        double  burgMag;
        int     seed;
        int     type, nThreads;
        bool    help;

        vector<string>  inpFiles, outFiles;
        vector<string>  auxFiles;
    
        vector<Var_t>     priVars;
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
    OPT_PRIVATEVALS,
    OPT_AUXFILE,
    OPT_THREADS,
	OPT_MAX
}OPT_t;

typedef	enum{
	FTYPE_AVERAGE_LINES = 0,
    FTYPE_PROC_EXTEND_DIS_DDD,
    FTYPE_PROC_EXTEND_DIS_MD,
    FTYPE_SPECIFY_EQUATIONS,
    FTYPE_GENERATE_DISLOCATION,
	FTYPE_MAX
}FTYPE_t;

typedef enum{
    INT_DATA = 0,
    DOUBLE_DATA,
    STRING_DATA
}DATA_TYPE;
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
        int     type;
        string  name;
        string  val;
}Variable_t;

typedef struct {
    vector<string>          variables;
    vector<Variable_t>      auxData;
    map<string, string>     aux;
    vector<vector<double> > data;
} Table_t;

typedef struct {
    vector<string>          variables;
    vector<vector<double> > data;
}LineList_t;


#endif 
