/****************************************************************************
 *
 *      Module:         Main.c 
 *      Description:    Contains the main routine for controlling
 *
 *      Included functions:
 *              main()
 *              CheckArgSanity()
 *              GetInArgs()
 *              InitDefaultValues()
 *              PrintHelp()
 *              Usage()
 *
 *      
 ***************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <cstring>

#include "Home.h"
#include "DDD.h"
#include "MD.h"
#include "Math.h"
#include "Parse.h"
#include "Util.h"

using namespace std;

/*
 *      Define and initialize an array of structures containing
 *      all the possible command line arguments, and some info
 *      determining how it is treated.
 *
 *      option    option  #characters     1 unless option
 *      type      name    in unique       has no associated
 *                        abbreviation    value
 */

Option_t optList[OPT_MAX] = {
        {OPT_HELP,          "help",     	1, 0},
        {OPT_INPFILE,       "inputfile",  	1, 1},
        {OPT_OUTFILE,       "outfile",  	1, 1},
        {OPT_SEED,          "seed",     	3, 1},
        {OPT_TYPE,          "type",     	1, 1},
        {OPT_PRIVATEVALS,   "d",            1, 1},
        {OPT_AUXFILE,       "auxfile",      3, 1},     
        {OPT_THREADS,       "nthreads",     1, 1}     
};


/*---------------------------------------------------------------------------
 *
 *      Function:       Usage
 *      Description:    Print out a brief message indicating the possible
 *                      command line options and terminated.
 *
 *      Arguments:
 *          program    name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void Usage(char *program)
{
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       PrintHelp
 *      Description:    Print out a detailed description of the available
 *                      program options, what they represent, how they
 *                      relate to each other, which are interdependent, or
 *                      mutually exclusive, default values, etc.
 *
 *      Arguments:
 *          program    name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void PrintHelp(char *program)
{
    Usage(program);

    printf("    Options may be abbreviated to the shortest non-ambiguous\n");
    printf("\n");

}

/*---------------------------------------------------------------------------
 *
 *      Function:      InitDefaultValues
 *      Description:   Set default values for all items that may be
 *                     specified via command line arguments.
 *
 *-------------------------------------------------------------------------*/
static void InitDefaultValues(InArgs_t *inArgs)
{
        inArgs->outFiles.push_back((char*)"output");
        inArgs->seed       = time(0) + getpid();
        inArgs->type       = 0;
        inArgs->help       = 0;
#ifdef _OPENMP
        inArgs->nThreads   = 4;
#else
        inArgs->nThreads   = 1;
#endif        	
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       CheckArgSanity
 *      Description:    Check the current values for the input arguments
 *                      (possibly a combination of default and user supplied
 *                      values), and verify that the values are reasonable.
 *
 *                      If there are any problems or inconsistencies, just
 *                      exit with an error message.
 *
 *-------------------------------------------------------------------------*/
static void CheckArgSanity(InArgs_t *inArgs)
{

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       PrintArgs
 *      Description:    Dump out the arguments (command line or default
 *                      values) that are being used to generate the problem.
 *                      Only those arguments appropriate to the type of
 *                      dislocation configuration being created will be
 *                      displayed.
 *
 *-------------------------------------------------------------------------*/
static void PrintArgs(InArgs_t *inArgs)
{
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       GetInArgs
 *      Description:    Parse the command line arguments.  Verify all options
 *                      are valid keywords and that all options requiring
 *                      associated values have them, and some simple sanity
 *                      checks.  All user supplied values will be stored in
 *                      a the inArgs structure, and a more in-depth sanity
 *                      check on those values is done elsewhere.
 *      Arguments
 *              argc    count of command line args
 *              argv    command line args
 *              inArgs  structure to hold user supplied values.  This
 *                      structure should have been populated with any
 *                      appropriate default values before this function
 *                      is called.
 *
 *-------------------------------------------------------------------------*/
static void GetInArgs(int argc, char *argv[], InArgs_t *inArgs)
{
        int     i, j, k;
        char    *argName;
        char    *token;
        Var_t   var;

        string                  argValue, arg;
        std::vector<string>     argValues,bakInps;
        vector<vector<string> > strs;

        for (i = 1; i < argc; i++) {
/*
 *          If the option doesn't begin with a '-' something
 *          is wrong, so notify the user and terminate.
 */
            if (argv[i][0] != '-' && atof(argv[i]) == 0) {
                Usage(argv[0]);
                exit(1);
            }

            argName = &argv[i][1];

/*
 *          Scan the array of valid options for the user supplied
 *          option name.  (This may be any unique abbreviation of
 *          one of any of the options.  If we don't find the option
 *          notify the user and terminate.
 */
            for (j = 0; j < OPT_MAX; j++) {
                if (!strncmp(argName, optList[j].optName,
                    optList[j].optMinAbbrev)) {
                    break;
                }
            }

            if (j == OPT_MAX) {
                Usage(argv[0]);
                exit(1);
            }

/*
 *          Verify that there is an associated value if the specified
 *          option is supposed to be paired with values.
 */
            if (optList[j].optPaired) {
                if (i+1 >= argc) {
                    Usage(argv[0]);
                    exit(1);
                } else {
                    i++;
                    vector<string>().swap(argValues);
                    while(i < argc && (argv[i][0] != '-' || 
                          argv[i][0] == '-' && atof(argv[i]) != 0)){
                        argValue = argv[i++];
                        argValues.push_back(argValue);
                        if(i >= argc)break;
                        if(argv[i][0] == '-' && atof(argv[i]) == 0)break;
                    }
                    i--;
                }
            }
/*
 *          Do any option-specific processing...
 */
            switch (j)  {
                case OPT_HELP:
                    PrintHelp(argv[0]);
                    inArgs->help = 1;
                    break;
                case OPT_INPFILE:
                    swap(inArgs->inpFiles, argValues);
                    break;
                case OPT_OUTFILE:
                    swap(inArgs->outFiles, argValues);
                    break;
                case OPT_SEED:
                    inArgs->seed = atoi(argValue.c_str());
                    break;
                case OPT_TYPE:
                    inArgs->type = atoi(argValue.c_str());
                    break;
                case OPT_PRIVATEVALS:
                    var.name = argName;
                    var.name.erase(0,1);
                    swap(var.vals, argValues);
                    inArgs->priVars.push_back(var);
                    break; 
                case OPT_AUXFILE:
                    swap(inArgs->auxFiles, argValues);
                    break;
                case OPT_THREADS:
                    inArgs->nThreads = atoi(argValue.c_str());
                    break;
                default:
                    break;
            }

        }
        
#if 0
        if(inArgs->priVars.size()>0){
            printf("Private vars for type %d:\n", inArgs->type);
            for(i=0; i<inArgs->priVars.size(); i++){
                printf("%s: ", inArgs->priVars[i].name.c_str());
                for(j=0; j<inArgs->priVars[i].vals.size(); j++){
                    printf("%s ", inArgs->priVars[i].vals[j].c_str());
                }
                printf("\n");
            }
        }
#endif
    
        strs.resize(inArgs->inpFiles.size());
        for(i=0; i<inArgs->inpFiles.size(); i++){
            strs[i] = GetFiles(inArgs->inpFiles[i]);
            for(j=0; j<strs[i].size(); j++){
                bakInps.push_back(strs[i][j]);
            }
        }
        swap(bakInps, inArgs->inpFiles);
        
        sort(inArgs->inpFiles.begin(), inArgs->inpFiles.end());
        if(inArgs->auxFiles.size() > 0){
            strs.resize(inArgs->auxFiles.size());
            vector<string>().swap(bakInps);
            for(i=0; i<inArgs->auxFiles.size(); i++){
                vector<string>().swap(strs[i]);
                strs[i] = GetFiles(inArgs->auxFiles[i]);
                for(j=0; j<strs[i].size(); j++){
                    bakInps.push_back(strs[i][j]);
                }
            }
            swap(bakInps, inArgs->auxFiles);
        }
        sort(inArgs->auxFiles.begin(), inArgs->auxFiles.end());
#if 1
        if(inArgs->inpFiles.size()<5){
            printf("Input files are:\n");
            for(i=0; i<inArgs->inpFiles.size(); i++){
                printf("%s \n", inArgs->inpFiles[i].c_str());
            }
        }
    
        if(inArgs->auxFiles.size() < 5){
            printf("Auxiliary file(s) :\n");
            for(i=0; i<inArgs->auxFiles.size(); i++){
                printf("%s \n", inArgs->auxFiles[i].c_str());
            }
        }
#endif
        return;
}



int main(int argc, char *argv[])
{
        int             i, j;
        InArgs_t        inArgs;

        InitDefaultValues(&inArgs);
        GetInArgs(argc, argv, &inArgs);
        CheckArgSanity(&inArgs);
        PrintArgs(&inArgs); 

#ifdef _OPENMP
        omp_set_num_threads(inArgs.nThreads); 

        printf("Threads count: %d \n", omp_get_max_threads());
#endif
/*
 *      All that's left is to invoke the proper function to
 *      handle data 
 */
#if 0
        
        real8 p1[3] = {0,0,0};
        real8 p2[3] = {0,0,0};
        real8 nDir[3] = {0,0,1.0};
        real8 point[3] = {1E10,1E80,0};
        real8 t;
        SegmentPlaneIntersection(p1, p2, nDir, point, t);

        int test = 0;
        test<<1;
        printf("test<<1: %d\n", test);

        test =1;
        test<<=1;
        printf("test<<=1: %d\n", test);

        i=10; test = (i<<1)+1;
        printf("(i<<1)+1 = %d\n", test);
            
        test = 0;
        i=3;
        test = 1+i<<1;
        i=3; printf("1+3<<1: %d\n", test);

        Fatal("t is %e",t);
#else
        
        switch (inArgs.type) {
            case FTYPE_AVERAGE_LINES:
                AverageLines(&inArgs);
                break;
            case FTYPE_PROC_EXTEND_DIS_DDD:
                HandleExtendedDislocation_DDD(&inArgs);
                break;
            case FTYPE_PROC_EXTEND_DIS_MD:
                HandleExtendedDislocation_MD(&inArgs);
                break;
            case FTYPE_SPECIFY_EQUATIONS:
                SpecifyEquations_PLTDATA(&inArgs);
                break;
            case FTYPE_GENERATE_DISLOCATION:
                GenerateDislocation(&inArgs);
                break;
            case FTYPE_DATA_HANDLING:
                HandleTecplotData(&inArgs);
                break;
            case FTYPE_ANIMATE_AUXDATA:
                AnimateAuxData(&inArgs);
            case FTYPE_ANIMATE_CURVE:
                AnimateCurve(&inArgs);
                break;
        }
#endif
        exit(0);
}


