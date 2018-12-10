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
        {OPT_HELP,      "help",     	1, 0},
        {OPT_INPFILE,   "inputfile",  	1, 1},
        {OPT_OUTFILE,   "outfile",  	1, 1},
        {OPT_SEED,      "seed",     	3, 1},
        {OPT_TYPE,      "type",     	1, 1}
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
        printf("\n");
        printf("Usage:\n");
        printf("\n");
        printf("\n");

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

    printf("\n");
 
    exit(0);
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
        inArgs->outFiles.push_back((char*)"output.plt");
        inArgs->seed       = time(0) + getpid();
        inArgs->type       = 0;
        	
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
        printf("\n");
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

        string              argValue;
        std::vector<string> argValues;

        for (i = 1; i < argc; i++) {
/*
 *          If the option doesn't begin with a '-' something
 *          is wrong, so notify the user and terminate.
 */
            if (argv[i][0] != '-') {
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
                    while(i < argc){
                        argValue = argv[++i];
                        argValues.push_back(argValue);
                        if(argv[i][0] == '-')break;
                    }
                }
            }

/*
 *          Do any option-specific processing...
 */
            switch (j)  {
                case OPT_HELP:
                    PrintHelp(argv[0]);
                    exit(0);
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
            }

        }
        return;
}



main(int argc, char *argv[])
{
        InArgs_t        inArgs;


        InitDefaultValues(&inArgs);
        GetInArgs(argc, argv, &inArgs);
        CheckArgSanity(&inArgs);
        PrintArgs(&inArgs); 

/*
 *      All that's left is to invoke the proper function to
 *      handle data 
 */
        switch (inArgs.type) {
            case FTYPE_AVERAGE_LINES:

                break;
            case FTYPE_PROC_EXTEND_DIS:

                break;
        }

        exit(0);
}


