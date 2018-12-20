#include <sstream>
#include "Home.h"
#include "Util.h"

#include "ProDataIO.h"
#include "Parse.h"
#include "Math.h"

using namespace std;

void AverageLines(InArgs_t *inArgs)
{
    int     index, i, j, k, colX, colY;
    int     readState;
    bool    firstFile = 1;
    real8   rsize = 20, min, max, effNums = 0;
    string  rsizeName("rsize"), varsName("vars");
    string  str("Ave_"), secLine; 

    LineList_t  list;
    Curve_t     curve;

    vector<int>         varID;
    vector<real8>       seq;
    vector<Table_t>     tables(inArgs->inpFiles.size());
    vector<string>      s1, s2;

    if((index = GetValID(inArgs->priVars, rsizeName)) < inArgs->priVars.size()){
        rsize = atof(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The remesh size (rsize) is %f\n", rsize);

    if((index = GetValID(inArgs->priVars, varsName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() < 2){
            Fatal("variables is not enough for aveage line");
        }

        list.variables.resize(inArgs->priVars[index].vals.size());
        varID.resize(inArgs->priVars[index].vals.size());

        printf("The variables (vars) are ");
        for(i=0; i<inArgs->priVars[index].vals.size(); i++){
            list.variables[i] = inArgs->priVars[index].vals[i];
            printf("%s ", list.variables[i].c_str());
        }
        printf("\n");
    }

    if(inArgs->help)return;

    if(inArgs->inpFiles.size() == 0){
        Fatal("There is no input file.");
    }

    for(i=0; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState)continue;

        if(tables[i].data.size() < 2){
            Fatal("The size of file %s is wrong", inArgs->inpFiles[i].c_str());
        }
        
        if(firstFile){
            if(list.variables.size()==0){
                list.variables.resize(tables[i].variables.size());
                varID.resize(list.variables.size());
            
                printf("The average variables (vars) are: ");
                for(j=0; j<tables[i].variables.size(); j++){
                    list.variables[j] = tables[i].variables[j];
                    printf("%s ",list.variables[j].c_str());
                }
                printf("\n");
            }
            min = tables[i].data[0][colX];
            max = tables[i].data[tables[i].data.size()-1][colX];
            firstFile = 0;
        }else{
            if(min < tables[i].data[0][colX]) min = tables[i].data[0][colX];
            if(max > tables[i].data[tables[i].data.size()-1][colX]){
                max = tables[i].data[tables[i].data.size()-1][colX];
            }
        }
    }

    printf("The effective range of %s is [%f,%f]\n", list.variables[0].c_str(), min, max);
    seq = GenerateSequence(min, max, rsize);
    list.data.resize(list.variables.size());
    for(i=0; i<list.data.size(); i++){
        list.data[i].resize(seq.size());
        for(j=0; j<list.data[i].size(); j++)list.data[i][j] = 0.0;
    }

    for(i=0; i<seq.size(); i++){
        list.data[0][i] = seq[i];
    }

    for(i=0; i<tables.size(); i++){
        if(tables[i].data.size()==0)continue;
        effNums++;

        for(j=0; j<varID.size(); j++){
            varID[j] = GetColIDFromTable(tables[i], list.variables[j]);
            if(varID[j] == tables[i].variables.size()){
                Fatal("there is no %s in the file %s ", list.variables[j].c_str(), 
                        inArgs->inpFiles[i].c_str());
            }

            if(j==0){
                curve.ax.resize(tables[i].data.size());
                curve.ay.resize(tables[i].data.size());
                for(k=0; k<curve.ax.size(); k++){
                    curve.ax[k] =  tables[i].data[k][varID[j]];
                }
            }else{
                for(k=0; k<curve.ax.size(); k++){
                    curve.ay[k] =  tables[i].data[k][varID[j]];
                }

                for(k=0; k<seq.size(); k++){
                    list.data[j][k] += LinearInterpolation(curve, seq[k], min, max);
                }
            }

        }
    }

    for(i=1; i<list.data.size(); i++){
        for(j=0; j<list.data[i].size(); j++){
            list.data[i][j] /= effNums;
        }
    }

    WriteTecplotNormalData(list, inArgs->outFiles[0], 10);

    return;
}
