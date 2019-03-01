#include <sstream>
#include "Home.h"
#include "Util.h"

#include "ProDataIO.h"
#include "Parse.h"
#include "Math.h"

using namespace std;

void AverageLines(InArgs_t *inArgs)
{
    int     index, i, j, k, colX = -1, colY;
    int     readState;
    bool    firstFile = 1, calTau = 0, specifyEqu = 0;
    real8   rsize = 20, min, max, effNums = 0, value;
    string  rsizeName("rsize"), varsName("vars"), overName("over"), specifyEquName("spe");
    string  str("Ave_"), secLine, tauName("tau"), overVar;

    LineList_t  list;
    Curve_t     curve;

    vector<int>         varID;
    vector<real8>       seq;
    vector<Table_t>     tables(inArgs->inpFiles.size());
    vector<string>      s1, s2;
    vector<vector<vector<real8> > > array;

    if((index = GetValID(inArgs->priVars, rsizeName)) < inArgs->priVars.size()){
        rsize = atof(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The remesh size (rsize) is %f\n", rsize);

    if((index = GetValID(inArgs->priVars, tauName)) < inArgs->priVars.size()){
        calTau = atoi(inArgs->priVars[index].vals[0].c_str());
    }
    if(calTau){
        printf("The variance (tau) will been caculated\n");
    }else{
        printf("The variance (tau) will not been caculated\n");
    }

    if((index = GetValID(inArgs->priVars, specifyEquName)) < inArgs->priVars.size()){
        specifyEqu = atoi(inArgs->priVars[index].vals[0].c_str());
    }
    if(specifyEqu){
        printf("Specifing equations (spe) is on.\n");
    }else{
        printf("Specifing equations (spe) is off.\n");
    }

    if((index = GetValID(inArgs->priVars, varsName)) < inArgs->priVars.size()){
        int index2;        

        if(inArgs->priVars[index].vals.size() < 2){
            Fatal("variables is not enough for aveage line");
        }

        list.variables.resize(inArgs->priVars[index].vals.size());
        varID.resize(inArgs->priVars[index].vals.size());

        printf("The variables (vars) are ");
        if((index2 = GetValID(inArgs->priVars, overName)) < inArgs->priVars.size()){
            overVar = inArgs->priVars[index2].vals[0];
        }else{
            overVar = inArgs->priVars[index].vals[0];
        }

        for(i=0; i<inArgs->priVars[index].vals.size(); i++){
            list.variables[i] = inArgs->priVars[index].vals[i];
            printf("%s ", list.variables[i].c_str());
            if(overVar == list.variables[i]){
                colX = i;
            }
        }
        printf(", averaged over through %s (over).\n", overVar.c_str());

        if(colX < 0){
            Fatal("no %s in the variables.", overVar.c_str());
        }

        if(colX != 0){
            swap(list.variables[0], list.variables[colX]);
        }
    }

    if(inArgs->help)return;

    if(inArgs->inpFiles.size() == 0){
        Fatal("There is no input file.");
    }

    printf("Ranges of files:\n ");
    for(i=0; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState)continue;

        if(specifyEqu){
            SpecifyEquations(tables[i]);
        }

        effNums++;
        if(tables[i].data.size() < 2){
            Fatal("The size of file %s is wrong", inArgs->inpFiles[i].c_str());
        }
        
        if(firstFile){
            if(list.variables.size()==0){
                list.variables.resize(tables[i].variables.size());
                varID.resize(list.variables.size());
            
                printf("The average variables (vars) are: ");

                if((index = GetValID(inArgs->priVars, overName)) < inArgs->priVars.size()){
                    overVar = inArgs->priVars[index].vals[0];   
                }else{
                    overVar = tables[i].variables[0];
                }

                for(j=0; j<tables[i].variables.size(); j++){
                    list.variables[j] = tables[i].variables[j];

                    if(overVar == list.variables[j]){
                        colX = j;
                    }
                    printf("%s ",list.variables[j].c_str());
                }
                printf(", averaged over through %s (over).\n", overVar.c_str());

                if(colX < 0){
                    Fatal("no %s in the variables.", overVar.c_str());
                }

                if(colX != 0){
                    swap(list.variables[0], list.variables[colX]);
                }
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

        printf("[%f,%f]\n", tables[i].data[0][colX], tables[i].data[tables[i].data.size()-1][colX]);
    }
    printf("The effective range of %s is [%f,%f]\n", overVar.c_str(), min, max);
    printf("%d files was read\n", (int)effNums);

    seq = GenerateSequence(min, max, rsize);
    if(calTau){
        array.resize((int)(effNums));
        for(i=0; i<array.size(); i++){
            array[i].resize(list.variables.size()-1);
            for(j=0; j<array[i].size(); j++){
                array[i][j].resize(seq.size());
            }
        }
    }

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
                    curve.ay[k] = tables[i].data[k][varID[j]];
                }

                for(k=0; k<seq.size(); k++){
                    value = LinearInterpolation(curve, seq[k], min, max);
                    list.data[j][k] += (value/effNums);
                    if(calTau){
                        array[i][j-1][k] = value; 
                    }
                }
            }
        }
    }

    if(calTau){
        int     oriVals;
        string  head("D_");
    
        oriVals = list.variables.size();
        list.variables.resize(2*list.variables.size()-1);
        list.data.resize(list.variables.size());

        for(i=oriVals; i<list.variables.size(); i++){
            list.data[i].resize(seq.size());
            list.variables[i] = head + list.variables[i-oriVals+1];
            for(j=0; j<list.data[i].size(); j++)list.data[i][j] = 0.0;
        }

        for(i=0; i<array.size(); i++){
            for(j=0; j<array[i].size(); j++){
                for(k=0; k<array[i][j].size(); k++){
                    list.data[j+oriVals][k] += (pow((array[i][j][k] - list.data[j+1][k]),2)/((real8)array.size()));
                }
            }
        }

        for(i=oriVals; i<list.data.size(); i++){
            for(j=0; j<list.data[i].size(); j++){
                list.data[i][j] = sqrt(list.data[i][j]);
            }
        }
    }

    WriteTecplotNormalData(list, inArgs->outFiles[0], 10);

    return;
}
