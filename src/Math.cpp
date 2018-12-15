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
    int     record = 0;
    real8   rsize = 20, min, max;
    string  rsizeName("rsize"), varsName("vars"), recordName("record");
    string  xName("X"), yName("Y"), str("Ave_"), secLine; 

    LineList_t  list;
    Curve_t     curve;

    vector<real8>       seq, yVals, aveVals;
    vector<Table_t>     tables(inArgs->inpFiles.size());
    vector<string>      s1, s2;

    if((index = GetValID(inArgs->priVars, rsizeName)) < inArgs->priVars.size()){
        rsize = atof(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The remesh size (rsize) is %f\n", rsize);

    if((index = GetValID(inArgs->priVars, varsName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 2){
            xName = inArgs->priVars[index].vals[0];
            yName = inArgs->priVars[index].vals[1];
        }
    }
    str += yName; 
    printf("The variables (vars) are %s %s\n", xName.c_str(), yName.c_str());

    if((index = GetValID(inArgs->priVars, recordName)) < inArgs->priVars.size()){
        record = atoi(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The Record (record) state is %d\n", record);

    for(i=0; i<inArgs->inpFiles.size(); i++){
        ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(i == 0){
            colX = GetColIDFromTable(tables[i], xName);
            colY = GetColIDFromTable(tables[i], yName);
            min = tables[i].data[0][colX];
            max = tables[i].data[tables[i].data.size()-1][colX];
      }else{
            if(colX != GetColIDFromTable(tables[i], xName) ||
               colY != GetColIDFromTable(tables[i], yName)){
                Fatal("The format of file %s is not same as the first one",
                      inArgs->inpFiles[i].c_str());
            }

            if(min < tables[i].data[0][colX]) min = tables[i].data[0][colX];
            if(max > tables[i].data[tables[i].data.size()-1][colX]){
                max = tables[i].data[tables[i].data.size()-1][colX];
            }
        }
        if(tables[i].data.size() < 2){
            Fatal("The size of file %s is wrong", inArgs->inpFiles[i].c_str());
        }
    }

    seq = GenerateSequence(min, max, rsize);
    yVals.resize(seq.size());
    aveVals.resize(seq.size());

    for(i=0; i<seq.size(); i++){
        aveVals[i] = 0.0;
    }

    list.variables.push_back(xName);
    list.data.push_back(seq);

    for(i=0; i<inArgs->inpFiles.size(); i++){
        curve.ax.resize(tables[i].data.size());
        curve.ay.resize(tables[i].data.size());

        for(j=0; j<curve.ax.size(); j++){
            curve.ax[j] = tables[i].data[j][colX];
            curve.ay[j] = tables[i].data[j][colY];
        }
        
        for(j=0; j<seq.size(); j++){
            yVals[j] = LinearInterpolation(curve, seq[j], min, max);
            aveVals[j] += yVals[j];
        }
        list.data.push_back(yVals);
         
        if(record){
            s1 = split(inArgs->inpFiles[i], "/");
            s2 = split(s1.back(), ".");
            list.variables.push_back(s2.front());
            list.data.push_back(yVals);
        }
    }

    list.variables.push_back(str);
    for(j=0; j<seq.size(); j++){
        aveVals[j] /= ((double)inArgs->inpFiles.size());
    }    
    list.data.push_back(aveVals);

    WriteTecplotNormalData(list, inArgs->outFiles[0], 10);

    return;
}
