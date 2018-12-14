#include <algorithm>
#include <functional>

#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"
#include "DDD.h"

using namespace std;


bool cmp(Point_t p, Point_t q)
{
    return(p.x < q.x);
}

void HandleExtendedDislocation(InArgs_t *inArgs)
{
    Table_t         table;
    real8           cubel = 20000, boundMin[3], boundMax[3];
    real8           remeshSize = 20.0, burgID1 = 7, burgID2 = 12;
    real8           rd, yi1, yj1, yk2, yk3, s;
    int             index,  i, j, k;
    int             colX, colY, colBurgID;
    string          cubelName = "cubel", burgIDName = "burgID", remeshSizeName = "rsize";  
    string          fileName;
    real8           separation = 0.0;
    Point_t         p;
    vector<Point_t> points1, points2;
    Curve_t         curve1, curve2;
    LineList_t      list;
    
    vector<real8>           seq;
    vector<vector<double> > data(7);  

    list.variables.push_back("x");
    list.variables.push_back("y1");
    list.variables.push_back("y2");
    list.variables.push_back("y");
    list.variables.push_back("Rd1");
    list.variables.push_back("Rd2");
    list.variables.push_back("Rd");

    if((index = GetValID(inArgs->priVars, cubelName)) < inArgs->priVars.size()){
        cubel = atof(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The cubel size is %f\n", cubel);
    
    boundMin[0] = -0.5*cubel;
    boundMin[1] = -0.5*cubel;
    boundMin[2] = -0.5*cubel;
    
    boundMax[0] =  0.5*cubel;
    boundMax[1] =  0.5*cubel;
    boundMax[2] =  0.5*cubel;

    if((index = GetValID(inArgs->priVars, burgIDName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 2){
            burgID1 = atoi(inArgs->priVars[index].vals[0].c_str());
            burgID2 = atoi(inArgs->priVars[index].vals[1].c_str());
        }
    }
    printf("The burgID of two partials are %2.0f %2.0f\n", burgID1, burgID2);

    if((index = GetValID(inArgs->priVars, remeshSizeName)) < inArgs->priVars.size()){
        remeshSize = atof(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The remesh size is %f\n", remeshSize);

    ReadTecplotNormalData(inArgs->inpFiles[0], table);

    if(table.data.size() < 2)Fatal("the data size is %d", table.data.size());

    colX = GetColIDFromTable(table, "X");
    colY = GetColIDFromTable(table, "Y");
    colBurgID = GetColIDFromTable(table, "burgID");
        
    p.x = 1.0E10; p.y = 1.0E10;
    for(i=0; i<table.data.size(); i++){
        if(p.x == table.data[i][colX] &&
           p.y == table.data[i][colY])continue;

        p.x = table.data[i][colX];
        p.y = table.data[i][colY];

        if(burgID1 == table.data[i][colBurgID]){
            points1.push_back(p);
        }
        if(burgID2 == table.data[i][colBurgID]){
            points2.push_back(p);
        }
    }

    sort(points1.begin(), points1.end(), cmp); 
    sort(points2.begin(), points2.end(), cmp); 

    curve1.ax.resize(points1.size());
    curve1.ay.resize(points1.size());
    curve2.ax.resize(points2.size());
    curve2.ay.resize(points2.size());

    for(i=0; i<points1.size(); i++){
        curve1.ax[i] = points1[i].x;
        curve1.ay[i] = points1[i].y;
    }

    for(i=0; i<points2.size(); i++){
        curve2.ax[i] = points2[i].x;
        curve2.ay[i] = points2[i].y;
    }

    seq = GenerateSequence(boundMin[0], boundMax[0], remeshSize);

    data[0].assign(seq.begin(), seq.end());
    for(i=1; i<data.size(); i++) data[i].resize(data[0].size());
    
    for(i=0; i<seq.size(); i++){
        data[1][i] = LinearInterpolation(curve1, seq[i], boundMin[0], boundMax[0]);
        data[2][i] = LinearInterpolation(curve2, seq[i], boundMin[0], boundMax[0]);
        data[3][i] = 0.5 * (data[1][i] + data[2][i]);
        separation += fabs(data[1][i] - data[2][i]);  
    }
    separation /= ((double)seq.size());

    for(i=0; i<seq.size(); i++){
        data[4][i] = 0.0;
        data[5][i] = 0.0;
        data[6][i] = 0.0;
    }

    for(i=0; i<seq.size(); i++){
        for(j=i; j<i+seq.size(); j++){
            k = (j < seq.size()) ? j : j-seq.size();
            data[4][j-i] += fabs(data[1][k] - data[1][i]);
            data[5][j-i] += fabs(data[2][k] - data[2][i]);
            data[6][j-i] += fabs(data[3][k] - data[3][i]);
        }
    }

    for(i=0; i<seq.size(); i++){
        data[0][i] += 0.5*cubel;
        data[4][i] /= ((double)seq.size());
        data[5][i] /= ((double)seq.size());
        data[6][i] /= ((double)seq.size());
    }

    swap(list.data, data);
    fileName = inArgs->outFiles[0];
    WriteTecplotNormalData(list, fileName, 10);

    return;
}





















