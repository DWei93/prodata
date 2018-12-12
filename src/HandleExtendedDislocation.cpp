#include <algorithm>
#include <functional>

#include "Home.h"
#include "Util.h"
#include "ReadData.h"

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
    int             index,  i;
    int             colX, colY, colBurgID;
    string          cubelName = "cubel", burgIDName = "burgID", remeshSizeName = "rsize";  

    Point_t         p;
    vector<Point_t> points1, points2;
    Curve_t         curve1, curve2;
    
    vector<real8>           seq;
    vector<vector<double> > data(3);  

    if((index = GetValID(inArgs->priVars, cubelName)) < inArgs->priVars.size()){
        cubel = stof(inArgs->priVars[index].vals[0]);
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
            burgID1 = stoi(inArgs->priVars[index].vals[0]);
            burgID2 = stoi(inArgs->priVars[index].vals[1]);
        }
    }
    printf("The burgID of two partials are %2.0f %2.0f\n", burgID1, burgID2);

    if((index = GetValID(inArgs->priVars, remeshSizeName)) < inArgs->priVars.size()){
        remeshSize = stof(inArgs->priVars[index].vals[0]);
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
        }else if(burgID2 == table.data[i][colBurgID]){
            points2.push_back(p);
        }else{
            continue;
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
    
    for(i=0; i<seq.size(); i++){
        data[1].push_back(LinearInterpolation(curve1, seq[i], boundMin[0], boundMax[0]));
        data[2].push_back(LinearInterpolation(curve2, seq[i], boundMin[0], boundMax[0]));
    }

    printf("variables = x, y1, y2\n");
    for(i=0; i<seq.size(); i++){
        printf("%f %f %f \n", data[0][i], data[1][i], data[2][i]);
    }

    return;
}





















