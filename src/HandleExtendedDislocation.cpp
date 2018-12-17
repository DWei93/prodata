#include <algorithm>
#include <functional>

#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"
#include "DDD.h"
#include "MD.h"

using namespace std;


bool cmp(Point_t p, Point_t q)
{
    return(p.x < q.x);
}

void HandleExtendedDislocation_DDD(InArgs_t *inArgs)
{
    Table_t         table;
    real8           cubel = 20000, boundMin[3], boundMax[3];
    real8           remeshSize = 20.0, burgID1 = 7, burgID2 = 12;
    real8           rd, yi1, yj1, yk2, yk3, s;
    int             index,  i, j, k, file, iPos;
    int             colX, colY, colBurgID;
    string          cubelName = "cubel", burgIDName = "burgID", remeshSizeName = "rsize";  
    string          fileName, secLine, fdir, curDir("./"), fname, aveFile;
    real8           separation = 0.0, position = 0.0;
    Point_t         p;
    vector<Point_t> points1, points2;
    Curve_t         curve1, curve2;
    LineList_t      list;
    
    vector<string>          words;
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

    if(inArgs->outFiles.size() < 2){
        aveFile = "Ave.plt";
    }else{
        aveFile = inArgs->outFiles[1];
    }
    ofstream out;
    out.open(aveFile.c_str(), ios::out);
    out << "variables = timeNow, separation, position" << endl;
    
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

    for(file=0; file<inArgs->inpFiles.size(); file++){

        vector<string>().swap(table.variables);

        ReadTecplotNormalData(inArgs->inpFiles[file], table, secLine);
        vector<string>().swap(words);
        words = split(secLine, "\"");
 
        if(table.data.size() < 2)Fatal("the data size is %d", table.data.size());
 
        colX = GetColIDFromTable(table, "X");
        colY = GetColIDFromTable(table, "Y");
        colBurgID = GetColIDFromTable(table, "burgID");
            
        vector<Point_t>().swap(points1);
        vector<Point_t>().swap(points2);
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
 
        vector<double>().swap(seq);
        seq = GenerateSequence(boundMin[0], boundMax[0], remeshSize);
 
        data.resize(7);
        for(i=0; i<data.size(); i++) data[i].resize(seq.size());
        for(i=0; i<seq.size(); i++){
            data[0][i] = seq[i];
        }
        
        separation = 0.0;
        position = 0.0;
        for(i=0; i<seq.size(); i++){
            data[1][i] = LinearInterpolation(curve1, seq[i], boundMin[0], boundMax[0]);
            data[2][i] = LinearInterpolation(curve2, seq[i], boundMin[0], boundMax[0]);
            data[3][i] = 0.5 * (data[1][i] + data[2][i]);
            separation += fabs(data[1][i] - data[2][i]);  
            position += data[3][i];
        }
        separation /= ((double)seq.size());
        position /= ((double)seq.size());
 
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

        iPos = inArgs->inpFiles[file].find_last_of('/');
        if(iPos == string::npos){
            fdir = curDir;
            fname = inArgs->inpFiles[file];
        }else{
            fdir = curDir + inArgs->inpFiles[file].substr(0,iPos);
            fname = inArgs->inpFiles[file].substr(iPos+1);
        }

        fileName += "-";
        fileName += fname;

        WriteTecplotNormalData(list, fileName, 10);
        vector<vector<double> >().swap(list.data);

        out << words[1] << " "; 
        out <<  setprecision(10)  << separation << " ";
        out <<  setprecision(10)  << position << endl;
    }

    out.close();
    return;
}


bool com(Atom_t p, Atom_t q){
    return(p.y < q.y);
}
void HandleExtendedDislocation_MD(InArgs_t *inArgs)
{
    int         i, j, file, index;
    double      cutofflen = 2.556;
    MgData_t    mg;
    LineList_t  list;
    string      cutoffName("cutoff"), dvarName("dvar"), dvar("c_vcna");

    if((index = GetValID(inArgs->priVars, cutoffName)) < inArgs->priVars.size()){
        cutofflen = atof(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The cut-off (cutoff) length is %f\n", cutofflen);

    if((index = GetValID(inArgs->priVars, dvarName)) < inArgs->priVars.size()){
        dvar = inArgs->priVars[index].vals[0];
    }
    printf("The determined var (dvar) is %s\n", dvar.c_str());

    ReadDataFromMDLogFile(inArgs->auxFiles, list);

    for(file=0; file<inArgs->inpFiles.size(); file++){
        ReadMGDataFile(inArgs->inpFiles[file], mg);
        sort(mg.atom.begin(), mg.atom.end(), com);
    
#if 0
        printf("Timestep %d, Atoms %d, Bouds %s %s %s, ", 
               mg.timestep, mg.atoms, mg.bounds[0].c_str(), 
               mg.bounds[1].c_str(), mg.bounds[2].c_str());
        printf("Box %e %e %e %e %e %e\n", mg.box[0], mg.box[1],
               mg.box[2], mg.box[3], mg.box[4], mg.box[5]);

        printf("variables are ");
        for(i=0; i<mg.variables.size(); i++){
            printf("%s(%d) ", mg.variables[i].c_str(), (int)mg.atom.size());
        }
        printf("\n");
#endif

    }
    return;
}


















