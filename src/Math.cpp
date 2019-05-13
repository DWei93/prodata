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
    real8   rsize = 0, min, max, effNums = 0, value;
    string  rsizeName("rsize"), varsName("vars"), overName("over"), specifyEquName("spe");
    string  str("Ave_"), secLine, tauName("tau"), overVar, weighName("weigh"), weightCoeff;

    LineList_t  list;
    Curve_t     curve;

    vector<int>         varID;
    vector<real8>       seq;
    vector<Table_t>     tables(inArgs->inpFiles.size());
    vector<string>      s1, s2;
    vector<vector<vector<real8> > > array;

    bool weigh = false;
    if((index = GetValID(inArgs->priVars, weighName)) < inArgs->priVars.size()){
        weigh = true;
        weightCoeff = inArgs->priVars[index].vals[0];
    }
    if(weigh){
        printf("Weight coffeicient (weigh) %s is defined\n", weightCoeff.c_str());
    }else{
        printf("No weight coffeicient (weigh) is defied\n");
    }
    
    if((index = GetValID(inArgs->priVars, rsizeName)) < inArgs->priVars.size()){
        rsize = atof(inArgs->priVars[index].vals[0].c_str());
    }
    if(rsize ==  0){
        printf("Not change the original remesh size (rsize)\n", rsize);
    }else{
        printf("The remesh size (rsize) is %f\n", rsize);
    }
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

    map<string, string>::iterator   iter;
    vector<real8>   weightList(inArgs->inpFiles.size());
    real8   totalWeight = 0.0;       
    int     firstReadFile = -1;        
    for(i=0; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState){
            weightList[i] = 0.0;
            continue;
        }else{
            if(firstReadFile < 0){
                firstReadFile = i;
            }
        }

        if(weigh){
//            for(iter=tables[i].aux.begin(); iter != tables[i].aux.end(); iter++){
//                printf("%s = %s\n", iter->first.c_str(), iter->second.c_str());
//            }
//            Fatal("T");
    
            if((iter = tables[i].aux.find(weightCoeff)) != tables[i].aux.end()){
                weightList[i] = atof(iter->second.c_str());
            }else{
                printf("waring: no weight variale %s in file %s\n", weightCoeff.c_str(), inArgs->inpFiles[i].c_str());
                weightList[i] = 1.0;
            }
        }else{
            weightList[i] = 1.0;
        }
        totalWeight += weightList[i];

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

        printf("File %s(%d): Range is [%f,%f], %d points\n", 
                inArgs->inpFiles[i].c_str(), i, tables[i].data[0][colX], tables[i].data[tables[i].data.size()-1][colX],
                (int)(tables[i].data.size()));
    }
    printf("The effective range of %s is [%f,%f]\n", overVar.c_str(), min, max);
    printf("%d files have been read, the total weight is %f\n", (int)effNums, totalWeight);

    if(rsize == 0){
        seq.resize(tables[firstReadFile].data.size());
        for(i=0; i<seq.size(); i++){
            seq[i] = tables[firstReadFile].data[i][varID[colX]];
        }
    }else{
        seq = GenerateSequence(min, max, rsize);
    }
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
//                    printf("%d %e %e\n", k, curve.ay[k], tables[i].data[k][varID[j]]);
                }

                if(rsize == 0){
                    if(seq.size() != curve.ay.size()){
                        Fatal("the %d table has a differnt size (%d) from the basic one (%d)",
                               i, (int)seq.size(), (int)curve.ay.size());
                    }
                    for(k=0; k<seq.size(); k++){
                        list.data[j][k] += (curve.ay[k]*weightList[i]/totalWeight);
                    }
//                    printf("%e %e %e %e %d\n", list.data[j][k],curve.ay[k],weightList[i],totalWeight, (int)curve.ay.size()); 
                    if(calTau)array[i][j-1][k] = curve.ay[k];

                }else{
                    for(k=0; k<seq.size(); k++){
                        value = LinearInterpolation(curve, seq[k], min, max);
                        list.data[j][k] += (value*weightList[i]/totalWeight);
                        if(calTau){
                            array[i][j-1][k] = value; 
                        }
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
                    list.data[j+oriVals][k] += (pow((array[i][j][k] - list.data[j+1][k]),2)*weightList[i]/totalWeight);
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


int PlanePlaneIntersection(double *n1, double *p1, double *n2, double*p2,
                           double *vec1, double *vec2)
{
    double  n1x = n1[0], n1y = n1[1], n1z = n1[2];
    double  n2x = n2[0], n2y = n2[1], n2z = n2[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  dir[3], dir1[3], dir2[3], dot, dis, vec3[3], vec4[4], vec5[3], vec6[3];
    int     type;

    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);

    if(sqrt(SQUARE(fabs(-(n1y*n2x) + n1x*n2y)) + 
       SQUARE(fabs(n1z*n2x - n1x*n2z)) + 
       SQUARE(fabs(-(n1z*n2y) + n1y*n2z))) < EPS1){
        if(fabs(n1x*(p1x - p2x) + n1y*(p1y - p2y) + n1z*(p1z - p2z)) < EPS2){
            VECTOR_COPY(vec1,n2);
            VECTOR_COPY(vec2,p2);
            return(PLANE_INTERSECTION);
        }

        VECTOR_COPY(vec1,n2);
        VECTOR_COPY(vec2,p2);
        return(NO_INTERSECTION);
    } 

    vec1[0] = -(n1z*n2y) + n1y*n2z;
    vec1[1] = n1z*n2x - n1x*n2z;
    vec1[2] = -(n1y*n2x) + n1x*n2y;
    NormalizeVec(vec1);

    if(sqrt(SQUARE(fabs(p1x - p2x)) + SQUARE(fabs(p1y - p2y)) + 
       SQUARE(fabs(p1z - p2z))) < EPS2){
        VECTOR_COPY(vec2,p2);
        return(LINE_INTERSECTION);
    }

    cross(vec1, n1, vec3);
    if(Normal(vec3)< 1.0E-4){
        Fatal("PlanePlaneIntersection in %s at %d\n", __FILE__, __LINE__);
    }
    NormalizeVec(vec3);

    type = LinePlaneIntersection(vec3, p1, n2, p2, vec4, vec2);
    if(type == NO_INTERSECTION){
        Fatal("PlanePlaneIntersection in %s at %d\n", __FILE__, __LINE__);
    }
#if 0
    dir[0] = p2[0] - p1[0];
    dir[1] = p2[1] - p1[1];
    dir[2] = p2[2] - p1[2];
    if(Normal(dir) < 1.0E-2){
        VECTOR_COPY(vec2, p2);
        return(LINE_INTERSECTION);
    }
    NormalizeVec(dir);

    dot = DotProduct(dir, n1);
    dir1[0] = dir[0] - dot*n1[0];
    dir1[1] = dir[1] - dot*n1[1];
    dir1[2] = dir[2] - dot*n1[2];
    if(Normal(dir1) < 1.0E-4){
        VECTOR_COPY(vec2, p1);
        return(LINE_INTERSECTION);
    }

    NormalizeVec(dir1);

    dot = DotProduct(dir, n2);
    dir2[0] = dir[0] - dot*n2[0];
    dir2[1] = dir[1] - dot*n2[1];
    dir2[2] = dir[2] - dot*n2[2];
    if(Normal(dir2) < 1.0E-4){
        VECTOR_COPY(vec2, p2);
        return(LINE_INTERSECTION);
    }

    NormalizeVec(dir2);

    type = LineLineIntersection(dir1, p1, dir2, p2, vec5, vec2, &dis);

    if(fabs(dis) > 2){
        FormatVector(vec3, "p1");
        FormatVector(vec4, "p2");
        printf("Warning: PlanePlaneIntersection dis is %f\n", dis);
    }
#endif
    return(LINE_INTERSECTION);
} 

int LinePlaneIntersection(double *d, double *p1, double *n, double *p2, 
                          double *vec1, double *vec2)
{
    double  dx = d[0], dy = d[1], dz = d[2];
    double  nx = n[0], ny = n[1], nz = n[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  dis;

    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);

    if(fabs(dx*nx + dy*ny + dz*nz) < EPS1){
        VECTOR_COPY(vec1, d);
        PointPlaneIntersection(p1, n, p2, vec2, &dis);
        if(fabs(dis) < EPS2){
            VECTOR_COPY(vec2, p1);
            return(LINE_INTERSECTION);
        }else{ 
            return(NO_INTERSECTION);
        }
    }

    double  t = (-(nx*p1x) - ny*p1y - nz*p1z + nx*p2x + ny*p2y + nz*p2z)/
                (dx*nx + dy*ny + dz*nz);

    vec2[0] = p1x + t*dx;    
    vec2[1] = p1y + t*dy;    
    vec2[2] = p1z + t*dz;

    if(isnan(vec2[0]) || isinf(vec2[0]) ||
       isnan(vec2[1]) || isinf(vec2[1]) ||
       isnan(vec2[2]) || isinf(vec2[2])){
        printf(" LinePlaneIntersection: %f,%f,%f\n", 
                vec2[0], vec2[1], vec2[2]);
    }

    if(isnan(vec2[0]))vec2[0] = 0.0;
    if(isnan(vec2[1]))vec2[1] = 0.0;
    if(isnan(vec2[2]))vec2[2] = 0.0;
    
    if(isinf(vec2[0]) || isinf(vec2[1]) || isinf(vec2[2]) ){
        printf("d: {%f,%f,%f}\n", d[0], d[1], d[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("n: {%f,%f,%f}\n", n[0], n[1], n[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        exit(0);
    }
#if 0
    if(sqrt(SQUARE(vec2[2]) + SQUARE(vec2[2]) + SQUARE(vec2[2])) > 10000){
        printf("vec2: {%f,%f,%f}\n", vec2[0], vec2[1], vec2[2]);
        printf("d: {%f,%f,%f}\n", d[0], d[1], d[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("n: {%f,%f,%f}\n", n[0], n[1], n[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        printf("PointPlaneIntersection\n");
    }
#endif

    return(POINT_INTERSECTION);    
}

int PointLineIntersection(double *p1, double *d, double *p2,
                         double *vec, double *dis){

    double  p1x = p1[0], p1y=p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y=p2[1], p2z = p2[2];
    double  dx = d[0], dy=d[1], dz = d[2];
    double  norm, norm2, invNorm, t;

    VECTOR_ZERO(vec);
    norm = sqrt(SQUARE(dx) + SQUARE(dy) + SQUARE(dz));
    invNorm = 1.0/norm;

    if(norm < EPS1){
        Fatal("PointLineIntersection: line direction is zero\n");
        return(NO_INTERSECTION);
    }

    norm2 = sqrt(SQUARE(p1x-p2x) + SQUARE(p1y-p2y) + SQUARE(p1z-p2z));
    if(norm2 < EPS1){
        vec[0] = p2[0];
        vec[1] = p2[1];
        vec[2] = p2[2];
        *dis = 0.0;
        printf("0");
        return(POINT_INTERSECTION);
    }

    if(fabs(dx*(p1x-p2x) + dy*(p1y-p2y) + dz*(p1z-p2z)) > norm*norm2*(1-EPS1)){
        vec[0] = p1[0];
        vec[1] = p1[1];
        vec[2] = p1[2];
        *dis = 0.0;
        printf("0");
        return(POINT_INTERSECTION);
    }

    t = (dx*(p1x - p2x) + dy*(p1y - p2y) + dz*(p1z - p2z))*(invNorm*invNorm);
    vec[0] = p2x + dx*t;
    vec[1] = p2y + dy*t;
    vec[2] = p2z + dz*t;

    *dis = sqrt(SQUARE(vec[0]-p1[0]) + SQUARE(vec[1]-p1[1]) + SQUARE(vec[2]-p1[2]) );

    if(isnan(*dis) || isinf(*dis)){
        printf(" PointLineIntersection: %f {%f,%f,%f}\n", 
               *dis, vec[0], vec[1], vec[2]);
        vec[0] = p1[0];
        vec[1] = p1[1];
        vec[2] = p1[2];
    }

    if(isnan(*dis))*dis = 0.0;
    if(isinf(*dis)) *dis = 0.0;
#if 0
    if(*dis < EPS2 && sqrt(SQUARE(vec[0]-p1[0]) + SQUARE(vec[0]-p1[0]) + SQUARE(vec[0]-p1[0])) > 10000){
        printf("vec: {%f,%f,%f}\n", vec[0], vec[1], vec[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("d: {%f,%f,%f}\n", d[0], d[1], d[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        printf("PointLineIntersection");
    }
#endif
    if(*dis < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}

int LineLineIntersection(double *d1, double *p1, double *d2, double *p2,
                         double *vec1, double *vec2, double *dis){
    double  d1x = d1[0], d1y = d1[1], d1z = d1[2];
    double  d2x = d2[0], d2y = d2[1], d2z = d2[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];

    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);

    if(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z) < EPS1*EPS1){
        Fatal("LineLineIntersection: line direction is zero L1\n");
        return(NO_INTERSECTION);
    }
    if(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z) < EPS1*EPS1){
        Fatal("LineLineIntersection: line direction is zero L2\n");
        return(NO_INTERSECTION);
    }

    if(fabs(1.0-fabs(d1x*d2x+d1y*d2y+d1z*d2z)) < EPS1){
        if(PointLineIntersection(p1, d2, p2, vec2, dis) == POINT_INTERSECTION){
            VECTOR_COPY(vec1, d1);
            *dis = fabs(*dis);
            return(LINE_INTERSECTION);
        }else{
            *dis = fabs(*dis);
            VECTOR_COPY(vec1, d2);
            VECTOR_COPY(vec2, p2);
            return(NO_INTERSECTION);
        }
    }

/*
 *  return the closed point on the first line
 */
    vec1[0] = p1x - (d1x*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/
              (4*SQUARE(d1x*d2x + d1y*d2y + d1z*d2z) - 
              4*(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z))*
              (SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z)));        
    vec1[1] = p1y - (d1y*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/
              (4*SQUARE(d1x*d2x + d1y*d2y + d1z*d2z) - 
              4*(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z))*
              (SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z)));
    vec1[2] = p1z - (d1z*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/
              (4*SQUARE(d1x*d2x + d1y*d2y + d1z*d2z) - 
              4*(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z))*
              (SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))); 

/*
 *  return the closed point on the second line
 */
    vec2[0] = p2x + (d2x*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/
              (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
              2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
              SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
              SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z)));

    vec2[1] = p2y + (d2y*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/
              (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
              2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
              SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
              SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z)));

    vec2[2] = p2z + (d2z*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/
              (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
              2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
              SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
              SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z)));

    *dis = sqrt(SQUARE(d1z*(d2y*(p1x - p2x) + d2x*(-p1y + p2y)) + 
           d1y*(-(d2z*p1x) + d2x*p1z + d2z*p2x - d2x*p2z) + 
           d1x*(d2z*p1y - d2y*p1z - d2z*p2y + d2y*p2z))/
           (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
           2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
           SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
           SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z))));


    if(isnan(*dis) || isinf(*dis)){
        printf("LineLineIntersection: %f\n", *dis);
    }
    if(isnan(*dis)) *dis = 0.0;
    if(isinf(*dis)) *dis = 0.0;

    if(fabs(*dis) < EPS2 && sqrt(SQUARE(vec2[0]-p1[0]) + SQUARE(vec2[0]-p1[0]) + SQUARE(vec2[0]-p1[0])) > 10000){
        printf("vec2: {%f,%f,%f}\n", vec2[0], vec2[1], vec2[2]);
        printf("d1: {%f,%f,%f}\n", d1[0], d1[1], d1[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("d2: {%f,%f,%f}\n", d2[0], d2[1], d2[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        printf("LineLineIntersection\n");
    }

    if(*dis < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}


int PointPlaneIntersection(double *p1, double *n, double *p2, double *vec, double *dis)
{
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  nx = n[0], ny = n[1], nz = n[2];
    double  temVec[3];

    if(nx*nx + ny*ny + nz*nz < EPS1*EPS1){
        printf("PointPlaneIntersection: normal is close to zero.\n");
        return(0);
    }

    *dis = sqrt(SQUARE(nx*p1x + ny*p1y + nz*p1z - nx*p2x - ny*p2y - nz*p2z)/
           (SQUARE(nx) + SQUARE(ny) + SQUARE(nz)));

    if(isnan(*dis) || isinf(*dis)){
        printf(" PointPlaneIntersection: %f\n", *dis);
    }

    if(isnan(*dis)) *dis = 0.0;
    if(isinf(*dis)) *dis = 0.0;

    vec[0] = (SQUARE(ny)*p1x + SQUARE(nz)*p1x + SQUARE(nx)*p2x + 
             nx*ny*(-p1y + p2y) + nx*nz*(-p1z + p2z))/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));
    vec[1] = (SQUARE(nx)*p1y + SQUARE(nz)*p1y + nx*ny*(-p1x + p2x) + 
             SQUARE(ny)*p2y + ny*nz*(-p1z + p2z))/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));
    vec[2] = (SQUARE(nx)*p1z + SQUARE(ny)*p1z + nx*nz*(-p1x + p2x) + 
             ny*nz*(-p1y + p2y) + SQUARE(nz)*p2z)/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));

    temVec[0] = p1[0] - vec[0];
    temVec[1] = p1[1] - vec[1];
    temVec[2] = p1[2] - vec[2];
    NormalizeVec(temVec);
    if(DotProduct(temVec, n) < 0){
        *dis *= -1.0;
    }
    
    if(fabs(*dis) < EPS2 && sqrt(SQUARE(vec[0]-p1[0]) + SQUARE(vec[0]-p1[0]) + SQUARE(vec[0]-p1[0])) > 10000){
        printf("vec: {%f,%f,%f}\n", vec[0], vec[1], vec[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("n: {%f,%f,%f}\n", n[0], n[1], n[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        printf("PointPlaneIntersection\n");
    }
    if(fabs(*dis) < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}

int PointPointIntersection(double *p1, double *p2, double *vec, double *dis){
    vec[0] = p2[0] - p1[0];
    vec[1] = p2[1] - p1[1];
    vec[2] = p2[2] - p1[2];

    *dis = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if(isnan(*dis) || isinf(*dis)){
        printf(" PointPointIntersection: %f\n", *dis);
    }
    if(isnan(*dis)) *dis = 0.0;
    if(isinf(*dis)) *dis = 0.0;
    
    if(*dis < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}


void cross(real8 a[3], real8 b[3], real8 c[3])
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

real8 Normal(real8 a[3])
{
        return( sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) );
}
