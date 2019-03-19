#include <algorithm>
#include <functional>

#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"
#include "DDD.h"
#include "MD.h"
#include "Math.h"

using namespace std;

void GenerateScrewDislocation(InArgs_t *inArgs)
{
    int         index, i;    
    MgData_t    mg;
    string      posName("pos"), lineName("line"), gDirName("gdir"), burgMagName("burgMag");
    

    real8       pos[3] = {0.0, 0.0, 0.0};
    real8       line[3] = {1.0, 0.0, 0.0};
    real8       gDir[3] = {0.0, 1.0, 0.0};
    real8       normal[3];
    real8       burgMag = 2.556, boundMin[3], boundMax[3], argument;

    if((index = GetValID(inArgs->priVars, posName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 3){
            pos[0] = atof(inArgs->priVars[index].vals[0].c_str());
            pos[1] = atof(inArgs->priVars[index].vals[1].c_str());
            pos[2] = atof(inArgs->priVars[index].vals[2].c_str());
        }
    }
    printf("The position of dislocation (pos): {%f, %f, %f}\n", pos[0], pos[1], pos[2]);

    if((index = GetValID(inArgs->priVars, lineName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 3){
            line[0] = atof(inArgs->priVars[index].vals[0].c_str());
            line[1] = atof(inArgs->priVars[index].vals[1].c_str());
            line[2] = atof(inArgs->priVars[index].vals[2].c_str());
        }
    }
    printf("The line direction of dislocation (line): {%f, %f, %f}\n", line[0], line[1], line[2]);

    if((index = GetValID(inArgs->priVars, gDirName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 3){
            gDir[0] = atof(inArgs->priVars[index].vals[0].c_str());
            gDir[1] = atof(inArgs->priVars[index].vals[1].c_str());
            gDir[2] = atof(inArgs->priVars[index].vals[2].c_str());
        }
    }
    printf("The glide direction of dislocation (gDir): {%f, %f, %f}\n", gDir[0], gDir[1], gDir[2]);

    if((index = GetValID(inArgs->priVars, burgMagName)) < inArgs->priVars.size()){
         burgMag = atof(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The magnitude of dislocation (burgMag): %f\n", burgMag);


    NormalizeVec(gDir);
    NormalizeVec(line);
    cross(gDir, line, normal);
    NormalizeVec(normal);

    if(fabs(DotProduct(line, gDir)) > 1.0E-5)Fatal("the glide direction should be prependicular to the line direction");

    if(inArgs->help)return;

    ReadMGDataFile(inArgs->inpFiles[0], mg);

    for(i=0; i<3; i++){
        boundMin[i] = mg.box[i][0];
        boundMax[i] = mg.box[i][1];
    } 

#pragma omp parallel private(i) shared(mg)
    {
        real8  point[3], dvec[3], dis, vec[3], a, b, shift;

        for(i=0; i<mg.atom.size(); i++){
            point[0] = mg.atom[i].x;
            point[1] = mg.atom[i].y;
            point[2] = mg.atom[i].z;
//            printf("%f %f %f %d\n", mg.atom[i].x, mg.atom[i].y, mg.atom[i].z, mg.atom[i].id);

            if(POINT_INTERSECTION == PointLineIntersection(point, line, pos, dvec, &dis)){
                continue;
            }

            vec[0] = dvec[0] - point[0];
            vec[1] = dvec[1] - point[1];
            vec[2] = dvec[2] - point[2];

            a = DotProduct(vec, gDir);
            b = DotProduct(vec, normal); 

            shift = (0.5*burgMag*atan2(a, b)/M_PI);

            mg.atom[i].x += (shift*line[0]);
            mg.atom[i].y += (shift*line[1]);
            mg.atom[i].z += (shift*line[2]);

            FoldBox(boundMin, boundMax, &(mg.atom[i].x), &(mg.atom[i].y), &(mg.atom[i].z));
        }
    }

//    WriteMGDataFile(inArgs->outFiles[0], mg); 
    MGToLMPDataFile(inArgs->outFiles[0], mg);
    return;
}


void GenerateDislocation(InArgs_t *inArgs){
    int     index;
    string  disTypeName("type");
    int     type = 0;
    
    if((index = GetValID(inArgs->priVars, disTypeName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 1){
            type = atoi(inArgs->priVars[index].vals[0].c_str());
        }
    }

    switch(type){
        case 0:
            GenerateScrewDislocation(inArgs);
            break;

        default:
            GenerateScrewDislocation(inArgs);
            break;
    }

    return;
}
