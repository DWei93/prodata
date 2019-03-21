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
    Dump_t    dum;
    string      posName("pos"), lineName("line"), gDirName("gdir"), burgMagName("burgMag");
    string      screwTypeName("stype"), multiDisName("multi");
    int         screwType = 0, nPairs = 1;
    double      rsize = 0, delta;

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

    if((index = GetValID(inArgs->priVars, screwTypeName)) < inArgs->priVars.size()){
         screwType = atoi(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The type of screw dislocation (stype): %d\n", screwType);

    if((index = GetValID(inArgs->priVars,  multiDisName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 3){
            nPairs = atoi(inArgs->priVars[index].vals[0].c_str());
            rsize = atof(inArgs->priVars[index].vals[1].c_str());
            delta = atof(inArgs->priVars[index].vals[2].c_str()); 
        }
    }
    printf("Multi-pairs (multi): %d pairs with a mesh size %f\n", nPairs, rsize);

    NormalizeVec(gDir);
    NormalizeVec(line);
    cross(gDir, line, normal);
    NormalizeVec(normal);
    FormatVector(normal, "normal");

    if(fabs(DotProduct(line, gDir)) > 1.0E-5)Fatal("the glide direction should be prependicular to the line direction");

    if(inArgs->help)return;

    ReadDumpFile(inArgs->inpFiles[0], dum);

    for(i=0; i<3; i++){
        boundMin[i] = dum.box[i][0];
        boundMax[i] = dum.box[i][1];
    } 

    real8  point[3], dvec[3], dis, vec[3], a, b, shift;

    while(nPairs > 0){

#pragma omp parallel for private(point, dvec, dis, vec, a, b, shift) shared(dum)
        for(i=0; i<dum.atom.size(); i++){
            point[0] = dum.atom[i].x;
            point[1] = dum.atom[i].y;
            point[2] = dum.atom[i].z;
//            printf("%f %f %f %d\n", dum.atom[i].x, dum.atom[i].y, dum.atom[i].z, dum.atom[i].id);
        
            if(POINT_INTERSECTION == PointLineIntersection(point, line, pos, dvec, &dis)){
//                printf("%d: %f %f %f - %f %f %f  @@ %f\n", i, pos[0], pos[1], pos[2], point[0], point[1], point[2], dis);
                continue;
            }
        
            vec[0] = dvec[0] - point[0];
            vec[1] = dvec[1] - point[1];
            vec[2] = dvec[2] - point[2];
        
            a = DotProduct(vec, gDir);
            b = DotProduct(vec, normal); 
        
            if(screwType == 0){
                shift = (0.5*burgMag*atan2(b, a)/M_PI);
            }else{
                shift = (0.5*burgMag*(M_PI + atan2(b, a)/M_PI));
            }
        
            dum.atom[i].x += (shift*line[0]);
            dum.atom[i].y += (shift*line[1]);
            dum.atom[i].z += (shift*line[2]);
        
//            FoldBox(boundMin, boundMax, &(dum.atom[i].x), &(dum.atom[i].y), &(dum.atom[i].z));
        }

        printf("pair %d: burgMag %f, position {%f,%f,%f}, rsize %f\n", nPairs, burgMag, pos[0], pos[1], pos[2], rsize);
        pos[0] += (rsize*gDir[0]);
        pos[1] += (rsize*gDir[1]);
        pos[2] += (rsize*gDir[2]);

        burgMag = -burgMag;

        rsize = delta;
        nPairs--;
        
    }

    WriteDumpFile(inArgs->outFiles[0], dum); 
    MGToLMPDataFile(inArgs->outFiles[0], dum);
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
