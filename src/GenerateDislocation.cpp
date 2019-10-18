#include <algorithm>
#include <functional>
#include <regex>
#include "math.h"

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
    Dump_t      dum;
    string      posName("pos"), lineName("line"), gDirName("gdir"), burgMagName("burgMag");
    string      screwTypeName("stype"), multiDisName("multi");
    int         screwType = 0, nPairs = 0;
    double      rsize = 0, delta;
    real8       pos[3] = {0.0, 0.0, 0.0};
    real8       line[3] = {1.0, 0.0, 0.0};
    real8       gDir[3] = {0.0, 1.0, 0.0};
    real8       normal[3];
    real8       burgMag = 2.556, boundMin[3], boundMax[3], argument;

    double      k[3] = {-1,-1,-1};
    double      shf[3] = {0,0,0};
    if((index = GetValID(inArgs->priVars, posName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 3){
            char    *token, str[64];
            printf("The dislocation position (pos) is: ");
            for(i=0; i<3; i++){
                snprintf(str, sizeof(str), "%s", inArgs->priVars[index].vals[i].c_str());
                if(strstr(str, "b") != NULL){
                    if((token = strtok(str, "b")) != NULL){
                        k[i] = atof(token);
                        if((token = strtok(NULL, "b"))!=NULL)shf[i] = atof(token);
                        printf(" (%f*box,%f) ", k[i], shf[i]);
                    }
                }else{
                    pos[i] = atof(inArgs->priVars[index].vals[i].c_str());
                    printf(" %f ", pos[i]);
                }
            }
            printf("\n");
        }
    }else{
        k[0] = 0.5;
        k[1] = 0.5;
        k[2] = 0.5;
        printf("The dislcoation (pos) is in the center of the box\n");
    }

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
    cross(line, gDir, normal);
    NormalizeVec(normal);
    FormatVector(normal, "normal");

    if(fabs(DotProduct(line, gDir)) > 1.0E-5)Fatal("the glide direction should be prependicular to the line direction");

    if(inArgs->help)return;

    if(strstr(inArgs->inpFiles[0].c_str(), ".lmp") != NULL){
        ReadLMPFile(inArgs->inpFiles[0], dum);
    }else{
        ReadDumpFile(inArgs->inpFiles[0], dum);
    }

    for(i=0; i<3; i++){
        boundMin[i] = dum.box[i][0];
        boundMax[i] = dum.box[i][1];
        if(k[i] != -1){
            pos[i] = boundMin[i] + k[i]*(boundMax[i]-boundMin[i])+shf[i];
        }
    } 

    real8  point[3], dvec[3], dis, vec[3];
    real8   a, b, shift;

    if(nPairs > 0){
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
                shift = (0.5*burgMag*std::atan2(b, a)/M_PI);
            }else{
                shift = (0.5*burgMag*(M_PI + std::atan2(b, a)/M_PI));
            }
        
            dum.atom[i].x += (shift*line[0]);
            dum.atom[i].y += (shift*line[1]);
            dum.atom[i].z += (shift*line[2]);
        
            FoldBox(7, boundMin, boundMax, &(dum.atom[i].x), &(dum.atom[i].y), &(dum.atom[i].z));
        }

        printf("pair %d: burgMag %f, position {%f,%f,%f}, rsize %f\n", nPairs, burgMag, pos[0], pos[1], pos[2], rsize);
        pos[0] += (rsize*gDir[0]);
        pos[1] += (rsize*gDir[1]);
        pos[2] += (rsize*gDir[2]);

        burgMag = -burgMag;

        rsize = delta;
        nPairs--;
    }    
    }else{
        for(i=0; i<dum.atom.size(); i++){
            FoldBox(7, boundMin, boundMax, &(dum.atom[i].x), &(dum.atom[i].y), &(dum.atom[i].z));
        }
    }

//    WriteDumpFile(inArgs->outFiles[0], dum); 
    MGToLMPDataFile(inArgs->outFiles[0], dum);
    return;
}


#define DisDisPlacement(r) (sqrt(M_PI)*(r)*(15.0+10.0*M_PI*(r)*(r) + 2.0*M_PI*M_PI*(r)*(r)*(r)*(r)) \
                             / (4.0*pow((2.0 + M_PI*(r)*(r)), 2.5)) )

void GenerateExtendedDislocation(InArgs_t *inArgs)
{
    int     index, i;
    string  posName("pos"), lineName("line"), separationName("separation"), gDirName("gdir");
    string  burgName("burg"), pbcName("bound"), debugName("debug");
    real8   boundMax[3], boundMin[3], pos[3];
    real8   line[3] = {1,0,0}, separation = 10.0, gDir[3] = {0,1,0}, normal[3];
    vector<vector<real8> > burgList(2), posList(2);
    int     pbc = 0;

    bool debug=false;
    double      k[3] = {-1,-1,-1};
    double      shf[3] = {0,0,0};
    double  mfactor = 10;
    if((index = GetValID(inArgs->priVars, debugName)) < inArgs->priVars.size()){
        debug = true;
        if(inArgs->priVars[index].vals.size() == 1){
            mfactor = atof(inArgs->priVars[index].vals[0].c_str());
        }
        printf("Debug Model on, ampplification %f times\n", mfactor);
    }

    if((index = GetValID(inArgs->priVars, posName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 3){
            char    *token, str[64];
            printf("The first dislocation position (pos) is: ");
            for(i=0; i<3; i++){
                snprintf(str, sizeof(str), "%s", inArgs->priVars[index].vals[i].c_str());
                if(strstr(str, "b") != NULL){
                    if((token = strtok(str, "b")) != NULL){
                        k[i] = atof(token);
                        if((token = strtok(NULL, "b"))!=NULL)shf[i] = atof(token);
                        printf(" (%f*box,%f) ", k[i], shf[i]);
                    }
                }else{
                    pos[i] = atof(inArgs->priVars[index].vals[i].c_str());
                    printf(" %f ", pos[i]);
                }
            }
            printf("\n");
        }
    }else{
        k[0] = 0.5;
        k[1] = 0.5;
        k[2] = 0.5;
        printf("The first dislcoation (pos) is in the center of the box\n");
    }

    printf("Reading atoms file %s ...\n", inArgs->inpFiles[0].c_str());

    Dump_t    dum;
    if(strstr(inArgs->inpFiles[0].c_str(), ".lmp") != NULL){
        ReadLMPFile(inArgs->inpFiles[0], dum);
    }else{
        ReadDumpFile(inArgs->inpFiles[0], dum);
    }

    for(i=0; i<3; i++){
        boundMin[i] = dum.box[i][0];
        boundMax[i] = dum.box[i][1];
        if(k[i] != -1){
            pos[i] = boundMin[i] + k[i]*(boundMax[i]-boundMin[i])+shf[i];
        }
    } 
    printf("%d atoms have been read\n", (int)dum.atom.size());
    printf("The centor of extended is {%f,%f,%f}\n", pos[0], pos[1], pos[2]);

    if((index = GetValID(inArgs->priVars, separationName)) < inArgs->priVars.size()){
        separation = atof(inArgs->priVars[index].vals[0].c_str());
    }
    printf("The separation of dislocations (separation): %f\n", separation);

    if((index = GetValID(inArgs->priVars, lineName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 3){
            line[0] = atof(inArgs->priVars[index].vals[0].c_str());
            line[1] = atof(inArgs->priVars[index].vals[1].c_str());
            line[2] = atof(inArgs->priVars[index].vals[2].c_str());
        }
        NormalizeVec(line);
    }
    printf("The line direction of dislocation (line): {%f, %f, %f}\n", line[0], line[1], line[2]);

    if((index = GetValID(inArgs->priVars, gDirName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 3){
            gDir[0] = atof(inArgs->priVars[index].vals[0].c_str());
            gDir[1] = atof(inArgs->priVars[index].vals[1].c_str());
            gDir[2] = atof(inArgs->priVars[index].vals[2].c_str());
        }
        NormalizeVec(gDir);
    }
    printf("The glide direction of dislocation (gdir): {%f, %f, %f}\n", gDir[0], gDir[1], gDir[2]);

    if((index = GetValID(inArgs->priVars, burgName)) < inArgs->priVars.size()){

        printf("two partial dislocations: \n");
        for(i=0; i<2; i++){
            burgList[i].resize(3);
            posList[i].resize(3);
            burgList[i][0] = atof(inArgs->priVars[index].vals[3*i].c_str());
            burgList[i][1] = atof(inArgs->priVars[index].vals[3*i+1].c_str());
            burgList[i][2] = atof(inArgs->priVars[index].vals[3*i+2].c_str());

            posList[i][0] = pos[0] + gDir[0]*separation*((real8)i-0.5);
            posList[i][1] = pos[1] + gDir[1]*separation*((real8)i-0.5);
            posList[i][2] = pos[2] + gDir[2]*separation*((real8)i-0.5);
            
            printf("%d: position {%f,%f,%f}, burg {%f,%f,%f}\n", i, posList[i][0], posList[i][1], posList[i][2],
                    burgList[i][0], burgList[i][1], burgList[i][2]);
        }
    }else{
        Fatal("use -dburg to define dislocation lists");
    }

    if((index = GetValID(inArgs->priVars, pbcName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() != 1){
            pbc = 7;
        }else{
            pbc = atof(inArgs->priVars[index].vals[0].c_str());    
        }
        printf("Boundary type (bound): ");
        if((pbc & 0x01) > 0)printf("p ");
        else printf("s ");

        if((pbc & 0x02) > 0)printf("p ");
        else printf("s ");

        if((pbc & 0x04) > 0)printf("p\n");
        else printf("s\n");

    }else{
        printf("Boundary type (bound): free surfaces along 3 dimensions have been applied\n");
    }

    if(inArgs->help)return;

    if(fabs(DotProduct(line, gDir)) > 1.0E-2){
        Fatal("glide direction should be perpendicular to line direction.");
    }

    cross(line, gDir, normal);
    NormalizeVec(normal);
    printf("The glide plane is {%f,%f,%f}\n", normal[0], normal[1], normal[2]);

    int     j;
    real8   point[3], vec1[3], vec2[3], displace, vec[3], dis1, dis2;
    real8   a, b, angle;

#pragma omp parallel for shared(dum) private (j,point,vec,vec1,vec2,dis1,dis2,displace,a,b, angle)
    for(i=0; i<dum.atom.size(); i++){

        point[0] = dum.atom[i].x; if((pbc & 0x01) > 0)point[0] -= pos[0];
        point[1] = dum.atom[i].y; if((pbc & 0x02) > 0)point[1] -= pos[1];
        point[2] = dum.atom[i].z; if((pbc & 0x04) > 0)point[2] -= pos[2];
        
        for(j=0; j<2; j++){
            vec[0] = posList[j][0]; if((pbc & 0x01) > 0)vec[0] -= pos[0];
            vec[1] = posList[j][1]; if((pbc & 0x02) > 0)vec[1] -= pos[1];
            vec[2] = posList[j][2]; if((pbc & 0x04) > 0)vec[2] -= pos[2];

            PointLineIntersection(point, line, vec, vec1, &dis1);
            vec1[0] = point[0] - vec1[0];
            vec1[1] = point[1] - vec1[1];
            vec1[2] = point[2] - vec1[2];

            NormalizeVec(vec1);

            a = DotProduct(vec1, gDir);
            b = DotProduct(vec1, normal); 

            angle = std::atan2(b, a)/M_PI;

//            if(angle < 0)angle += 1;
//            angle /= 2.0;
        
            displace = angle*DisDisPlacement(dis1);
            if(std::isnan(displace)){
                displace = 0;
                printf("nan: %f %f %f %f %f\n", a, b, angle, dis1, DisDisPlacement(dis1));
            }
            if(debug)displace *= mfactor;
            dum.atom[i].x += burgList[j][0]*displace;
            dum.atom[i].y += burgList[j][1]*displace;
            dum.atom[i].z += burgList[j][2]*displace;
            
        }
        if(!debug)FoldBox(pbc, boundMin, boundMax, &(dum.atom[i].x), &(dum.atom[i].y), &(dum.atom[i].z));
    }

    for(i=0; i<2; i++){
        if((pbc & 0x01) == 0){
            dum.box[0][0] -= (burgList[i][0]/2.0);
            dum.box[0][1] += (burgList[i][0]/2.0);
        }
        if((pbc & 0x02) == 0){
            dum.box[1][0] -= (burgList[i][1]/2.0);
            dum.box[1][1] += (burgList[i][1]/2.0);
        }
        if((pbc & 0x04) == 0){
            dum.box[2][0] -= (burgList[i][2]/2.0);
            dum.box[2][1] += (burgList[i][2]/2.0);
        }
    }

    MGToLMPDataFile(inArgs->outFiles[0], dum);

    if(inArgs->outFiles.size() == 2){
        WriteDumpFile(inArgs->outFiles[1], dum); 
    }
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
            GenerateExtendedDislocation(inArgs);
            break;

        case 1:
            GenerateScrewDislocation(inArgs);
            break;

        default:
            GenerateExtendedDislocation(inArgs);
            break;
    }

    return;
}
