#include <algorithm>
#include <functional>

#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"
#include "DDD.h"
#include "MD.h"

using namespace std;


bool cmp(Point_t &p, Point_t &q)
{
    return(p.x < q.x);
}

bool cmpvec(vector<double> &p, vector<double> &q)
{
    return(p[1] < q[1]);
}

void HandleExtendedDislocation_DDD(InArgs_t *inArgs)
{
    Table_t         table;
    real8           cubel = 20000, boundMin[3], boundMax[3];
    real8           remeshSize = 20.0, burgID1 = 7, burgID2 = 12;
    real8           rd, yi1, yj1, yk2, yk3, s, absBakPos, bakPos, cyc;
    int             index,  i, j, k, file, iPos;
    int             colX, colY, colBurgID;
    bool            plu, min, logFile=0;
    string          cubelName = "cubel", burgIDName = "burgID", remeshSizeName = "rsize";  
    string          fileName, secLine, fdir, curDir("./"), fname, aveFile, auxFile;
    real8           separation = 0.0, position = 0.0, position1, position2, timenow, dt;
    real8           pos1[3] = {0,0,0}, pos2[3]={0,0,0}, disp[3] = {0,0,0};
    Point_t         p;
    vector<Point_t> points1, points2;
    Curve_t         curve1, curve2;
    LineList_t      list;
    
    vector<string>          words;
    vector<real8>           seq, vec;
    vector<vector<double> > data(7), outs;  

    Table_t         auxTable;
    vector<Table_t> auxTables;

    if(inArgs->help){
        printf("Function:           Handle extended dislocation (2D) in DDD simulation\n");
        printf("    -dcubel:            cubel size of the cell\n");
        printf("    -dburgID:           two burg indexs of the extended dislocation\n");
        printf("    -dremesh:           remesh size of dislocation\n");
        printf("    -auxfile:           stress-strain curve files\n");
        printf("    -outfile:           <standardized file> <structure file> <stress-strain curve file>\n");
        return;
    }

    vec.resize(12);

    InitList(list);
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

    if(inArgs->outFiles.size() < 3){
        auxFile = "ss.plt";
    }else{
        auxFile = inArgs->outFiles[2];
    }

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

    if(inArgs->auxFiles.size()>0){
        auxTables.resize(inArgs->auxFiles.size());
        for(file=0; file<inArgs->auxFiles.size(); file++){
            ReadTecplotNormalData(inArgs->auxFiles[file], auxTables[file], secLine);
        }
        StitchTecplotData(auxTables, auxTable, 0);
        printf("auxTable size %d\n", (int)auxTable.data.size());
        if(auxTable.data.size() == 0){
            Fatal("the aux table size is zero.");
        }

        logFile = 1;
        auxTable.i=int(auxTable.data.size());
        WriteTecplotNormalData(auxTable, auxFile, 10);
    }
    if(inArgs->help)return;

    ofstream out;
    out.open(aveFile.c_str(), ios::out);
    out << "variables = file, timenow, separation, p1, p2, p, disp,";
    out << " velocity, v_p1, v_p2, v_p, v_separation";

    if(logFile){
        vec.resize(11+auxTable.variables.size());        
        for(i=1; i<auxTable.variables.size(); i++){
            out << ", " << auxTable.variables[i];
        }
    }
    out << endl;    

    bakPos = 0.0;
    for(file=0; file<inArgs->inpFiles.size(); file++){

        vector<string>().swap(table.variables);

        if(!ReadTecplotNormalData(inArgs->inpFiles[file], table, secLine))continue;

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

//        if(points1.size() == 0 || points2.size() == 0)continue;
 
        sort(points1.begin(), points1.end(), cmp); 
        sort(points2.begin(), points2.end(), cmp); 
 
        curve1.ax.resize(points1.size());
        curve1.ay.resize(points1.size());
        curve2.ax.resize(points2.size());
        curve2.ay.resize(points2.size());
 
        plu = 0; min = 0;
        for(i=0; i<points1.size(); i++){
            curve1.ax[i] = points1[i].x;
            curve1.ay[i] = points1[i].y;
            if(curve1.ay[i] > 0.5*boundMax[1])plu = 1;
            if(curve1.ay[i] < 0.5*boundMin[1])min = 1;
        }
 
        for(i=0; i<points2.size(); i++){
            curve2.ax[i] = points2[i].x;
            curve2.ay[i] = points2[i].y;
            if(curve2.ay[i] > 0.5*boundMax[1])plu = 1;
            if(curve2.ay[i] < 0.5*boundMin[1])min = 1;
        }

        if(plu && min){
            for(i=0; i<curve1.ax.size(); i++){
                if(curve1.ay[i] > 0.5*boundMax[1]){
                    curve1.ay[i] -= cubel;
                }
            }

            for(i=0; i<curve2.ax.size(); i++){
                if(curve2.ay[i] > 0.5*boundMax[1]){
                    curve2.ay[i] -= cubel;
                }
            }
        }
 
        vector<double>().swap(seq);
        seq = GenerateSequence(boundMin[0], boundMax[0], remeshSize);
 
        data.resize(7);
        for(i=0; i<data.size(); i++) data[i].resize(seq.size());
        for(i=0; i<seq.size(); i++){
            data[0][i] = seq[i];
        }
        
        position1 = 0.0;
        position2 = 0.0;
        for(i=0; i<seq.size(); i++){
            data[1][i] = LinearInterpolation(curve1, seq[i], boundMin[0], boundMax[0]);
            if(curve2.ax.size()<2) data[2][i]=data[1][i];
            else data[2][i] = LinearInterpolation(curve2, seq[i], boundMin[0], boundMax[0]);
            data[3][i] = 0.5*(data[1][i] + data[2][i]);
            if(i<seq.size()-1){
                position1 += data[1][i];
                position2 += data[2][i];
            }
        }

        position1 /= ((double)(seq.size()-1));
        position2 /= ((double)(seq.size()-1));
        position = 0.5*(position1 + position2);
        separation = fabs(position1 - position2);

        if(position < boundMin[1]){
            position  += cubel;
            position1 += cubel;
            position2 += cubel;
            for(i=0; i<seq.size(); i++){
                data[1][i] += cubel;
                data[2][i] += cubel;
                data[3][i] += cubel;
            }
        }

        absBakPos = (bakPos/cubel - floor(bakPos/cubel))*cubel;
        cyc = floor(bakPos/cubel);        

        if(fabs(absBakPos - position) > 0.5*cubel)cyc++;

        position += (cyc*cubel);
        position1 += (cyc*cubel);
        position2 += (cyc*cubel);

        for(i=0; i<seq.size(); i++){
            data[1][i] += (cubel*cyc);
            data[2][i] += (cubel*cyc);
            data[3][i] += (cubel*cyc);
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
 
        vec[8] = 0.0;
        vec[9] = 0.0;
        vec[10] = 0.0;
        vec[11] = 0.0;

        for(i=0; i<seq.size()-1; i++){
            vec[8] += ((data[1][i] - position1)*(data[1][i] - position1));
            vec[9] += ((data[2][i] - position2)*(data[2][i] - position2));
            vec[10] += ((data[3][i] - position)*(data[3][i] - position));
            vec[11] += ((fabs(data[1][i]-data[2][i]) - separation)*
                       (fabs(data[1][i]-data[2][i]) - separation));
        }

        vec[8] /= (double(seq.size()-1));
        vec[9] /= (double(seq.size()-1));
        vec[10] /= (double(seq.size()-1));
        vec[11] /= (double(seq.size()-1));

        vec[8] = sqrt(vec[8]);
        vec[9] = sqrt(vec[9]);
        vec[10] = sqrt(vec[10]);
        vec[11] = sqrt(vec[11]);

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

        list.i=list.data[0].size();
        WriteTecplotNormalData(list, fileName, 10);
        vector<vector<double> >().swap(list.data);

        vector<string>().swap(words);
        words = split(secLine, "\"");

        timenow = atof(words[1].c_str());
//      printf("T %s\n",table.T.c_str());
//      timenow = atof(table.T.c_str());

        vec[0] = (real8)file;
        vec[1] = timenow;
        vec[2] = separation;
        vec[3] = position1;
        vec[4] = position2;
        vec[5] = position;
        vec[6] = 0.0;
        vec[7] = 0.0;

        bakPos = position;

        outs.push_back(vec);

    }

    sort(outs.begin(), outs.end(), cmpvec); 

    if(logFile){
        int     startID = 0;
        for(i=0; i<outs.size(); i++){
            
            for(j=startID; j<auxTable.data.size(); j++){
                if(outs[i][1] <= auxTable.data[j][0]*1.0E9)break;
            }

            if(j == auxTable.data.size()){
                j--;
            }

            startID = j;

            k = 11;

            for(j=1; j< auxTable.variables.size(); j++){
                k++;
                outs[i][k] = auxTable.data[startID][j];
            }
        }
    }

    for(i=0; i<outs.size(); i++){
        if(i>0){
            pos1[1] = outs[i-1][5]; 
        }else{
            pos1[1] = 0.0;
        }

        dt = (i>0) ? (outs[i][1]-outs[i-1][1]) : outs[i][1];

        pos2[1] = outs[i][5]; 
        disp[1] = pos2[1] - pos1[1];

        if(fabs(disp[1]) > 0.5*cubel){
            for(j=i; j<outs.size(); j++){
                outs[j][5] += cubel;
            }
        }
        
        ZImage(7, boundMin, boundMax, disp, disp+1, disp+2);
        outs[i][6] = disp[1];
        outs[i][7] = disp[1]/dt;

        for(j=0; j<outs[i].size(); j++){
            out <<  setprecision(10)  << outs[i][j] << " ";
        }
        out << endl;
    }
    out.close();
    return;
}

typedef struct {
    real8               x,y,z;
    real8               width,scatter;
    vector<Atom_t *>    nbr;
}Probe_t;

typedef struct {
    long int timestep;
    double x,y,z;
    double absy, absz;
    double vy,vz;
    double separation;
    vector<double> vars;
}EDState_t;

bool com(Atom_t p, Atom_t q){
    if(fabs(p.y - q.y) < 1.0E-20){
        if(fabs(p.z - q.z)< 1.0E-20){
            return(p.x<q.x);
        }else{
            return(p.z < q.z);
        }
    }else{
        return(p.y < q.y);
    }
}

bool com2(Probe_t p, Probe_t q)
{
    return(p.y<q.y);
}
bool com3(EDState_t p, EDState_t q)
{
    return(p.timestep<q.timestep);
}
bool com4(Probe_t p, Probe_t q)
{
    return(p.nbr.size() > q.nbr.size());
}

void HandleExtendedDislocation_MD(InArgs_t *inArgs)
{
    int         i, j, file, index, firstNbr, indexVar;
    int         nums[2] = {18, 24}, logfile=0;
    int         stepID;
    double      cutofflen = 2.556, alpha = 0.1, beta = 1.0;
    double      position[3];
    Dump_t    dum;
    LineList_t  list;
    double      dir[3] = {0, 1, 0}, p0[3], dval=5, effeSepRange[2] = {1.0,40.0}; 
    bool        noP0 = 1;
    string      paraName("para"), dvarName("dvar"), dvar("c_cna"), p0Name("p0"), dirName("dir"), rangeName("range");
    string      numsName("nums"), ppstring("pp"), separationName("separation"), auxName("auxiliary"), marginName("margin");
    int     lastIndex = 0;
    double  dis, minDis, margin_x=5, margin_y=5, margin_z, range[2]={0,1};

    vector<Atom_t>::iterator     it; 
    vector<EDState_t>   states;

    if(inArgs->help){
        printf("Function:           Handle extended dislocation (2D) in MD simulation\n");
        printf("    -dmargin:           margin of surface atoms (disorder atoms)\n");
        printf("    -dpara:             <cut-off length> <resolution>\n");
        printf("    -ddvar:             control variable\n");
        printf("    -p0:                initial position of the extended dislocation\n");
        printf("    -ddir:              dislocation direction\n");
        printf("    -dnums:             range number of partial dislocation core atoms\n");
        printf("    -dseparation:       range separation of the extended dislocations\n");
        printf("    -auxfile:           log file\n");
        printf("    -input:             dump file\n");

        return;
    }
    Table_t table; InitTable(table);

    if((index = GetValID(inArgs->priVars, auxName)) < inArgs->priVars.size()){
        int numV = inArgs->priVars[index].vals.size()/2;
        for(i=0; i<numV; i++)table.aux[inArgs->priVars[index].vals[2*i]] = inArgs->priVars[index].vals[2*i+1];
    }
    if((index = GetValID(inArgs->priVars, marginName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size()==1){
            margin_x = atof(inArgs->priVars[index].vals[0].c_str());
            margin_y = atof(inArgs->priVars[index].vals[0].c_str());
            margin_z = atof(inArgs->priVars[index].vals[0].c_str());
        }else if(inArgs->priVars[index].vals.size() ==3){
            margin_x = atof(inArgs->priVars[index].vals[0].c_str());
            margin_y = atof(inArgs->priVars[index].vals[1].c_str());
            margin_z = atof(inArgs->priVars[index].vals[2].c_str());
        }
    }
    printf("margin (%f,%f,%f)\n", margin_x, margin_y, margin_z);

    if((index = GetValID(inArgs->priVars, paraName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 2){
            cutofflen = atof(inArgs->priVars[index].vals[0].c_str());
            alpha = atof(inArgs->priVars[index].vals[1].c_str());
        }
    }
    printf("The parameters (para): cut-off length %f and resolution value %f\n", cutofflen, alpha);

    if((index = GetValID(inArgs->priVars, dvarName)) < inArgs->priVars.size()){
        dvar = inArgs->priVars[index].vals[0];
        if(inArgs->priVars[index].vals.size() == 2){
            dval = atof(inArgs->priVars[index].vals[1].c_str());
        }
    }
    printf("The determined var (dvar) is %s, value is %f\n", dvar.c_str(), dval);

    if((index = GetValID(inArgs->priVars, p0Name)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() != 3){
            Fatal("You have to determine the initial point by  -dp0 x y z");
        }
        p0[0] = atof(inArgs->priVars[index].vals[0].c_str());
        p0[1] = atof(inArgs->priVars[index].vals[1].c_str());
        p0[2] = atof(inArgs->priVars[index].vals[2].c_str());
        noP0 = 0;
        printf("The initial probe point (p0) is %f %f %f\n", p0[0], p0[1], p0[2]);
    }

    if((index = GetValID(inArgs->priVars, dirName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() != 3){
            Fatal("at least 3 vals for %s", dirName.c_str());
        }
        dir[0] = atof(inArgs->priVars[index].vals[0].c_str());
        dir[1] = atof(inArgs->priVars[index].vals[1].c_str());
        dir[2] = atof(inArgs->priVars[index].vals[2].c_str());
        NormalizeVec(dir);
    }
    printf("The moving direction of probe (dir) is along %f, %f, %f\n", dir[0], dir[1], dir[2]);

    if((index = GetValID(inArgs->priVars, numsName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 2){
            nums[0] = atoi(inArgs->priVars[index].vals[0].c_str());
            nums[1] = atoi(inArgs->priVars[index].vals[1].c_str());
        }
    }
    printf("The range of effective number of a paritial (nums) is [%d,%d]\n", nums[0], nums[1]);

    if((index = GetValID(inArgs->priVars, separationName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 2){
            effeSepRange[0] = atof(inArgs->priVars[index].vals[0].c_str());
            effeSepRange[1] = atof(inArgs->priVars[index].vals[1].c_str());
        }
    }
    printf("The range of effective separtion of the extended dislocation (separation) is [%f,%f]\n", effeSepRange[0], effeSepRange[1]);

    if((index = GetValID(inArgs->priVars, rangeName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 2){
            range[0] = atof(inArgs->priVars[index].vals[0].c_str());
            range[1] = atof(inArgs->priVars[index].vals[1].c_str());
        }
    }
    printf("The determined var (dvar) is %s, value is %f\n", dvar.c_str(), dval);

    if(inArgs->auxFiles.size()>0){
        logfile = ReadDataFromMDLogFile(inArgs->auxFiles, list);
        list.i=(int)list.data[0].size();
        if(inArgs->outFiles.size()>1)WriteTecplotNormalData(list, inArgs->outFiles[1], 10);
    }
    if(inArgs->help)return;

    states.resize(inArgs->inpFiles.size());
    for(file=0; file<inArgs->inpFiles.size(); file++){
        states[file].timestep = 1E10;
    }

    double boundMin[3], boundMax[3];
#pragma omp parallel for private(dum, it, i, j, indexVar, firstNbr, lastIndex, dis, minDis, boundMin, boundMax) shared(states) 
    for(file=0; file<inArgs->inpFiles.size(); file++){
        double separation, path[3];
        #pragma omp critical
        {
            ReadDumpFile(inArgs->inpFiles[file], dum);
            printf("Finish reading %s, %d atoms\n",inArgs->inpFiles[file].c_str(), (int)dum.atom.size());
            ModifyDumpFile(inArgs->inpFiles[file], dum);
            printf("After modify %d\n", (int)dum.atom.size());
            WriteDumpFile(inArgs->outFiles[file], dum); 
        }

            continue;
        for(i=0; i<3; i++){
            boundMin[i] = dum.box[i][0];
            boundMax[i] = dum.box[i][1];
        }
        for(i=0; i<dum.variables.size(); i++){
            if(dvar == dum.variables[i])break;
        }
        if((indexVar = i) == dum.variables.size()){
            Fatal("There is no %s in the dum file %s", dvar.c_str(), 
                    inArgs->inpFiles[file].c_str());
        }
        
        for(i=0; i<dum.atom.size(); ){
            if(dum.atom[i].x < dum.box[0][0] + margin_x ||
               dum.atom[i].y < dum.box[1][0] + margin_y ||
               dum.atom[i].z < dum.box[2][0] + margin_z ||
               dum.atom[i].x > dum.box[0][1] - margin_x ||
               dum.atom[i].y > dum.box[1][1] - margin_y ||
               dum.atom[i].z > dum.box[2][1] - margin_z){
                dum.atom.erase(dum.atom.begin()+i);
                continue;
            }
            if(dum.atom[i].vars[indexVar] != dval){
                dum.atom.erase(dum.atom.begin()+i);
                continue;
            }

            i++;
        }

        if(dum.atom.size() == 0){
            continue;
        }
        
        if(noP0 ==1 && file == 0){
            Fatal("The initial probe point (p0) is %f %f %f\n", p0[0], p0[1], p0[2]);
        }

        sort(dum.atom.begin(), dum.atom.end(), com);

        vector<Probe_t> probes;
        Probe_t         probe;

        probe.y = p0[1];
        probe.z = p0[2];

        lastIndex = 0;  
        while(1){
            vector<Atom_t *>().swap(probe.nbr);

            firstNbr = 1;
            minDis = 1E10;
            for(i=lastIndex; i<dum.atom.size(); i++){
                path[0] = 0;
                path[1] = dum.atom[i].y-probe.y;
                path[2] = dum.atom[i].z-probe.z;
                ZImage(7, boundMin, boundMax, path, path+1, path+2);
                dis = sqrt(path[0]*path[0]+path[1]*path[1]+path[2]*path[2]);

                if(dis < minDis){
                    minDis = dis;
                }

                if(dis < cutofflen){
                    if(firstNbr){
                        lastIndex = i;
                        firstNbr = 0;
                    }
                    probe.nbr.push_back(&(dum.atom[i]));
                }
            }
            if(probe.nbr.size() > nums[0] && probe.nbr.size()<nums[1]){
                probes.push_back(probe);
                vector<Atom_t *>().swap(probe.nbr);
            }

            if(minDis < alpha){
                probe.x += (dir[0]*alpha);
                probe.y += (dir[1]*alpha);
                probe.z += (dir[2]*alpha);
            }else{
                probe.x += (dir[0]*minDis);
                probe.y += (dir[1]*minDis);
                probe.z += (dir[2]*minDis);
            }
            if(probe.y > dum.atom.back().y+alpha)break;
            FoldBox(7, boundMin, boundMax, &(probe.x), &(probe.y), &(probe.z));
        }

        for(i=0; i<probes.size(); i++){
            double bakPos[3] = {probes[i].x,probes[i].y,probes[i].z};
            probes[i].x=0;  probes[i].y=0;  probes[i].z=0;
            for(j=0; j<probes[i].nbr.size(); j++){
                path[0] = probes[i].nbr[j]->x - bakPos[0];
                path[1] = probes[i].nbr[j]->y - bakPos[1];
                path[2] = probes[i].nbr[j]->z - bakPos[2];
                ZImage(7, boundMin, boundMax, path, path+1, path+2);
                probes[i].x += path[0]; probes[i].y += path[1]; probes[i].z += path[2];
            }
            probes[i].x /= (double)probes[i].nbr.size();
            probes[i].y /= (double)probes[i].nbr.size();
            probes[i].z /= (double)probes[i].nbr.size();
            probes[i].x += bakPos[0]; probes[i].y += bakPos[1]; probes[i].z += bakPos[2];
            FoldBox(7, boundMin, boundMax, &(probes[i].x), &(probes[i].y), &(probes[i].z));

            probes[i].width=0;
            for(j=0; j<probes[i].nbr.size(); j++){
                path[0] = probes[i].nbr[j]->x - probes[i].x;
                path[1] = probes[i].nbr[j]->y - probes[i].y;
                path[2] = probes[i].nbr[j]->z - probes[i].z;

                ZImage(7, boundMin, boundMax, path, path+1, path+2);
                probes[i].width += sqrt(path[1]*path[1]+path[2]*path[2]);
            }
            probes[i].width /= double(probes[i].nbr.size());

            probes[i].scatter = 0;
            for(j=0; j<probes[i].nbr.size(); j++){
                path[0] = probes[i].nbr[j]->x - probes[i].x;
                path[1] = probes[i].nbr[j]->y - probes[i].y;
                path[2] = probes[i].nbr[j]->z - probes[i].z;
                ZImage(7, boundMin, boundMax, path, path+1, path+2);
                probes[i].scatter += pow(probes[i].width-sqrt(path[1]*path[1]+path[2]*path[2]),2);
            }
            probes[i].scatter = sqrt(probes[i].scatter)/double(probes[i].nbr.size());
        }
        
        separation = 0;

        if(probes.size() > 0){
            
            printf("File %s: ", inArgs->inpFiles[file].c_str());
            path[0] = probes.back().x-probes[0].x;
            path[1] = probes.back().y-probes[0].y;
            path[2] = probes.back().z-probes[0].z;
            ZImage(7, boundMin, boundMax, path, path+1, path+2);
            separation = sqrt(path[0]*path[0]+path[1]*path[1]+path[2]*path[2]);
            if(effeSepRange[0] > 0 && effeSepRange[1] >0 ){
                 sort(probes.begin(), probes.end(), com2);
                if(separation >  effeSepRange[0] && separation <  effeSepRange[1] ){
                    states[file].x = probes[0].x + path[0]*0.5;
                    states[file].y = probes[0].y + path[1]*0.5;
                    states[file].z = probes[0].z + path[2]*0.5;
                    FoldBox(7, boundMin, boundMax, &(states[file].x), &(states[file].y), &(states[file].z));
                    states[file].timestep = dum.timestep;
                    states[file].separation = separation;
                    
                    path[0] =0;
                    NormalizeVec(path);
                    printf("(%f,%f,#%d,%f,%f)--%f--LD(0,%f,%f)-->(%f,%f,#%d,%f,%f)\n", probes[0].y,probes[0].z,int(probes[0].nbr.size()),probes[0].width,
                            probes[0].scatter, separation,
                            path[1], path[2],
                            probes.back().y,probes.back().z,int(probes.back().nbr.size()), probes.back().width,probes.back().scatter);
 
                    if(logfile ==1){
                        states[file].vars.resize(list.variables.size()-1);
                        for(i=0; i<list.data[0].size(); i++){
                            if(states[file].timestep == ((int)list.data[0][i]))break;
                        }
                        for(j=0; j<states[file].vars.size(); j++){
                            states[file].vars[j] = list.data[j+1][i];
                        }
                    }
                }else printf("%d probes, 1: #%d width %f, scatter %f; 2: #%d width %f, scatter %f, separation %f, Skip!\n",
                            int(probes.size()),int(probes[0].nbr.size()),probes[0].width,probes[0].scatter,
                                           int(probes.back().nbr.size()),probes.back().width, probes.back().scatter, separation);
            }else{
                sort(probes.begin(), probes.end(), com4);
                
                states[file].x = probes[0].x;
                states[file].y = probes[0].y;
                states[file].z = probes[0].z;
                states[file].timestep = dum.timestep;
                states[file].separation = 0;
                if(logfile ==1){
                    states[file].vars.resize(list.variables.size()-1);
                    for(i=0; i<list.data[0].size(); i++){
                        if(states[file].timestep == ((int)list.data[0][i]))break;
                    }
                    for(j=0; j<states[file].vars.size(); j++){
                        states[file].vars[j] = list.data[j+1][i];
                    }
                }
                printf("(%f,%f,#%d,%f,%f)\n", probes[0].y,probes[0].z,int(probes[0].nbr.size()),probes[0].width, probes[0].scatter);
            }
        }
        
    }

    sort(states.begin(), states.end(), com3);

    table.variables.push_back("timestep");
    table.variables.push_back("separation");
    table.variables.push_back("y");
    table.variables.push_back("absy");
    table.variables.push_back("vy");
    table.variables.push_back("z");
    table.variables.push_back("absz");
    table.variables.push_back("vz");
    if(logfile == 1){
        for(i=1; i<list.variables.size(); i++){
            table.variables.push_back(list.variables[i]);
        }
    }

    real8    lastAbsy, vy, lasty, lastTimestep;
    real8    lastAbsz, vz, lastz;
    ReadDumpFile(inArgs->inpFiles[0], dum);
    real8    ybox = dum.box[1][1]-dum.box[1][0], zbox = dum.box[2][1] -dum.box[2][0];
    if(states.size()==0)Fatal("can not find any dislocation");
    int initStep = states[0].timestep;
    for(i=0; i<states.size(); i++){
        if(states[i].timestep > 1E9)break;}
    table.data.resize(i);

    printf("effective states %d, boxy is %f x %f\n", (int)table.data.size(), ybox, zbox);
    for(i=0; i<table.data.size(); i++){
        j=0;
        states[i].timestep -= initStep;
        table.data[i].resize(table.variables.size());
        lastAbsy = ((i==0) ? states[i].y : states[i-1].absy);
        lastAbsz = ((i==0) ? states[i].z : states[i-1].absz);

        lasty = ((i==0) ? states[i].y : states[i-1].y);
        lastz = ((i==0) ? states[i].z : states[i-1].z);

        lastTimestep = ((i==0) ? states[i].timestep : states[i-1].timestep);

        table.data[i][j++] = states[i].timestep;
        table.data[i][j++] = states[i].separation;

        table.data[i][j++] = states[i].y;
        vy = states[i].y - lasty;
        vy -= rint(vy * 1.0/ybox) * ybox;
        states[i].absy = lastAbsy+vy;
        table.data[i][j++] =  states[i].absy;
        states[i].vy = ((i==0) ? 0 : (vy)/(states[i].timestep - lastTimestep));
        table.data[i][j++] = states[i].vy;

        table.data[i][j++] = states[i].z;
        vz = states[i].z - lastz;
        vz -= rint(vz * 1.0/zbox) * zbox;
        states[i].absz = lastAbsz+vz;
        table.data[i][j++] =  states[i].absz;
        states[i].vz = ((i==0) ? 0 : (vz)/(states[i].timestep - lastTimestep));
        table.data[i][j++] = states[i].vz;

        for(auto k=0; k<states[i].vars.size(); k++){
            table.data[i][j+k] = states[i].vars[k];
        }
    }


    int initial=int(double(table.data.size())*range[0]-1); if(initial < 0)initial=0;
    int last=int(double(table.data.size())*range[1]-1); if(last >= table.data.size())last= int(table.data.size())-1;
    double velocity_y = (states[last].absy-states[initial].absy)/(states[last].timestep-states[initial].timestep);
    double velocity_z = (states[last].absz-states[initial].absz)/(states[last].timestep-states[initial].timestep);

    char val[50];
    snprintf(val, sizeof(val), "%e", sqrt(velocity_z*velocity_z+velocity_y*velocity_y));
    table.aux[std::string("velocity")] = val;

    double width=0;
    for(i=initial; i<last+1; i++){
        width += states[i].separation;
    }  
    width /= double(last+1-initial);
    snprintf(val, sizeof(val), "%e", width);
    table.aux[std::string("separation")] = val;
    table.i=int(table.data.size());

    WriteTecplotNormalData(table,  inArgs->outFiles[0], 10);
 
    return;
}


















