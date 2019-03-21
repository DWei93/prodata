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

bool cmpvec(vector<double> p, vector<double> q)
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

    vec.resize(12);

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
        if(auxTable.data.size() == 0){
            Fatal("the aux table size is zero.");
        }

        logFile = 1;
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

        vector<string>().swap(words);
        words = split(secLine, "\"");

        if(words.size() == 0){
            Fatal("can not find the timesolution in file %s", inArgs->inpFiles[file].c_str());
        }
 
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

        if(points1.size() == 0 || points2.size() == 0)continue;
 
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
            data[2][i] = LinearInterpolation(curve2, seq[i], boundMin[0], boundMax[0]);
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

        WriteTecplotNormalData(list, fileName, 10);
        vector<vector<double> >().swap(list.data);

        timenow = atof(words[1].c_str());

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
        
        ZImage(boundMin, boundMax, disp, disp+1, disp+2);
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
    vector<Atom_t *>    nbr;
}Probe_t;

typedef struct {
    long int timestep;
    double x,y,z;
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
    return(p.y<p.y);
}
bool com3(EDState_t p, EDState_t q)
{
    return(p.timestep<q.timestep);
}

void HandleExtendedDislocation_MD(InArgs_t *inArgs)
{
    int         i, j, file, index, firstNbr, indexVar;
    int         nums[2] = {18, 24}, logfile;
    int         stepID;
    double      cutofflen = 2.556, alpha = 0.1, beta = 1.0;
    double      position[3];
    Dump_t    dum;
    LineList_t  list;
    double      dir[3] = {0, 1, 0}, p0[3], dval=5, effeSepRange[2] = {1.0,40.0}; 
    bool        noP0 = 1;
    string      paraName("para"), dvarName("dvar"), dvar("c_vcna"), p0Name("p0"), dirName("dir");
    string      numsName("nums"), ppstring("pp"), separationName("separation");
    int     lastIndex = 0;
    double  dis, minDis;

    vector<Atom_t>::iterator     it; 
    vector<EDState_t>   states;

    ofstream out;
    out.open(inArgs->outFiles[0].c_str(), ios::out);

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

    if(inArgs->help)return;
    logfile = ReadDataFromMDLogFile(inArgs->auxFiles, list);

    states.resize(inArgs->inpFiles.size());
    for(file=0; file<inArgs->inpFiles.size(); file++){
        states[file].timestep = 1E10;
    }

#pragma omp parallel for private(dum, it, i, j, indexVar, firstNbr, lastIndex, dis, minDis) shared(states) 
    for(file=0; file<inArgs->inpFiles.size(); file++){
        #pragma omp critical
        {
            ReadDumpFile(inArgs->inpFiles[file], dum);
        }

        for(i=0; i<dum.variables.size(); i++){
            if(dvar == dum.variables[i])break;
        }
        if((indexVar = i) == dum.variables.size()){
            Fatal("There is no %s in the dum file %s", dvar.c_str(), 
                    inArgs->inpFiles[file].c_str());
        }
        
        for(i=0; i<dum.atom.size(); ){
            if(dum.atom[i].vars[indexVar] != dval){
                dum.atom.erase(dum.atom.begin()+i);
                continue;
            }

            if(dum.atom[i].x < dum.box[0][0] + 2 ||
               dum.atom[i].y < dum.box[1][0] + 2 ||
               dum.atom[i].z < dum.box[2][0] + 2 ||
               dum.atom[i].x > dum.box[0][1] - 2 ||
               dum.atom[i].y > dum.box[1][1] - 2 ||
               dum.atom[i].z > dum.box[2][1] - 2){
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

        probe.y = dum.atom[0].y;
        probe.z = p0[2];

        lastIndex = 0;  
        while(1){
            vector<Atom_t *>().swap(probe.nbr);

            firstNbr = 1;
            minDis = 1E10;
            for(i=lastIndex; i<dum.atom.size(); i++){
                dis = sqrt(pow((dum.atom[i].y-probe.y), 2) +
                      pow((dum.atom[i].z-probe.z), 2) );
                if(dis < minDis){
                    minDis = dis;
                }

                if(dis < cutofflen){
                    if(firstNbr){
                        lastIndex = i;
                        firstNbr = 0;
                    }
                    probe.nbr.push_back(&dum.atom[i]);
                }
            }
            if(probe.nbr.size() > nums[0] && probe.nbr.size()<nums[1]){
                probes.push_back(probe);
                vector<Atom_t *>().swap(probe.nbr);
            }
            if(probe.y > dum.atom.back().y)break;

            if(minDis < alpha){
                probe.x += (dir[0]*alpha);
                probe.y += (dir[1]*alpha);
                probe.z += (dir[2]*alpha);
            }else{
                probe.x += (dir[0]*minDis);
                probe.y += (dir[1]*minDis);
                probe.z += (dir[2]*minDis);
            }
        }

        for(i=0; i<probes.size(); i++){
            probes[i].x = 0.0;
            probes[i].y = 0.0;
            probes[i].z = 0.0;
            for(j=0; j<probes[i].nbr.size(); j++){
                probes[i].x += probes[i].nbr[j]->x;
                probes[i].y += probes[i].nbr[j]->y;
                probes[i].z += probes[i].nbr[j]->z;
            }
            probes[i].x /= (double)probes[i].nbr.size();
            probes[i].y /= (double)probes[i].nbr.size();
            probes[i].z /= (double)probes[i].nbr.size();

        }
        sort(probes.begin(), probes.end(), com2);

        double separation = 0;

        if(probes.size() > 0){
            separation = fabs(probes[0].y-probes.back().y);
            if(separation >  effeSepRange[0] && separation <  effeSepRange[1]){
                states[file].x = (probes[0].x+probes.back().x)*0.5;
                states[file].y = (probes[0].y+probes.back().y)*0.5;
                states[file].z = (probes[0].z+probes.back().z)*0.5;
                states[file].timestep = dum.timestep;
                states[file].separation = separation;

                if(logfile ==1){
                    states[file].vars.resize(list.variables.size()-1);
                    for(i=0; i<list.data[0].size(); i++){
                        if(states[file].timestep == ((int)list.data[0][i]))break;
                    }
                    for(j=0; j<states[file].vars.size(); j++){
                        states[file].vars[j] = list.data[j+1][i];
                    }
                }
            }
        }
        
    }

    sort(states.begin(), states.end(), com3);

    out << "variables = timestep, separation, y";
    if(logfile == 1){
        for(i=1; i<list.variables.size(); i++){
            out << ", " << list.variables[i];
        }
    }
    out << endl;

    for(i=0; i<states.size(); i++){
        if(states[i].timestep > 1E9)break;
        
        out << states[i].timestep << " ";
        out << states[i].separation << " ";
        out << states[i].y  << " ";
        for(j=0; j< states[i].vars.size(); j++){
            out << setprecision(10) << states[i].vars[j] << " "; 
        }
        out << endl;

//        printf("%d %f %f %d\n", states[i].timestep, states[i].separation, 
//                states[i].y, ((int)(states[i].vars.size())));
    }

    out.close();
    return;
}


















