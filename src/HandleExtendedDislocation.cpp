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
            data[1][i] = linearinterpolation(curve1, seq[i], boundmin[0], boundmax[0]) - cubel;
            data[2][i] = linearinterpolation(curve2, seq[i], boundmin[0], boundmax[0]) - cubel;
            data[3][i] = 0.5 * (data[1][i] + data[2][i]);
            if(i<seq.size()-2){
                separation += fabs(data[1][i] - data[2][i]);  
                position += data[3][i];
            }
        }
        separation /= ((double)(seq.size()-1));
        position /= ((double)(seq.size()-1));
        position += cubel;

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
            data[1][i] += cubel;
            data[2][i] += cubel;
            data[3][i] += cubel;
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

typedef struct {
    real8               x,y,z;
    vector<Atom_t *>    nbr;
}Probe_t;

typedef struct {
    int timestep;
    double x,y,z;
    double separation;
    vector<double> vars;
}EDState_t;

bool com(Atom_t p, Atom_t q){
    if(fabs(p.y - q.y) < 0.1){
        if(fabs(p.z - q.z)< 0.1){
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
    int         i, j, file, index, lastIndex, firstNbr, indexVar;
    int         nums[2] = {18, 24}, logfile;
    int         stepID;
    double      cutofflen = 2.556, alpha = 0.005, beta = 1.0;
    double      position[3], separation;
    MgData_t    mg;
    LineList_t  list;
    double      dir[3] = {0, 1, 0}, p0[3], dval=5; 
    bool        noP0 = 1;
    string      paraName("para"), dvarName("dvar"), dvar("c_vcna"), p0Name("p0"), dirName("dir");
    string      numsName("nums");

    vector<Atom_t>::iterator     it; 
    vector<EDState_t>   states;
    EDState_t   state;

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
    printf("The moving direction of probe is along %f, %f, %f\n", dir[0], dir[1], dir[2]);

    if((index = GetValID(inArgs->priVars, numsName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals.size() == 2){
            nums[0] = atoi(inArgs->priVars[index].vals[0].c_str());
            nums[1] = atoi(inArgs->priVars[index].vals[1].c_str());
        }
    }
    printf("The parameters (para): cut-off length %f and resolution value %f\n", cutofflen, alpha);
    logfile = ReadDataFromMDLogFile(inArgs->auxFiles, list);

    for(file=0; file<inArgs->inpFiles.size(); file++){
        ReadMGDataFile(inArgs->inpFiles[file], mg);

        for(i=0; i<mg.variables.size(); i++){
            if(dvar == mg.variables[i])break;
        }
        if((i=indexVar) == mg.variables.size()){
            Fatal("There is no %s in the mg file %s", dvar.c_str(), 
                    inArgs->inpFiles[file].c_str());
        }
        
        for(it = mg.atom.end(), i=mg.atom.size()-1; it != mg.atom.begin(); it--, i--){
            if((*it).x < mg.box[0] + (mg.box[1]-mg.box[0])*0.05 ||
               (*it).y < mg.box[2] + (mg.box[3]-mg.box[2])*0.05 ||
               (*it).z < mg.box[4] + (mg.box[5]-mg.box[4])*0.05 ||
               (*it).x > mg.box[0] + (mg.box[1]-mg.box[0])*0.95 ||
               (*it).y > mg.box[2] + (mg.box[3]-mg.box[2])*0.95 ||
               (*it).z > mg.box[4] + (mg.box[5]-mg.box[4])*0.95){
        
                if(mg.atom[i].vars[indexVar] == dval){
//                    vector<double>().swap(mg.atom[i].vars);
                    mg.atom.erase(it); 
                }
            }
        }

        if(noP0 ==1 && file == 0){
            j = 0;
            for(i=0; i<mg.atom.size(); i++){
                 if(mg.atom[i].vars[indexVar] != dval)continue;
                 j++;
                p0[0] += mg.atom[i].x;
                p0[2] += mg.atom[i].z;
            }
            if(j==0){
                Fatal("The");
            }
            p0[1]  = mg.box[2];
            p0[0] /= ((double)j);
            p0[2] /= ((double)j);
            printf("The initial probe point (p0) is %f %f %f\n", p0[0], p0[1], p0[2]);
        }

        sort(mg.atom.begin(), mg.atom.end(), com);

        vector<Probe_t> probes;
        Probe_t         probe;
        double          d = 0.0;

        probe.y = p0[1];
        probe.z = p0[2];

        lastIndex = 0;
        while(1){
            vector<Atom_t *>().swap(probe.nbr);
            probe.x += (dir[0]*d*cutofflen*alpha);
            probe.y += (dir[1]*d*cutofflen*alpha);
            probe.z += (dir[2]*d*cutofflen*alpha);
            
            firstNbr = 1;
            for(i=lastIndex; i<mg.atom.size(); i++){
                if((pow((mg.atom[i].y-probe.y), 2) + 
                   pow((mg.atom[i].z-probe.z), 2) < cutofflen*cutofflen)
                   && mg.atom[i].vars[indexVar] == dval){
                    if(firstNbr){
                        lastIndex = i;
                        firstNbr = 0;
                    }
                    probe.nbr.push_back(&mg.atom[i]);
                }
            }
            if(probe.nbr.size() > nums[0] && probe.nbr.size()<nums[1]){
                probes.push_back(probe);
                vector<Atom_t *>().swap(probe.nbr);
            }
            if(probe.y > mg.atom.back().y ||
               probe.z > mg.atom.back().z)break;
            d++;
        }

#if 0
        printf("Timestep %d, Atoms %d, Bouds %s %s %s, ", 
               mg.timestep, (int)mg.atom.size(), mg.bounds[0].c_str(), 
               mg.bounds[1].c_str(), mg.bounds[2].c_str());
//        printf("Box %e %e %e %e %e %e\n", mg.box[0], mg.box[1],
//               mg.box[2], mg.box[3], mg.box[4], mg.box[5]);

        printf("variables are ");
        for(i=0; i<mg.variables.size(); i++){
            printf("%s(%d) ", mg.variables[i].c_str(), (int)mg.atom.size());
        }
        printf("\n");
#endif
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

        separation = 0;

        if(probes.size() > 0){
            separation = fabs(probes[0].y-probes.back().y);
            if(separation > 0.1 && separation < 20){
                state.x = (probes[0].x+probes.back().x)*0.5;
                state.y = (probes[0].y+probes.back().y)*0.5;
                state.z = (probes[0].z+probes.back().z)*0.5;
                state.timestep = mg.timestep;
                state.separation = separation;
                if(logfile ==1){
                    state.vars.resize(list.variables.size()-1);
                    for(i=0; i<list.data[0].size(); i++){
                        if(state.timestep == ((int)list.data[0][i]))break;
                    }
                    for(j=0; j<state.vars.size(); j++){
                        state.vars[j] = list.data[j+1][i];
                    }
                }
                states.push_back(state);
            }
        }
        
    }

    if(states.size() == 0){
        Fatal("No out data can be dumped, please check parameters");
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


















