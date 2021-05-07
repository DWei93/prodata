
#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"


void Fatal(const char *format, ...) 
{
        char    msg[512];
        va_list args;

        va_start(args, format);
        vsnprintf(msg, sizeof(msg)-1, format, args);
        msg[sizeof(msg)-1] = 0;
        va_end(args);
        printf("ERROR: %s\n", msg);

        exit(1);
}

int DataType(const string &str)
{
    int     i, c = 0, n = 0, type;
    char    s[512];

    str.copy(s, str.length(), 0);
    *(s+str.length()) = '\0';

    for(i=0; i<str.length(); i++){
        if(s[i] >= '0' && s[i] >= '9' && n==0){
            n = 1;
            if(c)break;
        }

        if((s[i] >= 'a' && s[i] <= 'z' ||
           s[i] >= 'A' && s[i] <= 'Z') && c==0){
            c = 1;
            if(n)break;
        }
    }

    type |= ((n==1) ? NUMBER_DATA : 0);
    type |= ((c==1) ? CHAR_DATA : 0);
    printf("%d %d, cd\n", n, c);

    return(type);
}


vector<string> split(const string& str, const string& delim) 
{
    vector<string> res;
    if("" == str) return res;
    char * strs = new char[str.length() + 1] ; 
    strcpy(strs, str.c_str()); 
    char * d = new char[delim.length() + 1]; 
    strcpy(d, delim.c_str());
    char *p = strtok(strs, d); 
    while(p) {
        string s = p;  
        res.push_back(s); 
        p = strtok(NULL, d); 
    }   
    return res;
}

void WashString(string &str)
{
    string::iterator it;
    for (it =str.begin(); it != str.end(); ++it){
        if ( *it == ' ')str.erase(it);
    }
    return;
}

int  GetValID(const vector<Var_t> &vals, const string &name)
{
    int             i;

    for(i=0; i<vals.size(); i++){
        if(vals[i].name == name)break;
    }

    return(i);    
}

int GetColIDFromTable(const Table_t &table, const string &name)
{
    int  i, j=-1;
    
    for(i=0; i<table.variables.size(); i++){
        if(name == table.variables[i]){j=i; break;}
    } 

    return(j);
}

vector<double>  GenerateSequence(double from, double to, double meshSize)
{
    double  nums, sign;
    int     i, n;
    
    nums = (to-from)/fabs(meshSize);
    sign = (nums >= 0) ? 1.0:-1.0;
    nums = fabs(nums);

    if(nums - floor(nums) != 0){
        n = (int)(floor(nums) + 1);
    }else{
        n = (int)nums;
    }

    meshSize = (to-from)/((double)n);

    n++;
    vector<double>  seq(n);
    
    for(i=0; i<n; i++){
        seq[i] = from + (double(i))*fabs(meshSize)*sign;
    }
    return(seq);
}

void FoldBox(int pbc, real8 boundMin[3], real8 boundMax[3], real8 *x, real8 *y, real8 *z)
{
        real8   xc, yc, zc;

        real8   invLx, invLy, invLz, Lx, Ly, Lz;

        if((pbc & 0x01) > 0){
            xc = (boundMin[0] + boundMax[0]) * 0.5;
            Lx = boundMax[0] - boundMin[0];
            invLx = 1.0/Lx;
            *x -= rint((*x-xc)*invLx) * Lx;
        }

        if((pbc & 0x02) > 0){
            yc = (boundMin[1] + boundMax[1]) * 0.5;
            Ly = boundMax[1] - boundMin[1];
            invLy = 1.0/Ly;
            *y -= rint((*y-yc)*invLy) * Ly;
        }

        if((pbc & 0x04) > 0){
            zc = (boundMin[2] + boundMax[2]) * 0.5;
            Lz = boundMax[2] - boundMin[2];
            invLz = 1.0/Lz;
            *z -= rint((*z-zc)*invLz) * Lz;
        }

        return;
}

void ZImage(int pbc, real8 boundMin[3], real8 boundMax[3], real8 *x, real8 *y, real8 *z)
{
        real8   invLx, invLy, invLz, Lx, Ly, Lz;

/*
 *      If periodic boundaries are not in use, the provided position
 *      of (x,y,z) will not be adjusted since there are no other
 *      images available.
 */
        if((pbc & 0x01) > 0){
             Lx = boundMax[0] - boundMin[0];
             invLx = 1.0/Lx;
             *x -= rint(*x * invLx) * Lx;
        }
         
        if((pbc & 0x02) > 0){
             Ly = boundMax[1] - boundMin[1];
             invLy = 1.0/Ly;
             *y -= rint(*y * invLy) * Ly;
        }

        if((pbc & 0x04) > 0){
             Lz = boundMax[2] - boundMin[2];
             invLz = 1.0/Lz;
             *z -= rint(*z * invLz) * Lz;
        }    
        return;
}

real8 LinearInterpolation(const Curve_t &curve, real8 x, real8 min, real8 max)
{
    real8   length, xc;
    real8   x0, y0, x1, y1, t, boundVal;
    int     i;    

    if(curve.ax.size()<2)return 0.0;
    if(curve.ax[0] == curve.ax.back())
        Fatal("the range of line is zero in LineInterpolation");

    if(min == max){
        min = curve.ax[0];
        max = curve.ax.back();
    }
 
    length = max - min;
    if(length == 0)Fatal("the range of line is zero in LineInterpolation");

    xc = 0.5*(min + max);    
    x -= rint((x-xc)/length) * length;
    

    if(x>=curve.ax[0] && x<=curve.ax.back()){
        for(i=1; i<curve.ax.size(); i++){
            if(x<=curve.ax[i])break;
        }
        x0 = curve.ax[i-1];
        y0 = curve.ay[i-1];
        x1 = curve.ax[i];
        y1 = curve.ay[i];
    }else{
        x0 = curve.ax.back() - length; 
        y0 = curve.ay.back();

        x1 = curve.ax[0]; 
        y1 = curve.ay[0];

        boundVal = y0 + (y1-y0)*(x1-min)/(x1-x0); 

        if(x < curve.ax[0]){
            x0 = min;
            y0 = boundVal;
            
            x1 = curve.ax[0]; 
            y1 = curve.ay[0]; 
        }else{
            x0 = curve.ax.back();
            y0 = curve.ay.back();
            
            x1 = max;
            y1 = boundVal;
        }
    }

    if(x<x0 || x>x1)
        Fatal("something wrong with linearInterpolation %f, %f, %f %f %f",
              x, x0, x1, min, max);

    if(x1-x0 == 1.0E-10){
        return(y0);
    }else{
        return(y0+(y1-y0)*(x-x0)/(x1-x0));
    }
} 

void SwapTable(Table_t &table){
    vector<string>().swap(table.variables);
    vector<vector<double> >().swap(table.data);
    return;
}

void SwapLineList(LineList_t &list){
    vector<string>().swap(list.variables);
    vector<vector<double> >().swap(list.data);
    return;
}


void CleanDump(Dump_t &dum)
{
    int i;
    dum.timestep = 0;

    vector<vector<double> >().swap(dum.box);
    vector<string>().swap(dum.bounds);
    vector<string>().swap(dum.variables);
    
    for(i=0; i<dum.atom.size(); i++){
        vector<double>().swap(dum.atom[i].vars);
    }
    vector<Atom_t>().swap(dum.atom);

    return;
}


void NormalizeVec(real8 vec[3])
{
        real8 a2, a;

        a2 = (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        if (a2 > 0.0) {
            a = sqrt(a2);
            vec[0] /= a;
            vec[1] /= a;
            vec[2] /= a;
        }

        return;
}

void StitchTecplotData(vector<Table_t> &tables, Table_t &table, int eigenID)
{
    int         i, j, startID = 0;

    InitTable(table);

    table.variables.resize(tables[0].variables.size());
    for(i=0; i<tables[0].variables.size(); i++){
        table.variables[i] = tables[0].variables[i];
    }

    for(i=0; i<tables[0].data.size(); i++){
        table.data.push_back(tables[0].data[i]);
    }    

    if(tables.size() == 1)return;

    for(i=1; i<tables.size(); i++){
        if(tables[i].variables.size() != tables[0].variables.size())continue;

        for(j=startID; j<table.data.size(); j++){
            if(table.data[j][eigenID] >= tables[i].data[0][eigenID]){
                break;
            }
        }
        startID = j;

        if(j < table.data.size()){
            table.data.erase(table.data.begin()+j, table.data.end());
        }
       
        for(j=0; j<tables[i].data.size(); j++){
            table.data.push_back(tables[i].data[j]);
        }
    }

    table.i = (int)table.data.size();
    return;            
}



void SpecifyEquations_PLTDATA(InArgs_t *inArgs)
{
    int     i, index, pointID=0;
    bool    readState;
    string  secLine, noBackupName = ("nobackup"), backupFile, sufBack = (".bak");
    char    name[256];
    bool    backup = 1;

    vector<Table_t> tables;

    if((index = GetValID(inArgs->priVars, noBackupName)) < inArgs->priVars.size()){
        backup = 0;
        printf("Backup: NO! (nobackup)\n");
    }

    if(inArgs->help)return;

    if(inArgs->inpFiles.size() == 0){
        Fatal("There is no input file.");
    }
    tables.resize(inArgs->inpFiles.size());

    for(i=0; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState)continue;
        if(backup){
            backupFile = inArgs->inpFiles[i] + sufBack;
            WriteTecplotNormalData(tables[i], backupFile, 10, secLine); 
        }
        
        SpecifyEquations(tables[i]);
        if(inArgs->outFiles.size() == inArgs->inpFiles.size()){
            WriteTecplotNormalData(tables[i], inArgs->outFiles[i], 10, secLine); 
        }else{
            WriteTecplotNormalData(tables[i], inArgs->inpFiles[i], 10, secLine); 
        }
    }
    return;
}


void AnimateAuxData(InArgs_t *inArgs){
    int index;
    string xName=("x"), yName=("y"), x("strain"), y("stress"), secLine, fileName, auxName("plt");

    if((index = GetValID(inArgs->priVars, xName)) < inArgs->priVars.size()){
        x=inArgs->priVars[index].vals[0];
    }   
    if((index = GetValID(inArgs->priVars, yName)) < inArgs->priVars.size()){
        y=inArgs->priVars[index].vals[0];
    }   


    bool    readState;
    Table_t auxTable;
    auxTable.variables.push_back(x.c_str());
    auxTable.variables.push_back(y.c_str());
    auxTable.variables.push_back("v1");
    auxTable.variables.push_back("v2");
    auxTable.i=1; auxTable.j=1; auxTable.k=1;

    map<string, string>::iterator   iter;
    vector<double> point(4); point[2]=0; point[3]=0;

    vector<Table_t> tables;
    tables.resize(inArgs->inpFiles.size());
    double xlast, ylast;
    int numData=0;
    bool first=true;
    for(int i=0; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState)continue;

        if((iter = tables[i].aux.find(x)) != tables[i].aux.end()){
                 point[0] = atof(iter->second.c_str());
        }else Fatal("can not find variable x");
        if((iter = tables[i].aux.find(y)) != tables[i].aux.end()){
                 point[1] = atof(iter->second.c_str());
        }else Fatal("can not find variable x");
        if(first){numData++; auxTable.data.push_back(point); first=false; xlast=point[0]; ylast=point[1];}
        else{point[2]=point[0]-xlast; point[3]=point[1]-ylast; xlast=point[0]; ylast=point[1];
            auxTable.data.push_back(point);
            numData++;
        }

        auxTable.T=tables[i].T; auxTable.i = int(auxTable.data.size()); auxTable.F="Point"; auxTable.solutionTime=tables[i].solutionTime;
        fileName = y + "-" + x + string(auxTable.T) +  auxName;
        WriteTecplotNormalData(auxTable, fileName,  10, secLine); 
    }
}

void HandleTecplotData(InArgs_t *inArgs)
{
    int     i, index, pointID=0;
    bool    readState;
    string  secLine, noBackupName = ("nobackup"), backupFile, sufBack = (".bak");
    string  pointIDName=("pointID"), sfName("sf");
    char    name[256];
    bool    backup = 1;
    double  sf=-1.0;

    vector<Table_t> tables;

    if((index = GetValID(inArgs->priVars, noBackupName)) < inArgs->priVars.size()){
        backup = 0;
        printf("Backup: NO! (nobackup)\n");
    }

    if(inArgs->help)return;

    if(inArgs->inpFiles.size() == 0){
        Fatal("There is no input file.");
    }

    if((index = GetValID(inArgs->priVars, sfName)) < inArgs->priVars.size()){
        sf = atof(inArgs->priVars[index].vals[0].c_str());
        printf("Schmid Factor = %f \n", sf);
    }

    if((index = GetValID(inArgs->priVars, pointIDName)) < inArgs->priVars.size()){
    }else index =-1;

    tables.resize(inArgs->inpFiles.size());
#ifndef GSL
    real8 crss;
    vector<real8> crsses;
        
    for(i=0; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState)continue;
        if(backup){
            backupFile = inArgs->inpFiles[i] + sufBack;
            WriteTecplotNormalData(tables[i], backupFile, 10, secLine); 
        }
        if(true==Analysis(tables[i], crss)){
            crsses.push_back(crss);
        }
    }

    crss=0;
    real8 scatter=0;
    for(i=0;i<crsses.size();i++)crss+=crsses[i];
    crss /= double(crsses.size());
    for(i=0;i<crsses.size(); i++)scatter += pow(crsses[i]-crss,2);
    
    printf("DATA: %g %g\n",crss,sqrt(scatter)/double(crsses.size()));

    return;
#else
    vector<vector<real8> > dl;
    vector<real8> d(10);
    real8 sigma, twindef, hard, thard, crss;
    for(i=0; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState)continue;
        if(backup){
            backupFile = inArgs->inpFiles[i] + sufBack;
            WriteTecplotNormalData(tables[i], backupFile, 10, secLine); 
        }
        
        if(index>-1 && i/3<inArgs->priVars[index].vals.size())pointID = atoi(inArgs->priVars[index].vals[i/3].c_str());
        else if (index < 0) pointID=-1;
        else pointID=atoi(inArgs->priVars[index].vals[0].c_str());

        if(tables[i].aux.empty() == false){
            auto a = tables[i].aux.find("schmid");
            if(a!=tables[i].aux.end()){
                sf = atof(a->second.c_str());
                printf("Schmid factor from aux data %f\n",sf);
            }
        }
        
        if(true==Analysis(pointID, sf, tables[i], sigma, hard, thard, twindef, crss)){
            d[0]=sigma; d[1]=hard; d[2]=thard; d[3]=twindef; d[4]=0; d[5]=0; d[6]=0; d[7]=0; d[8]=crss; d[9]=0;
            dl.push_back(d);
        }else{
            printf("FAIELD: %s \n", inArgs->inpFiles[i].c_str());
        }
//      SpecifyEquations(tables[i]);
//  
        WriteTecplotNormalData(tables[i], inArgs->inpFiles[i], 10, secLine); 
    }

    for(int j=0; j<10; j++)d[j]=0; real8 num = double(dl.size());
    for(int j=0; j<dl.size(); j++){
        d[0] += (dl[j][0]/num);
        d[1] += (dl[j][1]/num);
        d[2] += (dl[j][2]/num);
        d[3] += (dl[j][3]/num);
        d[8] += (dl[j][8]/num);
    }

    if(num>1){
        for(int j=0;j<dl.size();j++){
            d[4]+=pow(dl[j][0]-d[0],2)/(num-1.0);
            d[5]+=pow(dl[j][1]-d[1],2)/(num-1.0);
            d[6]+=pow(dl[j][2]-d[2],2)/(num-1.0);
            d[7]+=pow(dl[j][3]-d[3],2)/(num-1.0);
            d[9]+=pow(dl[j][8]-d[8],2)/(num-1.0);
        }
    }

    printf("DATA: %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e\n",
        d[0],d[1],d[2],d[3],sqrt(d[4]),sqrt(d[5]),sqrt(d[6]),sqrt(d[7]), d[8], sqrt(d[9]));
    return;
#endif
}

void FormatVector(real8 vec[3], const char *msg){
	printf("%s vector:\n",msg);
	printf("{%.15f,%.15f,%.15f}\n", vec[0], vec[1], vec[2]);
}

void InitList(LineList_t &list){
    vector<string>().swap(list.variables);
    vector<vector<real8> >().swap(list.data);
    list.aux.clear();
    list.i=1;
    list.j=1;
    list.k=1;
    list.solutionTime = -1;
    list.T="";
    list.F = "Point";
    return;
}

void InitTable(Table_t &table){
    vector<string>().swap(table.variables);
    vector<vector<real8> >().swap(table.data);
    table.aux.clear();
    table.i=1;
    table.j=1;
    table.k=1;
    table.solutionTime = -1;
    table.T="";
    table.F = "Point";
    return;
}

bool cmp(vector<double> &p, vector<double> &q)
{
    return p[0]<q[0];
}
void SortTable(Table_t &table, int sortColID)
{
    swap(table.variables[0], table.variables[sortColID]);
    for(auto &a : table.data){
        swap(a[0],a[sortColID]);
    }
    sort(table.data.begin(), table.data.end(),cmp);
}
bool FindLinearPart(real8 (*line)[3], const int nums, int range[2])
{
    if(nums < 3)return 0;

    int i, j;

    return 0;
} 



void AnimateCurve(InArgs_t *inArgs){
    int index, fps=100,i;
    string xName=("x"), yName=("y"), x("strain"), y("stress"), nPointsName("fps"), secLine, fileName("aux.plt"), auxName("plt");

    if((index = GetValID(inArgs->priVars, xName)) < inArgs->priVars.size()){
        x=inArgs->priVars[index].vals[0];
    }   
    if((index = GetValID(inArgs->priVars, yName)) < inArgs->priVars.size()){
        y=inArgs->priVars[index].vals[0];
    }   

    if((index = GetValID(inArgs->priVars, nPointsName)) < inArgs->priVars.size()){
        fps=atoi(inArgs->priVars[index].vals[0].c_str());
    }   

    if(inArgs->outFiles.size()>0)fileName=inArgs->outFiles[0];

    Table_t table;
    bool    readState;
    if(inArgs->inpFiles.size() == 0)Fatal("no inpurt file");
    readState = ReadTecplotNormalData(inArgs->inpFiles[0], table, secLine);
    if(!readState)Fatal("can not read file %s", inArgs->inpFiles[0].c_str());

    Table_t auxTable;
    int cx = GetColIDFromTable(table, x);
    int cy = GetColIDFromTable(table, y);
    auxTable.variables.resize(2); 
    auxTable.variables[0]=x;
    auxTable.variables[1]=y;
    auxTable.i=1; auxTable.j=1; auxTable.k=1;

    int nPoints =int(table.data.size());
    if(nPoints<fps)Fatal("fps(%d) > nPoints(%d), please reset fps through -dfps",fps,nPoints);

    int freq=nPoints/fps;
    bool addLast=false;
    if((nPoints%fps)!=0){
        freq = nPoints/(fps-1); addLast=true;
    }

    vector<double> point(2);
    if(cx<0 || cy<0)Fatal("can not find %s or %s",x.c_str(),y.c_str());

    double x_range[2]={table.data[0][cx],table.data[0][cx]}, y_range[2]={table.data[0][cy],table.data[0][cy]};

    int nZones=0;
    for(i=0; i<table.data.size(); i++){
        point[0]=table.data[i][cx];
        point[1]=table.data[i][cy];
        if(point[1]>1E10)point[1]/=1E12;
        if(point[0]<x_range[0])x_range[0]=point[0];
        if(point[0]>x_range[1])x_range[1]=point[0];
        if(point[1]<y_range[0])y_range[0]=point[1];
        if(point[1]>y_range[1])y_range[1]=point[1];
        auxTable.data.push_back(point);
        if((i%freq)==0){
            auxTable.solutionTime = table.data[i][0];
            auxTable.T=to_string(table.data[i][0]); auxTable.i = i+1; auxTable.F="Point"; 
            nZones++;
            if(i==0){
                WriteTecplotNormalData(auxTable, fileName,  10, secLine, std::ios::out); 
            }else{
                WriteTecplotNormalData(auxTable, fileName,  10, secLine, std::ios::app); 
            }
        }
    }
    if(x_range[0]<1E-2)x_range[0]=0;
    x_range[1]=floor(x_range[1]/0.1)*0.1+0.1;
    if(y_range[0]<1)y_range[0]=0;
    y_range[1] = floor(y_range[1]/50)*50+50;

    if(addLast){
        auxTable.solutionTime = table.data[i-1][0];
        auxTable.T=to_string(table.data[i-1][0]); auxTable.i = i; auxTable.F="Point"; 
        WriteTecplotNormalData(auxTable, fileName,  10, secLine, std::ios::app); 
         nZones++;
    }
    printf("Suggestion range: %g %g %g %g\n",x_range[0],x_range[1],y_range[0],y_range[1]);
    printf("nZones: %d\n",nZones);
    return;
}

