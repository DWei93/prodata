
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
    int  i;
    
    for(i=0; i<table.variables.size(); i++){
        if(name == table.variables[i])break;

    } 
    return(i);
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

void FoldBox(real8 boundMin[3], real8 boundMax[3], real8 *x, real8 *y, real8 *z)
{
        real8   xc, yc, zc;

        real8   invLx, invLy, invLz, Lx, Ly, Lz;

        xc = (boundMin[0] + boundMax[0]) * 0.5;
        yc = (boundMin[1] + boundMax[1]) * 0.5;
        zc = (boundMin[2] + boundMax[2]) * 0.5;

        Lx = boundMax[0] - boundMin[0];
        invLx = 1.0/Lx;
        Ly = boundMax[1] - boundMin[1];
        invLy = 1.0/Ly;
        Lz = boundMax[2] - boundMin[2];
        invLz = 1.0/Lz;
    
        *x -= rint((*x-xc)*invLx) * Lx;

        *y -= rint((*y-yc)*invLy) * Ly;

        *z -= rint((*z-zc)*invLz) * Lz;
    
        return;
}

void ZImage(real8 boundMin[3], real8 boundMax[3], real8 *x, real8 *y, real8 *z)
{
        real8   xc, yc, zc;

        real8   invLx, invLy, invLz, Lx, Ly, Lz;

        xc = (boundMin[0] + boundMax[0]) * 0.5;
        yc = (boundMin[1] + boundMax[1]) * 0.5;
        zc = (boundMin[2] + boundMax[2]) * 0.5;

        Lx = boundMax[0] - boundMin[0];
        invLx = 1.0/Lx;
        Ly = boundMax[1] - boundMin[1];
        invLy = 1.0/Ly;
        Lz = boundMax[2] - boundMin[2];
        invLz = 1.0/Lz;
/*
 *      If periodic boundaries are not in use, the provided position
 *      of (x,y,z) will not be adjusted since there are no other
 *      images available.
 */
        *x -= rint(*x * invLx) * Lx;

        *y -= rint(*y * invLy) * Ly;

        *z -= rint(*z * invLz) * Lz;
    
        return;
}

real8 LinearInterpolation(const Curve_t &curve, real8 x, real8 min, real8 max)
{
    real8   length, xc;
    real8   x0, y0, x1, y1, t, boundVal;
    int     i;    

    if(curve.ax.size()<2)
        Fatal("there is no enough data for interpolation (%d)", (int)curve.ax.size());
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


    vector<string>().swap(table.variables);
    vector<vector<double> >().swap(table.data);

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

    return;            
}


void SpecifyEquations_PLTDATA(InArgs_t *inArgs)
{
    int     i, index;
    bool    readState;
    string  secLine, noBackupName = ("nobackup"), backupFile, sufBack = ("bak.");
    char    name[256];
    bool    backup = 1;

    vector<Table_t> tables;

    if((index = GetValID(inArgs->priVars, noBackupName)) < inArgs->priVars.size()){
        backup = 0;
    }
    if(backup){
        printf("Backup: Yes.\n");
    }else{
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
            backupFile = sufBack + inArgs->inpFiles[i];
            WriteTecplotNormalData(tables[i], backupFile, 10, secLine); 
        }
        SpecifyEquations(tables[i]);
    
        WriteTecplotNormalData(tables[i], inArgs->inpFiles[i], 10, secLine); 
    }

    return;

}

void FormatVector(real8 vec[3], const char *msg){
	printf("%s vector:\n",msg);
	printf("{%.15f,%.15f,%.15f}\n", vec[0], vec[1], vec[2]);
}
