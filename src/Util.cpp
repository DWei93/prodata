
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


void CleanMgData(MgData_t &mg)
{
    int i;
    mg.timestep = 0;
    for(i=0; i<6; i++){
        mg.box[i] = 0.0;
    }

    vector<string>().swap(mg.bounds);
    vector<string>().swap(mg.variables);
    
    for(i=0; i<mg.atom.size(); i++){
        vector<double>().swap(mg.atom[i].vars);
    }
    vector<Atom_t>().swap(mg.atom);

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

void SpecifyEquations(Table_t &table)
{
#if 1
    int     i, j;

    int colBurgID =  GetColIDFromTable(table, "burgID");
    int colIndex =  GetColIDFromTable(table, "index");

    for(i=0; i<table.data.size(); i++){
        real8 &burgID = table.data[i][colBurgID];
        real8 &index = table.data[i][colIndex];

        if(burgID > 6.5){
            if(burgID == 7 || burgID == 12 || burgID == 17){
                burgID = -1;
            }else{
                burgID = 7;
            }
        }
    }

#else
    int     i, j;
    bool    changed;
    real8   tauAve, diff;
    
    vector<real8>   vec(3);
    table.variables.push_back("tau_averaged_single");
    table.variables.push_back("diff");

    int colSigmaS1 = GetColIDFromTable(table, "sigma_sin1");
    int colSS1 = GetColIDFromTable(table, "s_sin1");
    int colTauS1 = GetColIDFromTable(table, "tau_sin1");

    int colSigmaS2 = GetColIDFromTable(table, "sigma_sin2");
    int colSS2 = GetColIDFromTable(table, "s_sin2");
    int colTauS2 = GetColIDFromTable(table, "tau_sin2");

    int colTauBic = GetColIDFromTable(table, "tau_bic");

    for(i=0; i<table.data.size(); i++){

        real8 &sigmaS1 = table.data[i][colSigmaS1];
        real8 &sS1 = table.data[i][colSS1];
        real8 &tauS1 = table.data[i][colTauS1];

        real8 &sigmaS2 = table.data[i][colSigmaS2];
        real8 &sS2 = table.data[i][colSS2];
        real8 &tauS2 = table.data[i][colTauS2];
        real8 tauBic = table.data[i][colTauBic];

        changed = ((sigmaS1 > sigmaS2) ? 1 : 0);
        if(changed){
            vec[0] = sigmaS1;
            vec[1] = sS1;
            vec[2] = tauS1;

            sigmaS1 = sigmaS2;
            sS1 = sS2;
            tauS1 = tauS2;

            sigmaS2 = vec[0];
            sS2 = vec[1];
            tauS2 = vec[2];
        }

        if(sS1 == sS2){
            tauAve = 0.5*(tauS1 + tauS2);
        }else{
            if(sS1 > sS2){
                tauAve = tauS1;
            }else{
                tauAve = tauS2;
            }
        }

        diff = tauBic - tauAve; 

        table.data[i].push_back(tauAve);
        table.data[i].push_back(diff);
    }
#endif
    return;
}


void SpecifyEquations_PLTDATA(InArgs_t *inArgs)
{
    int     i;
    bool    readState;
    string  secLine;

    vector<Table_t> tables;

    if(inArgs->help)return;

    if(inArgs->inpFiles.size() == 0){
        Fatal("There is no input file.");
    }

    tables.resize(inArgs->inpFiles.size());
    for(i=0; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState)continue;

        SpecifyEquations(tables[i]);
    
        WriteTecplotNormalData(tables[i], inArgs->outFiles[i], 10, secLine); 
    }

    return;

}
