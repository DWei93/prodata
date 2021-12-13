
#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"


void SpecifyEquations1(Table_t &table)
{
        double vy = (table.data.back()[3] - table.data[0][3])/
                    (table.data.back()[0] - table.data[0][0]);
        char val[50];
        snprintf(val, sizeof(val), "%e", vy);
        map<string, string>::iterator   iter;
        if((iter=table.aux.find(std::string("velocity")))!=table.aux.end()){
            iter->second = val;
        }else Fatal("something wrong");

        return;
}


void SpecifyEquations2(Table_t &table)
{
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
 
//            changed = ((sigmaS1 > sigmaS2) ? 1 : 0);
            changed = 0;
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
        return;
}


void SpecifyEquations3(Table_t &table)
{
        int     i, j;
        real8   vec[3];

        int colBurgID = GetColIDFromTable(table, "burgID");
        int colIndex = GetColIDFromTable(table, "index");

        for(i=0; i<table.data.size(); i++){
            real8 &burgID = table.data[i][colBurgID];
            real8 &index = table.data[i][colIndex];

            if(index>3.5){
                index++;
            }else if (index>2.5){
                if(burgID > 6.5)index++;
            }
        }

        return;
}

bool cmp1(vector<double> p, vector<double> q){

    return(p[5]>q[5]);

    if(p[0] == q[0]){
        return(p[3]<q[3]);
    }else{
        return(p[0]>q[0]);
    }
}

void SpecifyEquations4(Table_t &table)
{
        int     i;
        bool    active;

        sort(table.data.begin(), table.data.end(), cmp1);
        printf("maxvtau %e\n", table.data[0][5]);

        for(i=table.data.size()-1; i>=0; i--){
            if(table.data[i][5] != 0){
                printf("min vtau %e\n", table.data[i][5]);
                break;
            }
        }
        return;

        real8   thread_vsqua_a = 0.0, thread_vtau_a = 0.0;
        real8   thread_vsqua_ia = 0.0, thread_vtau_ia = 0.0, na = 0, nia = 0;
#pragma omp parallel for reduction(+:thread_vsqua_a) reduction(+:thread_vsqua_ia) reduction(+:thread_vtau_a) reduction(+:thread_vtau_ia) reduction(+:na) reduction(+:nia)
        for(i=0; i<table.data.size(); i++){
            if(table.data[i][0] == 1){
                na += 1.0;
                thread_vsqua_a += table.data[i][3]*table.data[i][3];
                thread_vtau_a += table.data[i][5];   
            }else{
                nia += 1.0;
                thread_vsqua_ia += table.data[i][3]*table.data[i][3];
                thread_vtau_ia += table.data[i][5];   
            }
        }

        printf("%f %f %f %f %f %f \n", thread_vsqua_a/na,
                                 thread_vsqua_ia/nia,
                                 thread_vtau_a/na,
                                 thread_vtau_ia/nia, na, nia);

        return;
}

void SpecifyEquations5(Table_t &table)
{
        int     i;
        real8   t0, lastAbsy;
        real8   ybox = 1.5368600500205644e+03*3.615E-10;
        real8   v1, v2;
        
        int colTimestep =  GetColIDFromTable(table, "timestep");
        int colY =  GetColIDFromTable(table, "y");

        table.variables.push_back("absy");
        int colAbsy =  GetColIDFromTable(table, "absy");

        t0 = table.data[0][colTimestep];

        real8 x0 = -1, y0 = 0;
        for(i=0; i<table.data.size(); i++){

            real8  &timestep = table.data[i][colTimestep];
            real8  &y = table.data[i][colY];

            timestep = (timestep-t0)*1E-15;

            y *= 3.615E-10;

            v1 = (i==0) ? ybox/2 : table.data[i-1][colY];
            lastAbsy = (i==0) ? ybox/2 : table.data[i-1][colAbsy];
            v2 = y - v1;
            v2 -= rint(v2 * 1.0/ybox) * ybox;

            table.data[i].push_back(lastAbsy+v2);
            if(table.data[i][colTimestep] > 2.0E-11 && x0 == -1){
                x0 = table.data[i][colTimestep];
                y0 = lastAbsy+v2; 
            }
        }

        printf("Slop %e \n", 
               (table.data.back()[colAbsy] - y0)/(table.data.back()[colTimestep] - x0));
        return;
}

void SpecifyEquations6(Table_t &table)
{
    double pla_step = 20000;

    int col_step = GetColIDFromTable(table, "step");
    int col_time = GetColIDFromTable(table, "cumulativeTime");
    int col_vsqua_act = GetColIDFromTable(table, "vsqua_active");
    int col_vsqua_inact = GetColIDFromTable(table, "vsqua_inactive");
    int col_diss_act = GetColIDFromTable(table, "tauv_active");
    int col_diss_inact = GetColIDFromTable(table, "tauv_inactive");

    real8 total_time=0.0, ave_vsqua_act=0.0, ave_vsqua_inact=0.0, ave_diss_act=0.0, ave_diss_inact=0.0;

#pragma omp parallel for reduction(+:total_time, ave_vsqua_act, ave_vsqua_inact, ave_diss_act, ave_diss_inact)
    for(int i=0; i<table.data.size(); i++){
        if(table.data[i][col_step] <= pla_step)continue;
        total_time += table.data[i][col_time];
        ave_vsqua_act += table.data[i][col_vsqua_act]*table.data[i][col_time];
        ave_vsqua_inact += table.data[i][col_vsqua_inact]*table.data[i][col_time];
        ave_diss_act += table.data[i][col_diss_act]*table.data[i][col_time];
        ave_diss_inact += table.data[i][col_diss_inact]*table.data[i][col_time];
    }
    ave_vsqua_act /=  total_time;
    ave_vsqua_inact /=  total_time;
    ave_diss_act /=  total_time;
    ave_diss_inact /=  total_time;

    char value[1024];

    snprintf(value, sizeof(value), "%f", total_time);
    table.aux["total_time"] = value;

    snprintf(value, sizeof(value), "%f", ave_vsqua_act);
    table.aux["ave_vsqua_act"] = value;

    snprintf(value, sizeof(value), "%f", ave_vsqua_inact);
    table.aux["ave_vsqua_ainct"] = value;

    snprintf(value, sizeof(value), "%f", ave_diss_act);
    table.aux["ave_diss_act"] = value;

    snprintf(value, sizeof(value), "%f", ave_diss_inact);
    table.aux["ave_diss_ainct"] = value;

    printf("the total time of the plastic stage %f ns.\n", total_time);
    printf("<v^2>_active=%e, <v^2>_ianctive=%e, <v\\tau>_active=%e, <v\\tau>_inactive=%e\n", ave_vsqua_act, ave_vsqua_inact, ave_diss_act, ave_diss_inact);
    return;
}



void SpecifyEquations8(Table_t &table)
{

    table.variables.push_back("p_neg_active");
    table.variables.push_back("p_pos_active");

    real8   total_neg = 0.0, total_pos = 0.0;
#pragma omp parallel for reduction(+: total_neg, total_pos)
    for(int i=0; i<200; i++){
        table.data[i].resize(5);
        if(i<100){
            table.data[i][3] = 0.0;
            table.data[i][4] = 0.0;
        }else{
            table.data[i][3] = table.data[100-(i-100)][2];
            table.data[i][4] = table.data[i][2];
            total_neg += (table.data[i][1]*table.data[100-(i-100)][2]);
            total_pos += (table.data[i][1]*table.data[i][2]);
        }
    }

#pragma omp parallel for 
    for(int i=100; i<200; i++){
        table.data[i][3] /= total_neg;
        table.data[i][4] /= total_pos;
    }

    return;
}

void SpecifyEquations9(Table_t &table)
{

    Table_t newTable;
    std::vector<double> vec(3);

    newTable.variables.push_back("x");
    newTable.variables.push_back("y");
    newTable.variables.push_back("tau");
    for(int i=0; i<table.data.size(); i++){
        if(fabs(table.data[i][1]) <300){
            vec[0] = table.data[i][0]; vec[1]=table.data[i][1]; vec[2]=table.data[i][2];
            newTable.data.push_back(vec);
        }
    }

    swap(newTable, table);
    table.i = (int)table.data.size();
    table.j = 1; table.k=1; table.F="Point";

    return;
}
void SpecifyEquations10(Table_t &table)
{
    table.variables.push_back("summation");

    real8 volume=pow(4000,3),burgMag=2.556E-10;
    for(auto &i:table.data){
        i[3] *= 266/volume; i[5] *=266/volume; i[7] *= 266/volume;
        i[1] *= ((i[5]+i[7])*1.0E6)/burgMag;
        i[2] *= i[3]*1.0E6/burgMag;
        i[4] *= i[5]*1.0E6/burgMag;  
        i[6] *= i[7]*1.0E6/burgMag; 
    }
    for(auto i=0; i<table.data.size(); i++){
        table.data[i].push_back(table.data[i][1]+table.data[i][2]+table.data[i][4]+table.data[i][6]);
    }

    int nID;
    Table_t newTable; InitTable(newTable);
    for(auto &i : table.variables)newTable.variables.push_back(i);
    newTable.variables.push_back("numPoints"); nID = (int)newTable.variables.size()-1;
    for(auto &i : table.variables){
        char name[50]; snprintf(name, sizeof(name),"tau_%s", i.c_str());
        newTable.variables.push_back(("tau_"+i));
    }
    for(auto i=0; i<table.data.size(); i++){
        int j=0; for(j=0; j<newTable.data.size(); j++){
            printf("%d %d: %e %e %e\n",(int)newTable.data.size(),j,table.data[i][0],newTable.data[j][0], fabs(newTable.data[j][0]-table.data[i][0]) );
            if(fabs(newTable.data[j][0]-table.data[i][0]) < 100)break;
        }
        if(j==newTable.data.size()){
            newTable.data.resize(j+1); newTable.data[j].resize(newTable.variables.size());
            for(auto k=0; k<newTable.data[j].size(); k++)newTable.data[j][k]=0;
            for(auto k=0; k<table.data[i].size(); k++){
               newTable.data[j][k]=table.data[i][k];
                
            }
            newTable.data[j][nID] = 1;
        }else{
            newTable.data[j][nID]++;
            for(auto k=1; k<table.data[i].size(); k++)newTable.data[j][k]+=table.data[i][k];
        }
    }

    for(auto i=0; i<newTable.data.size(); i++)for(auto j=1; j<table.data[i].size(); j++)newTable.data[i][j]/=newTable.data[i][nID];

    for(auto i=0; i<table.data.size(); i++){
        int j=0; for(j=0; j<nID; j++){
            if(fabs(newTable.data[j][0]-table.data[i][0]) < 100)break;
        }
        if(j==newTable.data.size()){
            Fatal("someting wrong");
        }else{
            newTable.data[j][nID];
            for(auto k=0; k<table.data[i].size(); k++)newTable.data[j][k+nID+1] += pow((table.data[i][k]-newTable.data[j][k]),2);
        }
    }

    for(auto i=0; i<newTable.data.size(); i++)for(auto j=nID+2; j<newTable.variables.size(); j++)
    {
        newTable.data[i][j] = sqrt(newTable.data[i][j]/newTable.data[i][nID]);
    }

    newTable.i = (int)newTable.data.size();
    swap(newTable, table);
    return;
}

void SpecifyEquations11(Table_t &table)
{
    int cIndex = GetColIDFromTable(table, "index"); // units %
    int cBurgID = GetColIDFromTable(table, "burgID"); // units %
    int cPlaneID = GetColIDFromTable(table, "normalID"); // units %
    int cInGrain = GetColIDFromTable(table, "inGrain"); // units %
    int cOnGB = GetColIDFromTable(table, "onGB"); // units %
    int ctau_g = GetColIDFromTable(table, "tau_g"); // MPa
    int cifMag = GetColIDFromTable(table, "ifMag"); // MPa
    int cconstraint=GetColIDFromTable(table, "constraint");
    double  tau_g=0.0, ntau_g=0.0, ifMag=0.0, dtau_g=0;
    for (auto i=0; i<table.data.size(); i++){
        if((table.data[i][cPlaneID]<3.9 || table.data[i][cPlaneID]>4.1) && table.data[i][cconstraint] <1){tau_g += table.data[i][ctau_g]; ntau_g += 1.0; ifMag+=table.data[i][cifMag];}
    }
    tau_g/=ntau_g;
    for (auto i=0; i<table.data.size(); i++){
        if((table.data[i][cPlaneID]<3.9 || table.data[i][cPlaneID]>4.1) && table.data[i][cconstraint] <1){dtau_g += pow(table.data[i][ctau_g]-tau_g,2); }
    }
    printf("DATA:%f %f %f\n",tau_g, ifMag/ntau_g, sqrt(dtau_g/(ntau_g-1.0)));
    return;
}

void SpecifyEquations13(Table_t &table)
{
    int cIndex = GetColIDFromTable(table, "index"); // units %
    int cBurgID = GetColIDFromTable(table, "burgID"); // units %
    int cPlaneID = GetColIDFromTable(table, "normalID"); // units %
    int cInGrain = GetColIDFromTable(table, "inGrain"); // units %
    int cOnGB = GetColIDFromTable(table, "onGB"); // units %
    int ctau_g = GetColIDFromTable(table, "tau_g"); // MPa
    int cifMag = GetColIDFromTable(table, "ifMag"); // MPa
    int cconstraint=GetColIDFromTable(table, "constraint");
    for (auto i=0; i<table.data.size(); ){
        if (!(table.data[i][cBurgID] > 4.5 && table.data[i][cBurgID] < 5.5 && table.data[i][cInGrain] < 0.5 && table.data[i][cInGrain] > -0.5)){
            table.data.erase(table.data.begin()+i, table.data.begin()+i+1);            
            printf("erase %d\n",i);
        }else{i++;printf("++ %d\n",i);}
    }
    printf("i=%d\n",table.i);
    table.i = table.data.size();
    printf("i=%d\n",table.i);
}
void SpecifyEquations12(Table_t &table)
{
    int cX = GetColIDFromTable(table, "Z"); // units %
    int cY = GetColIDFromTable(table, "Y"); // units %
    int cZ = GetColIDFromTable(table, "Z"); // units %
    int cIndex = GetColIDFromTable(table, "index"); // units %
    int cConstraint = GetColIDFromTable(table, "constraint"); // units %
    int cBurgID = GetColIDFromTable(table, "burgID"); // units %
    double p[3]={0,0,0};
    std::vector<int> ids;
    int constraint=-1;

    for (auto i=0; i<table.data.size(); i++){
        double &x=table.data[i][cX];
        double &y=table.data[i][cY];
        double &z=table.data[i][cZ];
        if((x-p[0])*(x-p[0])+(y-p[1])*(y-p[1])+(z-p[2])*(z-p[2]) > 1.0 || i==table.data.size()-1){
            for(auto j=0;j<ids.size();j++)table.data[ids[j]][cConstraint]=constraint;
            ids.resize(1); ids[0]=i; constraint=-1;
        }else{
            ids.push_back(i);
        }
        if(fabs(table.data[i][cIndex]-5)<0.1 && table.data[i][cConstraint]>3.5)table.data[i][cConstraint]=6;
        constraint = (table.data[i][cConstraint]>constraint)?table.data[i][cConstraint]:constraint;
        if(table.data[i][cIndex]>4){
            if(table.data[i][cBurgID]<6) table.data[i][cIndex] = 5;
            else if(table.data[i][cBurgID]>17 && table.data[i][cBurgID]<21)table.data[i][cIndex] = 6;
            else table.data[i][cIndex] = 7;
        }
    }
    return;
}

void SpecifyEquations14(Table_t &table)
{
    int cep = GetColIDFromTable(table, "plaStn"); // units %
    int cz = GetColIDFromTable(table, "Z"); // units %

    int nxy = table.i*table.j;
    printf("I=%d, J=%d, K=%d, nxy=%d\n",table.i, table.j, table.k, nxy);


    Table_t auxTable;
    auxTable.variables.resize(2);
    auxTable.variables[0] = string("Z");
    auxTable.variables[1] = string("ep");

    auxTable.data.resize(table.k);
    auxTable.i=table.k;
    auxTable.j=1; auxTable.k=1; auxTable.F="point";
    auxTable.solutionTime = table.solutionTime;
    auxTable.T =  table.T;

    for(auto i=0; i<auxTable.data.size(); i++){
        auxTable.data[i].resize(2);
        auxTable.data[i][0] = table.data[i*nxy][cz];
        auxTable.data[i][1] = 0.0;
        for(auto j=0;j<nxy;j++){
            auxTable.data[i][1] += table.data[i*nxy+j][cep];
        }
    }
#if 1
    auxTable.data[0][1]  =(auxTable.data[1][1] +auxTable.data[0][1]+auxTable.data[40][1] )/3.0;
    auxTable.data[41][1] =(auxTable.data[1][1] +auxTable.data[0][1]+auxTable.data[40][1] )/3.0;
    auxTable.data[17][1] =(auxTable.data[16][1]+auxTable.data[17][1]+auxTable.data[18][1])/3.0;
    auxTable.data[24][1] =(auxTable.data[23][1]+auxTable.data[24][1]+auxTable.data[25][1])/3.0;

    auxTable.data.erase(auxTable.data.begin()+40);
    auxTable.data.erase(auxTable.data.begin()+25);
    auxTable.data.erase(auxTable.data.begin()+23);
    auxTable.data.erase(auxTable.data.begin()+18);
    auxTable.data.erase(auxTable.data.begin()+16);
    auxTable.data.erase(auxTable.data.begin()+1);

#else

    auxTable.data[17][1] = (auxTable.data[16][1]+auxTable.data[17][1]+auxTable.data[18][1])/3.0;
    auxTable.data.erase(auxTable.data.begin()+18);
    auxTable.data.erase(auxTable.data.begin()+16);

    auxTable.data[5][1] = auxTable.data[4][1]+auxTable.data[5][1]+auxTable.data[6][1];
    auxTable.data[36][1] = auxTable.data[35][1]+auxTable.data[36][1]+auxTable.data[37][1];

    auxTable.data.erase(auxTable.data.begin()+37);
    auxTable.data.erase(auxTable.data.begin()+35);
    auxTable.data.erase(auxTable.data.begin()+6);
    auxTable.data.erase(auxTable.data.begin()+4);
#endif

    auxTable.i = int(auxTable.data.size());
    string secLine, fileName, baseName("ep-"), pltName(".plt");
    auto a = table.aux.find("id");
    if(a!=table.aux.end()){
        fileName = baseName+ a->second +pltName;
        WriteTecplotNormalData(auxTable, fileName,  10, secLine, std::ios::out); 
    }
    return;
}

void SpecifyEquations15(Table_t &table)
{
        int     i, j;
        real8   vec[3];

        int colStress = GetColIDFromTable(table, "stress");

        for(i=0; i<table.data.size(); i++){
            table.data[i][colStress] *= 0.8791208791;
        }

        return;
}

void SpecifyEquations16(Table_t &table)
{
        int     i, j;
        real8   vec[3];

        int colmd_fs = GetColIDFromTable(table, "md_fs");
        int colsfi = GetColIDFromTable(table, "sfi");
        
        double md=0.0, msf=0.408, bmag=0.2556;

        for(i=0; i<table.data.size(); i++){
            if(fabs(table.data[i][colsfi] - msf)<0.01){
                if(md < table.data[i][colmd_fs]*0.2556){
                    md =  table.data[i][colmd_fs]*bmag;
                }
            }
        }
        printf("DATA: %.2f\n",md);

        return;
}

void SpecifyEquations17(Table_t &table)
{
        int     i, j;
        real8   vec[3];

        int cPst = GetColIDFromTable(table, "plastic_strain");
        int cStress = GetColIDFromTable(table, "stress");

        double minActStress = 1E5;
        for(i=0; i<table.data.size(); i++){
            if(table.data[i][cPst] > 0.05){
                if(minActStress  > table.data[i][cStress]){
                    minActStress = table.data[i][cStress];
                }
            }
        }
        printf("DATA: %.2f\n",minActStress);

        return;
}

void SpecifyEquations18(Table_t &table)
{
        int     i, j;
        real8   vec[3];

        int ctime = GetColIDFromTable(table, "timenow");
        int cseparation = GetColIDFromTable(table, "separation");
        int cfile = GetColIDFromTable(table, "file");
        int cp = GetColIDFromTable(table, "p");

        int cDseparation = GetColIDFromTable(table, "D_separation");
        int cDfile = GetColIDFromTable(table, "D_file");
        int cDp = GetColIDFromTable(table, "D_p");

        int cstr_zy = GetColIDFromTable(table, "str_zy");

        Table_t auxTable;
        InitTable(auxTable);
        auxTable.variables.push_back("timenow_ns");
        auxTable.variables.push_back("appStress_MPa");
        auxTable.variables.push_back("displacement_b");
        auxTable.variables.push_back("stdev_displacement_b");
        
        auxTable.variables.push_back("separation_b");
        auxTable.variables.push_back("stdev_separation_b");

        vector<double> p(auxTable.variables.size());
        for(i=0; i<table.data.size(); i++){
            p[0]=table.data[i][ctime];
            p[1]=table.data[i][cstr_zy]/1E6;
            p[2]=table.data[i][cp];
            p[3]=table.data[i][cDp];
            p[4]=table.data[i][cseparation];
            p[5]=table.data[i][cDseparation];
            auxTable.data.push_back(p);
        }
        auxTable.i = auxTable.data.size();

        string secLine("");
        WriteTecplotNormalData(auxTable, "loading-unloading-simplied.plt",  10, secLine, std::ios::out); 

        return;
}

void SpecifyEquations19(Table_t &table)
{
        int     i, j;
        real8   vec[3];

        int ctime = GetColIDFromTable(table, "timenow");
        int cseparation = GetColIDFromTable(table, "separation");
        int cfile = GetColIDFromTable(table, "file");
        int cp = GetColIDFromTable(table, "p");

        int cDseparation = GetColIDFromTable(table, "D_separation");
        int cDfile = GetColIDFromTable(table, "D_file");
        int cDp = GetColIDFromTable(table, "D_p");

        int cstr_zy = GetColIDFromTable(table, "str_zy");

        Table_t auxTable;
        InitTable(auxTable);


        auxTable.variables.push_back("appliedStress_MPa");
        
        auxTable.variables.push_back("separation_b");
        auxTable.variables.push_back("stdev_separation_b");

        vector<double> p(auxTable.variables.size());
        for(i=0; i<table.data.size(); i++){
            p[0]=table.data[i][cfile]*10;
            p[1]=table.data[i][cseparation];
            p[2]=table.data[i][cDseparation];
            auxTable.data.push_back(p);
        }
        auxTable.i = auxTable.data.size();

        string secLine("");
        WriteTecplotNormalData(auxTable, "simplified.plt",  10, secLine, std::ios::out); 

        return;
}
#define MAX(a,b) ((a)>(b)?(a):(b))
void SpecifyEquations20(Table_t &table)
{
        int     i, j;
        int cPst = GetColIDFromTable(table, "plaStn");
        int cdisDen = GetColIDFromTable(table, "disDen");
        int cX = GetColIDFromTable(table, "X");
        int cY = GetColIDFromTable(table, "Y");
        int cZ = GetColIDFromTable(table, "Z");

        double totPstn=0.0, box[3]={0,0,0}, totDisDen=0;
        for(i=0; i<table.data.size(); i++){
            totPstn += table.data[i][cPst];
            totDisDen += table.data[i][cdisDen];
            box[0] = MAX(2.0*table.data[i][cX], box[0]);
            box[1] = MAX(2.0*table.data[i][cY], box[1]);
            box[2] = MAX(2.0*table.data[i][cZ], box[2]);
        }
        double bMag=2.556E-10, vol=box[0]*box[1]*box[2];

        printf("DATA: %g %g m^-2\n", totPstn/box[0]/box[1]/box[2]*100.0, totDisDen/vol/bMag/bMag);

        return;
}
void SpecifyEquations(Table_t &table)
{
        int     i, j;

        SpecifyEquations20(table);
//      ClearRepetition(table);
    
        return;
}

#ifdef GSL
#include <gsl/gsl_fit.h>
bool Analysis(int pointID, double sf, Table_t &table, real8 &sigma, real8 &hard, real8 &thard, real8 &twindef, real8 &crss, bool &aveHard){
    int cStrain = GetColIDFromTable(table, "strain"); // units %
    int cStress = GetColIDFromTable(table, "stress"); // units MPa
    int cPst = GetColIDFromTable(table, "plastic_strain");
    int cTpst = GetColIDFromTable(table, "twinning_deformation"); 
    int cDen = GetColIDFromTable(table, "disDen");
    int cFit = GetColIDFromTable(table, "fit_stress");
    
    static const double schmid[28]={
                                2.7216553e-01,
                                4.0409609e-01,
                                4.4647817e-01,
                                5.0000000e-01,
                                4.3301270e-01,
                                4.5306962e-01,
                                4.0824829e-01,
                                4.2060056e-01,
                                4.8629953e-01,
                                4.9279928e-01,
                                4.9240388e-01,
                                4.9240388e-01,
                                5.0000000e-01,
                                4.5810976e-01,
                                4.7384710e-01,
                                4.7303654e-01,
                                4.5465119e-01,
                                4.8571785e-01,
                                4.9934573e-01,
                                4.9317787e-01,
                                4.9400820e-01,
                                4.8766921e-01,
                                4.4297535e-01,
                                4.6424283e-01,
                                4.7140452e-01,
                                4.6424283e-01,
                                4.4297535e-01,
                                4.0824829e-01};


    double strain, stress, pst, tpst, den;
    int last=0;
    int i;
    for(i=table.data.size()-1; i>0; i--){
        if(table.data[i][cStrain] < 0.5){last=i; break;};
    }
    last = table.data.size()-1;
    double l_strain=table.data[last][cStrain];
    double l_stress=table.data[last][cStress];
    double l_pst=table.data[last][cPst];
    double l_tpst=(cTpst>-1)?table.data[last][cTpst]:0;
    double l_den=table.data[last][cDen];
    bool fit = (cFit<0)?true:false;
    printf("cFit = %d\n",cFit);

    double critPst= 0.199, youngs0 = 123.275, youngs,critStn=0.201;
    int  critI,elaI=-1;

#if 0
    critPst=0.002;
    critStn=0.1;
#else
    if(l_pst<0.05){
        printf("warning: can not find yield stress, last (%d) pst=%.2f, cPst=%d \n", last, l_pst,cPst);
        sigma = 0.0;
        return false;
    }
#endif
    aveHard = (l_pst > 2.0*critStn)?true:false;
    for(i=0; i<last+1; i++){
        pst=table.data[i][cPst];
        if(pst>1.0E-3 && elaI==-1){elaI=i;}
        if(pst > critPst){
            sigma = table.data[i][cStress];
            critStn=table.data[i][cStrain];
            critI = i;
            break;
        }
    }
    twindef = 100.0*table.data[i][cTpst]/table.data[i][cPst];
    if(last+1==i){
        printf("Warning: plastic strain is small, but you can use it to define yield stress\n");
        sigma = l_stress;
        crss = sf*sigma;
        return true;
    }else{
        printf("i=%d\n sigma=%.2f MPa\n",i, sigma);
    }

    for(i=0; i<last; i++){
        pst=table.data[i][cPst];
        if(pst > critPst){
            break;
        }
    }
    
    int nPoints = last-i+1;
    double *stnList, *strList, *tpstList, *eStnList, *eStrList; 
    stnList = (double *)malloc(nPoints*sizeof(double));
    strList = (double *)malloc(nPoints*sizeof(double));
    tpstList = (double *)malloc(nPoints*sizeof(double));
    eStrList = (double *)malloc((elaI+1)*sizeof(double));
    eStnList = (double *)malloc((elaI+1)*sizeof(double));

    int j=0;
    double c0, c1, cov00, cov01, cov11, chisq, c00, c11;
    for(j=0; j<elaI+1; j++){
        eStnList[j] = table.data[j][cStrain];
        eStrList[j] = table.data[j][cStress];
    } 
    gsl_fit_mul(eStnList,1,eStrList,1,elaI+1,&c1,&cov11,&chisq);
    youngs = c1/10.0; 

    free(eStrList);
    free(eStnList);

    j=0;
    for(/*skip*/; i<last+1; i++){
        stnList[j] = table.data[i][cStrain];
        strList[j] = table.data[i][cStress];
        tpstList[j] = (cTpst>-1)?-youngs0*table.data[i][cTpst]:0.0;
        j++;
    }

    gsl_fit_linear(stnList,1,strList,1,nPoints,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
    hard = 0.1*c1; c00=c0; c11=c1;

//  if(l_pst > 0.5){
//      sigma = c0+c1*critStn;
//  }
    
    if(cTpst>-1){
        gsl_fit_linear(stnList,1,tpstList,1,nPoints,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
        thard = 0.1*c1;
    }else thard = 0.0;

    if(sf<0){
        pointID -=1;
        if(pointID>-1){
            crss = sigma*schmid[pointID];
            printf("point %d: %7.4e=%7.4e*%7.4e\n",pointID+1, crss, sigma, schmid[pointID]);
        }else crss = 0.0;
    }else{
        crss = sf*sigma;
        printf("crss=%f*%f=%f\n",sf,sigma,crss);
    }
    
#if 0
    printf ("# best fit: Y = %g + %g X\n", c0, c1);
    printf ("# covariance matrix:\n");
    printf ("# [ %g, %g\n#   %g, %g]\n",
          cov00, cov01, cov01, cov11);
    printf ("# chisq = %g\n", chisq);
#endif

    free(stnList);
    free(strList);
    free(tpstList);

    double inp = c00/(10.0*youngs-c11);
    if(fit){
        printf("Plot Fitting curve ... intersection point (%5.2f,%5.2f)\n",inp, 10.0*youngs*inp);
        table.variables.push_back("fit_stress");
        for(i=0;i<table.data.size();i++){
            strain = table.data[i][cStrain];
            if(strain<inp){
                table.data[i].push_back(youngs*10.0*strain);
            }else{
                table.data[i].push_back(c11*strain+c00);
            }
        }   
    }

    return true;
}
#else
bool Analysis(Table_t &table, real8 &crss){
    crss = -1;
    real8 max=0;
    for(int i=0; i<table.data.size(); i++){
        if(table.data[i][10]>max){
            crss = table.data[i][0];
            max = table.data[i][10];
        } 
    }
    return true;
}
#endif



void InitialStructure(InArgs_t *inArgs){
    typedef enum{
        cfrlen = 0,
        cdpb,
        cdfs1,
        cdtb1,
        cdtb2,
        csi,
        cso,
        cstd1,
        cstd2,
        careai,
        careao,
        cbx,
        cby,
        cbz,
        cnx,
        cny,
        cnz,
        ccbetafr,
        cstb,
        cdfs2,
        ccbetafs1,
        ccbetafs2,
        ctbdir,
        ccbetatb1,
        ccbetatb2,
        ctauPin1,
        ctauPin2,
        cmathematica,
        ctau0,
        ctauc,
        csigmac,
        cpeneType,
        crc,
        cnPileUps,
        cpileLength,
        ctauBack,
        csigmay,
        clinter,
        clpo,
        cmax
    }InitaiVals;

    if(inArgs->inpFiles.size() == 0){
        Fatal("There is no input file.");
    }
    vector<Table_t> tables(inArgs->inpFiles.size());
    int i;
    Table_t auxTable;
    auxTable.variables.push_back("index");
    auxTable.variables.push_back("tauc");
    auxTable.variables.push_back("frlen");
    auxTable.variables.push_back("dfs");
    auxTable.variables.push_back("peneType");
    auxTable.variables.push_back("nPileUps");
    auxTable.variables.push_back("pileLength");
    auxTable.variables.push_back("tauBack");
    auxTable.variables.push_back("nActiveSrcs");
    auxTable.variables.push_back("tauPin");
    
    vector<double> dPoint(auxTable.variables.size());
    vector<double> ave(dPoint.size()), sigma(dPoint.size());
    double num=0.0, sfMax=-1, bMag=0.2556;
    int index;

    if((index = GetValID(inArgs->priVars, "sf")) < inArgs->priVars.size()){
        sfMax = atof(inArgs->priVars[index].vals[0].c_str());
    }
    for(i=0;i<dPoint.size();i++){
        ave[i]=0.0; sigma[i]=0.0;
    }

//#pragma omp parallel for 
    for(i=0; i<inArgs->inpFiles.size(); i++){
        string secLine;
        bool readState;
//        #pragma omp critical
        {
            readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        }
        
        if(readState){
            vector<vector<double> > &data = tables[i].data; 
            dPoint[0] = double(i);
            if(sfMax<=0){
                for(int j=0; j<data.size(); j++){
                    sfMax = sfMax<data[j][csi]?data[j][csi]:sfMax;
                }
            }
            dPoint[1] = data[0][ctauc];
            dPoint[2] = data[0][cfrlen]*bMag; 
            dPoint[3] = 0;
            dPoint[4] = data[0][cpeneType];
            dPoint[5] = data[0][cnPileUps];
            dPoint[6] = data[0][cpileLength]*bMag;
            dPoint[7] = data[0][ctauBack];
            dPoint[8] = 0;
            
            // Find the most easily activated colinear reactions
            dPoint[9] = (data[0][ctauPin1]<data[0][ctauPin2])?data[0][ctauPin1]:data[0][ctauPin2];
#if 1
            // If src is tructated, the do not consider colinear reaction
            if(data[0][cmathematica]>0)dPoint[9] = dPoint[1];
            if(data[0][ctauc] - data[0][ctau0] >0.1){dPoint[9]  = dPoint[1];}
#endif
            // If pinning stress is too small or higher than crss, skip it
            if(dPoint[9] < 0.1 || dPoint[9]>dPoint[1])dPoint[9] = dPoint[1];

            for(int j=0; j<data.size(); j++){
                if(fabs(data[j][csi] -sfMax) < 0.02) {
                    dPoint[8]++;
                    if(data[j][cdfs1] > dPoint[3])dPoint[3]=data[j][cdfs1];
                    if(data[j][cdfs2] > dPoint[3])dPoint[3]=data[j][cdfs2];
                }
            }
            dPoint[3] *= bMag;
            auxTable.data.push_back(dPoint);
            for(int j=0;j<dPoint.size();j++){
                ave[j]+=dPoint[j];
                printf("%g ",dPoint[j]);
            }
            printf("\n");
            num++;
        }
    }
    auxTable.i = auxTable.data.size();
   
    for(int j=0;j<dPoint.size();j++)ave[j]/=num;
    for(i=0;i<auxTable.data.size();i++){
        for(int j=0; j<dPoint.size(); j++){
            sigma[j] += pow(auxTable.data[i][j]-ave[j], 2);
        }
    }
    for(int j=0; j<dPoint.size(); j++){
        if(num>1.1){
            sigma[j] = sqrt(sigma[j]/(num-1));
        }else{
            sigma[j] = 0.0;
        }
    }
    printf("DAT:     number tau_c(MPa), frlen(nm), dfs(nm), peneType nPileUps plieLength(nm) tauBack(MPa) nActiveSrcs tauPin\n");
    printf("AVE: ");
    for(int j=0; j<dPoint.size(); j++)printf(" %9.2f", ave[j]);
    printf("\nVAR: ");
    for(int j=0; j<dPoint.size(); j++)printf(" %9.2f", sigma[j]);
    printf("\n");

    return;
}

void GNDAnalysis(InArgs_t *inArgs){
    int i, j, k;
    string secLine(""), outfile;
    bool readState = false;
    std::size_t iPos;
    std::vector<Table_t> tables(inArgs->inpFiles.size());

    if(inArgs->inpFiles.size()>0){
        readState = ReadTecplotNormalData(inArgs->inpFiles[0], tables[0], secLine);
    }
    if(readState==false){
        Fatal("can not find first gnd files.");
    }

    for(i=1; i<inArgs->inpFiles.size(); i++){
        readState = ReadTecplotNormalData(inArgs->inpFiles[i], tables[i], secLine);
        if(!readState)continue;

        if(tables[i].variables.size()!= tables[0].variables.size()
           || tables[i].variables.size()<=2 
           || tables[i].data[0].size() != tables[0].data[0].size()){
            Fatal("File %s: varible sizes is wrong", inArgs->inpFiles[i].c_str());
        }

        for(j=0; j<tables[i].data.size(); j++){
            for(k=2; k<tables[i].data[j].size(); k++){
                tables[i].data[j][k] -= tables[0].data[j][k];
            }
        }

        iPos = inArgs->inpFiles[i].find_last_of("/\\");
        outfile = "delta-" + inArgs->inpFiles[i].substr(iPos+1);
        printf("Output delta gnd file %s:\n", outfile.c_str());
        WriteTecplotNormalData(tables[i], outfile, 10, secLine, std::ios::out);
    }
    return;
}

void CustomHandleTecplotData(InArgs_t *inArgs){
    int     index;
    string  functionName("func");

    if((index = GetValID(inArgs->priVars, functionName)) < inArgs->priVars.size()){
        if(inArgs->priVars[index].vals[0].find("InitialStructure") != string::npos){
            InitialStructure(inArgs);
            return;
        }else if(inArgs->priVars[index].vals[0].find("GNDAnalysis") != string::npos){
            GNDAnalysis(inArgs);
        }else{
            Fatal("Can not support Analysis function %s", inArgs->priVars[index].vals[0].c_str());
        }
    }

    return;
}

