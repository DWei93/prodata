
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
        int colPlaneID = GetColIDFromTable(table, "planeID");

        printf("size %d %d: \n", (int)table.data.size(), (int)table.variables.size());
        
        table.variables.push_back("bx");
        table.variables.push_back("by");
        table.variables.push_back("bz");

        table.variables.push_back("nx");
        table.variables.push_back("ny");
        table.variables.push_back("nz");

        for(i=0; i<table.data.size(); i++){
            real8 &burgID = table.data[i][colBurgID];

            int intBurgID = (int)burgID;
            switch(intBurgID){
                case 0:
                    vec[0] = 0.0; vec[1] = 1.0; vec[2] = 1.0;
                    break;
                case 1:    
                    vec[0] = 0.0; vec[1] = 1.0; vec[2] = -1.0;
                    break;
                case 2:
                    vec[0] = -1.0; vec[1] = 0.0; vec[2] = 1.0;
                    break;
                case 3:
                    vec[0] = 1.0; vec[1] = 0.0; vec[2] = 1.0;
                    break;
                case 4:
                    vec[0] = -1.0; vec[1] = 1.0; vec[2] = 0.0;
                    break;
                case 5:
                    vec[0] = 1.0; vec[1] = 1.0; vec[2] = 0.0;
                    break;
                default:
                    vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0;
                    break;
            }

            NormalizeVec(vec);
            table.data[i].push_back(vec[0]);
            table.data[i].push_back(vec[1]);
            table.data[i].push_back(vec[2]);
            

            real8 &planeID = table.data[i][colPlaneID];
            int   intPlaneID = (int)planeID;
            switch(intPlaneID){
                case -2:
                    vec[0] = 1.0; vec[1] = 0.0; vec[2] = 0.0;
                    break;
                case -1:
                    vec[0] = 1.0; vec[1] = 1.0; vec[2] = 0.0;
                    break;
                case 0:
                    vec[0] = 1.0; vec[1] = 1.0; vec[2] = 1.0;
                    break;
                case 1:
                    vec[0] = 1.0; vec[1] = 1.0; vec[2] = -1.0;
                    break;
                case 2:
                    vec[0] = -1.0; vec[1] = 1.0; vec[2] = 1.0;
                    break;
                case 3:
                    vec[0] = 1.0; vec[1] = -1.0; vec[2] = 1.0;
                    break;
                default:
                    vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0;
                    break;
            }
            NormalizeVec(vec);
            table.data[i].push_back(vec[0]);
            table.data[i].push_back(vec[1]);
            table.data[i].push_back(vec[2]);

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
void SpecifyEquations(Table_t &table)
{
        int     i, j;

        SpecifyEquations1(table);
    
        return;
}

bool Analysis(Table_t &table, real8 &sigma, real8 &hard, real8 &twindef){
    int last=table.data.size()-1; real8 strain;
    for(int i=0; i<last; i++){
        if(table.data[i][3] >0.1){
            if(table.data[last][3]-table.data[i][3]<0.05)return false;
            sigma = table.data[i][2];
            strain = table.data[i][1];
            hard = (table.data[last][2]-sigma)/(table.data[last][1]-strain);
            twindef = 1.0E2*fabs(table.data[i][4]/table.data[i][3]);
            return true;
        }
    }
    return false;
}
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

