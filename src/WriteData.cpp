#include <ctime>

#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"
#include "MD.h"
#include "DDD.h"
#include "Math.h"

using namespace std;

void WriteTecplotNormalData(const LineList_t &list, const string &file, double precision) 
{
    int     i, j;
    string  fn(file), plt(".plt");
    if(fn.length() > 4){
        int loc = fn.find(plt, fn.length()-5); 
        if(loc == string::npos){
            fn += ".plt";
        }
    }else{
        fn += ".plt";
    }

    ofstream out;
    out.open(fn.c_str(), ios::out);

    out << "variables = "; 
    for(i=0; i<list.variables.size(); i++){
        if(i<list.variables.size()-1){
            out << list.variables[i] << ", ";
        }else{
            out << list.variables[i] << endl;
        }
    }

    if(((int)list.data[0].size()) != list.i*list.j*list.k){
        Fatal("somthing wrong with the list size %d != %dx%dx%d",
              (int)list.data[0].size(), list.i, list.j, list.k);
    }

    out << "zone i = " << list.i << ", j = " << list.j << ", k = " << list.k;
    if(list.T.length() != 0) out << ", T = \"" << list.T << "\"";
    if(list.solutionTime >=0) out << ", SOLUTIONTIME = " << list.solutionTime;
    out << ", F = " << list.F << "\n";

    if(!list.aux.empty()){
        for(const auto &pair : list.aux){
            out << "AUXDATA " << pair.first << " = \"" << pair.second << "\"\n";
        }
    }

    for(i=0; i<list.data[0].size(); i++){
        for(j=0; j<list.variables.size(); j++){
            out << setprecision(precision) <<list.data[j][i] << " ";
        }
        out << endl;
    }
     
    out.close();
//    printf("Finsish writing output file %s\n", fn.c_str());    
    return;
}

void WriteTecplotNormalData(const Table_t &table, const string &file, double precision, string secLine) 
{
    int     i, j;
    string  fn(file), plt(".plt");

    if(fn.length() > 4){
        int loc = fn.find(plt, fn.length()-5); 
        if(loc == string::npos){
            fn += ".plt";
        }
    }else{
        fn += ".plt";
    }

    ofstream out;
    out.open(fn.c_str(), ios::out);

    out << "variables = "; 
    for(i=0; i<table.variables.size(); i++){
        if(i<table.variables.size()-1){
            out << table.variables[i] << ", ";
        }else{
            out << table.variables[i] << endl;
        }
    }

    if(((int)table.data.size()) != table.i*table.j*table.k){
        Fatal("somthing wrong with the table size %d != %dx%dx%d",
              (int)table.data.size(), table.i, table.j, table.k);
    }
    out << "zone i = " << table.i << ", j = " << table.j << ", k = " << table.k;
    if(table.T.length() != 0) out << ", T = \"" << table.T << "\"";
    if(table.solutionTime >=0) out << ", SOLUTIONTIME = " << table.solutionTime;
    out << ", F = " << table.F << "\n";

    if(!table.aux.empty()){
        for(const auto &pair : table.aux){
            out << "AUXDATA " << pair.first << " = \"" << pair.second << "\"\n";
        }
    }

    for(i=0; i<table.data.size(); i++){
        for(j=0; j<table.variables.size(); j++){
            out << setprecision(precision) << table.data[i][j] << " ";
        }
        out << endl;
    }
     
    out.close();
//    printf("Finsish writing output file %s\n", fn.c_str());    
    return;
}


int WriteDumpFile(const string &file, Dump_t &dum, int precision) 
{
    int     i, j;
    string  fn(file), plt(".dum");

    if(fn.length() > 4){
        int loc = fn.find(plt, fn.length()-5); 
        if(loc == string::npos){
            fn += ".dum";
        }
    }else{
        fn += ".dum";
    }

    clock_t     timestart = clock();
    ofstream out;
    out.open(fn.c_str(), ios::out);
    
    out << "ITEM: TIMESTEP" << endl;
    out << dum.timestep << endl;

    out << "ITEM: NUMBER OF ATOMS" << endl;
    out << dum.atom.size() << endl;
    
    out << "ITEM: BOX BOUNDS";
    for(i=0; i<dum.bounds.size(); i++){
        out << " " << dum.bounds[i];
    }
    out << endl;

    for(i=0; i<dum.box.size(); i++){
        for(j=0; j<dum.box[i].size(); j++){
            out <<  setprecision(precision) << dum.box[i][j] << " ";
        }
        out << endl;
    }

    out << "ITEM: ATOMS id type x y z";
    for(i=0; i<dum.variables.size(); i++){
        out << " " << dum.variables[i];
    }
    out << endl;

    for(i=0; i<dum.atom.size(); i++){
        out << dum.atom[i].id << " ";
        out << dum.atom[i].type << " ";
        out << setprecision(precision) << dum.atom[i].x << " ";
        out << setprecision(precision) << dum.atom[i].y << " ";
        out << setprecision(precision) << dum.atom[i].z << " ";
        
        for(j=0; j<dum.atom[i].vars.size(); j++){
            out << setprecision(precision) << dum.atom[i].vars[j] << " "; 
        }
        out << endl;
    }

    out.close();
    printf("The output to file %s was completed (%fs)\n",file.c_str(), ((real8)(clock()-timestart))/CLOCKS_PER_SEC);
    return(0);
}



int MGToLMPDataFile(const string &file, Dump_t &dum, int precision)
{
    int     i, j, nTypes = 0;
    string  fn(file), plt(".lmp");
    clock_t     timestart = clock();

    if(fn.length() > 4){
        int loc = fn.find(plt, fn.length()-5); 
        if(loc == string::npos){
            fn += ".lmp";
        }
    }else{
        fn += ".lmp";
    }

    ofstream out;
    out.open(fn.c_str(), ios::out);

    time_t  now = time(0);
    char* dt = ctime(&now);

    out << "LMMPS: " << dt << endl;

    out << dum.atom.size() << " atoms" << endl << endl;

    for(i=0; i<dum.atom.size(); i++){
        nTypes = MAX(nTypes, dum.atom[i].type);
    }
    out << nTypes << " atom types" << endl << endl;

    out << setprecision(precision) << dum.box[0][0] << " " << dum.box[0][1] << " xlo xhi" << endl;
    out << setprecision(precision) << dum.box[1][0] << " " << dum.box[1][1] << " ylo yhi" << endl;
    out << setprecision(precision) << dum.box[2][0] << " " << dum.box[2][1] << " zlo zhi" << endl;

    out << endl << "Atoms" << endl << endl;

    for(i=0; i<dum.atom.size(); i++){
        out << dum.atom[i].id << " "<< dum.atom[i].type << " ";
        out << setprecision(precision) << dum.atom[i].x << " ";
        out << setprecision(precision) << dum.atom[i].y << " ";
        out << setprecision(precision) << dum.atom[i].z << " ";
        out << setprecision(precision) << dum.atom[i].vx << " ";
        out << setprecision(precision) << dum.atom[i].vy << " ";
        out << setprecision(precision) << dum.atom[i].vz << endl;
    }

    out.close();

    printf("The output to file %s was completed (%fs)\n",file.c_str(), ((real8)(clock()-timestart))/CLOCKS_PER_SEC);
    return 0;
}
