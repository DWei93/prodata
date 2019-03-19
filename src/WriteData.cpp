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

    if(secLine != ""){
        out << secLine << endl;
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


int WriteMGDataFile(const string &file, MgData_t &mg, int precision) 
{
    int     i, j;
    string  fn(file), plt(".mg");

    if(fn.length() > 3){
        int loc = fn.find(plt, fn.length()-4); 
        if(loc == string::npos){
            fn += ".mg";
        }
    }else{
        fn += ".mg";
    }

    ofstream out;
    out.open(fn.c_str(), ios::out);
    
    out << "ITEM: TIMESTEP" << endl;
    out << mg.timestep << endl;

    out << "ITEM: NUMBER OF ATOMS" << endl;
    out << mg.atom.size() << endl;
    
    out << "ITEM: BOX BOUNDS";
    for(i=0; i<mg.bounds.size(); i++){
        out << " " << mg.bounds[i];
    }
    out << endl;

    for(i=0; i<mg.box.size(); i++){
        for(j=0; j<mg.box[i].size(); j++){
            out <<  setprecision(precision) << mg.box[i][j] << " ";
        }
        out << endl;
    }

    out << "ITEM: ATOMS id type x y z";
    for(i=0; i<mg.variables.size(); i++){
        out << " " << mg.variables[i];
    }
    out << endl;

    for(i=0; i<mg.atom.size(); i++){
        out << mg.atom[i].id << " ";
        out << mg.atom[i].type << " ";
        out << setprecision(precision) << mg.atom[i].x << " ";
        out << setprecision(precision) << mg.atom[i].y << " ";
        out << setprecision(precision) << mg.atom[i].z << " ";
        
        for(j=0; j<mg.atom[i].vars.size(); j++){
            out << setprecision(precision) << mg.atom[i].vars[j] << " "; 
        }
        out << endl;
    }

    out.close();
    return(0);
}



int MGToLMPDataFile(const string &file, MgData_t &mg, int precision)
{
    int     i, j, nTypes = 0;
    string  fn(file), plt(".lmp");

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

    out << mg.atom.size() << " atoms" << endl << endl;

    for(i=0; i<mg.atom.size(); i++){
        nTypes = MAX(nTypes, mg.atom[i].type);
    }
    out << nTypes << " atom types" << endl << endl;

    out << setprecision(precision) << mg.box[0][0] << " " << mg.box[0][1] << " xlo xhi" << endl;
    out << setprecision(precision) << mg.box[1][0] << " " << mg.box[1][1] << " ylo yhi" << endl;
    out << setprecision(precision) << mg.box[2][0] << " " << mg.box[2][1] << " zlo zhi" << endl;

    out << endl << "Atoms" << endl << endl;

    for(i=0; i<mg.atom.size(); i++){
        out << mg.atom[i].id << " "<< mg.atom[i].type << " ";
        out << setprecision(precision) << mg.atom[i].x << " ";
        out << setprecision(precision) << mg.atom[i].y << " ";
        out << setprecision(precision) << mg.atom[i].z << endl;
    }

    out.close();
    return 0;
}
