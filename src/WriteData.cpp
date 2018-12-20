
#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"
#include "MD.h"
#include "DDD.h"

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
