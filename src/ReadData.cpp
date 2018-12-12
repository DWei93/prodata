
#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"

using namespace std;

int ReadTecplotNormalData(string &file, Table_t &table)
{
    int         lineNo = 1, i, j;
    string      str;
    ifstream    infile;

    vector<string>  line;
    vector<double>  col;
    infile.open(file.c_str());

    if(!infile){
        Fatal("cant not the file %s", file.c_str());
    }

    vector<string>().swap(table.variables);
    vector<vector<real8> >().swap(table.data);
/*
 *  Read the first line.
 */
    getline(infile, str);
    if(str.empty())Fatal("there is noting in the file %s", file.c_str());
    line = split(str, "=");
    table.variables = split(line[1], ",");
    for(i=0; i<table.variables.size(); i++)WashString(table.variables[i]);
    col.resize(table.variables.size());

    j = 0;
    while(getline(infile,str))
    {   
        lineNo++;
        line = split(str, " ");
       
        if(line.size() != table.variables.size())continue; 
        for(i=0; i<col.size(); i++){
            col[i] = atof(line[i].c_str());
        }
        table.data.push_back(col);
        j++;
    }
    infile.close();
    printf("Finsish reading input file %s\n", file.c_str());    
    return 1;            
}
