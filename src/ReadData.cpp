
#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"

#include "DDD.h"
#include "MD.h"

using namespace std;

int ReadTecplotNormalData(string &file, Table_t &table, string &secLine)
{
    int         lineNo = 1, i, j, type;
    string      str, v("v");
    ifstream    infile;

    vector<string>  line;
    vector<double>  col;
    infile.open(file.c_str());

    if(!infile){
        printf("Warning: cant not the file %s", file.c_str());
        return (0);
    }

    vector<string>().swap(table.variables);
    vector<vector<real8> >().swap(table.data);
/*
 *  Read the first line.
 */
    getline(infile, str);
    if(str.empty()){
        printf("Warning: there is noting in the file %s", file.c_str());
        return (0);
    }

    line = split(str, " ");
    for(i=0; i<line.size(); i++){
        if(atof(line[i].c_str()) == 0)break;
    }
    vector<string>().swap(line);

    if(i == 0){
        line = split(str, "=");
        table.variables = split(line[1], ",");
        for(i=0; i<table.variables.size(); i++){
            WashString(table.variables[i]);
        }
        vector<string>().swap(line);
        col.resize(table.variables.size());

    }else{
        line = split(str, " ");
        table.variables.resize(line.size());
        col.resize(table.variables.size());
        for(i=0;i<line.size();i++){
            table.variables[i] = v + (char)i;
            col[i] = atof(line[i].c_str());
        }
        table.data.push_back(col);
    }


    j = 0;
    while(getline(infile,str))
    {   
        lineNo++;
        line = split(str, " ");
        if(lineNo == 2)secLine = str;
       
        if(line.size() != table.variables.size()){
            continue; 
        }
        for(i=0; i<col.size(); i++){
            col[i] = atof(line[i].c_str());
        }
        table.data.push_back(col);
        j++;
        vector<string>().swap(line);
    }
    infile.close();

//    printf("Finsish reading input file %s\n", file.c_str());    
    return 1;            
}


int ReadMGDataFile(const string &file, MgData_t &mgdata)
{
    int         i, j, atoms;
    int         idCol = 0, typeCol = 0;
    int         xCol = 0, yCol = 0, zCol = 0;
    ifstream    infile;
    string      str;

    CleanMgData(mgdata);

    infile.open(file.c_str());

    vector<string>  words, subwords;
    vector<int>     varCol;

    if(!infile){
        printf("Warning: cant not the file %s", file.c_str());
        return (0);
    }

    while(getline(infile,str))
    {
        vector<string>().swap(words);
        words = split(str, " ");
        if(words.size() == 0)continue;

        if(words[0] == "ITEM:"){
            if(words[1] == "TIMESTEP"){
                getline(infile,str);
                mgdata.timestep = atoi(str.c_str());
            }else if(words[1] == "NUMBER"){
                if(words.size()==4){
                    if(words[2] == "OF" && words[3] == "ATOMS"){
                        getline(infile,str);
                        atoms = atoi(str.c_str());
                    }
                }
            }else if(words[1] == "BOX"){
                if(words.size()==6){
                    if(words[2] == "BOUNDS"){
                        mgdata.bounds.resize(3);
                        mgdata.bounds[0] = words[3];
                        mgdata.bounds[1] = words[4];
                        mgdata.bounds[2] = words[5];
                        for(i=0; i<3; i++){
                            vector<string>().swap(subwords);
                            getline(infile,str);
                            subwords = split(str, " ");
                            if(subwords.size() != 2){
                                Fatal("in file %s, can not read %s", file.c_str(), str.c_str());
                            }
                            mgdata.box[i*2] = atof(subwords[0].c_str());
                            mgdata.box[1+i*2] = atof(subwords[1].c_str());
                        }
                    }
                } 
             }else if(words[1] == "ATOMS"){
                if(words.size()<=7){
                    printf("Waring: in file %s, can not read %s", file.c_str(), str.c_str());
                    return(0);
                }

                varCol.resize(words.size()-7);
                mgdata.variables.resize(words.size()-7);
                j = 0;
                for(i=2; i<words.size(); i++){
                    if(words[i] == "x"){xCol = i-2;        continue;}
                    if(words[i] == "y"){yCol = i-2;        continue;}
                    if(words[i] == "z"){zCol = i-2;        continue;}
                    if(words[i] == "id"){idCol = i-2;      continue;}
                    if(words[i] == "type"){typeCol = i-2;  continue;}
                    mgdata.variables[j] = words[i];
                    varCol[j] = i-2;
                    j++;
                } 

                mgdata.atom.resize(atoms);
                for(i=0; i<mgdata.atom.size(); i++){
                    vector<string>().swap(subwords);
                    getline(infile,str);
                    subwords = split(str, " ");
                    if(subwords.size() != mgdata.variables.size()+5){
                        printf("Warning: in file %s, can not read %s", file.c_str(), str.c_str());
                        return (0);
                    }

                    mgdata.atom[i].x = atof(subwords[xCol].c_str());
                    mgdata.atom[i].y = atof(subwords[yCol].c_str());
                    mgdata.atom[i].z = atof(subwords[zCol].c_str());
                    mgdata.atom[i].id = atoi(subwords[idCol].c_str());
                    mgdata.atom[i].type = atoi(subwords[typeCol].c_str());

                    mgdata.atom[i].vars.resize(mgdata.variables.size());
                    for(j=0; j<mgdata.variables.size(); j++){
                        mgdata.atom[i].vars[j] = atof(subwords[varCol[j]].c_str());
                    }
                }
            }else{
                printf("Warning: not support the ITEM %s", words[1].c_str());
                return (0);
            }
        } /*End of if */
        
    }/* End of while(getline(file, str))*/
    
    return(1);
}


int ReadDataFromMDLogFile(const vector<string> &files, LineList_t &list)
{
    int         i, j, firstTime;
    ifstream    infile;
    string      str;

    vector<string>  words, subwords;
    SwapLineList(list);
   
    firstTime = 1;
    for(i=0; i<files.size(); i++){
        infile.open(files[i].c_str());
        
        while(getline(infile,str)){
            vector<string>().swap(words);
            words = split(str, " ");
            if(words.size() < 1)continue;
            if(words[0] != "Step")continue;
            if(firstTime){
                list.variables.assign(words.begin(), words.end());
                list.data.resize(list.variables.size());
                firstTime = 0;
            }
            while(getline(infile, str)){
                vector<string>().swap(subwords);
                subwords = split(str, " ");
                if(subwords.size() != list.variables.size())break;
                for(j=0; j<list.variables.size(); j++){
                    list.data[j].push_back(atof(subwords[j].c_str()));
                }
            }
        }
        infile.close();
    }
#if 0
    printf("List: ");
    for(j=0; j<list.variables.size(); j++){
        printf("%s ", list.variables[j].c_str());
    }
    for(i=0; i<list.data[0].size(); i++){
        for(j=0; j<list.variables.size(); j++){
            printf("%f ", list.data[j][i]);
        }
        printf("\n");
    }
#endif
    if(list.data.size()==0){
        return(0);
    }
    return (1);
}


