
#include <strings.h>
#include <regex>
#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"

#include "DDD.h"
#include "MD.h"

#define     MAXLINELENGTH 512

using namespace std;
const char* strstri(const char* str, const char* subStr)
{
    int len = strlen(subStr);
    if(len == 0)
    {
        return NULL;
    }
    while(*str)
    {
        if(strncasecmp(str, subStr, len) == 0)
        {
            return str;
        }
        ++str;
    }
    return NULL;
}

int ReadTecplotNormalData(string &file, Table_t &table, string &secLine)
{
    int  i, numPoints = 1, j;
    char str[MAXLINELENGTH], delim[] = " \t\",\n", *p, *token;
    char quotation[]="\"\n";
    char name[256], *rt;
    FILE *fp;
    if((fp = fopen(file.c_str(), "r+")) == NULL){
        printf("Warning: can not open fie %s\n", file.c_str());
        return(0);
    }else{
//        printf("Reading tecplot file %s ...\n", file.c_str());
    }

    secLine = "";
    vector<string>().swap(table.variables);
    vector<vector<real8> >().swap(table.data);
    table.aux.clear();
    table.i=1;
    table.j=1;
    table.k=1;
    table.solutionTime = -1;
    table.T="";

    std::regex equation("\\S+\\s*=\\s*[-\"\'+.a-zA-Z0-9]+");
    std::regex aux_equation("AUXDATA\\s+\\S+\\s*=\\s*[-\"\'+.a-zA-Z0-9]+");
    std::regex var_name("[a-zA-Z]+[a-zA-z0-9_\\-]*");
    std::regex value("[a-zA-z0-9_\\.\\+\\-]+");

    std::smatch equation_match;
    std::smatch v_match;

    std::string buff0, buff1, buff2;

    Variable_t  auxVar;

    bool    firstPoint = 1;
    int     currSize;
    char    *p2;
    while(1){
        rt = fgets(str, MAXLINELENGTH, fp);
        if(feof(fp))break;
        if(str == "\n" || str == NULL)continue;

        if(firstPoint){
            if(strstr(str, "variables") != NULL){
                token = strtok(str, "=");
                while(1){
                    if((token = strtok(NULL, delim)) != NULL){
                        table.variables.push_back(token);
                    }else{
                        break;
                    }
//                    table.variables.push_back(strtok(NULL, delim));
                }
//                printf("variables: ");
//                for(i=0;i<table.variables.size();i++)printf("%s ", table.variables[i].c_str());
//                printf("\n");
                continue;
            }

            buff0 = str;
            if(std::regex_search(buff0, equation_match, aux_equation)){
                buff1 = equation_match[0].str();
                if(std::regex_search(buff1, v_match, var_name)){
                    buff2 = v_match.suffix().str();
                    if(std::regex_search(buff2, v_match, var_name)){
                        buff1 = v_match[0].str();
                        buff2 = v_match.suffix().str();
                        if(std::regex_search(buff2, v_match, value)){
                            table.aux[buff1] = v_match[0].str();
                        }
                    }
                }
                continue;
            }else if(std::regex_search(buff0,equation_match, equation)){
                 buff1 = buff0;
                 while(std::regex_search(buff1,equation_match, equation)){
                    buff2 = equation_match[0].str();
                    if(std::regex_search(buff2,v_match, var_name)){
                        buff2 = v_match[0].str();
                        buff1 = v_match.suffix().str();
                        if(std::regex_search(buff1, v_match, value)){
                            if(buff2 == "T"){
                                table.T = v_match[0].str(); 
                            }else if (buff2 == "I" || buff2 == "i"){
                                table.i = atoi(v_match[0].str().c_str());
                            }else if (buff2 == "J" || buff2 == "j"){
                                table.j = atoi(v_match[0].str().c_str());
                            }else if (buff2 == "K" || buff2 == "k"){
                                table.k = atoi(v_match[0].str().c_str());
                            }else if (buff2 == "SOLUTIONTIME" || buff2 == "solutiontime"){
                                table.solutionTime = atof(v_match[0].str().c_str());
                            }else{
                                printf("Warning: can not identify %s", buff2.c_str());
                            }
                        }
                    }
                    buff1 = equation_match.suffix().str();
                }
                continue;
            }

            firstPoint = 0;
            if(table.variables.size() == 0){
                char strBak[MAXLINELENGTH];
                sprintf(strBak, "%s", str);
                i = 0;
                token = strtok(str, " \n");
                if(token != NULL){
                    snprintf(name, sizeof(name), "v%d", i); 
                    table.variables.push_back(name);
                    
                    while((token = strtok(NULL, " \n")) != NULL){
                        i++;
                        snprintf(name, sizeof(name), "v%d", i); 
                        table.variables.push_back(name);
                    } 
                    printf("Undifined variables: ");
                    for(i=0;i<table.variables.size();i++)printf("%s ", table.variables[i].c_str());
                    printf("\n");
                }
                sprintf(str, "%s", strBak);
            }
            numPoints = table.i*table.j*table.k;
        }


//        if(strtok(str, " ") == NULL)continue;
        if(numPoints > 1){
            currSize = (int)table.data.size();
            table.data.resize(currSize+numPoints);
            for(i=currSize; i<table.data.size() && !feof(fp); i++){
                table.data[i].resize(table.variables.size());
                table.data[i][0] = atof(strtok(str, " ")); 
                for(j=1; j<table.variables.size(); j++){
                    table.data[i][j] = atof(strtok(NULL, " \n"));
                }
                rt = fgets(str, MAXLINELENGTH, fp);
            }
        }else{
            currSize = (int)table.data.size();
            table.data.resize(currSize+numPoints);
            table.data[currSize].resize(table.variables.size());
            table.data[currSize][0] = atof(strtok(str, " ")); 
            for(j=1; j<table.variables.size(); j++){
                table.data[currSize][j] = atof(strtok(NULL, " \n"));
            }
//            printf("%d: ", currSize);
//            for(i=0;i<table.data[currSize].size(); i++)printf("%f ", table.data[currSize][i]);
//            printf("\n");
        }
    }
    fclose(fp);

    if((int)table.data.size() != table.i*table.j*table.k)table.i = (int)table.data.size();
//    printf("Finish reading input file %s, %d points\n", file.c_str(), currSize+1);    
    return 1;            
}


int ReadDumpFile(const string &file, Dump_t &dum)
{
    int         i, j, atoms;
    int         idCol = 0, typeCol = 0;
    int         xCol = 0, yCol = 0, zCol = 0;
    ifstream    infile;
    string      str;

    CleanDump(dum);

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
                dum.timestep = atoi(str.c_str());
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
                        dum.bounds.resize(3);
                        dum.bounds[0] = words[3];
                        dum.bounds[1] = words[4];
                        dum.bounds[2] = words[5];
                        
                        dum.box.resize(3);
                        for(i=0; i<3; i++){
                            vector<string>().swap(subwords);
                            getline(infile,str);
                            subwords = split(str, " ");
                            if(subwords.size() != 2){
                                Fatal("in file %s, can not read %s", file.c_str(), str.c_str());
                            }
                            dum.box[i].resize(2);
                            dum.box[i][0] = atof(subwords[0].c_str());
                            dum.box[i][1] = atof(subwords[1].c_str());
                        }
                    }
                }else if(words.size() == 9){
                    dum.bounds.resize(6);
                    for(i=0; i<6; i++){
                        dum.bounds[i] = words[3+i];
                    }
                    dum.box.resize(3);

                    for(i=0; i<3; i++){
                        vector<string>().swap(subwords);
                        getline(infile,str);
                        subwords = split(str, " ");
                        if(subwords.size() != 3){
                            Fatal("in file %s, can not read %s", file.c_str(), str.c_str());
                        }
                        dum.box[i].resize(3);
                        dum.box[i][0] = atof(subwords[0].c_str());
                        dum.box[i][1] = atof(subwords[1].c_str());
                        dum.box[i][2] = atof(subwords[2].c_str());

                    }
                } 
             }else if(words[1] == "ATOMS"){
                if(words.size()<7){
                    printf("Waring: in file %s, can not read %s - 1\n", file.c_str(), str.c_str());
                    return(0);
                }

                varCol.resize(words.size()-7);
                dum.variables.resize(words.size()-7);
                j = 0;
                for(i=2; i<words.size(); i++){
                    if(words[i] == "x"){xCol = i-2;        continue;}
                    if(words[i] == "y"){yCol = i-2;        continue;}
                    if(words[i] == "z"){zCol = i-2;        continue;}
                    if(words[i] == "id"){idCol = i-2;      continue;}
                    if(words[i] == "type"){typeCol = i-2;  continue;}
                    dum.variables[j] = words[i];
                    varCol[j] = i-2;
                    j++;
                } 

                dum.atom.resize(atoms);
                for(i=0; i<dum.atom.size(); i++){
                    vector<string>().swap(subwords);
                    getline(infile,str);
                    subwords = split(str, " ");
                    if(subwords.size() != dum.variables.size()+5){
                        printf("Warning: in file %s, can not read %s - 2\n", file.c_str(), str.c_str());
                        return (0);
                    }

                    dum.atom[i].x = atof(subwords[xCol].c_str());
                    dum.atom[i].y = atof(subwords[yCol].c_str());
                    dum.atom[i].z = atof(subwords[zCol].c_str());
                    dum.atom[i].id = atoi(subwords[idCol].c_str());
                    dum.atom[i].type = atoi(subwords[typeCol].c_str());

                    dum.atom[i].vars.resize(dum.variables.size());
                    for(j=0; j<dum.variables.size(); j++){
                        dum.atom[i].vars[j] = atof(subwords[varCol[j]].c_str());
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


bool ReadLMPFile(const string file, Dump_t &dum)
{
    CleanDump(dum);

    FILE *fp;
    if((fp = fopen(file.c_str(), "r+")) == NULL){
        return(0);
    }

    char str[MAXLINELENGTH], delim[] = " \t", *p, *rt;
    dum.box.resize(3);
    dum.bounds.resize(3);
    for(int i = 0; i<3; i++)dum.box[i].resize(2);

    while(!feof(fp)){
        rt = fgets(str, MAXLINELENGTH, fp);
        
        if((p = strstr(str, "atoms")) != NULL){
            dum.atom.resize(atoi(str));
            printf("atoms: %d\n",(int)dum.atom.size());
            continue;
        }

        if((p = strstr(str, "xlo xhi")) != NULL){
            dum.box[0][0] = atof(strtok(str, delim));
            dum.box[0][1] = atof(strtok(NULL, delim));
            printf("box x: [%f, %f]\n", dum.box[0][0], dum.box[0][1]);
            dum.bounds[0] = "xx";
            continue;
        }

        if((p = strstr(str, "ylo yhi")) != NULL){
            dum.box[1][0] = atof(strtok(str, delim));
            dum.box[1][1] = atof(strtok(NULL, delim));
            dum.bounds[1] = "yy";
            printf("box y: [%f, %f]\n", dum.box[1][0], dum.box[1][1]);
            continue;
        }

        if((p = strstr(str, "zlo zhi")) != NULL){
            dum.box[2][0] = atof(strtok(str, delim));
            dum.box[2][1] = atof(strtok(NULL, delim));
            dum.bounds[2] = "zz";
            printf("box z: [%f, %f]\n", dum.box[2][0], dum.box[2][1]);
            continue;
        }

        if((p = strstr(str, "xy xz yz")) != NULL){
            dum.box[0].resize(3);
            dum.box[0].resize(3);
            dum.box[0].resize(3);
            dum.bounds.resize(6);

            dum.box[0][2] = atof(strtok(str, delim));
            dum.box[1][2] = atof(strtok(NULL, delim));
            dum.box[2][2] = atof(strtok(NULL, delim));
            dum.bounds[3] = "xy";
            dum.bounds[4] = "xz";
            dum.bounds[5] = "yz";

            printf("triclinic box: [%f, %f %f]\n", dum.box[0][2], dum.box[1][2], dum.box[2][2]);
            continue;
        }

        if((p = strstr(str, "Atoms")) != NULL){
            rt = fgets(str, MAXLINELENGTH, fp);
            for(int i=0; i<dum.atom.size(); i++){
                rt = fgets(str, MAXLINELENGTH, fp);
                dum.atom[i].id = atoi(strtok(str, delim));
                dum.atom[i].type = atoi(strtok(NULL, delim));
                dum.atom[i].x = atof(strtok(NULL, delim));
                dum.atom[i].y = atof(strtok(NULL, delim));
                dum.atom[i].z = atof(strtok(NULL, delim));
            }
        }
    }

    fclose(fp);
    return 1;
}
