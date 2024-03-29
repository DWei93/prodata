
#include <strings.h>
#include <regex>
#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"

#include "DDD.h"
#include "MD.h"
#include "Math.h"

#define     MAXLINELENGTH 1024

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

    InitTable(table);
    secLine = "";

//  std::regex aux_equation("AUXDATA\\s+\\S+\\s*=\\s*[-\"\'+.a-zA-Z0-9]+");
    std::regex aux_equation("AUXDATA\\s+\\S+\\s*=");
    std::regex equation("\\S+\\s*=\\s*[-\"\'+.a-zA-Z0-9]+");
    std::regex var_name("[a-zA-Z]+[a-zA-z0-9_\\-]*");
    std::regex value("[a-zA-z0-9_\\.\\+\\-]+");
    std::smatch equation_match;
    std::smatch v_match;

    std::string buff0, buff1, buff2, buff3;

    Variable_t  auxVar;

    bool    firstPoint = 1;
    int     currSize, rowNow=0;
    char    *p2;
    int     line=0;
    while(1){
        rt = fgets(str, MAXLINELENGTH, fp);
        if(line==0 && str[0] == '\0')return 0;
        line++;
        if(line==2)secLine=str;
        if(feof(fp))break;
        if(str == "\n" || str == NULL)continue;

        if(firstPoint){
            if(strstri(str, "variables") != NULL){
                token = strtok(str, "=");
                while(1){
                    if((token = strtok(NULL, delim)) != NULL){
                        table.variables.push_back(token);
                    }else{
                        break;
                    }
//                    table.variables.push_back(strtok(NULL, delim));
                }
                continue;
            }

            buff0 = str;
            if(std::regex_search(buff0, equation_match, aux_equation)){
                buff1 = equation_match[0].str();
                buff3 = equation_match.suffix().str();
//              printf("buf1=%s; buff3=%s\n",buff1.c_str(),buff3.c_str());
                if(std::regex_search(buff1, v_match, var_name)){
                    buff2 = v_match.suffix().str();
                    if(std::regex_search(buff2, v_match, var_name)){
                        buff1 = v_match[0].str();
                        if(std::regex_search(buff3, v_match, value)){
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
                            }else if (buff2 == "F"){
                                table.F = v_match[0].str(); 
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
                token = strtok(str, " \n\t");
                if(token != NULL){
                    snprintf(name, sizeof(name), "v%d", i); 
                    table.variables.push_back(name);
                    
                    while((token = strtok(NULL, " \n\t")) != NULL){
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
                table.data[i][0] = atof(strtok(str, " \t")); 
                for(j=1; j<table.variables.size(); j++){
                    table.data[i][j] = atof(strtok(NULL, " \n\t"));
                }
                rt = fgets(str, MAXLINELENGTH, fp);
                rowNow++;
            }
        }else{
            currSize = (int)table.data.size();
            if(currSize<=rowNow){
                currSize += 100;
                table.data.resize(currSize);
                for(int row=0;row<100;row++)table.data[rowNow+row].resize(table.variables.size());
            }
            p = strtok(str, " \t");
            if(p == nullptr)Fatal("Read error in %s at row %d",file.c_str(),currSize);
            table.data[rowNow][0] = atof(p); 
            for(j=1; j<table.variables.size(); j++){
                p = strtok(NULL, " \t");
                if(p == nullptr)Fatal("Read error in %s at row %d",file.c_str(),rowNow+1);
                table.data[rowNow][j] = atof(p);
                if(!Numerical(table.data[rowNow][j])){
                    if(isnan(table.data[rowNow][j]))table.data[rowNow][j]=0.0; 
                    else 
                    Fatal("Read error in %s at row %d, not numerical",file.c_str(),rowNow+1);
                }
            }
            rowNow++;
        }
    }
    fclose(fp);

    table.data.resize(rowNow);
    if((int)table.data.size() != table.i*table.j*table.k)table.i = (int)table.data.size();
    printf("Finish reading input file %s, %d points\n", file.c_str(), (int)table.data.size());    
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
    InitList(list);
   
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
    list.i = int(list.data[0].size());
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

            bool readVelocity=false;
            rt = fgets(str, MAXLINELENGTH, fp);
            rt=strtok(str, delim);if(rt!=(char*)NULL)dum.atom[0].id=atoi(rt);else Fatal("can not read atom id");
            rt=strtok(NULL, delim);if(rt!=(char*)NULL)dum.atom[0].type=atoi(rt);else Fatal("can not read atom type");
            rt=strtok(NULL, delim);if(rt!=(char*)NULL)dum.atom[0].x=atof(rt);else Fatal("can not read atom coordinate x");
            rt=strtok(NULL, delim);if(rt!=(char*)NULL)dum.atom[0].y=atof(rt);else Fatal("can not read atom coordinate y");
            rt=strtok(NULL, delim);if(rt!=(char*)NULL)dum.atom[0].z=atof(rt);else Fatal("can not read atom coordinate z");
            
            rt=strtok(NULL, delim);if(rt!=(char*)NULL){
                dum.atom[0].vx = atof(rt);
                rt=strtok(NULL, delim);if(rt!=(char*)NULL)dum.atom[0].vy=atof(rt);else Fatal("can not read atom velocity y");
                rt=strtok(NULL, delim);if(rt!=(char*)NULL)dum.atom[0].vz=atof(rt);else Fatal("can not read atom velocity z");
            }else{dum.atom[0].vx = 0.0;dum.atom[0].vy = 0.0;dum.atom[0].vz = 0.0;}

            for(int i=1; i<dum.atom.size(); i++){
                rt = fgets(str, MAXLINELENGTH, fp);
                dum.atom[i].id = atoi(strtok(str, delim));
                dum.atom[i].type = atoi(strtok(NULL, delim));
                dum.atom[i].x = atof(strtok(NULL, delim));
                dum.atom[i].y = atof(strtok(NULL, delim));
                dum.atom[i].z = atof(strtok(NULL, delim));
                dum.atom[i].vx = 0.0;dum.atom[i].vy = 0.0;dum.atom[i].vz = 0.0;
            }
        }
    }

    fclose(fp);
    return 1;
}

int ModifyDumpFile(const string &file, Dump_t &dum)
{
    double tot=double(dum.atom.size()), i=0;
    for (auto it=dum.atom.begin(); it!=dum.atom.end(); ){
        i++;
        if((*it).type==1) {
            it=dum.atom.erase(it); 
        } else {
            ++it;
        }
        
        if(int(i*10000/tot)%100==0)printf("%f\n",i*100/tot);
    }
    return 0;    
}
