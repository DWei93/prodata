
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>

#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"
#include "Parse.h"

using namespace std;

vector<string> GetFiles(const string &file) 
{
    int     iPos;
    DIR     *dir;
    char    basePath[256], *getcwdReturn;

    struct dirent   *ptr;
    string          f(file), fdir, fname, dname, bpath, curDir("./");
    vector<string>  strs;
    
    memset(basePath, '\0', sizeof(basePath));
    getcwdReturn = getcwd(basePath, sizeof(basePath));
//    bpath = basePath;

    iPos = f.find_last_of('/');
    if(iPos == string::npos){
        fdir = curDir;
        fname = f;
    }else{
        fdir = curDir + f.substr(0,iPos);
        fname = f.substr(iPos+1);
    }

    
    if((dir = opendir(fdir.c_str())) == NULL){
        Fatal("Can not open dir %s", fdir.c_str());
    }

    while((ptr = readdir(dir)) != NULL){
        if(strcmp(ptr->d_name, ".") == 0 || 
           strcmp(ptr->d_name, "..") == 0){     //current dir OR paraent dir
          continue;
        }else if(ptr->d_type == 8){   //file
            dname = ptr->d_name;
            int loc = dname.find(fname, 0); 
            if(loc != string::npos){
                strs.push_back(fdir + "/" + dname);
            }
        }else if(ptr->d_type == 10){    // link file

        }else if(ptr->d_type == 4){     //  dir

        }
    }

    return(strs);
}











