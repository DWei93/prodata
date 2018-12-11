#include <algorithm>
#include <functional>

#include "Home.h"
#include "Util.h"
#include "ReadData.h"

using namespace std;

void HandleExtendedDislocation(InArgs_t *inArgs)
{
    Table_t         table;
    real8           boundMin[3], boundMax[3], size = 10000;
    int             index;
    string          bound = "bound";  
    

    ReadTecplotNormalData(inArgs->inpFiles[0], table);

    if((index = GetValID(inArgs->priVars, bound)) < inArgs->priVars.size()){
        size = stof(inArgs->priVars[index].vals[0]);
        printf("size is %f\n", size);
    }
    

    
    
    return;
}
