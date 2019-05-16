
#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"


void SpecifyEquations1(Table_t &table)
{
        int     i;
        real8   tau;
 
        vector<vector<double> >::iterator   it;
        int colSeparation =  GetColIDFromTable(table, "separation");

        for(it = table.data.begin(); it != table.data.end(); ){

            if((*it)[colSeparation] < 25){
                table.data.erase(it);
                continue;
            }

            it++;
        }

        return;
}

void SpecifyEquations(Table_t &table)
{
        int     i, j;

        SpecifyEquations1(table);
    
        return;
}
