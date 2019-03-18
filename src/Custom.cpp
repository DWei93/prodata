#include "Home.h"
#include "Util.h"
#include "ProDataIO.h"

void SpecifyEquations1(Table_t &table)
{
        int     i, j;
 
        int colBurgID =  GetColIDFromTable(table, "burgID");
        int colIndex =  GetColIDFromTable(table, "index");
 
        for(i=0; i<table.data.size(); i++){
            real8 &burgID = table.data[i][colBurgID];
            real8 &index = table.data[i][colIndex];
 
            if(burgID > 6.5){
                if(burgID == 7 || burgID == 12 || burgID == 17){
                    burgID = -1;
                }else{
                    burgID = 7;
                }
            }
        }
        return;
}


void SpecifyEquations(Table_t &table)
{
        int     i, j;

        SpecifyEquations1(table);
    
        return;
}
