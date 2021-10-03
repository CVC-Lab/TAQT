#include <par/reg3data.h>
#include <stdio.h>

int main(int argc, char* argv[]) 
{
        Reg3Data *pData = new Reg3Data(argv[1]);
        printf("func min = %f, max = %f\n", pData->getFuncMin(), pData->getFuncMax()
);
        delete pData;
        return 0;
}

