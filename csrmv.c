#include <stdio.h>

void main(){
    float columnid[6] = {2,1,2,3,0,1};
    float csrrow[4] = {0,1,2,4};
    float value[6] = {3,7,3,7,4,1};
    float vector[4] = {0,1,2,4};
    float solution[4];
    int na = 4;
    int nb , nc;
    int size, size1;  
    int a1, c1 , b1;
    int i;
    
    size = sizeof(columnid);
    size1 = sizeof(columnid[0]);
    nb = size/size1;

    for(i=0;i<na;i++){
        solution[i] = 0 ;
    }

    c1=0;
    a1=-1;
    for(i=0;i<nb;i++){
        nc = csrrow[c1];
        if(nc == i){
            a1 = a1+1;
            c1 = c1+1;
        }
        b1=columnid[i];
        solution[a1] = solution[a1] + value[i] * vector[b1];
    }

    for(i=0;i<na;i++){
        printf("%f\n",solution[i]);
    }
}