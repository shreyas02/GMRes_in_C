#include <stdio.h>

void main(){
    float a[4][4] = {{0,0,3,0} , {0,7,0,0} , {0,0,3,7} , {4,1,0,0}};
    float columnid[100], csrrow[100], value[100];
    int n = 4;
    int i,j,k,l,m,n1,n2;
    int test ;

    n1 = 0;
    n2 = 0;
    test = 1;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(a[i][j] != 0){
                columnid[n1] = j;
                value[n1] = a[i][j];
                if(test ==1){
                    csrrow[n2] = n1;
                    n2 =n2 + 1;
                    test = 0;
                }
                n1 = n1+1; 
            }
            if(j == n-1){
                test =1;
            }
        }
    }

//printing values array
for(i =0; i<n1; i++){
    printf("%f\t",value[i]);
}
printf("\n");


//printing coulmn id arrray
for(i =0; i<n1; i++){
    printf("%f\t",columnid[i]);
}
printf("\n");


//printing row csr array
for(i =0; i<n2; i++){
    printf("%f\t",csrrow[i]);
}
printf("\n");

}