#include <stdio.h>
//#include<math.h>

//CSR IMPLEMENTATION 
//GMRES
float norm2d(int rows,int columns, float arr[100][100]){
    int i,j;
    float norm = 0;
    for(i = 0; i < rows; i++){
        for(j = 0; j<columns ; j++){
            norm = norm + arr[i][j] * arr[i][j];
        }
    }
    norm = pow(norm,0.5);
    return norm;
}

float norm1d(int rows, float arr[100]){
    int i;
    float norm = 0;
    for(i = 0; i < rows; i++){
            norm = norm + arr[i] * arr[i];
    }
    norm = pow(norm,0.5);
    return norm;
}


int main(){
    int n ; //no of equations 
    //float a[100][100];
    //float b[100];
    float v[1000][1000], h[1000][1000];
    float x[1000], sn[1000], cs[1000], r[1000] , e1[1000] , beta[1000], avj[1000] , y[1000];
    int m = 1000;
    float tollerance = 0.0000001;
    int iter =200000 , n1 ; //enter number of iteration here
    int i,j,k,ni,nj;
    float temp , rnorm , bnorm, error, t,sum;

    // matrix inputs - 
    float columnid[9] = {0,1,2,1,2,3,2,3,3};
    float csrrow[4] = {0,3,6,8};
    float value[9] = {1,2,3,1,2,3,1,2,1};
    float b[4] = {1,5,7,4};

    //new variables required for CSR 
    int nb, nc;
    int size, size1;  
    int a1, c1 , b1;

    //taking size of the vector
    size = sizeof(b);
    size1 = sizeof(b[0]);
    n = size/size1;

    
    //initial guess 
    for(i = 0; i <n; i ++){
        x[i]=0;
    }

    for(n1=0; n1<iter; n1++){

    //calculate residue 
    
    for(i=0;i<n;i++){
        r[i] = 0;
    }
    
    //new code here 

    size = sizeof(columnid);
    size1 = sizeof(columnid[0]);
    nb = size/size1;

    c1=0;
    a1=-1;
    for(i=0;i<nb;i++){
        nc = csrrow[c1];
        if(nc == i){
            a1 = a1+1;
            c1 = c1+1;
        }
        b1=columnid[i];
        r[a1] = r[a1] - value[i] * x[b1];
    }

    for(i=0;i<n;i++){
        r[i] = r[i] + b[i];
    }

    rnorm = norm1d(n,r);
    

    for(i=0;i<n;i++){
        v[i][0] = r[i]/rnorm ; 
        //printf("%f\n",v[i][0]);
    }

    bnorm = norm1d(n,b);
    //printf("%f\n",bnorm);
    error = rnorm / bnorm;
    //printf("%f\n",error);

    //zeroing e1 vector of size m+1
    for(i=0;i<m+1;i++){
        e1[i] = 0;
    }
    e1[0] = 1;

    //getting beta 
    for(i=0;i<m+1;i++){
        beta[i] = rnorm*e1[i];
        //printf("%f\n",beta[i]);
    }

    if(error < tollerance) {
            break;
        }

    //starting arnoldi iterations 
    for( nj=0;nj<m;nj++){

        for(i=0;i<n;i++){
            avj[i] = 0;
        }

        // new code here 
        size = sizeof(columnid);
        size1 = sizeof(columnid[0]);
        nb = size/size1;

        c1=0;
        a1=-1;
        for(i=0;i<nb;i++){
            nc = csrrow[c1];
            if(nc == i){
             a1 = a1+1;
             c1 = c1+1;
            }
            b1=columnid[i];
            avj[a1] = avj[a1] + value[i] * v[b1][nj];
        }
        
        // for i loop
        for(ni = 0; ni<=nj; ni++){
            //hij = avj * v
            temp = 0;
            for(i=0;i<n;i++){
                temp = temp + avj[i] * v[i][ni];               //no corellation between i and ni
            }                                                  //ni goes to m and i goes to m 
            h[ni][nj] = temp ;
            // avj = avj - hij*vi
            for(i=0;i<n;i++){
                avj[i] = avj[i] - h[ni][nj] * v[i][ni];
            }
        }
        temp = norm1d(n,avj);
        h[nj+1][nj] = temp;
        for(i=0;i<n;i++){
            v[i][nj+1] = avj[i] / temp ;
            //printf("%f\n",v[i][nj+1]); 
        }

        //applying rotation 
        for(ni = 0; ni <= nj-1 ; ni++){
            temp = cs[ni] * h[ni][nj] + sn[ni]*h[ni+1][nj];
            h[ni+1][nj] = - sn[ni] * h[ni][nj] + cs[ni]*h[ni+1][nj]; 
            h[ni][nj] = temp;
        }

        //update sin cos values for rotation 

        t = h[nj][nj] * h[nj][nj] + h[nj+1][nj] * h[nj+1][nj] ; 
        t = pow(t,0.5);
        //printf("%f\n",t);
        cs[nj] = h[nj][nj] / t;
        sn[nj] = h[nj+1][nj]/t ; 
        //printf("%f\t%f\n",cs[nj],sn[nj]);

        //eliminate hj+1i
        h[nj][nj] = cs[nj] * h[nj][nj] + sn[nj]*h[nj+1][nj];
        h[nj+1][nj] = 0;

        beta[nj+1] = -sn[nj] *beta[nj];
        beta[nj] = cs[nj] * beta[nj];
        error = fabs(beta[nj+1] / bnorm) ;
        //printf("%f\n",error); 

        if(error < tollerance) {
            break;
        }
    }
    // calculate the result 
    nj = nj-1;
    y[nj] = beta[nj]/h[nj][nj];
    //printf("%f\n",y[nj]);
    for(i = nj-1; i >= 0; i--){
        sum = beta[i];
        //printf("%f\n",sum);
        for(k = i+1; k <= nj; k++){
            sum = sum - h[i][k]*y[k];
            //printf("%f\n",h[i][k]);
        }
        y[i] = sum / h[i][i] ;
    }

    }

    for(i = 0;i<n; i++){
       for(j = 0; j<= nj; j++){
        x[i] = x[i] + v[i][j] * y[j];
       }
    }

    //printing results
    printf("X is - \n");
    printf("\n");

    for(i = 0;i<n; i++){
       printf("%f\n",x[i]);
    }
    printf("\n");

    return 0;
}