#include <stdio.h>
#include<math.h>
#include<omp.h> 

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
    int n = 4; //no of equations 
    //float a[100][100];
    //float b[100];
    float v[100][100], h[100][100];
    float x[100], sn[100], cs[100], r[100] , e1[100] , beta[100], avj[100] , y[100];
    int m = 100;
    float tollerance = 0.00000001;
    int iter =20 , n1 ; //enter number of iteration here
    int i,j,k,ni,nj;
    float temp , rnorm , bnorm, error, t,sum;

    // matrix inputs - 
    float a[4][4] = {{1,2,3,5} , {-1,7,8,0} , {2,15,3,7} , {4,1,7,9}};
    float b[4] = {1,5,7,4};

    
    //initial guess 
    for(i = 0; i <n; i ++){
        x[i]=0;
    }

    for(i = 0; i <m; i ++){
        //sn[i]=0;
        //cs[i]=0;
    }

    for(n1=0; n1<iter; n1++){

    //calculate residue 
    
    for(i=0;i<n;i++){
        r[i] = 0;
    }

    //can be done in a single loop 
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            r[i] = r[i] - a[i][j]*x[j];
        }
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
        //printf("%d\n",nj);
        //avj = a * v 

        for(i=0;i<n;i++){
            avj[i] = 0;
        }

        for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            avj[i] = avj[i] + a[i][j] * v[j][nj];
        }
        //printf("%f\n",avj[i]);
        }
        
        // for i loop
        for(ni = 0; ni<=nj; ni++){
            //hij = avj * v
            temp = 0;
            for(i=0;i<n;i++){
                temp = temp + avj[i] * v[i][ni];
            }
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


    for(i = 0;i<n; i++){
       for(j = 0; j<= nj; j++){
        x[i] = x[i] + v[i][j] * y[j];
       }
       printf("%f\n",x[i]);
    }
    printf("\n");
    
    //bug code
    //break;
    }

    return 0;
}