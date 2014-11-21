#include<stdio.h>
#include<math.h>
int main()
{
    int i, n; //size of matrix to be Gauss-Jordan'd
    float A[20][20], b[10]; //define matix of more indices than will likely ever be used in the problem

    //NOTE: program currently only working for 3x3 matrices 

    printf("\nEnter size of matrix: ");
    scanf("%d",&n);
    printf("\nEnter elements of matrix A by row (one element at a time):\n");
    for(i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            printf(" A[%d][%d]:", i,j);
            scanf("%f",&A[i][j]);
        }
    }
    printf("\nEnter elements of vector b (one element at a time):\n");
    for(i=0; i<n; i++)
    {
            printf(" b[%d]:", i);
            scanf("%f",&b[i]);
    }
    //see if user wants pivoting
    int pivot;
     printf(" Do you want pivoting? (1 = yes, 0 = no)");
     scanf("%d",&pivot);

   
    //combine matrix A and vector b so all operations will be done to both
    float Ab[20][20];
    for (int i = 0; i<n; i++){
        for (int j = 0; j< n+1; j++){
            if (j<n){
	        Ab[i][j] = A[i][j];
	    }
	    if (j==n){
	        Ab[i][j] = b[i];
	    }
	}
    }

  //print the matrix prior to sorting
  for (int i=0; i<n; i++){
        for( int j = 0; j<n+1; j++){
	   printf("%f ", Ab[i][j]);
	}
	printf("\n");
    }
    printf("\n");
	
    //if pivot = yes, sort the rows so largest element is at the top
    float temp[n+1];
    if (pivot == 1){
        for (int i = 0; i<n-1; i++){
	    if( fabs(Ab[i][0]) < fabs(Ab[i+1][0])){
		for (int j=0; j<n+1; j++){
	            temp[j] = Ab[i][j];
		    Ab[i][j] = Ab[i+1][j];
		    Ab[i+1][j] = temp[j];
	 	}
            }
            if( fabs(Ab[0][0]) < fabs(Ab[i][0])){
		for (int j=0; j<n+1; j++){
	            temp[j] = Ab[0][j];
		    Ab[0][j] = Ab[i][j];
		    Ab[i][j] = temp[j];
	 	}
	    }
	}
    }
    //print the matrix after sorting
   for (int i=0; i<n; i++){
        for( int j = 0; j<n+1; j++){
	   printf("%f ", Ab[i][j]);
	}
	printf("\n");
    }
    printf("\n");

   
    //divide row 0 by Ab[0][0]
    //float temp;
    float a00 = Ab[0][0];
    for (int j=0; j<n+1; j++){
        temp[j] = Ab[0][j];
	Ab[0][j] = temp[j]/a00;
    }
    //subtract stuff from row 1 to get A[1][0] = 0
    float a10 = Ab[1][0];
    for (int j=0; j<n+1; j++){
        temp[j] = Ab[1][j];
	Ab[1][j] = temp[j] - a10*Ab[0][j];
    } 
   //subtract stuff from row 2 to get A[2][0] = 0
    float a20 = Ab[2][0];
    for (int j=0; j<n+1; j++){
        temp[j] = Ab[2][j];
	Ab[2][j] = temp[j] - a20*Ab[0][j];
    } 
  

   /*
   for (int i=0; i<n; i++){
        for( int j = 0; j<n+1; j++){
	   printf("%f ", Ab[i][j]);
	}
	printf("\n");
    }
    printf("\n");*/

    
    //divide row 1 by A[1][1]
    float a11 = Ab[1][1];
    for (int j=0; j<n+1; j++){
        temp[j] = Ab[1][j];
	Ab[1][j] = temp[j]/a11;
    }
    //subtract stuff from row 2 to get A[2][1] = 0
    float a21 = Ab[2][1];
    for (int j=0; j<n+1; j++){
        temp[j] = Ab[2][j];
	Ab[2][j] = temp[j] - a21*Ab[1][j];
    } 

    //divide row 2 by A[2][2]
    float a22 = Ab[2][2];
    for (int j=0; j<n+1; j++){
        temp[j] = Ab[2][j];
	Ab[2][j] = temp[j]/a22;
    }

    //print the matrix after all operations done
    for (int i=0; i<n; i++){
        for( int j = 0; j<n+1; j++){
	   printf("%f ", Ab[i][j]);
	}
	printf("\n");
    }

    /*
    //subtract off the second and third term of row 0
    float a12 = Ab[1][2];
    for (int j=0; j<n+1; j++){
        temp[j] = Ab[1][j];
	Ab[1][j] = temp[j] - a12*A[2][j];
    }
    for (int i=0; i<n; i++){
        for( int j = 0; j<n+1; j++){
	   printf("%f ", Ab[i][j]);
	}
	printf("\n");
    }
	printf("\n");

    float a01 = Ab[0][1];
    float a02 = Ab[0][2];
    for (int j=0; j<n+1; j++){
        temp[j] = Ab[0][j];
	Ab[0][j] = temp[j] - a01*A[1][j];
        temp[j] = Ab[0][j];
	Ab[0][j] = temp[j] - a02*A[2][j];
    }*/

    //now we have upper triangular matrix, from here we can print out the x's
    float x3 = Ab[n-1][n];
    float x2 = Ab[n-2][n] - Ab[n-2][n-1]*x3;
    float x1 = Ab[n-3][n] - Ab[n-3][n-2]*x2 - Ab[n-3][n-1]*x3;
   
    printf("x1 = %f\n", x1);
    printf("x2 = %f\n", x2);
    printf("x3 = %f\n", x3);
   
    return(0);
}

