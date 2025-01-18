#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <sys/time.h>

#define MIN 0.0
#define MAX 10.0

int n;
double start_ttime=0.0, start_ptime=0.0, end_ptime=0.0, end_ttime=0.0, tot_ptime=0.0, tot_time=0.0;

//function to check the size of the matrix
int checkPower(int num){
  if (num<=0) return 0;
  while(num%2==0){
	  num/=2;
  }
  return (num==1);
}

//function to print the matrix
void print(float** mat){
  for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			printf("%0.8f	", mat[i][j]);
		}
		printf("\n"); 
	}
}

//function to check if the transpose is correct
int check_correctness(float** T1, float** T2){
  int check=0;
  for(int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (T1[i][j]!=T2[i][j]){
        check=1;
        break;
      }
    }
  }
  return check;
}

//function to check the symmetry sequentially
int checkSym(float** m){
  start_ttime=MPI_Wtime();
  //assume the matrix is symmetric
  int check=1;
  for(int i=0; i<n; i++){
    for (int j=0; j<n; j++){//check if there are non-symmetric values
      if (m[i][j]!=m[j][i]){
        check=0;
      }
    }
  }
  end_ttime=MPI_Wtime();
  tot_time=end_ttime-start_ttime;
  
  printf("  - Wall-clock time needed for the sequential symmetric check: %0.8f ms\n", tot_time*1000);
  return check;
}

//function to transpose the matrix sequentially
void matTranspose(float** m, float** t){
  start_ttime=MPI_Wtime();
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      t[i][j]=m[j][i];
    }
  }
  end_ttime=MPI_Wtime();
  tot_time=end_ttime-start_ttime;
  
  printf("  - Wall-clock time needed for the sequential matrix transposition: %0.8f ms\n", tot_time*1000);
}

//function to check the symmetry with OpenMPI
int checkSymMPI(float** m, int rank, int s){
  int check=1;
   //assume the matrix is symmetric
  if (rank!=0){
    //allocate the matrix in each process
    m=(float**) malloc(n*sizeof(float*));
    for (int i=0; i<n; i++){
      m[i]=(float*) malloc(n*sizeof(float));
    }
    //control if the allocation is done successfully
    if (m==NULL){
      fprintf(stderr, "Error: memory allocation failed\n"); 
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  //copy the matrix in each process
  for(int i=0; i<n; i++){
    MPI_Bcast(m[i], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
   }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    start_ttime=MPI_Wtime();
  }
  //each process work on a different portion of the matrix
  int chunk=n/s;
   
  for(int i=rank*chunk; i<chunk+rank*chunk; i++){
    for (int j=i; j<n; j++){
      if (m[i][j]!=m[j][i] || i==j){
        check=0;
      } 
    }
  }
  int global_check;
  //reducing the result of every process into a global variable and pass it to all processes
  MPI_Reduce(&check, &global_check, 1, MPI_INT, MPI_PROD, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_check, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if(rank!=0){
    for (int i=0; i<n; i++){
      free(m[i]);
    }
    free(m);
  }else{
    end_ttime=MPI_Wtime();
    tot_time=end_ttime-start_ttime;
    printf("  - Wall-clock time needed for the MPI symmetric check: %0.8f ms\n", tot_time*1000);
  }
  return global_check;
}

//function to transpose the matrix with OpenMPI
void matTransposeMPI(float** m, float** t, int s, int rank){   
  if (rank==0){start_ptime=MPI_Wtime();}
  //divide the matrix into rows per process
  int rows=n/s;
  
  //allocate the local matrix, local transpose and flat matrices
  float* local_mat=(float*) malloc(rows*n*sizeof(float));
  float* local_t=(float*) malloc(rows*n*sizeof(float));
  float* flat_mat=NULL;  
  float* flat_t=NULL;
  
  //check if the allocation is done successfully
  if (local_mat==NULL || local_t==NULL){
    fprintf(stderr, "Error: memory allocation failed\n"); 
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  //allocate the matrix in a linear way
  if (rank==0){
    flat_mat=(float*) malloc(n*n*sizeof(float));
    flat_t=(float*) malloc(n*n*sizeof(float));
    if(flat_mat==NULL || flat_t==NULL){
      fprintf(stderr, "Error: memory allocation failed\n"); 
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++){
        flat_mat[n*i+j]=m[i][j];
      }
    }
  }
  //send to each process their rows
  MPI_Scatter(flat_mat, rows*n, MPI_FLOAT, local_mat, rows*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){start_ttime=MPI_Wtime();}
  
  //transpose locally the rows
  for (int i=0; i<rows; i++){
    for (int j=0; j<n; j++){
      local_t[rows*j+i]=local_mat[n*i+j];
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){end_ttime=MPI_Wtime();}
  
  //gathering the transposed rows into a flat matrix transposed
  MPI_Gather(local_t, rows*n, MPI_FLOAT, flat_t, rows*n, MPI_FLOAT, 0, MPI_COMM_WORLD);  
  
  //assign the correct value of the flat transposed to the square matrix transposed
  if (rank==0){ 
    for (int i=0; i<n; i++){
      for (int j=0;j<n; j+=rows){
        for (int k=0; k<rows; k++){
          t[i][j+k]=flat_t[i*rows+j*n+k];
        }
      }
    } 
    end_ptime=MPI_Wtime();
    tot_ptime=end_ptime-start_ptime;
    tot_time=end_ttime-start_ttime;
    printf("  - Wall-clock time needed only for the MPI matrix transposition: %0.8f ms\n", tot_time*1000);
    printf("  - Wall-clock time needed for the total MPI matrix transposition: %0.8f ms\n\n", tot_ptime*1000);
    free(flat_mat); 
    free(flat_t);
  }
  free(local_mat);
  free(local_t);  
  
}

int main(int argc, char *argv[]) {
  srand(time(0));
  int size, myrank;
  
  //initialize the MPI environment
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  if (myrank==0){ 
    //control the number of arguments
    if (argc != 2) { 
        fprintf(stderr, "Error: provide the matrix size\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    //record the size of the matrix from arguments
    n=atoi(argv[1]);
    //control if the sizee of the matrix is valid
    if (!checkPower(n)){
      fprintf(stderr, "Error: size of the matrix non valid (must be a power of 2)\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    //control if the number of processes is valid
    if (n%size!=0 || size>n){
      fprintf(stderr, "Error: number of processes non valid (must be a power of 2 and less than size of the matrix)\n");
      MPI_Abort(MPI_COMM_WORLD, 1); 
    }
  } 
  
  //send to all processes the size of the matrix
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  
  float** matrix=NULL; 
  float** T_serial=NULL;
  float** T_MPI=NULL;
  
  if (myrank==0){
    //allocate the matrices
    matrix=(float**) malloc(n*sizeof(float*));
    for (int i=0; i<n; i++){
      matrix[i]=(float*) malloc(n*sizeof(float));
    }
    T_serial=(float**) malloc(n*sizeof(float*));
    for (int i=0; i<n; i++){
      T_serial[i]=(float*) malloc(n*sizeof(float));
    }
    T_MPI=(float**) malloc(n*sizeof(float*));
    for (int i=0; i<n; i++){
      T_MPI[i]=(float*) malloc(n*sizeof(float));
    }
    //check if the allocation is done successfully
    if (matrix == NULL || T_serial==NULL || T_MPI==NULL){
      fprintf(stderr, "Error: memory allocation failed\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    //populate the matrices
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++){
        matrix[i][j]=MIN + (float)rand() / (float)(RAND_MAX / (MAX - MIN));
      }
    }
    //print(matrix);
    
    printf("\n  SEQUENTIAL EXECUTION _________________________________________________\n");
    //check if the matrix is symmetric
    if (checkSym(matrix)==1){
       fprintf(stderr, "Error: memory allocation failed\n");
       MPI_Abort(MPI_COMM_WORLD, 1);
    }else{ 
      //if the matrix is not symmetric transpose
      matTranspose(matrix, T_serial);
      //print(T_serial);
    }
    
    printf("\n  MPI PARALLEL EXECUTION _______________________________________________\n");
  }
  
  //check if the matrix is symmetric
  if (checkSymMPI(matrix, myrank, size)==1){
    if(myrank==0){
      fprintf(stderr, "Error: the matrix is symmetric\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }else{
    //if the matrix is not symmetric transpose
    matTransposeMPI(matrix, T_MPI, size, myrank);
    if (myrank==0){
      //control if the transposition gave the right result
      if (check_correctness(T_serial, T_MPI)==1){
        print(T_MPI);
        fprintf(stderr, "Error: wrong matrix transposition computation\n"); 
        MPI_Abort(MPI_COMM_WORLD, 1);
      }else{
        //print(T_MPI);
      }
    }
  }
  
  //deallocate the matrices
  if(myrank==0){
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
        free(T_serial[i]);
        free(T_MPI[i]);
    }
    free(matrix);
    free(T_MPI);
    free(T_serial);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
    
  MPI_Finalize();
  
 return 0;
 }
