/* Declarations */

/* Definitions */

# define ZERO 0
# define ONE 1
# define TWO 2
# define THREE 3
# define FOUR 4
# define SEVEN 7
# define EIGHT 8
# define LENGSOLN 3
# define LENGCODSTA 5
# define LENGDOMES 10
# define LENGPARA 25

# define TRAN 0.0001e0        /* 0.1 mm */
# define SCAL 1.5678559429e-2 /* 0.1 mm converted into ppb */              
# define ROTA 3.2339352256e-6 /* 0.1 mm converted into arc seconds */
# define TINY 1.e-20          /* Threshold for rootmat function (cf. matrix_test_networks) */

/* Structure for Genetically Modified Networks */

typedef struct 
{
	int nbrpar;      /* Number of parameters */
	int nbrsta;      /* Number of stations */
	int nbreop;      /* Number of EOP */

	double nbrobs;   /* Number of observations */
	double nbrunk;   /* Number of unknowns (including the reduced parameters) */
	double ssores;   /* Square sum of residuals */
	double sigma;    /* Unit variance factor of the least-square computation with minimum constraints */

	char **codsta;   /* Station codes */
	char **domes;    /* Station DOMES numbers */
	char **soln;     /* Station SOLN */
	char **parnam;   /* Parameter names */
  
	double *norvec;  /* Vector of the normal system without any constraint */
	double *solvec;  /* Solution computed with minimum constraints */
	double *solori;  /* Solution provided in the SINEX file */
	double **normat; /* Normal matrix without any constraint */
	double **stamat; /* Matrix for transformation of station positions */
	double **prtmat; /* Matrix for transformation of Xp and Yp */
	double **solmat; /* Variance-covariance matrix of solution computed with minimum constraints */
	double **covori; /* Variance-covariance matrix provided in the SINEX file */
	double **aprmat; /* A priori variance-covariance matrix provided in the SINEX file */
} gmn;

/* Functions */

/* auxiliary_test_networks */

/* Function to compute the minimum value of a double vector */
double minvect(double *vector,int n);
/* Function to compute the maximum value of a double vector */
double maxvect(double *vector,int n);
/* Function to compute the weighted mean value of a double vector */
double weightmean(double *vector,double *var,int n);
/* Function to compute the mean value of a double vector */
double meanvect(double *vector,int n);
/* Function to sort the vector in crescent order */
double *sortvect(double *vector,int n);
/* Function to compute the median value of a double vector */
double medvect(double *vector,int n);
/* Function to compute the standard deviation of the average of decorrelated variables */
double meanvar1(double **matcov,int n);
/* Function to compute the standard deviation of the average of correlated variables */
double meanvar2(double **matcov,int n);
/* Routine to compute the distribution of a given station network wrt the three directions X, Y, and Z */
void distribution3(gmn *pgmn,int *datum,double *distri);
/* Routine to compute the distribution of a given station network wrt the three directions X, Y, and Z */
void distribution8(gmn *pgmn,int *datum,double *distri);
/* Function to compute the histogram of the spherical errors of a given station position set */
double *histosta(gmn *pgmn,double **matrix);
/* Function to compute the WRMS of a vector of numerical values provided with their formal errors */
double wrms(double *value,double *error,int n);
/* Function to compute the mean presence of the stations in a given station network */
double mean_presence(gmn *pgmn,int *datum,char *file_presence);
/* Function to check if a parameter label matches a station */
int check_parameter(gmn *pgmn,int index_par,int index_sta);

/* matrix_test_networks */

/* Function to dynamically allocate memory for an integer vector with l rows */
int *iallovect(int l);
/* Routine to free the memory allocated by iallovect */
void ifreevect(int *vector);
/* Function to dynamically allocate memory for a double vector with l rows */
double *allovect(int l);
/* Routine fo free the memory allocated by allovect */
void freevect(double *vector);
/* Function to dynamically allocate memory for a double matrix with l rows and c columns */
double **allomat(int l,int c);
/* Routine fo free the memory allocated by allomat */
void freemat(double **matrix,int l);
/* Routine to initialize an integer vector of dimension l */
void initivect(int *vector,int l);
/* Routine to initialize a double vector of dimension l */
void initvect(double *vector,int l);
/* Routine to initialize a double matrix of dimension (l,c) */
void initmat(double **matrix,int l,int c);
/* Routine to scale the matrix with the double factor */
void scalemat(double **matrix,double factor,int l,int c);
/* Routine to scale the vector with the double factor */
void scalevect(double *vector,double factor,int l);
/* Routine to print the matrix on screen */
void printmat(double **matrix,int l,int c);
/* Routine to print the vector on screen */
void printvect(double *vector,int l);
/* Routine to add the two matrices mat1 and mat2 */
void addmat(double **mat1,double **mat2,int l,int c);
/* Routine to add the two vectors vect1 and vect2 */
void addvect(double *vect1,double *vect2,int l);
/* Routine to substract the matrix mat2 from mat1 */
void submat(double **mat1,double **mat2,int l,int c);
/* Routine to substract the vector vect2 from vect1 */
void subvect(double *vect1,double *vect2,int l);
/* Function to multiply the matrices mat1 and mat2 */
double **matmul(double **mat1,int l1,int c1,double **mat2,int l2,int c2);
/* Function to multiply matrix by vector */
double *matvectmul(double **matrix,int l,int c,double *vector,int n);
/* Function to transpose the matrix */
double **transpose(double **matrix,int l,int c);
/* Function to compute the infinite norm of the matrix */
double infnorm(double **matrix,int n);
/* Function to compute the root of the matrix by Cholevsky's method */
double **rootmat(double **matrix,int n);
/* Function to compute the inverse of the matrix by Cholevsky's method */
double **cholev(double **matrix,int n);
/* Function to compute the determinant of the matrix */
double det(double **matrix,int n);
/* Function to compute the condition number of the matrix */
double cond(double **matrix,int n);
/* Function to dynamically allocate memory for a char matrix with l rows and c columns */
char **allocharmat(int l,int c);
/* Routine to free the memory allocated by allocharmat */
void freecharmat(char **matrix,int l);
/* Routine to print the vector of strings on screen */
void printcharmat(char **charmat,int l);

/* genmodnet_test_networks */

/* Function to initialize the gmn structure by the provided files */
gmn *initgmn(char *file_para,char *file_tran);
/* Routine to free the gmn structure pointed by pgmn allocated with initgmn */
void freegmn(gmn *pgmn);
/* Function to compute a solution with minimum constraints applied on given datums */
int calsol(gmn *pgmn,int *datum,int *keys);
/* Function to project matrix with minimum constraints over a datum */
double **proj_matrix(gmn *pgmn,double **matrix,int *datum,int *keys);
/* Function to compute the stability matrix for a given datum for stations */
double **stability_matrix(gmn *pgmn,int *datum,int *keys);
/* Function to compute the reference system effect on the station (key==1) or source position (key==2) basis for pgmn->solmat */
double **sysreffect(gmn *pgmn,double **matrix,int key);
/* Function to compute the Helmert transformation parameters on the station (key==1) or source position (key==2) basis */
/* with (key_apr=1) or without (key_apr=0) the a priori variance-covariance matrix for station positions */
int transfo(gmn *pgmn,double *solution,double **variance,double *partra,double **vartra,int key,int key_apr);
/* Function to compute the Helmert transformation for a given datum of stations and all the transformed parameters */
/* with (key_apr=1) or without (key_apr=0) the a priori variance-covariance matrix for station positions */
int transfo_datum(gmn *pgmn,int *datum,double *solution,double **variance,double *partra,double **vartra,double *soltransfo,double **vartransfo,int key_apr);

