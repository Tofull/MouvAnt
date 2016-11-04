/* Bibliotheque de calcul pour l'utilisation de mouv'ant*/
/* Inspiré du */
/* Program for network testing for genetic algorithms */
/* Author: D. Coulot */
/* Last Modification: 2015/12/03 */
/* Toutes les fonctionnalités nouvelles apportées par Mouv'Ant se trouvent à la fin de ce fichier */

/* Inclusions */

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
#include <string.h>

# include "global.h"

/* Function to compute the minimum value of a double vector */
/* The vector is allocated outside this function */

double minvect(double *vector,int n)
{
	int i;

	double minval;

	minval=1.e+16;

	for (i=0;i<n;i++)
		if (vector[i]<minval)
			minval=vector[i];

	return minval;
}

/* Function to compute the maximum value of a double vector */
/* The vector is allocated outside this function */

double maxvect(double *vector,int n)
{
	int i;

	double maxval;

	maxval=-1.e+16;

	for (i=0;i<n;i++)
		if (vector[i]>maxval)
			maxval=vector[i];

	return maxval;
}

/* Function to compute the mean value of a double vector */
/* The vector is allocated outside this function */

double meanvect(double *vector,int n)
{
	int i;

	double mean;

	mean=0.;

	for (i=0;i<n;i++)
		mean+=vector[i];

	return mean/((double) n);
}

/* Function to compute the weighted mean value of a double vector */
/* The vector and the vector of the variances are allocated outside this function */

double weightmean(double *vector,double *var,int n)
{
	int i;

	double sum;
	double mean;

	sum=0.;
	mean=0.;

	for (i=0;i<n;i++)
		{
			sum+=1./var[i];
			mean+=vector[i]/var[i];
		}
		
	return mean/sum;
}

/* Function to sort the vector in crescent order */
/* The vector containing the result is allocated in this function */
/* This vector must be freed outside this function */

double *sortvect(double *vector,int n)
{
	int i,j;

	int flag;

	double temp;

	double *vectint;

	vectint=(double *) allovect(n);

	for (i=0;i<n;i++)
		vectint[i]=vector[i];

	flag=0;

	for (i=0;i<=n;i++)
		{
			flag=1;
			for (j=0;j<n-1-i-1;j++)
				if (vectint[j+1]<=vectint[j])
					{
						temp=vectint[j];
						vectint[j]=vectint[j+1];
						vectint[j+1]=temp;
						flag=0;
					}
				if (flag==1) 
					break;
		}

	return vectint;
}

/* Function to compute the median value of a double vector */
/* The vector is allocated outside this function */

double medvect(double *vector,int n)
{
	double *vectint;

	vectint=(double *) sortvect(vector,n);

	if (n%2==1)
		return vectint[(n+1)/2];
	else
		return 0.5*(vectint[n/2]+vectint[(n/2)+1]);
 }

/* Function to compute the standard deviation of the average of decorrelated variables */
/* The variance-covariance matrix of these variables is provided as input and is allocated outside the function */

double meanvar1(double **matcov,int n)
{
	int i;

	double var;

	var=0.;

	for (i=0;i<n;i++)
		var+=matcov[i][i];

	var/=((double) n*n);

	return sqrt(var);
}

/* Function to compute the standard deviation of the average of correlated variables */
/* The variance-covariance matrix of these variables is provided as input and is allocated outside the function */

double meanvar2(double **matcov,int n)
{
	int i,j;

	double var;

	var=0.;

	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			var+=matcov[i][j];

	var/=((double) n*n);

	return sqrt(var);
}



/* Routine to compute the distribution of a given station network wrt the three directions X, Y, and Z */

void distribution3(gmn *pgmn,int *datum,double *distri)
{
	int i;
	int sum;
	int sumx,sumy,sumz;

	double coordx,coordy,coordz; 

	sum=0;
	sumx=sumy=sumz=0;

	for (i=0;i<pgmn->nbrsta;i++)
		if (datum[i]==1)
			{
				sum+=1;

				coordx=(pgmn->stamat[3*i][3])*1.e9;
				coordy=(pgmn->stamat[3*i+1][3])*1.e9;
				coordz=(pgmn->stamat[3*i+2][3])*1.e9;

				if (coordx>0) 
					sumx+=1;
				if (coordy>0) 
					sumy+=1;
				if (coordz>0) 
					sumz+=1;
			}

	distri[0]=(((double) sumx)/((double) sum))*100.;
	distri[1]=(((double) sumy)/((double) sum))*100.;
	distri[2]=(((double) sumz)/((double) sum))*100.;
}

/* Routine to compute the distribution of a given station network wrt the three directions X, Y, and Z (%) */

void distribution8(gmn *pgmn,int *datum,double *distri)
{
	int i;
	int sum=0;
	int sumxpypzp=0;
	int sumxpypzm=0;
	int sumxpymzp=0;
	int sumxpymzm=0;
	int sumxmypzp=0;
	int sumxmypzm=0;
	int sumxmymzp=0;
	int sumxmymzm=0;

	double coordx,coordy,coordz; 

	for (i=0;i<pgmn->nbrsta;i++)
		if (datum[i]==1)
			{
				sum+=1;

				coordx=(pgmn->stamat[3*i][3])*1.e9;
				coordy=(pgmn->stamat[3*i+1][3])*1.e9;
				coordz=(pgmn->stamat[3*i+2][3])*1.e9;

				if ((coordx>0)&&(coordy>0)&&(coordz>0))
					sumxpypzp+=1;
				if ((coordx>0)&&(coordy>0)&&(coordz<0))
					sumxpypzm+=1;
				if ((coordx>0)&&(coordy<0)&&(coordz>0))
					sumxpymzp+=1;
				if ((coordx>0)&&(coordy<0)&&(coordz<0))
					sumxpymzm+=1;
				if ((coordx<0)&&(coordy>0)&&(coordz>0))
					sumxmypzp+=1;
				if ((coordx<0)&&(coordy>0)&&(coordz<0))
					sumxmypzm+=1;
				if ((coordx<0)&&(coordy<0)&&(coordz>0))
					sumxmymzp+=1;
				if ((coordx<0)&&(coordy<0)&&(coordz<0))
					sumxmymzm+=1;
			}

	distri[0]=(((double) sumxpypzp)/((double) sum))*100.;
	distri[1]=(((double) sumxpypzm)/((double) sum))*100.;
	distri[2]=(((double) sumxpymzp)/((double) sum))*100.;
	distri[3]=(((double) sumxpymzm)/((double) sum))*100.;
	distri[4]=(((double) sumxmypzp)/((double) sum))*100.;
	distri[5]=(((double) sumxmypzm)/((double) sum))*100.;
	distri[6]=(((double) sumxmymzp)/((double) sum))*100.;
	distri[7]=(((double) sumxmymzm)/((double) sum))*100.;
}

/* Function to compute the histogram of the spherical errors of a given station position set */
/* The vector containing the results is allocated in this function */
/* This vector must be freed outside this function */

double *histosta(gmn *pgmn,double **matrix)
{
	int i;
	int count;
  
	double *error;

	error=(double *) allovect(pgmn->nbrsta);

	count=-1;
	for (i=0;i<pgmn->nbrpar;i++)
		if ((pgmn->parnam[i][0]=='S')&&(pgmn->parnam[i][3]=='X'))
			{
				count++;
				error[count]=sqrt(matrix[i][i]+matrix[i+1][i+1]+matrix[i+2][i+2]+2.*matrix[i][i+1]+2.*matrix[i][i+2]+2.*matrix[i+1][i+2]);
			}

	return error;
}

/* Function to compute the WRMS of a vector of numerical values provided with their formal errors */

double wrms(double *value,double *error,int n)
{
	int i;

	double stat1,stat2;

	stat1=0;
	stat2=0;

	for (i=0;i<n;i++)
		{
			stat1+=(value[i]*value[i])/(error[i]*error[i]);
			stat2+=1./(error[i]*error[i]);
		}

	return sqrt(stat1/stat2);
}

/* Function to compute the mean presence of the stations in a given station network */

double mean_presence(gmn *pgmn,int *datum,char *file_presence)
{
	int i;
	int intpresence;
	int dimnet=0;

	double sum;
	double presence;

	char char_station[17];

	FILE *fp;

	sum=0.;

	/* Opening of the file */

	fp=fopen(file_presence,"r");

	/* Loop over the stations in the file */

	while ((fscanf(fp,"%s %lf %d\n",char_station,&(presence),&(intpresence))) != EOF)
		{	
			for (i=0;i<pgmn->nbrsta;i++)
				if (datum[i]==1)
					{
						if ((pgmn->codsta[i][0]==char_station[0])&&(pgmn->codsta[i][1]==char_station[1]))
							if ((pgmn->codsta[i][2]==char_station[2])&&(pgmn->codsta[i][3]==char_station[3]))
								if ((pgmn->codsta[i][6]==char_station[6])&&(pgmn->domes[i][0]==char_station[8]))
									if((pgmn->domes[i][1]==char_station[9])&&(pgmn->domes[i][2]==char_station[10])) 
										if((pgmn->domes[i][3]==char_station[11])&&(pgmn->domes[i][4]==char_station[12])) 
											if((pgmn->domes[i][5]==char_station[13])&&(pgmn->domes[i][6]==char_station[14])) 
												if((pgmn->domes[i][7]==char_station[15])&&(pgmn->domes[i][8]==char_station[16])) 
													{
														dimnet+=1;
														sum+=presence;
													}
					}
		}

	/* Closing of the file */
  
	fclose(fp);

	/* End */

	return sum/((double) dimnet);
}

/* Function to check if a parameter label matches a given station */

int check_parameter(gmn *pgmn,int index_par,int index_sta)
{
	int flag=0;
  
	if ((pgmn->parnam[index_par][5]==pgmn->codsta[index_sta][0])&&(pgmn->parnam[index_par][6]==pgmn->codsta[index_sta][1]))
		if ((pgmn->parnam[index_par][7]==pgmn->codsta[index_sta][2])&&(pgmn->parnam[index_par][8]==pgmn->codsta[index_sta][3]))
			if (pgmn->parnam[index_par][11]==pgmn->codsta[index_sta][6])
				if ((pgmn->parnam[index_par][13]==pgmn->domes[index_sta][0])&&(pgmn->parnam[index_par][14]==pgmn->domes[index_sta][1]))
					if ((pgmn->parnam[index_par][15]==pgmn->domes[index_sta][2])&&(pgmn->parnam[index_par][16]==pgmn->domes[index_sta][3]))
						if ((pgmn->parnam[index_par][17]==pgmn->domes[index_sta][4])&&(pgmn->parnam[index_par][18]==pgmn->domes[index_sta][5]))
							if ((pgmn->parnam[index_par][19]==pgmn->domes[index_sta][6])&&(pgmn->parnam[index_par][20]==pgmn->domes[index_sta][7]))
								if ((pgmn->parnam[index_par][21]==pgmn->domes[index_sta][8])&&(pgmn->parnam[index_par][22]==pgmn->soln[index_sta][0]))
									if (pgmn->parnam[index_par][23]==pgmn->soln[index_sta][1])
										flag=1;

	return flag;
}



/* Function to initialize the gmn structure by the provided files */
/* The pointed gmn structure is allocated in this function */
/* It must be freed outside this function */ 

gmn *initgmn(char *file_para,char *file_tran)
{
	int i,j;
	int ii,jj;
	
	char commen[80];

	FILE *fp1,*fp2;

	gmn *pgmn;

	/* Opening of both files */

	fp1=fopen(file_para,"r");
	fp2=fopen(file_tran,"r");

	/* Allocation of the pointer pgmn */

	pgmn=(gmn *) malloc((size_t) sizeof(gmn));

	/* Reading of the number of parameters */

	fscanf(fp1,"%s\n",commen);
	fscanf(fp1," %d\n",&(pgmn->nbrpar));

	/* Reading of the numbers of stations and EOP */

	fscanf(fp2,"%s\n",commen);
	fscanf(fp2," %d\n",&(pgmn->nbrsta));

	fscanf(fp2,"%s\n",commen);
	fscanf(fp2," %d\n",&(pgmn->nbreop));

	fscanf(fp2,"%s\n",commen);
	fscanf(fp2,"%s\n",commen);

	/* Allocations of the members of the pointed gmn structure */

	/* Allocation of codsta */

	pgmn->codsta=(char **) allocharmat(pgmn->nbrsta,LENGCODSTA);

	/* Allocation of domes */

	pgmn->domes=(char **) allocharmat(pgmn->nbrsta,LENGDOMES);

	/* Allocation of soln */

	pgmn->soln=(char **) allocharmat(pgmn->nbrsta,LENGSOLN);

	/* Allocation of parnam */

	pgmn->parnam=(char **) allocharmat(pgmn->nbrpar,LENGPARA);

	/* Allocation of norvec */

	pgmn->norvec=(double *) allovect(pgmn->nbrpar);

	/* Allocation of solvec */

	pgmn->solvec=(double *) allovect(pgmn->nbrpar);

	/* Allocation of solori */

	pgmn->solori=(double *) allovect(pgmn->nbrpar);

	/* Allocation of normat */

	pgmn->normat=(double **) allomat(pgmn->nbrpar,pgmn->nbrpar);

	/* Allocation of stamat */

	pgmn->stamat=(double **) allomat((3*pgmn->nbrsta),SEVEN);

	/* Allocation of prtmat */

	pgmn->prtmat=(double **) allomat(pgmn->nbreop,TWO);

	/* Allocation of solmat */

	pgmn->solmat=(double **) allomat(pgmn->nbrpar,pgmn->nbrpar);

	/* Allocation of covori */

	pgmn->covori=(double **) allomat(pgmn->nbrpar,pgmn->nbrpar);

	/* Allocation of aprmat */

	pgmn->aprmat=(double **) allomat(pgmn->nbrpar,pgmn->nbrpar);

	/* Reading of parnam */

	fscanf(fp1,"%s\n",commen);
	for (i=0;i<pgmn->nbrpar;i++)
		fscanf(fp1," %s\n",pgmn->parnam[i]);
    
	/* Reading of nbrobs */

	fscanf(fp1,"%s\n",commen);
	fscanf(fp1," %lf\n",&(pgmn->nbrobs));
    
	/* Reading of nbrunk */

	fscanf(fp1,"%s\n",commen);
	fscanf(fp1," %lf\n",&(pgmn->nbrunk));
    
	/* Reading of ssores */

	fscanf(fp1,"%s\n",commen);
	fscanf(fp1," %lf\n",&(pgmn->ssores));
    
	/* Reading of the estimated solution */

	fscanf(fp1,"%s\n",commen);
	for (i=0;i<pgmn->nbrpar;i++)
		fscanf(fp1," %d %lf\n",&ii,&(pgmn->solori[i]));

	/* Reading of the variance covariance matrix */

	fscanf(fp1,"%s\n",commen);
	for (i=0;i<pgmn->nbrpar;i++)
		{
			for (j=0;j<=i;j++)
				fscanf(fp1," %d %d %lf\n",&ii,&jj,&(pgmn->covori[i][j]));
			for (j=0;j<i;j++)
				pgmn->covori[j][i]=pgmn->covori[i][j];
		}

	/* Reading of the a priori variance covariance matrix */

	fscanf(fp1,"%s\n",commen);
	for (i=0;i<pgmn->nbrpar;i++)
		{
			for (j=0;j<=i;j++)
				fscanf(fp1," %d %d %lf\n",&ii,&jj,&(pgmn->aprmat[i][j]));
			for (j=0;j<i;j++)
				pgmn->aprmat[j][i]=pgmn->aprmat[i][j];
		}
		
	/* Reading of normat */

	fscanf(fp1,"%s\n",commen);
	for (i=0;i<pgmn->nbrpar;i++)
		{
			for (j=0;j<=i;j++)
				fscanf(fp1," %d %d %lf\n",&ii,&jj,&(pgmn->normat[i][j]));
			for (j=0;j<i;j++)
				pgmn->normat[j][i]=pgmn->normat[i][j];
		}

	/* Reading of norvec */

	fscanf(fp1,"%s\n",commen);
	for (i=0;i<pgmn->nbrpar;i++)
		fscanf(fp1," %d %lf\n",&ii,&(pgmn->norvec[i]));
 
	/* Reading of codsta */

	fscanf(fp2,"%s\n",commen);
	for (i=0;i<pgmn->nbrsta;i++)
		fscanf(fp2," %s\n",pgmn->codsta[i]);

	/* Reading of domes */

	fscanf(fp2,"%s\n",commen);
	for (i=0;i<pgmn->nbrsta;i++)
		fscanf(fp2," %s\n",pgmn->domes[i]);

	/* Reading of soln */

	fscanf(fp2,"%s\n",commen);
	for (i=0;i<pgmn->nbrsta;i++)
		fscanf(fp2," %s\n",pgmn->soln[i]);

	/* Reading of stamat */
    
	fscanf(fp2,"%s\n",commen);
	for (i=0;i<(3*pgmn->nbrsta);i++)
		for (j=0;j<SEVEN;j++)
			fscanf(fp2," %d %d %lf\n",&ii,&jj,&(pgmn->stamat[i][j]));

	/* Reading of prtmat */
    
	fscanf(fp2,"%s\n",commen);
	for (i=0;i<pgmn->nbreop;i++)
		for (j=0;j<TWO;j++)
			fscanf(fp2," %d %d %lf\n",&ii,&jj,&(pgmn->prtmat[i][j]));

	/* Closure of both files */

	fclose(fp1);
	fclose(fp2);

	/* End */

	return pgmn;
}

/* Routine to free the gmn structure pointed by pgmn allocated with initgmn */

void freegmn(gmn *pgmn) 
{
	freevect(pgmn->norvec);
	freevect(pgmn->solvec);
	freevect(pgmn->solori);
	
	freemat(pgmn->normat,pgmn->nbrpar);
	freemat(pgmn->stamat,3*(pgmn->nbrsta));
	freemat(pgmn->prtmat,pgmn->nbreop);
	freemat(pgmn->solmat,pgmn->nbrpar);
freemat(pgmn->covori,pgmn->nbrpar);
	freemat(pgmn->aprmat,pgmn->nbrpar);
	
	freecharmat(pgmn->codsta,pgmn->nbrsta);
	freecharmat(pgmn->domes,pgmn->nbrsta);
	freecharmat(pgmn->soln,pgmn->nbrsta);
	freecharmat(pgmn->parnam,pgmn->nbrpar);

	free(pgmn);
	pgmn=NULL;
}

/* Function to compute a solution with minimum constraints applied on a given datum for stations */
/* keys = array of keys for the minimum constraints (0=parameter to be excluded, 1=parameter to be considered) */
/* From (Sillard & Boucher, J Geodesy, 2001) */

int calsol(gmn *pgmn,int *datum,int *keys)
{
	int i,j,k;
	int size;
	int counter;

	double *solution;
	double *resmin;
  
	double **matd,**matb;
	double **tmatd,**tmatb;
	double **dtd,**dtdi;
	double **sigtheta;
	double **sb,**btsb;
	double **norinv;
	double **mat_temp;
  
	/* Computation of size (number of parameters to be considered) */

	size=0;
	for (i=0;i<SEVEN;i++)
		size+=keys[i];

	/* Computation of the matrix D related to the station positions of the network for all the transformation parameters */

	matd=(double **) allomat(pgmn->nbrpar,SEVEN);
	initmat(matd,pgmn->nbrpar,SEVEN);

	for (i=0;i<pgmn->nbrsta;i++)
		{
			if (datum[i]==1)
				{
					for (j=0;j<pgmn->nbrpar;j++)
						{
							if ((pgmn->parnam[j][0]=='S')&&(pgmn->parnam[j][3]=='X'))
								{
									if (check_parameter(pgmn,j,i)==1)
										{
											for (k=0;k<SEVEN;k++)
												{
													matd[j][k]=pgmn->stamat[3*i][k];
													matd[j+1][k]=pgmn->stamat[3*i+1][k];
													matd[j+2][k]=pgmn->stamat[3*i+2][k];
												}
										}
								}
						}
				}
		}
  
	/* Computation of D^T */

	tmatd=(double **) transpose(matd,pgmn->nbrpar,SEVEN);

	/* Computation of D^T.D */

	dtd=(double **) matmul(tmatd,SEVEN,pgmn->nbrpar,matd,pgmn->nbrpar,SEVEN);
	freemat(matd,pgmn->nbrpar);

	/* Computation of (D^T.D)^(-1) */

	dtdi=(double **) cholev(dtd,SEVEN);
	if (dtdi==NULL)
		{
			freemat(tmatd,SEVEN);
			freemat(dtd,SEVEN);
			return 1;
		}
	freemat(dtd,SEVEN);

	/* Computation of B=(D^T.D)^(-1).D^T */
	/* Here we follow (Wu et al., J Geophys Res, 2015) */
  
	matb=(double **) matmul(dtdi,SEVEN,SEVEN,tmatd,SEVEN,pgmn->nbrpar);
	freemat(tmatd,SEVEN);
	freemat(dtdi,SEVEN);

	mat_temp=(double **) allomat(size,pgmn->nbrpar);
	counter=-1;
	for (i=0;i<SEVEN;i++)
		{
			if (keys[i]==0)
				continue;
			counter++;
			for (j=0;j<pgmn->nbrpar;j++)
				mat_temp[counter][j]=matb[i][j];
		}
	freemat(matb,SEVEN);
  
	matb=(double **) allomat(size,pgmn->nbrpar);
	for (i=0;i<size;i++)
		for (j=0;j<pgmn->nbrpar;j++)
			matb[i][j]=mat_temp[i][j];
	freemat(mat_temp,size);
  
	/* Computation of B^T */
  
	tmatb=(double **) transpose(matb,size,pgmn->nbrpar);
  
	/* Computation of B^T.SIGMA_THETA^(-1).B */

	sigtheta=(double **) allomat(size,size);
	initmat(sigtheta,size,size);

	counter=-1;
	for (i=0;i<SEVEN;i++)
		{
			if (keys[i]==0)
				continue;
			counter++;
			switch(i)
				{
					case 0:
					case 1:
					case 2:
						sigtheta[counter][counter]=1./((TRAN)*(TRAN));
						break;
					case 3:
						sigtheta[counter][counter]=1./((SCAL)*(SCAL));
						break;
					default:
						sigtheta[counter][counter]=1./((ROTA)*(ROTA));
						break;
				}
		}

	sb=(double **) matmul(sigtheta,size,size,matb,size,pgmn->nbrpar);

	btsb=(double **) matmul(tmatb,pgmn->nbrpar,size,sb,size,pgmn->nbrpar);
	freemat(sb,size);
	freemat(tmatb,pgmn->nbrpar);

	addmat(pgmn->normat,btsb,pgmn->nbrpar,pgmn->nbrpar);

	/* Computation of the solution of the regularized normal system */

	norinv=cholev(btsb,pgmn->nbrpar);
	if (norinv==NULL)
		{
			freemat(btsb,pgmn->nbrpar);
			return 1;
		}

	solution=(double *) matvectmul(norinv,pgmn->nbrpar,pgmn->nbrpar,pgmn->norvec,pgmn->nbrpar);

	resmin=(double *) matvectmul(matb,size,pgmn->nbrpar,solution,pgmn->nbrpar);

	pgmn->sigma=pgmn->ssores;

	for (i=0;i<pgmn->nbrpar;i++)
		pgmn->sigma-=pgmn->norvec[i]*solution[i];

	for (i=0;i<size;i++)
		pgmn->sigma+=resmin[i]*resmin[i]*sigtheta[i][i];

	pgmn->sigma=(pgmn->sigma)/(pgmn->nbrobs+size-pgmn->nbrunk);

	for (i=0;i<pgmn->nbrpar;i++)
		{
			pgmn->solvec[i]=solution[i];
			for (j=0;j<pgmn->nbrpar;j++)
				pgmn->solmat[i][j]=(pgmn->sigma)*norinv[i][j];
		}

	/* Freedom */
  
	freemat(matb,size);
	freemat(sigtheta,size);
	freemat(btsb,pgmn->nbrpar);
	freemat(norinv,pgmn->nbrpar);
	freevect(solution);
  
	/* End */

	return 0;
}

/* Function to compute the stability matrix for a given datum for stations */
/* keys = array of keys for the minimum constraints (0=parameter to be excluded, 1=parameter to be considered) */
/* From (Kotsakis, J Geodesy, 2012) */

double **stability_matrix(gmn *pgmn,int *datum,int *keys)
{
	int i,j,k;
	int size;
	int counter;

	double **matd,**mate,**matb;
	double **tmatd,**tmate;
	double **dtd,**dtdi;
	double **matbtmate;
	double **mat_temp;
	double **mats;

  /* Computation of size (number of parameters to be considered) */

	size=0;
	for (i=0;i<SEVEN;i++)
		size+=keys[i];

	/* Computation of the matrix D related to the station positions for all the stations involved in the datum */
	/* Computation of the matrix E for all the stations */

	matd=(double **) allomat(3*pgmn->nbrsta,SEVEN);
	initmat(matd,3*pgmn->nbrsta,SEVEN);

	mate=(double **) allomat(size,3*pgmn->nbrsta);
	initmat(mate,size,3*pgmn->nbrsta);

	for (i=0;i<pgmn->nbrsta;i++)
		{
			counter=-1;
			for (k=0;k<SEVEN;k++)
				{
					if (keys[k]==0)
						continue;
					counter++;
					mate[counter][3*i]=pgmn->stamat[3*i][k];
					mate[counter][3*i+1]=pgmn->stamat[3*i+1][k];
					mate[counter][3*i+2]=pgmn->stamat[3*i+2][k];
				}		
			if (datum[i]==1)
				{
					for (k=0;k<SEVEN;k++)
						{
							matd[3*i][k]=pgmn->stamat[3*i][k];
							matd[3*i+1][k]=pgmn->stamat[3*i+1][k];
							matd[3*i+2][k]=pgmn->stamat[3*i+2][k];
						}
				}
		}

	/* Computation of E^T */

	tmate=(double **) transpose(mate,size,3*pgmn->nbrsta);
	freemat(mate,size);
  
	/* Computation of D^T */

	tmatd=(double **) transpose(matd,3*pgmn->nbrsta,SEVEN);

	/* Computation of D^T.D */

	dtd=(double **) matmul(tmatd,SEVEN,3*pgmn->nbrsta,matd,3*pgmn->nbrsta,SEVEN);
	freemat(matd,3*pgmn->nbrsta);

	/* Computation of (D^T.D)^(-1) */

	dtdi=(double **) cholev(dtd,SEVEN);
	if (dtdi==NULL)
		{
			freemat(tmate,3*pgmn->nbrsta);
			freemat(tmatd,SEVEN);
			freemat(dtd,SEVEN);
			return NULL;
		}
	freemat(dtd,SEVEN);

	/* Computation of B=(D^T.D)^(-1).D^T */
	/* Here we follow (Wu et al., J Geophys Res, 2015) */

	matb=(double **) matmul(dtdi,SEVEN,SEVEN,tmatd,SEVEN,3*pgmn->nbrsta);
	freemat(tmatd,SEVEN);
	freemat(dtdi,SEVEN);

	counter=-1;
	mat_temp=(double **) allomat(size,3*pgmn->nbrsta);
	for (i=0;i<SEVEN;i++)
		{
			if (keys[i]==0)
				continue;
			counter++;
			for (j=0;j<3*pgmn->nbrsta;j++)
				mat_temp[counter][j]=matb[i][j];
		}
  
	freemat(matb,SEVEN);
	matb=(double **) allomat(size,3*pgmn->nbrsta);
	for (i=0;i<size;i++)
		for (j=0;j<3*pgmn->nbrsta;j++)
			matb[i][j]=mat_temp[i][j];
	freemat(mat_temp,size);
  
	/* Computation of B.E^T */
  
	matbtmate=(double **) matmul(matb,size,3*pgmn->nbrsta,tmate,3*pgmn->nbrsta,size);
	freemat(tmate,3*pgmn->nbrsta);
  
	/* Computation of the stability matrix S=(B.E^T)^(-1) */

	mats=(double **) cholev(matbtmate,size);

	/* Freedom */
  
	freemat(matbtmate,size);

	/* End */
	
	return mats;
}

/* Function to compute the reference system effect on the station position (key=1) or EOP (key=2) basis for matrix */

double **sysreffect(gmn *pgmn,double **matrix,int key)
{
	int i,j;
	int counti,countj;
  
	int dline;
	int size;

	double **matd,**tmatd;
	double **matint,**siginv;
	double **sigd,**tdsigd;
	double **effect;

	/* Initializations */

	if (key==1)
		{
			dline=3*(pgmn->nbrsta);
			size=SEVEN;
		}
	else
		{
			dline=pgmn->nbreop;
			size=TWO;
		}
		
	/* Computation of the matrix D for the parameters involved */

	matd=(double **) allomat(dline,size);
	initmat(matd,dline,size);
      
	for (i=0;i<dline;i++)
		for (j=0;j<size;j++)
			switch(key)
				{
					case 1:
							matd[i][j]=pgmn->stamat[i][j];
							break;
					case 2:
							matd[i][j]=pgmn->prtmat[i][j];
							break;
				}
      
	/* Computation of D^T */
      
	tmatd=(double **) transpose(matd,dline,size);
      
	/* Inversion of the variance-covariance matrix of parameters of interest */ 
      
	counti=-1;
	matint=(double **) allomat(dline,dline);
      
	switch(key)
		{
			case 1:
				for (i=0;i<pgmn->nbrpar;i++)
					if (pgmn->parnam[i][0]=='S')
						{
							counti+=1;
							countj=-1;
							for (j=0;j<pgmn->nbrpar;j++)
								if (pgmn->parnam[j][0]=='S')
									{
										countj+=1;
										matint[counti][countj]=matrix[i][j];
									}
						}
				break;
			case 2:
				for (i=0;i<pgmn->nbrpar;i++)
					if ((pgmn->parnam[i][0]=='X')||(pgmn->parnam[i][0]=='Y'))
						{
							counti+=1;
							countj=-1;
							for (j=0;j<pgmn->nbrpar;j++)
								if ((pgmn->parnam[j][0]=='X')||(pgmn->parnam[j][0]=='Y'))
									{
										countj+=1;
										matint[counti][countj]=matrix[i][j];
									}
						}
					break;
		}
      
	siginv=(double **) cholev(matint,dline);
	if (siginv==NULL)
		{
			freemat(matd,dline);
			freemat(tmatd,size);
			freemat(matint,dline);
			return (double **) NULL;
		}
      
	freemat(matint,dline);

	/* Computation of D^T.SIGMA^(-1).D */
      
	sigd=(double **) matmul(siginv,dline,dline,matd,dline,size);
	freemat(matd,dline);
      
	tdsigd=(double **) matmul(tmatd,size,dline,sigd,dline,size);
	freemat(tmatd,size);
	freemat(siginv,dline);
	freemat(sigd,dline);
      
	/* Computation of the reference system effect */
      
	effect=(double **) cholev(tdsigd,size);
	if (effect==NULL)
		{
			freemat(tdsigd,size);
			return (double **) NULL;
		}
	freemat(tdsigd,size);

	/* End */

	return effect;
}

/* Function to project matrix with minimum constraints over a datum */
/* keys = array of keys for the minimum constraints (0=parameter to be excluded, 1=parameter to be considered) */
/* From (Sillard & Boucher, J Geodesy, 2001) */

double **proj_matrix(gmn *pgmn,double **matrix,int *datum,int *keys)
{
	int i,j,k;
	int size;
	int counter;

	double **matd,**matb;
	double **tmatd,**tmatb;
	double **dtd,**dtdi;
	double **sigtheta;
	double **ssolbt,**bssol;
	double **bssolbt,**bssolbti;
	double **mat_temp;
	double **mat_proj;
  
	/* Computation of size (number of parameters to be considered) */

	size=0;
	for (i=0;i<SEVEN;i++)
		size+=keys[i];

	/* Computation of the matrix D related to the station positions for all the parameters involved */

	matd=(double **) allomat(pgmn->nbrpar,SEVEN);
	initmat(matd,pgmn->nbrpar,SEVEN);

	for (i=0;i<pgmn->nbrsta;i++)
		{
			if (datum[i]==1)
				{
					for (j=0;j<pgmn->nbrpar;j++)
						{
							if ((pgmn->parnam[j][0]=='S')&&(pgmn->parnam[j][3]=='X'))
								if (check_parameter(pgmn,j,i)==1)
									{
										for (k=0;k<SEVEN;k++)
											{
												matd[j][k]=pgmn->stamat[3*i][k];
												matd[j+1][k]=pgmn->stamat[3*i+1][k];
												matd[j+2][k]=pgmn->stamat[3*i+2][k];
											}
									}
						}
				}
		}
  
	/* Computation of D^T */

	tmatd=(double **) transpose(matd,pgmn->nbrpar,SEVEN);

	/* Computation of D^T.D */

	dtd=(double **) matmul(tmatd,SEVEN,pgmn->nbrpar,matd,pgmn->nbrpar,SEVEN);
	freemat(matd,pgmn->nbrpar);

	/* Computation of (D^T.D)^(-1) */

	dtdi=(double **) cholev(dtd,SEVEN);
	if (dtdi==NULL)
		{
			freemat(tmatd,size);
			freemat(dtd,size);
			return NULL;
		}
	freemat(dtd,size);

	/* Computation of B=(D^T.D)^(-1).D^T */
	/* Here we follow (Wu et al., J Geophys Res, 2015) */

	matb=(double **) matmul(dtdi,SEVEN,SEVEN,tmatd,SEVEN,pgmn->nbrpar);
	freemat(tmatd,SEVEN);
	freemat(dtdi,SEVEN);

	mat_temp=(double **) allomat(size,pgmn->nbrpar);
	counter=-1;
	for (i=0;i<SEVEN;i++)
		{
			if (keys[i]==0)
				continue;
			counter++;
			for (j=0;j<pgmn->nbrpar;j++)
				mat_temp[counter][j]=matb[i][j];
		}
  
	freemat(matb,SEVEN);
  
	matb=(double **) allomat(size,pgmn->nbrpar);
	for (i=0;i<size;i++)
		for (j=0;j<pgmn->nbrpar;j++)
			matb[i][j]=mat_temp[i][j];
	freemat(mat_temp,size);
  
	/* Computation of B^T */
  
	tmatb=(double **) transpose(matb,size,pgmn->nbrpar);
  
	/* Projection of matrix (Sigma_sol) */
	/* Sigma_{proj}=Sigma_{sol}-Sigma_{sol}*B^T*(B*Sigma_{sol}*B^T+Sigma_{theta})^(-1)*B*Sigma_{sol} */

	/* Matrix Sigma_{theta}*/
  
	sigtheta=(double **) allomat(size,size);
	initmat(sigtheta,size,size);

	counter=-1;
	for (i=0;i<SEVEN;i++)
		{
			if (keys[i]==0)
				continue;
			counter++;
			switch(i)
				{
					case 0:
					case 1:
					case 2:
						sigtheta[counter][counter]=(TRAN)*(TRAN);
						break;
					case 3:
						sigtheta[counter][counter]=(SCAL)*(SCAL);
						break;
					default:
						sigtheta[counter][counter]=(ROTA)*(ROTA);
						break;
				}
		}

	/* Matrix B*Sigma_{sol}*B^T+Sigma_{theta} */

	ssolbt=(double **) matmul(matrix,pgmn->nbrpar,pgmn->nbrpar,tmatb,pgmn->nbrpar,size);
	bssolbt=(double **) matmul(matb,size,pgmn->nbrpar,ssolbt,pgmn->nbrpar,size);

	addmat(sigtheta,bssolbt,size,size);

	/* Matrix (B*Sigma_{sol}*B^T+Sigma_{theta})^(-1) */
  
	bssolbti=(double **) cholev(bssolbt,size);
	if (bssolbti==NULL)
		{
			freemat(matb,size);
			freemat(tmatb,pgmn->nbrpar);
			freemat(sigtheta,size);
			freemat(ssolbt,pgmn->nbrpar);
			freemat(bssolbt,size);
			return NULL;
		}
	freemat(tmatb,pgmn->nbrpar);
	freemat(sigtheta,size);
	freemat(bssolbt,size);

	/* Matrix (B*Sigma_{sol}*B^T+Sigma_{theta})^(-1)*B*Sigma_{sol} */

	bssol=(double **) matmul(matb,size,pgmn->nbrpar,matrix,pgmn->nbrpar,pgmn->nbrpar);
	freemat(matb,size);

	mat_temp=(double **) matmul(bssolbti,size,size,bssol,size,pgmn->nbrpar);
	freemat(bssolbti,size);
	freemat(bssol,size);

	/* Matrix Sigma_{sol}*B^T*(B*Sigma_{sol}*B^T+Sigma_{theta})^(-1)*B*Sigma_{sol} */

	mat_proj=(double **) matmul(ssolbt,pgmn->nbrpar,size,mat_temp,size,pgmn->nbrpar);
	freemat(ssolbt,pgmn->nbrpar);
	freemat(mat_temp,size);

	/* Projection of matrix */

	submat(matrix,mat_proj,pgmn->nbrpar,pgmn->nbrpar);

	/* End */
  
	return mat_proj;
}

/* Function to compute the Helmert transformation parameters on the station position (key=1) or EOP (key=2) basis */
/* with (key_apr=1) or without (key_apr=0) the a priori variance-covariance matrix */

int transfo(gmn *pgmn,double *solution,double **variance,double *partra,double **vartra,int key,int key_apr)
{
	int i,j;
	int counti,countj;
  
	int dline;
	int size;
	
	double sigma0;

	double *vecint,*vecint2,*vecint3,*vecint4;
	double *res;

	double **matd,**tmatd;
	double **matint,**siginv;
	double **sigd,**tdsigd;
	double **effect;

	/* Initializations */

	switch(key)
		{
			case 1:
				dline=3*(pgmn->nbrsta);
				size=SEVEN;
				break;
			case 2:
				dline=pgmn->nbreop;
				size=TWO;
				break;
		}
		
	initvect(partra,size);
	initmat(vartra,size,size);
      
	/* Computation of the matrix D for the transformation parameters involved */

	matd=(double **) allomat(dline,size);
	initmat(matd,dline,size);
      
	for (i=0;i<dline;i++)
		for (j=0;j<size;j++)
			switch(key)
				{
					case 1:
						matd[i][j]=pgmn->stamat[i][j];
						break;
					case 2:
						matd[i][j]=pgmn->prtmat[i][j];
						break;
				}
			
	/* Computation of D^T */
      
	tmatd=(double **) transpose(matd,dline,size);
      
	/* Inversion of the variance-covariance matrix for the parameters of interest */ 

	counti=-1;
	matint=(double **) allomat(dline,dline);
	vecint=(double *) allovect(dline);

	switch(key)
		{
			case 1:
				for (i=0;i<pgmn->nbrpar;i++)
					if (pgmn->parnam[i][0]=='S')
						{
							counti+=1;
							vecint[counti]=-solution[i];
							countj=-1;
							for (j=0;j<pgmn->nbrpar;j++)
								if (pgmn->parnam[j][0]=='S')
									{
										countj+=1;
										if (key_apr==0)
											matint[counti][countj]=variance[i][j];
										else
											matint[counti][countj]=variance[i][j]+pgmn->aprmat[i][j];
									}
						}
				break;
			case 2:
				for (i=0;i<pgmn->nbrpar;i++)
					if ((pgmn->parnam[i][0]=='X')||(pgmn->parnam[i][0]=='Y'))
						{
							counti+=1;
							vecint[counti]=-solution[i];
							countj=-1;
							for (j=0;j<pgmn->nbrpar;j++)
								if ((pgmn->parnam[j][0]=='X')||(pgmn->parnam[j][0]=='Y'))
									{
										countj+=1;
										matint[counti][countj]=variance[i][j];
										if (key_apr==0)
											matint[counti][countj]=variance[i][j];
										else
											matint[counti][countj]=variance[i][j]+pgmn->aprmat[i][j];
									}
						}
				break;
		}
					
	siginv=(double **) cholev(matint,dline);
	if (siginv==NULL)
		{
			freevect(vecint);
			freemat(matd,dline);
			freemat(tmatd,SEVEN);
			freemat(matint,dline);
			return 1;
		}
	freemat(matint,dline);

	/* Computation of D^T.SIGMA^(-1).SOL */

	vecint2=(double *) matvectmul(siginv,dline,dline,vecint,dline);

	vecint3=(double *) matvectmul(tmatd,size,dline,vecint2,dline);
	freevect(vecint2);

	/* Computation of D^T.SIGMA^(-1).D */
      
	sigd=(double **) matmul(siginv,dline,dline,matd,dline,size);
      
	tdsigd=(double **) matmul(tmatd,size,dline,sigd,dline,size);
	freemat(tmatd,size);
	freemat(sigd,dline);
      
	/* Computation of the reference system effect */
      
	effect=(double **) cholev(tdsigd,size);
	if (effect==NULL)
		{
			freevect(vecint3);
			freemat(tdsigd,size);
			return 1;
		}
	freemat(tdsigd,size);

	/* Computation of the parameters */

	vecint4=(double *) matvectmul(effect,size,size,vecint3,size);
	for (i=0;i<size;i++)
		partra[i]=vecint4[i];
		
	freevect(vecint3);
	freevect(vecint4);
			
	/* Computation of the residuals */
			
	res=(double *) matvectmul(matd,dline,size,partra,size);
	subvect(vecint,res,dline);
			
	/* Computation of the unit variance factor */
			
	vecint2=(double *) matvectmul(siginv,dline,dline,res,dline);
			
	sigma0=0.;
	for (i=0;i<dline;i++)
		sigma0+=res[i]*vecint2[i];
	sigma0/=((double) (dline-size));
			
	/* Computation of the variance-covariance matrix */
			
	scalemat(effect,sigma0,size,size);
			
	for (i=0;i<size;i++)
		for (j=0;j<size;j++)
			vartra[i][j]=effect[i][j];

	/* Freedom */

	freevect(vecint);
	freevect(vecint2);
	freevect(res);
	freemat(matd,dline);
	freemat(siginv,dline);
	freemat(effect,size);
	
	/* End */

	return 0;
}

/* Function to compute the Helmert transformation for a given datum of stations and all the transformed parameters */
/* with (key_apr=1) or without (key_apr=0) the a priori variance-covariance matrix for station positions */

int transfo_datum(gmn *pgmn,int *datum,double *solution,double **variance,double *partra,double **vartra,double *soltransfo,double **vartransfo,int key_apr)
{
	int i,j,k;
	int size_datum;
	int size_no_datum;
	int count1,count2,count3;
	int counti;

	int *ind_dX;
	int *ind_dZ;
	int *ind_dEOP;

	double *vec_temp;
	double *deltaX;
	double *deltaZ;
	double *deltaEOP;
	double *deltaX0;
	double *deltaZ0;
	double *deltaEOP0;
	double *PdX;
	double *GTPdX;

	double **mat_temp;
	double **mat_temp1;
	double **mat_temp2;
	double **G;
	double **Gtilde;
	double **GtildeT;
	double **Ftilde;
	double **FtildeT;
	double **SigmadX;
	double **SigmadZ;
	double **SigmadEOP;
	double **SigmadXdZ;
	double **SigmadXdEOP;
	double **P;
	double **Pinv;
	double **GT;
	double **GTP;
	double **PG;
	double **GTPG;
	double **SigthetaGT;
	double **SigmadX0;
	double **SigmadZ0;
	double **SigmadEOP0;
	double **SigmadXdZ0;
	double **hatM;
	double **hatMT;

	/* !!! Only the useful (for computation of reference system effects) parts !!! */
	/* !!! of the variance covariance matrix of the transformed parameters are computed !!! */

	/* Number of stations in the datum */

	size_datum=0;

	for (i=0;i<pgmn->nbrsta;i++)
		size_datum+=datum[i];

	size_no_datum=pgmn->nbrsta-size_datum;

	/* Indices */
	/* ind_dX = indices for stations in the datum */
	/* ind_dZ = indices for stations absent from the datum */
	/* ind_dEOP = indices for EOP */

	count1=-1;
	count2=-1;
	count3=-1;

	ind_dX=(int *) iallovect(3*size_datum);
	initivect(ind_dX,3*size_datum);

	if (size_no_datum>0)
		{
			ind_dZ=(int *) iallovect(3*size_no_datum);
			initivect(ind_dZ,3*size_no_datum);
		}

	ind_dEOP=(int *) iallovect(pgmn->nbreop);
	initivect(ind_dEOP,pgmn->nbreop);

	for (j=0;j<pgmn->nbrpar;j++)
		if ((pgmn->parnam[j][0]=='S')&&(pgmn->parnam[j][3]=='X'))
			{
				for (i=0;i<pgmn->nbrsta;i++)
					if (check_parameter(pgmn,j,i)==1)
						switch (datum[i])
							{
								case 0:
									count2++;
									ind_dZ[count2]=j;
									count2++;
									ind_dZ[count2]=j+1;
									count2++;
									ind_dZ[count2]=j+2;
									break;
								case 1:
									count1++;
									ind_dX[count1]=j;
									count1++;
									ind_dX[count1]=j+1;
									count1++;
									ind_dX[count1]=j+2;
									break;
							}
			}
		else if ((pgmn->parnam[j][0]=='X')||(pgmn->parnam[j][0]=='Y'))
			{
				count3++;
				ind_dEOP[count3]=j;
			}
  
	/* Estimation of the transformation parameters */
	
	/* Observation system */
  
	deltaX=(double *) allovect(3*size_datum);
	initvect(deltaX,3*size_datum);
  
	G=(double **) allomat(3*size_datum,SEVEN);
	initmat(G,3*size_datum,SEVEN);
  
	SigmadX=(double **) allomat(3*size_datum,3*size_datum);
	initmat(SigmadX,3*size_datum,3*size_datum);

	Pinv=(double **) allomat(3*size_datum,3*size_datum);
	initmat(Pinv,3*size_datum,3*size_datum);

	for (i=0;i<3*size_datum;i++)
		{
			deltaX[i]=-solution[ind_dX[i]];
			for (j=0;j<pgmn->nbrsta;j++)
				{
					if (check_parameter(pgmn,ind_dX[i],j)==0)
						continue;
					if (pgmn->parnam[ind_dX[i]][3]=='X')
						counti=3*j;
					if (pgmn->parnam[ind_dX[i]][3]=='Y')
						counti=3*j+1;
					if (pgmn->parnam[ind_dX[i]][3]=='Z')
						counti=3*j+2;
					for (k=0;k<SEVEN;k++)
						G[i][k]=pgmn->stamat[counti][k];
					break;
				}
			for (j=0;j<3*size_datum;j++)
				{
					SigmadX[i][j]=variance[ind_dX[i]][ind_dX[j]];
					if (key_apr==0)
						Pinv[i][j]=SigmadX[i][j];
					else
						Pinv[i][j]=SigmadX[i][j]+pgmn->aprmat[ind_dX[i]][ind_dX[j]];
				}
		}

	/* Weight matrix */

	P=(double **) cholev(Pinv,3*size_datum);
	freemat(Pinv,3*size_datum);
	if (P==NULL)
		{
			ifreevect(ind_dX);
			if (size_no_datum>0)
				ifreevect(ind_dZ);
			ifreevect(ind_dEOP);
			freevect(deltaX);
			freemat(G,3*size_datum);
			freemat(SigmadX,3*size_datum);
			return 1;
		}
		
	/* Least-squares */

	GT=(double **) transpose(G,3*size_datum,SEVEN);
	PG=(double **) matmul(P,3*size_datum,3*size_datum,G,3*size_datum,SEVEN);
	GTPG=(double **) matmul(GT,SEVEN,3*size_datum,PG,3*size_datum,SEVEN);
	GTP=(double **) matmul(GT,SEVEN,3*size_datum,P,3*size_datum,3*size_datum);
	PdX=(double *) matvectmul(P,3*size_datum,3*size_datum,deltaX,3*size_datum);
	GTPdX=(double *) matvectmul(GT,SEVEN,3*size_datum,PdX,3*size_datum);

	mat_temp=(double **) cholev(GTPG,SEVEN);
	if (mat_temp==NULL)
		{
			ifreevect(ind_dX);
			if (size_no_datum>0)
				ifreevect(ind_dZ);
			ifreevect(ind_dEOP);
			freevect(deltaX);
			freevect(PdX);
			freevect(GTPdX);
			freemat(G,3*size_datum);
			freemat(SigmadX,3*size_datum);
			freemat(P,3*size_datum);
			freemat(GT,SEVEN);
			freemat(PG,3*size_datum);
			freemat(GTPG,SEVEN);
			return 1;
		}

	hatM=(double **) matmul(mat_temp,SEVEN,SEVEN,GTP,SEVEN,3*size_datum);
	hatMT=(double **) transpose(hatM,SEVEN,3*size_datum);
	freemat(GTP,SEVEN);

	/* Vector solution */

	vec_temp=(double *) matvectmul(mat_temp,SEVEN,SEVEN,GTPdX,SEVEN);

	for (i=0;i<SEVEN;i++)
		partra[i]=vec_temp[i];

	freevect(vec_temp);

	/* Vector of transformed stations */

	deltaX0=(double *) matvectmul(G,3*size_datum,SEVEN,partra,SEVEN);
	subvect(deltaX,deltaX0,3*size_datum);
	scalevect(deltaX0,-1.,3*size_datum);

	/* Variance-covariance matrix of the solution */

	for (i=0;i<SEVEN;i++)
		for (j=0;j<SEVEN;j++)
			vartra[i][j]=mat_temp[i][j];

	freemat(mat_temp,SEVEN);
	
	/* Variance-covariance matrix of the transformed stations */
	
	SigthetaGT=(double **) matmul(vartra,SEVEN,SEVEN,GT,SEVEN,3*size_datum);
	SigmadX0=(double **) matmul(G,3*size_datum,SEVEN,SigthetaGT,SEVEN,3*size_datum);
	addmat(SigmadX,SigmadX0,3*size_datum,3*size_datum);
		
	mat_temp1=(double **) matmul(hatM,SEVEN,3*size_datum,SigmadX,3*size_datum,3*size_datum);
	
	mat_temp2=(double **) matmul(G,3*size_datum,SEVEN,mat_temp1,SEVEN,3*size_datum);
	freemat(mat_temp1,SEVEN);
	
	scalemat(mat_temp2,-1.,3*size_datum,3*size_datum);
	mat_temp1=(double **) transpose(mat_temp2,3*size_datum,3*size_datum);
	
	addmat(mat_temp1,SigmadX0,3*size_datum,3*size_datum);
	addmat(mat_temp2,SigmadX0,3*size_datum,3*size_datum);
	
	freemat(mat_temp1,3*size_datum);
	freemat(mat_temp2,3*size_datum);

	/* Transformed parameters */

	/* Stations absent from the datum */

	if (size_no_datum>0)
		{
			deltaZ=(double *) allovect(3*size_no_datum);
			initvect(deltaZ,3*size_no_datum);
  
			Gtilde=(double **) allomat(3*size_no_datum,SEVEN);
			initmat(Gtilde,3*size_no_datum,SEVEN);
  
			SigmadZ=(double **) allomat(3*size_no_datum,3*size_no_datum);
			initmat(SigmadZ,3*size_no_datum,3*size_no_datum);

			for (i=0;i<3*size_no_datum;i++)
				{
					deltaZ[i]=solution[ind_dZ[i]];
					for (j=0;j<pgmn->nbrsta;j++)
						{
							if (check_parameter(pgmn,ind_dZ[i],j)==0)
								continue;
							if (pgmn->parnam[ind_dZ[i]][3]=='X')
								counti=3*j;
							if (pgmn->parnam[ind_dZ[i]][3]=='Y')
								counti=3*j+1;
							if (pgmn->parnam[ind_dZ[i]][3]=='Z')
								counti=3*j+2;
							for (k=0;k<SEVEN;k++)
								Gtilde[i][k]=pgmn->stamat[counti][k];
							break;
						}
					for (j=0;j<3*size_no_datum;j++)
						SigmadZ[i][j]=variance[ind_dZ[i]][ind_dZ[j]];
				}

			GtildeT=(double **) transpose(Gtilde,3*size_no_datum,SEVEN);

			deltaZ0=(double *) matvectmul(Gtilde,3*size_no_datum,SEVEN,partra,SEVEN);
			addvect(deltaZ,deltaZ0,3*size_no_datum);

			SigmadXdZ=(double **) allomat(3*size_datum,3*size_no_datum);
			initmat(SigmadXdZ,3*size_datum,3*size_no_datum);
	
			for (i=0;i<3*size_datum;i++)
				for (j=0;j<3*size_no_datum;j++)
					SigmadXdZ[i][j]=variance[ind_dX[i]][ind_dZ[j]];

			mat_temp1=(double **) matmul(vartra,SEVEN,SEVEN,GtildeT,SEVEN,3*size_no_datum);
	
			SigmadZ0=(double **) matmul(Gtilde,3*size_no_datum,SEVEN,mat_temp1,SEVEN,3*size_no_datum);
			freemat(mat_temp1,SEVEN);
	
			addmat(SigmadZ,SigmadZ0,3*size_no_datum,3*size_no_datum);
	
			mat_temp1=(double **) matmul(hatM,SEVEN,3*size_datum,SigmadXdZ,3*size_datum,3*size_no_datum);
	
			mat_temp2=(double **) matmul(Gtilde,3*size_no_datum,SEVEN,mat_temp1,SEVEN,3*size_no_datum);
			freemat(mat_temp1,SEVEN);
			freemat(Gtilde,3*size_no_datum);
	
			scalemat(mat_temp2,-1.,3*size_no_datum,3*size_no_datum);
			mat_temp1=(double **) transpose(mat_temp2,3*size_no_datum,3*size_no_datum);
	
			addmat(mat_temp1,SigmadZ0,3*size_no_datum,3*size_no_datum);
			addmat(mat_temp2,SigmadZ0,3*size_no_datum,3*size_no_datum);
	
			freemat(mat_temp1,3*size_no_datum);
			freemat(mat_temp2,3*size_no_datum);

			mat_temp1=(double **) matmul(vartra,SEVEN,SEVEN,GtildeT,SEVEN,3*size_no_datum);
			SigmadXdZ0=(double **) matmul(G,3*size_datum,SEVEN,mat_temp1,SEVEN,3*size_no_datum);
			freemat(mat_temp1,SEVEN);
			
			addmat(SigmadXdZ,SigmadXdZ0,3*size_datum,3*size_no_datum);
			
			mat_temp1=(double **) matmul(hatM,SEVEN,3*size_datum,SigmadXdZ,3*size_datum,3*size_no_datum);
			freemat(SigmadXdZ,3*size_datum);

			mat_temp2=(double **) matmul(G,3*size_datum,SEVEN,mat_temp1,SEVEN,3*size_no_datum);
			freemat(mat_temp1,SEVEN);
			
			scalemat(mat_temp2,-1.,3*size_datum,3*size_no_datum);
			addmat(mat_temp2,SigmadXdZ0,3*size_datum,3*size_no_datum);
			freemat(mat_temp2,3*size_datum);
			
			mat_temp1=(double **) matmul(hatMT,3*size_datum,SEVEN,GtildeT,SEVEN,3*size_no_datum);
			freemat(GtildeT,SEVEN);

			mat_temp2=(double **) matmul(SigmadX,3*size_datum,3*size_datum,mat_temp1,3*size_datum,3*size_no_datum);
			scalemat(mat_temp2,-1.,3*size_datum,3*size_no_datum);
			addmat(mat_temp2,SigmadXdZ0,3*size_datum,3*size_no_datum);
			freemat(mat_temp2,3*size_datum);
		}

	/* EOP */

	deltaEOP=(double *) allovect(pgmn->nbreop);
	initvect(deltaEOP,pgmn->nbreop);
  
	Ftilde=(double **) allomat(pgmn->nbreop,SEVEN);
	initmat(Ftilde,pgmn->nbreop,SEVEN);
 
	SigmadEOP=(double **) allomat(pgmn->nbreop,pgmn->nbreop);
	initmat(SigmadEOP,pgmn->nbreop,pgmn->nbreop);

	for (i=0;i<pgmn->nbreop;i++)
		{
			deltaEOP[i]=solution[ind_dEOP[i]];
			for (j=0;j<pgmn->nbreop;j++)
				SigmadEOP[i][j]=variance[ind_dEOP[i]][ind_dEOP[j]];
		}

	counti=-1;
	for (i=0;i<pgmn->nbrpar;i++)
		{
			if ((pgmn->parnam[i][0]!='X')&&(pgmn->parnam[i][0]!='Y'))
				continue;
			counti++;
			Ftilde[counti][4]=pgmn->prtmat[i][0];
			Ftilde[counti][5]=pgmn->prtmat[i][1];
		}
	
	FtildeT=(double **) transpose(Ftilde,pgmn->nbreop,SEVEN);
  
	deltaEOP0=(double *) matvectmul(Ftilde,pgmn->nbreop,SEVEN,partra,SEVEN);
	addvect(deltaEOP,deltaEOP0,pgmn->nbreop);

	SigmadXdEOP=(double **) allomat(3*size_datum,pgmn->nbreop);
	initmat(SigmadXdEOP,3*size_datum,pgmn->nbreop);
	
	for (i=0;i<3*size_datum;i++)
		for (j=0;j<pgmn->nbreop;j++)
			SigmadXdEOP[i][j]=variance[ind_dX[i]][ind_dEOP[j]];

	mat_temp1=(double **) matmul(vartra,SEVEN,SEVEN,FtildeT,SEVEN,pgmn->nbreop);
	freemat(FtildeT,SEVEN);
	
	SigmadEOP0=(double **) matmul(Ftilde,pgmn->nbreop,SEVEN,mat_temp1,SEVEN,pgmn->nbreop);
	freemat(mat_temp1,SEVEN);
	
	addmat(SigmadEOP,SigmadEOP0,pgmn->nbreop,pgmn->nbreop);
	
	mat_temp1=(double **) matmul(hatM,SEVEN,3*size_datum,SigmadXdEOP,3*size_datum,pgmn->nbreop);
	freemat(SigmadXdEOP,3*size_datum);
	
	mat_temp2=(double **) matmul(Ftilde,pgmn->nbreop,SEVEN,mat_temp1,SEVEN,pgmn->nbreop);
	freemat(mat_temp1,SEVEN);
	freemat(Ftilde,pgmn->nbreop);
	
	scalemat(mat_temp2,-1.,pgmn->nbreop,pgmn->nbreop);
	mat_temp1=(double **) transpose(mat_temp2,pgmn->nbreop,pgmn->nbreop);
	
	addmat(mat_temp1,SigmadEOP0,pgmn->nbreop,pgmn->nbreop);
	addmat(mat_temp2,SigmadEOP0,pgmn->nbreop,pgmn->nbreop);
	
	freemat(mat_temp1,pgmn->nbreop);
	freemat(mat_temp2,pgmn->nbreop);
	
	/* All parameters */
	
	initvect(soltransfo,pgmn->nbrpar);
	initmat(vartransfo,pgmn->nbrpar,pgmn->nbrpar);

	for (i=0;i<3*size_datum;i++)
		{
			soltransfo[ind_dX[i]]=deltaX0[i];
			for (j=0;j<3*size_datum;j++)
				vartransfo[ind_dX[i]][ind_dX[j]]=SigmadX0[i][j];
		}
		
	if (size_no_datum>0)
		{
			for (i=0;i<3*size_no_datum;i++)
				{
					soltransfo[ind_dZ[i]]=deltaZ0[i];
					for (j=0;j<3*size_no_datum;j++)
						vartransfo[ind_dZ[i]][ind_dZ[j]]=SigmadZ0[i][j];
				}
			for (i=0;i<3*size_datum;i++)
				for (j=0;j<3*size_no_datum;j++)
					{
						vartransfo[ind_dX[i]][ind_dZ[j]]=SigmadXdZ0[i][j];
						vartransfo[ind_dZ[j]][ind_dX[i]]=vartransfo[ind_dX[i]][ind_dZ[j]];
					}
		}
			
	for (i=0;i<pgmn->nbreop;i++)
		{
			soltransfo[ind_dEOP[i]]=deltaEOP0[i];
			for (j=0;j<pgmn->nbreop;j++)
				vartransfo[ind_dEOP[i]][ind_dEOP[j]]=SigmadEOP0[i][j];
		}
		
	/* Freedom */

	ifreevect(ind_dX);
	if (size_no_datum>0)
		ifreevect(ind_dZ);
	ifreevect(ind_dEOP);

	freevect(deltaX);
	freevect(PdX);
	freevect(GTPdX);
	freevect(deltaX0);
	if (size_no_datum>0)
		{
			freevect(deltaZ);
			freevect(deltaZ0);
		}
	freevect(deltaEOP);	
	freevect(deltaEOP0);
		
	freemat(G,3*size_datum);
	freemat(SigmadX,3*size_datum);
	freemat(P,3*size_datum);
	freemat(GT,SEVEN);
	freemat(PG,3*size_datum);
	freemat(GTPG,SEVEN);
	freemat(SigthetaGT,SEVEN);
	freemat(SigmadX0,3*size_datum);
	freemat(hatM,SEVEN);
	freemat(hatMT,3*size_datum);
	if (size_no_datum>0)
		{
			freemat(SigmadZ,3*size_no_datum);
			freemat(SigmadZ0,3*size_no_datum);
			freemat(SigmadXdZ0,3*size_no_datum);
		}
	freemat(SigmadEOP,pgmn->nbreop);	
	freemat(SigmadEOP0,pgmn->nbreop);	
	
	return 0;
}
/* Function to dynamically allocate memory for an integer vector with l rows */
/* The vector is allocated by this function */
/* The vector must be freed ouside this function */

int *iallovect(int l)
{
	int *vector;

	/* Dynamic allocation of array */

	vector=(int *) malloc((size_t) l*sizeof(int));

	if (vector == NULL)
		{
			//printf("Allocation failure in iallovect()\n");
			exit(EXIT_FAILURE);
		}

	return vector;
}

/* Routine to free the memory allocated by iallovect */

void ifreevect(int *vector)
{
	free(vector);
	vector=NULL;
}

/* Function to dynamically allocate memory for a double vector with l rows */
/* The vector is allocated by this function */
/* The vector must be freed ouside this function */

double *allovect(int l)
{
	double *vector;

	/* Dynamic allocation of array */

	vector=(double *) malloc((size_t) l*sizeof(double));
	if (vector == NULL)
		{
			//printf("Allocation failure in allovect()\n");
			exit(EXIT_FAILURE);
		}

	return vector;
}

/* Routine to free the memory allocated by allovect */

void freevect(double *vector)
{
	free(vector);
	vector=NULL;
}

/* Function to dynamically allocate memory for a double matrix with l rows and c columns */
/* The matrix is allocated by this function */
/* The matrix must be freed ouside this function */

double **allomat(int l,int c)
{
	int i;

	double **matrix;

	/* Dynamic allocation of pointers to rows */

	matrix=(double **) malloc((size_t) (l*sizeof(double *)));

	if (matrix == NULL) 
		{
			//printf("Allocation failure in allomat()\n");
			exit(EXIT_FAILURE);
		}

	/* Allocation of rows and setting of pointers to them */

	for (i=0;i<l;i++)
		{
			matrix[i]=(double *) malloc((size_t) c*sizeof(double));
			if (matrix[i] == NULL)
				{
					//printf("Allocation failure in allomat()\n");
					exit(EXIT_FAILURE);
				}
		}
		
	return matrix;
}

/* Routine to free the memory allocated by allomat */

void freemat(double **matrix,int l)
{
	int i;

	for (i=0;i<l;i++) 
		{
			free(matrix[i]);
			matrix[i]=NULL;
		}

	free(matrix);
	matrix=NULL;
}

/* Routine to initialize an integer vector of dimension l */
/* The vector must be allocated outside this function */
/* The vector is directly modified */

void initivect(int *vector,int l)
{
	int i;

	for (i=0;i<l;i++)
		vector[i]=0;
}

/* Routine to initialize a double vector of dimension l */
/* The vector must be allocated outside this function */
/* The vector is directly modified */

void initvect(double *vector,int l) 
{
	int i;

	for (i=0;i<l;i++)
		vector[i]=0.;
}

/* Routine to initialize a double matrix of dimension (l,c) */
/* The matrix must be allocated outside this function */
/* The matrix is directly modified */

void initmat(double **matrix,int l,int c)
{
	int i,j;

	for (i=0;i<l;i++)
		for (j=0;j<c;j++)
			matrix[i][j]=0.;  
}
 
/* Routine to scale the matrix with the double factor */
/* The matrix must be allocated outside this function */
/* The matrix is directly modified */

void scalemat(double **matrix,double factor,int l,int c)
{
	int i,j;

	double temp;

	for (i=0;i<l;i++)
		for (j=0;j<c;j++)
			{
				temp=matrix[i][j];
				temp *= factor;
				matrix[i][j]=temp;
			}
}

/* Routine to scale the vector with the double factor */
/* The vector must be allocated outside this function */
/* The vector is directly modified */

void scalevect(double *vector,double factor,int l)
{
	int i;

	double temp;

	for (i=0;i<l;i++)
		{
			temp=vector[i];
			temp *= factor;
			vector[i]=temp;
		}
}

/* Routine to print the matrix on screen */
/* The matrix must be allocated outside this function */

void printmat(double **matrix,int l,int c)
{
	int i,j;

	for (i=0;i<l;i++)
		for (j=0;j<c;j++)
			printf("%d %d %f\n",i+1,j+1,matrix[i][j]);
  
	printf("\n");
}

/* Routine to print the vector on screen */
/* The vector must be allocated outside this function */

void printvect(double *vector,int l)
{
	int i;

	for (i=0;i<l;i++)
    printf("%d %f\n",i+1,vector[i]);
  
	printf("\n");
}

/* Routine to add the two matrices mat1 and mat2 */
/* The matrices must be allocated outside this function */
/* The matrix mat2 contains the result obtained */

void addmat(double **mat1,double **mat2,int l,int c)
{
	int i,j;

	double temp;

	for (i=0;i<l;i++)
		for (j=0;j<c;j++)
			{
				temp=mat1[i][j];
				temp+=mat2[i][j];
				mat2[i][j]=temp;
			}
}

/* Routine to add the two vectors vect1 and vect2 */
/* The vectors must be allocated outside this function */
/* The vector vect2 contains the result obtained */

void addvect(double *vect1,double *vect2,int l)
{
	int i;

	double temp;

	for (i=0;i<l;i++)
		{
			temp=vect1[i];
			temp+=vect2[i];
			vect2[i]=temp;
		}
}

/* Routine to substract the matrix mat2 from mat1 */
/* The matrices must be allocated outside this function */
/* The matrix mat2 contains the result obtained */

void submat(double **mat1,double **mat2,int l,int c)
{
	double factor=-1.;

	scalemat(mat2,factor,l,c);
	addmat(mat1,mat2,l,c);
}

/* Routine to substract the vector vect2 from vect1 */
/* The vectors must be allocated outside this function */
/* The vector vect2 contains the result obtained */

void subvect(double *vect1,double *vect2,int l)
{
	double factor=-1.;

	scalevect(vect2,factor,l);
	addvect(vect1,vect2,l);
}

/* Function to multiply the matrices mat1 and mat2 */
/* The matrices must be allocated outside this function */
/* The matrix containing the result is allocated in this function */
/* This matrix must be freed outside this function */

double **matmul(double **mat1,int l1,int c1,double **mat2,int l2,int c2)
{
	int i,j,k;

	double **matrix;

	if (c1!=l2) 
		{
			//printf("Dimension problem in matmul!\n");
			return (double **) NULL;
		}

	matrix=allomat(l1,c2);
	initmat(matrix,l1,c2);

	for (i=0;i<l1;i++)
		for (j=0;j<c2;j++)
			for (k=0;k<c1;k++)
				matrix[i][j]+=mat1[i][k]*mat2[k][j];

	return matrix;
}

/* Function to multiply matrix by vector */
/* The matrix and the vector must be allocated outside this function */
/* The vector containing the result is allocated in this function */
/* This vector must be freed outside this function */

double *matvectmul(double **matrix,int l,int c,double *vector,int n) 
{
	int i,k;

	double *vect;

	if (c!=n) 
		{
			//printf("Dimension problem in matvectmul!\n");
			return (double *) NULL;
		}

	vect=allovect(l);
	initvect(vect,l);

	for (i=0;i<l;i++)
		for (k=0;k<c;k++)
			vect[i]+=matrix[i][k]*vector[k];

	return vect;
}

/* Function to transpose the matrix */
/* The matrix must be allocated outside this function */
/* The matrix containing the result is allocated in this function */
/* This matrix must be freed outside this function */

double **transpose(double **matrix,int l,int c)
{
	int i,j;

	double **tmatrix;

	tmatrix=allomat(c,l);
	initmat(tmatrix,c,l);

	for (i=0;i<c;i++)
		for (j=0;j<l;j++)
			tmatrix[i][j]=matrix[j][i];

	return tmatrix;
}

/* Function to compute the infinite norm of the matrix */
/* The matrix must be allocated outside this function */

double infnorm(double **matrix,int n)
{
	int i,j;

	double sum;
	double norm;

	norm=0.;

	for (i=0;i<n;i++)
		{
			sum=0;
			for (j=0;j<n;j++)
				sum+=fabs(matrix[i][j]);
			if (sum>norm) 
				norm=sum;
		}
		
	return norm;
}

/* Function to compute the root of the matrix by Cholevsky's method */
/* The matrix must be allocated outside this function */
/* The matrix containing the result is allocated in this function */
/* This matrix must be freed outside this function */

double **rootmat(double **matrix,int n)
{
	int j,k,l;

	double temp;

	double **matl;

	matl=allomat(n,n);
	initmat(matl,n,n);

	/* Computation of matl such as matrix=matl*tmatl */

	for (l=0;l<n;l++)
		{
			for (j=0;j<=l-1;j++)
				{
					temp=matrix[l][j];
					for (k=0;k<=j-1;k++)
						temp-=matl[l][k]*matl[j][k];
					temp=temp/matl[j][j];
					matl[l][j]=temp;
				}

			temp=matrix[l][l];
			for (k=0;k<=l-1;k++)
				temp-=matl[l][k]*matl[l][k];
      
			if (temp<=TINY)
				{
					//printf("Matrix singular in rootmat()\n");
					freemat(matl,n);
					return (double **) NULL;
				}
      
			matl[l][l]=sqrt(temp);
		}

	return matl;
}

/* Function to compute the inverse of the matrix by Cholevsky's method */
/* The matrix must be allocated outside this function */
/* The matrix containing the result is allocated in this function */
/* This matrix must be freed outside this function */

double **cholev(double **matrix,int n)
{
	int i,j,k;

	double temp;

	double **matl,**matinv;

	/* Computation of matl such as matrix=matl*tmatl */

	matl=rootmat(matrix,n);
	if (matl==NULL) 
		{
			//printf("Matrix non invertible in cholev()\n");
			return (double **) NULL;
		}

	/* Computation of matl^(-1) */

	for (i=0;i<n;i++)
		{
			for (j=0;j<=i-1;j++)
				{
					temp=-matl[i][j]*matl[j][j];
					for (k=j+1;k<=i-1;k++)
						temp-=matl[i][k]*matl[k][j];
					temp=temp/matl[i][i];
					matl[i][j]=temp;
				}
			temp=matl[i][i];
			matl[i][i]=1./temp;
		}

	/* Computation of matrix^(-1) with matrix^(-1)=tmatl*matl */

	matinv=allomat(n,n);
	initmat(matinv,n,n);

	for (i=0;i<n;i++)
		{
			for (j=0;j<=i;j++)
				{
					for (k=i;k<n;k++)
						matinv[i][j]+=matl[k][i]*matl[k][j];
					matinv[j][i]=matinv[i][j];
				}
		}

	freemat(matl,n);
	
	return matinv;
}

/* Function to compute the determinant of the matrix */
/* The matrix must be allocated outside this function */

double det(double **matrix,int n)
{
	int i;

	double d;

	double **root;

	root=rootmat(matrix,n);

	if (root==NULL)
		{
			//printf("Singular matrix in det()\n");
			return 0.;
		}

	d=1.;

	for (i=0;i<n;i++)
		d*=root[i][i]*root[i][i];

	freemat(root,n);

	return d;
}

/* Function to compute the condition number of the matrix */
/* The matrix must be allocated outside this function */

double cond(double **matrix,int n)
{
	double condnumber;

	double **matinv;

	matinv=cholev(matrix,n);

	if (matinv==NULL)
		{
			//printf("Singular matrix in cond()\n");
			return 1.e16;
		}

	condnumber=infnorm(matrix,n)*infnorm(matinv,n);

	freemat(matinv,n);

	return condnumber;
}

/* Function to dynamically allocate memory for a char matrix with l rows and c columns */
/* The matrix is allocated by this function */
/* The matrix must be freed ouside this function */

char **allocharmat(int l,int c)
{
	int i;

	char **matrix;

	/* Dynamic allocation of pointers to rows */

	matrix=(char **) malloc((size_t) (l*sizeof(char *)));

	if (matrix == NULL) 
		{
			//printf("Allocation failure in allocharmat()\n");
			exit(EXIT_FAILURE);
		}

	/* Allocation of rows and setting of pointers to them */

	for (i=0;i<l;i++)
		{
			matrix[i]=(char *) malloc((size_t) c*sizeof(char));
			if (matrix[i] == NULL)
				{
					//printf("Allocation failure in allocharmat()\n");
					exit(EXIT_FAILURE);
				}
		}
		
	return matrix;
}

/* Routine to free the memory allocated by allocharmat */

void freecharmat(char **matrix,int l)
{
	int i;

	for (i=0;i<l;i++) 
		{
			free(matrix[i]);
			matrix[i]=NULL;
		}

	free(matrix);
	matrix=NULL;
}

/* Routine to print the vector of strings on screen */
/* The char matrix must be allocated outside this function */

void printcharmat(char **charmat,int l)
{
	int i;

	for (i=0;i<l;i++)
		printf("%d %s\n",i+1,charmat[i]);
  
	printf("\n");
}







/* Main program */


gmn* lectureFichier(char ** p)
{
    char* txt = malloc(1000 * sizeof(char));
    char list_para[100];
    char list_tran[100];

    txt[0] = '\0';
    strcat(txt, *p);
    
    
    
    
    strcat(txt, "_para.dat");
    strcpy(list_para, txt);
    
    txt[0] = '\0';
    strcat(txt, *p);
    strcat(txt, "_tran.dat");
    strcpy(list_tran, txt);
	
	gmn *pgmn;
    
	pgmn=(gmn *) initgmn(list_para,list_tran);
	

    free(txt);
	return pgmn;
}





int main(gmn* pgmn, char ** p)
{
	/*
	 * gcc -shared -Wl,-soname,mouvant -o mouvant.so -fPIC mouvant.c
	 * 
	 * 
	 */
    
    
    char* txt = malloc(1000 * sizeof(char));
    char nomFichierObjectif[100];
    char nomFichierDatum[100];
    
    txt[0] = '\0';
    strcat(txt, *p);
    strcat(txt, ".ecartType");
    //printf("%s\n", txt);
    strcpy(nomFichierObjectif, txt);
    
    txt[0] = '\0';
    strcat(txt, *p);
    strcat(txt, ".datum");
    //printf("%s\n", txt);
    strcpy(nomFichierDatum, txt);
    
    
    free(txt);
    
	
	int i;
	int check;

	
	int *keys;
	int *datum;
	
    int nAPasPlante;
    
	double **effect_sta;
	double **effect_eop;

	char chaine[100000] = "";

	FILE* fp1 = NULL;
	FILE* fichierObjectif = NULL;
	
	fp1=fopen(nomFichierDatum,"r");
	if (fp1 != NULL)
	{
		fichierObjectif = fopen(nomFichierObjectif, "w+");
		while(fgets(chaine, 100000, fp1) != NULL)
		{
			if (fichierObjectif != NULL)
			{
			    nAPasPlante=1;
				datum=(int *) iallovect(pgmn->nbrsta);
				for (i=0;i<pgmn->nbrsta;i++)
					{
						if (chaine[i]==49) //Truc bizarre parce que chaine[i] vaut soit 0 ou 1 (mais 48 ou 49 dans C)
						{
							/* REMPLIR AVEC LE CHOIX DES FOURMIS*/
							datum[i]=1;
						}
						else
						{
							datum[i]=0;
						}
					}

				keys=(int *) iallovect(SEVEN);
				initivect(keys,SEVEN);
				/* On met les trois dernieres clés à 1 -> sigma x, sigma y, sigma z */
				keys[4]=1;
				keys[5]=1;
				keys[6]=1;

				check=(int) calsol(pgmn,datum,keys);
				if (check==1)
				{
					
					nAPasPlante=0;
				}
                if (nAPasPlante==1)
                {
                    effect_sta=(double **) sysreffect(pgmn,pgmn->solmat,ONE);
				    if (effect_sta==NULL)
					    {
						    
						    nAPasPlante=0;
					    }
                }
				
                if (nAPasPlante==1)
                {
                    effect_eop=(double **) sysreffect(pgmn,pgmn->solmat,TWO);
				    if (effect_eop==NULL)
					    {
						    freemat(effect_sta,SEVEN);
				            
						    nAPasPlante=0;
					    }
                }
				


				if (nAPasPlante==1)
				{
				    fprintf(fichierObjectif,"%.10e ",sqrt(effect_eop[0][0]));
				    fprintf(fichierObjectif,"%.10e ",sqrt(effect_eop[1][1]));
				    fprintf(fichierObjectif,"%.10e \n",sqrt(effect_sta[6][6]));
				    freemat(effect_sta,SEVEN);
				    freemat(effect_eop,TWO);
				}
				else
				{
				    fprintf(fichierObjectif,"1000 ");
				    fprintf(fichierObjectif,"1000 ");
				    fprintf(fichierObjectif,"1000 \n");
				}
				
				


			}
		}
		fclose(fichierObjectif); // On ferme le fichier qui a été ouvert
		
	}
	
	
	fclose(fp1);// On ferme le fichier qui a été ouvert

	
	/* PENSER A LIBERER LE GMN */
	ifreevect(datum);
	//freegmn(pgmn); // Libéré dans le script python
	ifreevect(keys);
	
	
    
	return 0;
}

void liberergmn(gmn* pgmn)
{
	freegmn(pgmn);
}


