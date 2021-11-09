#	include <stdio.h>
#	include <stdlib.h>	
#	include	<ankhead.h>
#	include <math.h>
#	include <string.h>
#	include <omp.h>


//	---------------------------------------------------------------------------------------
//				Function to read a single baseline
//															AB	26 July 2018															
//	---------------------------------------------------------------------------------------

int	read_a_baseline (BaseParType basepar,	char *namestring){
	
	int		p,d,r,c;
	char	fname[100];
	
	FILE	*fp;
	sprintf(fname,"%s_%d_%d.array",namestring,basepar.anta,basepar.antb);
	
	if((fp	=	fopen(fname,"r"))==NULL){
		printf("File not found	%s\n",fname);
		return 1;
	}
		
	for (p = 0; p < basepar.polsize; p++)
		for (d = 0; d < basepar.datsize; d++)			
			for (r=0;r<basepar.recsize;r++)
				fread(basepar.data[p][d][r], sizeof(float), basepar.chansize, fp);
	//printf("Test value	%f\n",basepar.data[0][2][100][12]);	
	fclose(fp);
			
	return 0;
}

//	---------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to write a single baseline
//															AB	26 July 2018
//	---------------------------------------------------------------------------------------

int	write_a_baseline (BaseParType basepar,	char *namestring){
	
	int		p,d,r,c;
	char	fname[100];
	float	test = 0.0;
	
	FILE	*fp;
	sprintf(fname,"%s_%d_%d_f.array",namestring,basepar.anta,basepar.antb);
	
	if((fp	=	fopen(fname,"w"))==NULL){
		printf("Cannot write to	%s\n",fname);
		return 1;
	}
	
	for (p = 0; p < basepar.polsize; p++)
		for (d = 0; d < basepar.datsize; d++)			
			for (r=0;r<basepar.recsize;r++)
				fwrite(basepar.data[p][d][r], sizeof(float), basepar.chansize, fp);
	//printf("Test value	%f\n",basepar.data[0][2][100][12]);
	fclose(fp);
				
	return 0;
}

//	-----------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to read baseline status file
//															AB	26 July 2018
//	---------------------------------------------------------------------------------------

BaseParType	*read_base_status (int ants, char *statusfile, int *nbaselines){
	
	int			i, j, k;
	int			nbase;
	BaseParType	*baseparams;
	
	nbase		=	ants*(ants-1)/2;
	*nbaselines	=	nbase;
	
	baseparams	=	(BaseParType *) malloc ( nbase * sizeof(BaseParType));
			
	FILE		*fp;
	if((fp		=	fopen(statusfile,"r"))==NULL)
		printf("Error reading status file .........!!!");
	
	for (i = 0; i < nbase; i++){
		
		fscanf(fp,"%d",&baseparams[i].anta);
		fscanf(fp,"%d",&baseparams[i].antb);
		fscanf(fp,"%d",&baseparams[i].exists);
		fscanf(fp,"%d",&baseparams[i].recsize);
		fscanf(fp,"%d",&baseparams[i].chansize);
		baseparams[i].polsize			=	POLARRS;
		baseparams[i].datsize			=	DATARRS;
		baseparams[i].flfrac[0]			=	baseparams[i].flfrac[1] = 0.0;
		baseparams[i].flfrac_after[0]	=	baseparams[i].flfrac_after[1] = 0.0;
		if(baseparams[i].exists==0){
			baseparams[i].flfrac[0]			=	baseparams[i].flfrac[1] = 1.0;
			baseparams[i].flfrac_after[0]	=	baseparams[i].flfrac_after[1] = 1.0;
		}
	}
	
	fclose(fp);
	
	return (baseparams);
}

//	-----------------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to read uvbin status file
//															AB	6 August 2018
//	---------------------------------------------------------------------------------------

BaseParType	*read_uvbin_status (int uvgrids[2], char *statusfile, int *nbaselines){
	
	int			i, j, k;
	int			nbase;
	BaseParType	*baseparams;
	
	nbase		=	uvgrids[0]*uvgrids[1];
	*nbaselines	=	nbase;
	
	baseparams	=	(BaseParType *) malloc ( nbase * sizeof(BaseParType));
			
	FILE		*fp;
	if((fp		=	fopen(statusfile,"r"))==NULL)
		printf("Error reading status file .........!!!");
	
	for (i = 0; i < nbase; i++){
		
		fscanf(fp,"%d",&baseparams[i].anta);
		fscanf(fp,"%d",&baseparams[i].antb);
		fscanf(fp,"%d",&baseparams[i].exists);
		fscanf(fp,"%d",&baseparams[i].recsize);
		fscanf(fp,"%d",&baseparams[i].chansize);
		baseparams[i].polsize			=	POLARRS;
		baseparams[i].datsize			=	DATARRS;
		baseparams[i].flfrac[0]			=	baseparams[i].flfrac[1] = 0.0;
		baseparams[i].flfrac_after[0]	=	baseparams[i].flfrac_after[1] = 0.0;
		if(baseparams[i].exists==0){
			baseparams[i].flfrac[0]			=	baseparams[i].flfrac[1] = 1.0;
			baseparams[i].flfrac_after[0]	=	baseparams[i].flfrac_after[1] = 1.0;
		}
	}
	
	fclose(fp);
	
	return (baseparams);
}

//	-----------------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to process baselines 
//															AB	5 August 2018
//	---------------------------------------------------------------------------------------

BaseParType *processbaselines(int ants, char *statusfile, int *nbaselines, char *fnames, char *flagparfile, FILE *badbasefile) {
	
	int			i,p,d,r,c,t,j,k,curound;
	int			nbase,flagmode,nbadbase;
	int			maxrecsize=0,maxchansize=0;
	int			uvgrid[2];
	BaseParType	*baseparams;
	FlagParType	*flagparams;
	int			bl_doflag, flagrounds, threads;
	float		mean_tolerance, rms_tolerance, bl_minfrac, blockpow;
	FILE		*fpflagpar;
	
	fpflagpar	=	fopen (flagparfile,"r");
	flagmode	=	init_baseflag (fpflagpar, &bl_doflag, &mean_tolerance, &rms_tolerance, &bl_minfrac, &blockpow, &flagrounds, &threads, uvgrid);
	//printf("%d	%d	%d	%d	%.2f	%.2f	%.2f	%.2f	%d %d\n", flagmode,bl_doflag,threads,flagrounds,mean_tolerance,rms_tolerance,blockpow,bl_minfrac,uvgrid[0],uvgrid[1]);
	
	flagparams	=	(FlagParType *) malloc (flagrounds * sizeof(FlagParType));
	
	for (curound = 0; curound < flagrounds; curound++)	{
		init_flagpars (fpflagpar, curound, &flagparams[curound]);
		//printf("\n%d %d %d %d %.2f %.2f %d %d %d %d %.2f %.2f\n",flagparams[curound].whatflag,flagparams[curound].doflag,flagparams[curound].domean,
					//flagparams[curound].qtype,flagparams[curound].tolerance,flagparams[curound].min_fraction,flagparams[curound].fitorder,
					//flagparams[curound].ascending,flagparams[curound].chanblockfac,flagparams[curound].recblockfac,flagparams[curound].chanmaxfrac,flagparams[curound].recmaxfrac);
	}
	if (flagmode)
		baseparams	=	read_uvbin_status(uvgrid, statusfile, nbaselines);
	else
		baseparams	=	read_base_status(ants, statusfile, nbaselines);
		
	nbase		=	*nbaselines;
	
	for (i = 0; i < nbase; i++) {
		if (baseparams[i].chansize > maxchansize) 
			maxchansize = baseparams[i].chansize; 
		if (baseparams[i].recsize > maxrecsize) 
			maxrecsize = baseparams[i].recsize;
	} 
	
	printf("\nMaximum channels	= %d	Records	= %d\n\n",maxchansize,maxrecsize);
	
	//	---------------------------------------------------------------------------------------------------
	//								Allocate buffer arrays
	//	---------------------------------------------------------------------------------------------------
	
	float		***datarray, *goutarr, ***buffarrc, ***buffarrr, *****baserawdata, ***buffdev1d, ****buffdev2d, ***basestatsarray;
	double		**coeffs;	
	
	basestatsarray	=	(float ***) malloc (nbase * sizeof(float **));
	for (i = 0; i < nbase; i++) {
		basestatsarray[i]	=	(float **) malloc (POLARRS * sizeof(float *));
		for (p=0; p<POLARRS; p++) 
			basestatsarray[i][p]	=	(float *) malloc (5* sizeof(float));
	}
	
	goutarr		=	(float *) malloc (GOUTLEN * sizeof(float));
	if ( init_gout (goutarr))	{
		free (goutarr);
		printf("Hmm... May be you forgot it at the clinic !!!\n");
		return (baseparams);
	}	
	
	baserawdata	=	(float *****) malloc ( threads * sizeof(float ****) );
	datarray	=	(float ***) malloc ( threads * sizeof(float **) );
	buffarrc	=	(float ***) malloc ( threads * sizeof(float **) );
	buffarrr	=	(float ***) malloc ( threads * sizeof(float **) );	
	coeffs		=	(double **) malloc ( threads * sizeof(double *) );
	buffdev1d	=	(float ***) malloc ( threads * sizeof(float **) );
	buffdev2d	=	(float ****) malloc ( threads * sizeof(float ***) );
	
	for (t = 0; t < threads; t++) {
	
		buffarrc[t]		=	(float **) malloc (12 * sizeof(float *));
		buffarrr[t]		=	(float **) malloc (5 * sizeof(float *));
		datarray[t]		=	(float **) malloc ( maxrecsize * sizeof(float *) );
		
		for (r = 0; r < maxrecsize; r++)
			datarray[t][r]	=	(float *) malloc ( maxchansize * sizeof(float) );
			
		for (j = 0; j < 12; j++)
			buffarrc[t][j]	=	(float *) malloc (maxchansize * sizeof(float));
			
		for (j = 0; j < 5; j++)
			buffarrr[t][j]	=	(float *) malloc (maxrecsize * sizeof(float));
			
		coeffs[t]			=	(double *) malloc (10 * sizeof(double));
		
		baserawdata[t]		=	(float ****) malloc ( POLARRS * sizeof(float ***) );
		for (p=0; p<POLARRS; p++) {
			baserawdata[t][p]	=	(float ***)malloc( DATARRS *sizeof(float **));
			
			for (d=0; d<DATARRS; d++) {
				baserawdata[t][p][d]	=	(float **)malloc( maxrecsize *sizeof(float *));
				
				for (r=0; r<maxrecsize; r++)
					baserawdata[t][p][d][r]	=	(float *)malloc( maxchansize *sizeof(float));
			}
		}
		
		buffdev1d[t]	=	(float **) malloc ( 4 * sizeof(float *) );
		buffdev2d[t]	=	(float ***) malloc ( 6 * sizeof(float **) );
		
		for (j = 0; j < 4; j++)
			buffdev1d[t][j]	=	(float *) malloc ( maxrecsize*maxchansize * sizeof(float) );
		
		for (j = 0; j < 6; j++) {
			buffdev2d[t][j]	=	(float **) malloc ( maxrecsize * sizeof(float *) );
			
			for (r = 0; r < maxrecsize; r++)
				buffdev2d[t][j][r]	=	(float *) malloc ( maxchansize * sizeof(float) );
		}
	}	
	
	//	---------------------------------------------------------------------------------------------------
	//							Process baselines in parallel
	//	---------------------------------------------------------------------------------------------------
	
	#	pragma omp parallel for num_threads (threads) shared (baseparams,fnames,fpflagpar,bl_minfrac,blockpow,flagrounds,flagparams) private (i,p,d,r,t)
		for (i = 0; i < nbase; i++) {
	
			if (baseparams[i].exists==0) {
				printf("Baseline	%d %d		Flagged......\n",baseparams[i].anta,baseparams[i].antb);
				continue;
			}
			
			t		=	omp_get_thread_num();
			
			baseparams[i].data	=	baserawdata[t];		
		
			if(read_a_baseline (baseparams[i], fnames)) 
				printf("Baseline	%d %d		Reading Error...... !!!!!\n",baseparams[i].anta,baseparams[i].antb);
		
			//	----------------------	Actual flagging function	------------------------------------------------		
			
			printf("\nBaseline	%d %d		taken by thread	%d\n\n",baseparams[i].anta,baseparams[i].antb,t);
			
			flagbaseline (&baseparams[i], bl_minfrac, blockpow, flagrounds, flagparams, datarray[t], goutarr, buffarrc[t], buffarrr[t], coeffs[t], 
																												buffdev1d[t], buffdev2d[t], basestatsarray[i]);			
		
			//	----------------------	Flagging done	---------------------------------------------------
		
			if(write_a_baseline (baseparams[i], fnames)) 
				printf("Baseline	%d %d		Writing Error...... !!!!!\n",baseparams[i].anta,baseparams[i].antb);
		
			printf("\nBaseline	%d %d	Fractions	",baseparams[i].anta,baseparams[i].antb);					
			for (p=0; p<POLARRS; p++)
				printf("%.3f %.3f	",baseparams[i].flfrac[p],baseparams[i].flfrac_after[p]);
			printf("\n\n");
		}
	
	//	--------------------------------------------------------------------------------------------------------
	//							Flag bad baselines
	//	-------------------------------------------------------------------------------------------------------
	
	if (flagmode==0) {
		for (i = 0; i < nbase; i++) {
			if (baseparams[i].exists==0)
				for (p=0; p<POLARRS; p++)
					basestatsarray[i][p][4]	=	1.0;
			//printf("%f	%f	%f	%f\n",basestatsarray[i][0][0],basestatsarray[i][1][0],basestatsarray[i][0][1],basestatsarray[i][1][1]);
		}	
	
		nbadbase	=	findbadbase (baseparams, nbase, basestatsarray, bl_doflag, mean_tolerance, rms_tolerance, bl_minfrac, maxchansize*maxrecsize, goutarr, badbasefile);
	}
	
	//	-------------------------------------------------------------------------------------------------------
	
	free (basestatsarray);
	free (buffdev1d);
	free (buffdev2d);
	free (baserawdata);
	free (buffarrc);
	free (buffarrr);
	free (coeffs);
	free (datarray);
	free (goutarr);
	free (flagparams);
	
	fclose (fpflagpar);
	return (baseparams);
}

//	-------------------------------------------------------------------------------------------------------





































































