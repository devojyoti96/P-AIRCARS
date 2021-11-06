#	include <stdio.h>
#	include <stdlib.h>	
#	include	<ankhead.h>
#	include <math.h>
#	include <string.h>


//	-----------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	int				i,j,k;
	
	BaseParType		*baseparams;
	int				ants, flagmode;
	char			statusfile[100],binfiles[100],flagparfile[100],badbasename[100];
	int				nbase	=	0;
	POLARRS		=	atoi(argv[1]);	

	sprintf(flagparfile,"scratch/flagpars.pars");
	sprintf(badbasename,"scratch/badbase.list");
	
	FILE		*fp,*badbasefile;
	fp		=	fopen(flagparfile,"r");
	badbasefile	=	fopen(badbasename,"w");
	
	fscanf (fp,"%d	%d",&ants,&flagmode);
	
	fclose(fp);
	
	if (flagmode) {
		sprintf(statusfile,"scratch/uvbin_status.txt");
		sprintf(binfiles,"scratch/uvbin");	
	}
	else {
		sprintf(statusfile,"scratch/baseline_status.txt");
		sprintf(binfiles,"scratch/baseline");
	}
	
	printf("\nFlagmode	= %d	Antennas	= %d\n\n",flagmode,ants);
	
	baseparams	=	processbaselines (ants, statusfile, &nbase, binfiles, flagparfile, badbasefile);	
	
	printf("Total baselines / uvbins		%d\n",nbase);
	
	/*for (i = 0; i < nbase; i++)
		printf("Baseline / UVbin	%d %d	%d	%.3f %.3f	%.3f %.3f\n",baseparams[i].anta,baseparams[i].antb,baseparams[i].exists,
				baseparams[i].flfrac[0],baseparams[i].flfrac[1],baseparams[i].flfrac_after[0],baseparams[i].flfrac_after[1]);*/
		
	free (baseparams);
	fclose (badbasefile);
	
	return 0;
}

























