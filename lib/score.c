/*---------------------------------------------------------
   score:   Compiled program.
   rf:      File containing all sequences (any order).
   wf1:     Output file containing alignment results 
   wf2:     Output file containing difference matrix.
   wf3:     Output file containing alignment.
   pen:     gap penalty (default=8).
   rbase:   random base
  
				       12/20/86  df    
------------------------------------------------------------*/
/* JMK 3/20/97	cut read_seq and replaced with fasta format 
	compile: gcc score.c -o score -lm
	command: score rf wf1 wf2 wf3 [DOO || FTA] [pen=8] [rbase=77] 
		!!! note
			everything in [] must be added in order.		    
------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define NSEQ 600        /*Number of sequences*/
#define MAXS 3000       /*Maximum size for any sequence*/
#define MAXG 1000	/* Maximum number of gaps in an alignment */
#define MAXAS MAXS+MAXG  /*Maximum size for aligned sequences*/
#define	CODESZ 5

#define DOO 101
#define FTA 103

#define Mt 25
#define Ms 20
#define Mr 18
#define Mq 17
#define Mp 15
#define Mo 14
#define Mn 13
#define Mm 12
#define Ml 11
#define Mk 10
#define Mj 9
#define Mi 8 
#define Mh 7
#define Mg 6
#define Mf 5
#define Me 4
#define Md 3
#define Mc 2
#define Mb 1
#define Ma 0
#define XB 0


short int mptable[21][21] = {
{Ms,Mi,Mg,Mf,Mg,Mf,Me,Md,Md,Md,Mf,Me,Md,Md,Mg,Mc,Mg,Me,Mi,Ma,XB},
{Mi,Mk,Mj,Mj,Mj,Mj,Mj,Mi,Mi,Mh,Mh,Mi,Mi,Mg,Mh,Mf,Mh,Mf,Mf,Mg,XB},
{Mg,Mj,Ml,Mi,Mj,Mi,Mi,Mi,Mi,Mh,Mh,Mh,Mi,Mh,Mi,Mg,Mi,Mf,Mf,Md,XB},
{Mf,Mj,Mi,Mo,Mj,Mh,Mh,Mh,Mh,Mi,Mi,Mi,Mh,Mg,Mg,Mf,Mh,Md,Md,Mc,XB},
{Mg,Mj,Mj,Mj,Mk,Mj,Mi,Mi,Mi,Mi,Mh,Mg,Mh,Mh,Mh,Mg,Mi,Me,Mf,Mc,XB},
{Mf,Mj,Mi,Mh,Mj,Mn,Mi,Mj,Mi,Mh,Mg,Mf,Mg,Mf,Mf,Me,Mh,Md,Md,Mb,XB},
{Me,Mj,Mi,Mh,Mi,Mi,Mk,Mk,Mj,Mj,Mk,Mi,Mj,Mg,Mg,Mf,Mg,Me,Mg,Me,XB},
{Md,Mi,Mi,Mh,Mi,Mj,Mk,Mm,Ml,Mk,Mj,Mh,Mi,Mf,Mg,Me,Mg,Mc,Me,Mb,XB},
{Md,Mi,Mi,Mh,Mi,Mi,Mj,Ml,Mm,Mk,Mj,Mh,Mi,Mg,Mg,Mf,Mg,Md,Me,Mb,XB},
{Md,Mh,Mh,Mi,Mi,Mh,Mj,Mk,Mk,Mm,Ml,Mj,Mj,Mh,Mg,Mg,Mg,Md,Me,Md,XB},
{Mf,Mh,Mh,Mi,Mh,Mg,Mk,Mj,Mj,Ml,Mo,Mk,Mi,Mg,Mg,Mg,Mg,Mg,Mi,Mf,XB},
{Me,Mi,Mh,Mi,Mg,Mf,Mi,Mh,Mh,Mj,Mk,Mo,Ml,Mi,Mg,Mf,Mg,Me,Me,Mk,XB},
{Md,Mi,Mi,Mh,Mh,Mg,Mj,Mi,Mi,Mj,Mi,Ml,Mn,Mi,Mg,Mf,Mg,Md,Me,Mf,XB},
{Md,Mg,Mh,Mg,Mh,Mf,Mg,Mf,Mg,Mh,Mg,Mi,Mi,Mo,Mk,Mm,Mk,Mi,Mg,Me,XB},
{Mg,Mh,Mi,Mg,Mh,Mf,Mg,Mg,Mg,Mg,Mg,Mg,Mg,Mk,Mn,Mk,Mm,Mj,Mh,Md,XB},
{Mc,Mf,Mg,Mf,Mg,Me,Mf,Me,Mf,Mg,Mg,Mf,Mf,Mm,Mk,Mo,Mk,Mk,Mh,Mg,XB},
{Mg,Mh,Mi,Mh,Mi,Mh,Mg,Mg,Mg,Mg,Mg,Mg,Mg,Mk,Mm,Mk,Mm,Mh,Mg,Mc,XB},
{Me,Mf,Mf,Md,Me,Md,Me,Mc,Md,Md,Mg,Me,Md,Mi,Mj,Mk,Mh,Mq,Mp,Mi,XB},
{Mi,Mf,Mf,Md,Mf,Md,Mg,Me,Me,Me,Mi,Me,Me,Mg,Mh,Mh,Mg,Mp,Mr,Mi,XB},
{Ma,Mg,Md,Mc,Mc,Mb,Me,Mb,Mb,Md,Mf,Mk,Mf,Me,Md,Mg,Mc,Mi,Mi,Mt,XB},
{XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB,XB},
                      };


short int xindex[26]={  4, 7, 0, 7, 8,17, 5,10,
		14,21,12,15,13, 6,21, 3,
		 9,11, 1, 2,21,16,19,20,
		18, 8,
		 };

char aseq1[MAXAS], aseq2[MAXAS];
char seq1[MAXAS], seq2[MAXAS];
char seqt[NSEQ][MAXG], code[NSEQ][CODESZ+1];
char seq[NSEQ][MAXS];
char *strrchr();

 int L[MAXS][MAXS], R[MAXS][MAXS], M[MAXS][MAXS];
 int lseq[NSEQ], nseq;
 int pen;     /*None negative gap penalty*/
 int bord_flag, Rbase;
 int family[NSEQ][NSEQ];     /*Groupings of Related Sequences*/
 int seq_format = 0;

double pro[NSEQ][NSEQ];     /*Changing Matrix*/
double D[NSEQ][NSEQ];
double maxD = 0;	    /* to normalize distance: 1.0 >= D >= 0.0 */

main(ac,av)
int ac;
char *av[];
{
	FILE *fd1, *fd2, *fd3, *fd4, *fopen();   /*File description*/ 
	int i,j,k,si,sj;
	double Align();

   /* Check for correct command line */
	if(ac<5)
		{
		system("/usr/ucb/clear");
		fprintf(stderr,"Purpose: Binary comparisons of protein sequences - LOM matrix.\n");
		fprintf(stderr,"Usage: ");
		fprintf(stderr,"%s rf wf1 wf2 wf3 [DOO || FTA] [pen] [rbase]\n",av[0]);
   		fprintf(stderr,"\trf:  All sequences (any order)\n");
		fprintf(stderr,"\twf1: Alignment results\n");
		fprintf(stderr,"\twf2: Similarity matrix\n");
		fprintf(stderr,"\twf3: All binary alignments\n");
		fprintf(stderr,"\t[DOO || FTA]: Optional sequence format, Doolittle or Fasta, default = FTA\n");
		fprintf(stderr,"\t[pen]: Optional Integer > 0, gap penalty, default = 8\n");
		fprintf(stderr,"\t[rbase]: Optional Integer > 0, random base, default = 77\n");
		fprintf(stderr,"!!!Note: Options must be set INORDER,  i.e.\n");
		fprintf(stderr,"pen requires setting seq. format, rbase requires both pen and seq. format\n\n");
		fprintf(stderr,"Current limitations: Length =< %d(unaligned), %d(aligned), ",MAXS,MAXAS);
		fprintf(stderr,"NSEQ =< %d\n",NSEQ);
		fprintf(stderr,"Score programs: \n");
		fprintf(stderr,"\tscore:  wf2 contains lower left distance matrix.\n");
		exit(0);
		}

   /* Open fd1 */
	if( (fd1 = fopen( av[1],"r" )) == NULL )   /*error check*/
		exit(fprintf(stderr,"%s: can't open %s\n",av[1]) );
	if( (fd2 = fopen( av[2],"w" )) == NULL )   /*error check*/
		exit(fprintf(stderr,"%s: can't open %s\n",av[2]) );
	if( (fd3 = fopen( av[3],"w" )) == NULL )   /*error check*/
		exit(fprintf(stderr,"%s: can't open %s\n",av[3]) );
	if( (fd4 = fopen( av[4],"w" )) == NULL )   /*error check*/
		exit(fprintf(stderr,"%s: can't open %s\n",av[4]) );

   /* Set default parameters */
     if (ac > 5) {
	if ((!strcmp ("DOO", av[5])) || (!strcmp ("doo", av[5])))
		seq_format = DOO;
	else if ((!strcmp ("FTA", av[5])) || (!strcmp ("fta", av[5])))
		seq_format = FTA;
     }
     else seq_format = FTA;

     if (ac > 6) {
	if((pen = atoi(av[6]))==0)
		pen=8;
     }
     else pen=8;

     if (ac > 7) {
	if((Rbase = atoi(av[7]))==0)
		Rbase=77;
     }
     else Rbase=77;

     bord_flag = 0; /* Assume all D > 0 */

   /* Read sequences */
	switch (seq_format) {
		case DOO: Read_doo_seq(fd1);
			  break;
		case FTA: Read_fta_seq(fd1);
			  break;
		default : fprintf (stderr, " Incorrect Input Format Specified \n");
			  break;
	} /* end switch */
	

   /* Align all combinatorial pairs */
	fprintf(fd2,"\n\n            Binary comparisons between %d ",nseq);
	fprintf(fd2,"protein sequences\n");
	fprintf(fd2,"\t\t\t(gap penalty = %d)\n\n",pen);
	fprintf(fd2,"\t\t\t(random base = %d)\n\n",Rbase);
	for(i=0; i<nseq; ++i)
	    fprintf(fd2,"   %2d)%s:\t%s\n",i+1,code[i],seqt[i]);
	fprintf(fd2,"\n");
/*	fprintf(fd2,"    Comparisons         Sreal    NAS     %%id");
	fprintf(fd2,"    Siden   Srand     D\n");
*/
	fprintf(fd2,"    Comparisons           Sreal  NAS      %%id");
	fprintf(fd2,"     Siden   Srand    D\n");
	for(i=0; i<nseq-1; ++i)
		for(j=i+1; j<nseq; ++j)
			{
			fprintf(fd2,"%-6s (%d) %-6s (%d) ",code[i],lseq[i]-1,
				code[j],lseq[j]-1);
			fprintf(fd4,"%s: %s",code[i],seqt[i]);
			fprintf(fd4,"     (sequence length=%d)\n",lseq[i]-1);
			fprintf(fd4,"%s: %s",code[j],seqt[j]);
			fprintf(fd4,"     (sequence length=%d)\n",lseq[j]-1);
		       
			// printf("\n\n%s", seq[i]);
			//	printf("\n\n%s", seq[j]);
			strcpy(seq1,seq[i]);
			strcpy(seq2,seq[j]);
                       
			D[i+1][j+1] = Align(fd2,fd4,lseq[i],lseq[j]);
			
			if (maxD < D[i+1][j+1])
				maxD = D[i+1][j+1];

			D[j+1][i+1] = D[i+1][j+1];
			}
/*	if(bord_flag == 1)
		{
		fprintf(fd2,"\nDifference scores from binary alignments:\n");
		fprintf(fd2,"      ");
		for(i=0; i<nseq; ++i)
			fprintf(fd2,"  %s  ",code[i]);
		for(i=0; i<nseq; ++i)
			{
			fprintf(fd2,"\n%s ",code[i]);
			for(j=0; j<nseq; ++j)
				fprintf(fd2,"%7.2lf ",D[i+1][j+1]);
			}
		fprintf(fd2,"\n\n  *** Distant pairs are indicated by -1.00 ***\n\n");
		look_for_neg();
		}
   */
   /* Copy difference matrix into a file */		
/*	for(i=1; i<nseq+1; ++i)				
		{					
		for(j=1; j<i+1; ++j)		
			fprintf(fd3,"%7.2lf ",D[i][j]);	
		fprintf(fd3,"\n");			
		}					
*/

/* JMK 3/20/97: for gentool-cluster program			*/
/*	Copy lower left difference matrix to file fd3(wf2)	*/
	fprintf(fd3,"six lead in lines for gentool-cluster program \n");
	fprintf(fd3,"\nLower left Matrix \n");
	fprintf(fd3,"\nDiagonals are included, yes to cluster program \n\n");
	for(i=1; i<nseq+1; ++i)				 	
		{						
		for(j=1; j<i+1; ++j)
		{
			if(j == i)
				fprintf(fd3,"%1.6f ",0);
			else
				fprintf(fd3,"%1.6f ",(1.1-(D[i][j]/maxD)));
		}
		fprintf(fd3,"\n");				
		}						

   /* Branching order */
	bord(fd2);

	fprintf(fd2,"\n");
	for(i=0; i<nseq; ++i)
		fprintf(fd2,"%2d)%s: %s",i+1,code[i],seqt[i]);
	fprintf(fd2,"\n");
	fclose(fd4);
	fclose(fd3);
	fclose(fd2);
	fclose(fd1);
return 0;
}

Read_doo_seq(fd)
FILE *fd;
{
	int i, m;
	char s;

   /* Enter sequences */
	m=nseq=0;
	for( ; ; )
		{
		if(fgetc(fd)==EOF)
			return;
		fgets(code[m],5L,fd); /*Get the 1st 4 letters as code*/
		fgetc(fd);
		fgetc(fd);
		while((s=fgetc(fd))==' ')
			;
		ungetc(s,fd);
		fgets(seqt[m], sizeof seqt[m], fd);  /*Get sequence title*/

 		seq[m][0]=' ';
       		fseek(fd,11,1); /*Skip 11 characters*/
       		for(i=1; (s=fgetc(fd)) && (s != '*'); ++i)
			{
			if(s=='\n')
				{
				--i;
				fseek(fd,11,1);
				}
			else if(s==' ')
				--i;
			else
				seq[m][i]=s;
	       		}
		fgetc(fd);
		lseq[m]=i;
		seq[m][i]='\0';
		printf("\n\n%d) %s",m, seq[m]);
		m++;
		nseq++;
		}
}


/* JMK 3/20/97	 fasta format	*/
Read_fta_seq(fd)
FILE *fd;
{
	int i, m;
	char s;

   /* Enter sequences */
	m=nseq=0;
	for( ; ; )
		{
		if((fgetc(fd))==EOF)
			return;
	/*	fgets(code[m],5L,fd); Get the 1st 4 letters as code*/
		while((s=fgetc(fd))==' ') 			;
	        ungetc(s,fd);
		fgets(seqt[m], sizeof(seqt[m]), fd);  /*Get sequence title*/
		strncpy(code[m],seqt[m], CODESZ);
		for(i = CODESZ; (i > 1 && !isalnum(code[m][i])); i--)
		{
			code[m][i] = ' ';
		}
		code[m][i+1] = '\0';
 		seq[m][0]=' ';
       		for(i=1; ((s=fgetc(fd)) && (s != '>') && (s != EOF)); ++i)
			{
			if(s=='\n')
				--i;
			else if(s==' ')
				--i;
			else
				seq[m][i]=s;
	       		}
		if (s == '>')
			ungetc(s,fd);
		lseq[m]=i;
		seq[m][i]='\0';
		m++;
		nseq++;
		}
}

double Align(fd,fda,si,sj)
FILE *fd, *fda;
int si,sj;
{
	int Mmax;
	int dum1, dum2;
	int i, ii, j, k, jj, kk, m, n, t;
	int ir,mu,nc,nar,ssigma;
	int sigma;   /*Similarity constant between pairs of letters*/
	int Sim;   /*Aligned sequences optimal similarity score*/
/*	int nMa,nMb,nMc,nMd,nMe,nMf,nMg,nMh,nMi,nMj;
#	int nMk,nMl,nMm,nMn,nMo,nMp,nMq,nMr,nMs,nMt;
	int n1, n2, n3, mv, 
*/
	int c;
	int gaps;     /*Total number of gaps introduced*/
	int blanks, kgps1, kgps2, gps1[MAXG], gps2[MAXG];
	int mp, np, pdum, maxp, ascore;
	int oh[4];   /* Overhang lengths */
	int bks1max, bks2max;

	double fidty, fsim, Srand, sum, diff;
	double bscore;
	double fn1, fn2, fn1n2;
	double lsmall, NAS;

   /* Construct the first row and column of the 3 matrices */
	for( j=0; j<si; ++j )
		L[j][0] = R[j][0] = M[j][0] = 0;
	for( k=0; k<sj; ++k )
		L[0][k] = R[0][k] = M[0][k] = 0;

   /* Construct the remaining elements of the 3 matrices */
	for( j=1; j<si; ++j )
		{
		jj = j - 1;
		for( k=1; k<sj; ++k )
			{
			kk = k - 1;
			sigma = mptable[xindex[seq1[j]-'A']]
				       [xindex[seq2[k]-'A']];
			Mmax = M[jj][kk];     /*Determine M[j][k]*/
			dum1 = L[jj][kk] - pen;
			dum2 = R[jj][kk] - pen;
			if( Mmax < dum1 )
				Mmax = dum1;
			if( Mmax < dum2 )
				Mmax = dum2;
			M[j][k] = Mmax + sigma;
			L[j][k] = L[jj][k];   /*Determine L[j][k]*/
			if( L[jj][k] < M[jj][k] )
				L[j][k] = M[jj][k];
			R[j][k] = R[j][kk];    /*Determine R[j][k]*/
			if( R[j][kk] < M[j][kk] )
				R[j][k] = M[j][kk];
			}
		}

   /* Determine the optimal score of the aligned sequences */
	m = si-1;
	n = sj-1;
	Sim = ( R[m][n] > L[m][n] ) ? R[m][n] : L[m][n];
	Sim = ( Sim > M[m][n] ) ? Sim : M[m][n];
	fsim = Sim;
	fsim = fsim/10.0;
	fprintf(fda,"Optimal (similarity) score of alignment = %3.2f\n",fsim);

   /* Tracing out the aligned sequences */
	j = m;
	k = n;
	c = 0;
	if( Sim == M[j][k] )
		{
		aseq1[c] = seq1[j];
		aseq2[c] = seq2[k];
		c++;
		branchM(&j,&k,&c);
		}
	else if( Sim == L[j][k] )
		{
		aseq1[c] = seq1[j];
		aseq2[c] = ' ';
		c++;
		branchL(&j,&k,&c);
		}
	else
		{
		aseq1[c] = ' ';
		aseq2[c] = seq2[k];
		c++;
		branchR(&j,&k,&c);
		}

   /* Catch the N-terminal overhangs */
	if( j==0 )
		{
		while( k>-1 )
			{
			aseq1[c] = ' ';
			aseq2[c++] = seq2[--k];
			}
		}
	if( k==0 )
		{
		while( j>-1 )
			{
			aseq1[c] = seq1[--j];
			aseq2[c++] = ' ';
			}
		}
	if( c>MAXAS )   /*error check*/
		exit(printf("job aborted, please increase MAXAS"));

   /* Compute percent identity */
	lsmall = si;
	if( si > sj )
		lsmall = sj;

	nc = ir= nar = 0;
	for( j=0; j<c; ++j)
		{
		if( aseq1[j] != aseq2[j] )
			;
		else if( aseq1[j] == ' ')
			;
		else if( aseq1[j] == 'C' )
			{
			nc++;
			ir++;
			}
		else
			ir++;
		if( aseq1[j] >= 'A' && aseq1[j] <= 'Z' )
			{
			if( aseq2[j] >= 'A' && aseq2[j] <= 'Z' )
				nar++;
			else
				;
			}
		else
			;
		}
	
	fidty = ir;
	fidty = fidty / nar;
	fprintf(fda,"%d matches, ",ir);
	fprintf(fda,"%6.2f%%id, ",fidty*100);
	
   /* Compute and print NAS */
	NAS = fsim*100.0 / nar;
	fprintf(fda,"NAS=%4.2f, ",NAS);

   /* Compute the total number of internal gaps */
	ssigma = 0;
	for( j=0; j<c; ++j )
		{
		if( (aseq1[j] == ' ') || (aseq2[j] == ' ') )
			;
		else
			ssigma += mptable[xindex[aseq1[j]-'A']]
					 [xindex[aseq2[j]-'A']];
		}
	gaps = ( ssigma - Sim ) / pen;
	fprintf(fda,"%d internal gaps.\n",gaps);

   /* Compute the average similarity score when seq1 is aligned to
      itself and seq2 is aligned to itself */
	sum=0;
	for(i=1; i<si; ++i)
		sum += mptable[xindex[seq1[i]-'A']][xindex[seq1[i]-'A']];
	for(i=1; i<sj; ++i)
		sum += mptable[xindex[seq2[i]-'A']][xindex[seq2[i]-'A']];
	sum=sum/20.0;

   /* Random score: Srand=Rbase*lsmall/100 */
	Srand = Rbase*(lsmall-1.0)/100.0;

   /* Print summary in a file */
	fprintf(fd,"%7.2lf %7.2lf %7.2lf",fsim,NAS,fidty*100.0);
	fprintf(fd," %7.2lf %7.2lf ",sum,Srand);
	if(fsim <= Srand)
		{
		bord_flag = 1;  /* At least one diff < 0 */
		fprintf(fd,"  Too distant!\n");
		diff = -1.00;
		}
	else
		{
		diff = -log((fsim-Srand)/(sum-Srand))*100.0;
		fprintf(fd," %7.2lf\n",diff);
		}

   /* Compute the internal gap length in each sequence */
	blanks = 0;
	bks1max = 0;
	for( j=1; j<MAXG; ++j )
		gps1[j] = 0;

	for( j=c-3; j>-1; )
		{
		k = j;
		if( (aseq1[j] != ' ') && (aseq1[--j] == ' ') )
			{
			blanks++;
			while( aseq1[--j] == ' ' )
				blanks++;
			if( j>-1 )
				{
				gps1[blanks]++;
				bks1max = (bks1max>blanks ) ? bks1max : blanks;
				blanks=0;
				}
			}
		if( j==k )
			j--;
		}


	blanks = 0;
	bks2max = 0;
	for( j=1; j<MAXG; ++j )
		gps2[j] = 0;

	for( j=c-3; j>-1; )
		{
		k = j;
		if( (aseq2[j] != ' ') && (aseq2[--j] == ' ') )
			{
			blanks++;
			while( aseq2[--j] == ' ' )
				blanks++;
			if( j>-1 )
				{
				gps2[blanks]++;
				bks2max = ( bks2max>blanks ) ? bks2max : blanks;
				blanks=0;
				}
			}
		if( j==k )
			j--;
		}
	
   /* Compute lengths of overhangs */
	blanks = 0;
	j = 0;
	while( aseq1[j] == ' ' )
		{
		blanks++;
		j++;
		}
	oh[0] = blanks;
	
	blanks = 0;
	j = c-3;
	while( aseq1[j] == ' ' )
		{
		blanks++;
		j--;
		}
	oh[1] = blanks;
	
	blanks = 0;
	j = 0;
	while( aseq2[j] == ' ' )
		{
		blanks++;
		j++;
		}
	oh[2] = blanks;
	
	blanks = 0;
	j = c-3;
	while( aseq2[j] == ' ' )
		{
		blanks++;
		j--;
		}
	oh[3] = blanks;

   /* Print the internal gap lengths of each sequence */
	gaps = 0;
	fprintf(fda,"Internal gaps: freq(gap length)\n  Seq 1: ");
	for( j=1; j<=bks1max; j++)
		{
		if( gps1[j] != 0 )
			{
			fprintf(fda,"%d(%d), ",gps1[j],j);
			gaps += gps1[j];
			}
		}
	fprintf(fda,"\n  Seq 2: ");
	for( j=1; j<=bks2max; j++)
		{
		if( gps2[j] != 0 )
			{
			fprintf(fda,"%d(%d), ",gps2[j],j);
			gaps += gps2[j];
			}
		}	
	fprintf(fda,"\n\n");

   /* Print the aligned sequences */
	fprintf(fda,"      --------------- SmfLOM alignment ---------------\n\n");
	j = c-3;
	while( j>-1 )
		{
		k = 0;
		t = j;
		while( k<30 )
			{
			if( j>-1 )
				fprintf(fda,"%2c",aseq1[j--]);
			k++;
			}
		fprintf(fda,"\n ");
		k = 0;
		while( k<30 )
			{
			if( t>-1 )
				{
				fprintf(fda,"%c",aseq2[t]);
				if( aseq1[t] != aseq2[t] )
					fprintf(fda," ");
				else if( aseq1[t] == ' ' )
					fprintf(fda," ");
				else
					fprintf(fda,"_");
				}
			k++;
			t--;
			}
		j = t;
		fprintf(fda,"\n\n");
		}
	fprintf(fda,"      ------------------------------------------------\n\n");

return(diff);
}

branchM(pj,pk,pc)
int *pj,*pk,*pc;
{
	int Mmax;

	if( (*pj ==0) || (*pk==0) )
		return;
	(*pj)--;
	(*pk)--;

	Mmax = M[*pj][*pk];
	if( Mmax < L[*pj][*pk] - pen )
		Mmax = L[*pj][*pk] - pen;
	if( Mmax < R[*pj][*pk] - pen )
		Mmax = R[*pj][*pk] - pen;

	if( Mmax == M[*pj][*pk] )
		{
		aseq1[*pc] = seq1[*pj];
		aseq2[*pc] = seq2[*pk];
		(*pc)++;
		branchM(pj,pk,pc);
		return;
		}

	else if( Mmax == (L[*pj][*pk] - pen) )
		{
		aseq1[*pc] = seq1[*pj];
		aseq2[*pc] = ' ';
		(*pc)++;
		branchL(pj,pk,pc);
		return;
		}

	else
		{
		aseq1[*pc] = ' ';
		aseq2[*pc] = seq2[*pk];
		(*pc)++;
		branchR(pj,pk,pc);
		return;
		}
}

branchL(pj,pk,pc)
int *pj,*pk,*pc;
{
	if( *pj == 0 )
		return;

	(*pj)--;
	aseq1[*pc] = seq1[*pj];

	if( L[*pj][*pk] > M[*pj][*pk] )
		{
		aseq2[*pc] = ' ';
		(*pc)++;
		branchL(pj,pk,pc);
		return;
		}

	else
		{
		aseq2[*pc] = seq2[*pk];
		(*pc)++;
		branchM(pj,pk,pc);
		return;
		}
}

branchR(pj,pk,pc)
int *pj,*pk,*pc;
{
	if( *pk == 0 )
		return;

	(*pk)--;
	aseq2[*pc] = seq2[*pk];

	if( R[*pj][*pk] > M[*pj][*pk] )
		{
		aseq1[*pc] = ' ';
		(*pc)++;
		branchR(pj,pk,pc);
		return;
		}

	else
		{
		aseq1[*pc] = seq1[*pj];
		(*pc)++;
		branchM(pj,pk,pc);
		return;
		}
}

bord(fd)
FILE *fd;
{
	int i,ii,ii1,j,k,q,qz,side,ss1,v,vv,w,z,zt,zz,Z;
	int minj1,mink1,r,s,t;

	side = 0;
	i=qz=1;
	i=s=Z=0;
	for(j=1;j<nseq;++j)
		{
		for(k=j+1;k<=nseq;++k)
			pro[j][k]=pro[k][j]=D[j][k];
		}
	fprintf(fd,"\n\nDistance scores for phylogenetic tree\n\n");
	if(nseq<=16)
		{
		ii=1;
		fprintf(fd,"      ");
		for(t=0;t<nseq;++t)
			fprintf(fd,"  %s  ",code[t]);
		ii=ii1=w=1;
		for(t=0;t<nseq;++t)
			{
			fprintf(fd,"\n%s ",code[t]);
			for(q=0;q<nseq;++q)
				fprintf(fd,"%7.2lf ",pro[ii1][w++]);
			w=1;
			ii1++;
			}
		}
	FINDMIN(&minj1, &mink1);
	zt=1;
	ss1=nseq-1;
	pro[mink1][minj1]=pro[minj1][mink1]=1000;
	family[0][0]=minj1;
	family[0][1]=mink1;
	for(zz=1;zz<=ss1;++zz)
		{
		NEWMATRIX();
		FINDMIN(&minj1, &mink1);
		FAMFIND(&Z, &ss1, minj1, mink1);
		if(family[qz][0]!=0)
			{
			zz--;
			qz++;
			}
		if(zz== ss1-1 && zt==1)
			goto rat2;
		}
	rat2: 	fprintf(fd,"\n\n\n");
		fprintf(fd,"Branching order: ");
		for(i=0;i<=Z;++i)
			{
			if(side==1)
				fprintf(fd,"\n       Clusters: ");
			else if(side>1)
				fprintf(fd,"\n                 ");
			fprintf(fd,"(%s,%s)",code[family[i][0]-1],
					     code[family[i][1]-1]);
			for(v=2;family[i][v]!=0;)
				{
				fprintf(fd,"%s)",code[family[i][v]-1]);
				v++;
				if(v%11==0)
					fprintf(fd,"\n                  ");
				}
			side++;
			}
		fprintf(fd,"\n");
}


FINDMIN(pminj1,pmink1)
int *pminj1, *pmink1;
{
	int j,k;

	j=1;
	k=2;
	*pminj1 = j;
	*pmink1 = k;
	while(j<nseq)
		{
		while(k<nseq)
			{
			if(pro[*pminj1][*pmink1] > pro[j][++k])
				{
				*pminj1 = j;
				*pmink1 = k;
				}
			else 
				continue;
			}
		++j;
		k=j;
		}
	return;
}

NEWMATRIX()
{
	int i,j,k,q,qq1,qq2;

	j=1;
	k=2;
	while(j<nseq)
		{
		while(k<=nseq)
			{
			for(i=0;family[i][0]!=0;++i)
				{
				for(q=0;family[i][q]!=0;q++)
					{
					if(pro[j][k]==1000)
						goto NEWK;
					else if(j==family[i][q])
						{
						qq1=k;
						qq2=j;
						AV(qq1,qq2,i);
						goto NEWK;
						}
					else if(k==family[i][q])
						{
						qq1=j;
						qq2=k;
						AV(qq1,qq2,i);
						goto NEWK;
						}
					else 
						continue;
					}
				}
	 		NEWK: ++k;
			}
		j++;
		k=j+1;
		}
	return;
}

FAMFIND(pZ,pss1,minj1,mink1)
int *pZ, *pss1, minj1, mink1;
{
   /*Will detect whether the pair with the*/
   /*minimum value belongs to any other family*/
	int i,q,q1,W;

	for(i=0;family[i][0]!=0;++i)
		{
  		for(q=0;family[i][q]!=0;++q)
			{
    			if(family[i][q]==minj1)
				{
				q1=mink1;
				OTHER(q1,i);
				goto tree;
				}
    			if(family[i][q]==mink1)
				{
				q1=minj1;
				OTHER(q1,i);
				goto tree;
				}
    			}
		}
	(*pZ)++;
	if(W=1)
		W=0;
		/*A new pair to be placed into a new family*/

	pro[minj1][mink1]=pro[mink1][minj1]=1000;
	family[*pZ][0]=minj1;
	family[*pZ][1]=mink1;
	W++;
	if(W==1)
		(*pss1)--;
	tree: return;
}

OTHER(q1,i)
int q1,i;
{
		/*Checks to see if the other min. value belongs*/
		/*to another family*/
	int ii,ii1,v,vv;

	for(ii=0;family[ii][0]!=0;ii++)
		{
		for(ii1=0;family[ii][ii1]!=0;++ii1)
			{
			if(q1==family[ii][ii1])
				{
				ANOTHER(&ii,i);
				return;
				}
			}
		}

	/*Does not belong to a new family*/

	for(v=0;v<NSEQ;++v)
		{
		if(family[i][v]==0)
			{
			family[i][v]=q1;
			break;
			}
		}
	for(v=0;family[i][v+1]!=0;++v)
		{
		for(vv=1;family[i][vv]!=0;++vv)
			{
			if(family[i][v]!=family[i][vv])
				pro[family[i][v]][family[i][vv]]=
				pro[family[i][vv]][family[i][v]]=1000;
			}
		}
	return;
}

ANOTHER(pii,i)
int *pii, i;
{
	      /*Min. does belong to another family...add to*/
	      /*family of first min. value*/
	int v,vv;

	for(v=0;v<NSEQ;v++)
		{
		if(family[i][v]==0)
			{
/* ------ Bug. Never checking for array outof bounds, (v, vv). 
			for(vv=0;family[*pii][vv]!=0;++vv)
				family[i][v++]=family[*pii][vv];
*/
			for(vv=0;(family[*pii][vv]!=0)&&(vv<NSEQ);++vv)
				family[i][v++]=family[*pii][vv];
			break;
			}
		}
	for(v=0;family[i][v+1]!=0;v++)
		{
		for(vv=1;family[i][vv]!=0;++vv)
			{
			if(family[i][v]!=family[i][vv])
				pro[family[i][v]][family[i][vv]]=
				pro[family[i][vv]][family[i][v]]=1000;
			}
		}
	return;
}

AV(qq1,qq2,i)
int qq1,qq2,i;
{
	 /*Determines whether qq1 also belongs to a*/
	 /*family; if so, average throughout both families*/

	int ii,ii1,m,n,q,v;

	m=n=v=0;
	pro[m][n]=0;
	for(ii=0;family[ii][0]!=0;++ii)
		{
		for(ii1=0;family[ii][ii1]!=0;++ii1)
			{
			if(qq1==family[ii][ii1])
				goto heman;
			else 
				continue;
			}
		}
	for(q=0;family[i][q]!=0;++q)
		{
		++v;
		pro[m][n] += D[family[i][q]][qq1];
		}
	pro[m][n]=pro[m][n]/v;
	pro[qq1][qq2]=pro[qq2][qq1]=pro[m][n];
	goto heman2;

	    /*qq1 belongs to a separate family average*/
	    /*thoughout both families*/

	heman:	for(ii1=0;family[ii][ii1]!=0;ii1++)
			{
			for(q=0;family[i][q]!=0;++q)
				{
				v++;
				pro[m][n] += D[family[i][q]]
						 [family[ii][ii1]];
				}
			}
		pro[m][n]=pro[m][n]/v;
		for(ii1=0;family[ii][ii1]!=0;++ii1)
			{
			for(q=0;family[i][q]!=0;++q)
				pro[family[ii][ii1]][family[i][q]]=
				pro[family[i][q]][family[ii][ii1]]=pro[m][n];
			}

	heman2: return;
}

look_for_neg()
{
	int i,ii,j,temp_j;

	for(i=1; i<nseq; i++)
		{
		for(j=i+1; j<nseq+1; j++)
			{
			if(D[i][j]<0 && j!=nseq+1)
				{
				temp_j = j-1;
				while(temp_j<nseq)
					{
					for(ii=temp_j; ii<nseq; ++ii)
						strcpy(code[temp_j],code[temp_j+1]);
					temp_j++;
					}
				temp_j = j-1;
				while(temp_j<nseq)
					{
					for(ii=temp_j; ii<nseq; ++ii)
						strcpy(seqt[temp_j],seqt[temp_j+1]);
					temp_j++;
					}
				temp_j = j;
				while(temp_j<nseq+1)
					{
					for(ii=1; ii<nseq+1; ++ii)
						D[ii][temp_j] = D[ii][temp_j+1];
					temp_j++;
					}
				temp_j = j;
				while(temp_j<nseq+1)
					{
					for(ii=1; ii<nseq+1; ++ii)
						D[temp_j][ii] = D[temp_j+1][ii];
					temp_j++;
					}
				nseq--;
				look_for_neg();
				}
			}
		}
	return;
}
