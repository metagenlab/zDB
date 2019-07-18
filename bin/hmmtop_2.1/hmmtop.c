/******************************************************************************/
/* To complile the program execute:                                           */
/*      cc hmmtop.c -lm -o hmmtop                                             */
/*                                                                            */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*       HH     HH  M       M  M       M  TTTTTTTT    OOOOO    PPPPPP         */
/*       HH     HH  MM     MM  MM     MM  TTTTTTTT   OOO OOO   PP    PP       */
/*       HHHHHHHHH  MMM   MMM  MMM   MMM     TT     OO     OO  PP    PP       */
/*       HHHHHHHHH  MM MMM MM  MM MMM MM     TT     OO     OO  PPPPPP         */
/*       HH     HH  MM  M  MM  MM  M  MM     TT      OOO OOO   PP             */
/*       HH     HH  MM     MM  MM     MM     TT       OOOOO    PP             */
/*                                                                            */
/*                                                                            */
/*     Prediction of transmembrane helices and topology for transmembrane     */
/*                   proteins based on Hidden Markov Model                    */
/*                                                                            */
/*                                                                            */
/*                                 2.1 version                                */
/*                                  2000,2001                                 */
/*                                                                            */
/*                                      by                                    */
/*                     Gabor E. Tusnady and Istvan Simon                      */
/*             Institute of Enzymology, Biological Research Center            */
/*                       Hungarian Academy of Sciences                        */
/*                                                                            */
/*                  H-1113 Budapest Karolina ut 29. HUNGARY                   */
/*                     H-1518 Budapest P.O.Box 7 HUNGARY                      */
/*    Telephone : (36-1) 466-9276, (36-1) 466-5633/162, (36-1) 466-5633/132   */
/*                         Telefax: (36-1) 466-5465                           */
/*                           E-mail: tusi@enzmi.hu                            */
/*                    Internet: http://www.enzim.hu/hmmtop                    */
/*                                                                            */
/*                 Copyright (C) 2000, 2001, Gabor E. Tusnady                 */
/*                                                                            */
/******************************************************************************/

#define HMMTOP_VERSION_NUMBER 2
#define HMMTOP_MINOR_NUMBER 1

/*Changelog: 2.1-2.0
  - Bug fixing: name of the first sequence in a fasta file contained the '>'.
	  Now, it is corrected. Thanks for Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>
		reporting this bug.
		May 03, 2006
  - Bug fixing: in hmmtop_optim() the limit of reporting Serious error was lowered
    from zero to -1e-6. Thanks for Michal Jaffee <mikij@volcani.agri.gov.il> 
		reporting this bug.
		Sept 11, 2003
  - Bug fixing: hmmtop_addlocate() contained a wrong memory allocation (sizeof(char) 
    instead of sizeof(char*)). Thanks for Ville Karhu <Ville.Karhu@ktl.fi> reporting 
    this bug.
    March 17, 2003
  - The setenv is not defined in sgi, therefore it is changed to putdef
    according to the applied architecture
*/
/*ChangeLog: 1.0-2.0
  - Final revision od the source and program manual.
    April 18, 2001.
  - The architecture were systematically varied to reach the best prediction 
    accuracy on single sequences. The results are the followings:
    on the 158 protein (reported on our JMB article):
      original     698     716  689  ( 97, 99%)     149  ( 94%)     135  ( 85%)
      fast (noit)  698     685  670  ( 97, 96%)     133  ( 84%)     119  ( 75%)
    on the 148 protein reported by Moller et al (Bioinformatics):
      original     711     740  691  ( 95, 97%)     133  ( 90%)     107  ( 73%)
      fast (noit)  711     681  653  ( 94, 92%)     107  ( 73%)      89  ( 61%)
    just for comparison on the 148 sequences:
      Name     Pred OK                  sum OK tm       sum OK top
      das      720  647  ( 90, 91%)     110  ( 74%)       0  (  0%)
      predtmr  677  644  ( 93, 91%)     104  ( 70%)       0  (  0%)
      sosui    677  634  ( 91, 89%)     108  ( 73%)       0  (  0%)
      tmap     693  636  ( 91, 89%)     101  ( 68%)      61  ( 41%)
      tmpred   717  655  ( 92, 92%)     108  ( 73%)      69  ( 47%)
      phdhtm   687  646  ( 92, 91%)     111  ( 75%)      87  ( 59%)
      memsat   692  641  ( 91, 90%)     119  ( 80%)      95  ( 64%)
      memsat2  673  642  ( 93, 90%)     113  ( 77%)      88  ( 60%)
      toppred  772  685  ( 92, 96%)     128  ( 86%)      96  ( 65%)
      tmhmm    690  658  ( 94, 93%)     113  ( 77%)      94  ( 64%)   !!!!!
      hmm1.0   733  681  ( 94, 96%)     128  ( 86%)     100  ( 68%)
    February 20, 2001.
  - one switch --print_probabilities, -pp is changed to print all probabilities
    in the architecture.
    February 19, 2001.
  - The old architecture was applied in a new frame. The whole source was
    rewritten to allow any particular architecture given in a file pointed 
    by the HMMTOP_ARCH environment variable, see the syntax there.
    February 6, 2001.
  - added a new switch: -locf (--loc_file). This have to point to a file
    containing the locates for each sequence (line by line). The syntax of
    the locates are the same as in -loc= switch.
  - Added a new function: there is an oppotunity to locate any part of the
    protein to a structural part, with the argument -loc=bpos-epos-spart 
    or --locate=bpos-epos-spart where bpos and epos are the begin and end 
    position of the location, respectively; spart is the structural part in 
    which locate the sequence piece and may be I,i,H,o,O. If the end position
    is 'E' then it means the last amino acid in the sequence.
    June 27, 2000
  - Added a new switch: HMMTOP_ITERATION. Default value is 1, but can be switch
    off by the -noit or --noiteration switch. In this case the program do the 
    same as memsat or tmhmm, i.e. calculate the best topology if the model is
    given.
    Januar 11, 2001
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define LIMIT 5e-4
#define IFUT 20
#define MPRED 1
#define SPRED 0
#define RANDOM 0
#define PSEUDO 1
#define ATA 5
#define HCUR 0
#define HLOG 1
#define HNEW 2
#define HSAV 3
#define HINI 4

typedef struct {
	char *code;
	int t,e,*f,*bt,*bm,nf,nb;
	char *name,let;
	double ini,end;
} STATE;

int hmmtop_readfasta(FILE *);
int hmmtop_readpir(FILE *);
int hmmtop_readswissprot(FILE *);
int  hmmtop_readseq(FILE *, char*);
void hmmtop_usage(char *);
void hmmtop_shortusage(char *);
void hmmtop_main(char*,char*,char*,char*,char*);
void hmmtop_optim(void);
void hmmtop_bestpath(int);
void hmmtop_output(FILE *);
void hmmtop_getpsv(void);
void hmmtop_error(char*);
void hmmtop_info(char*);
void hmmtop_setmem(void);
void hmmtop_newsequence(void);
void hmmtop_newaminoacid(char);
void hmmtop_getname(FILE *);
void hmmtop_onenorm(double*,int);
void hmmtop_onelog(double*,double*,int);
void hmmtop_norm(void);
void hmmtop_copy(int,int);
void hmmtop_zero(int);
void hmmtop_backward(int);
void hmmtop_addlocate(char *);
void hmmtop_setlocate(void);
void hmmtop_readarchitecture(void);
void hmmtop_readlocates();
double hmmtop_log(double);
double hmmtop_forward(int);
double hmmtop_core(void);


char 	AS[]="ACDEFGHIKLMNPQRSTVWY";
char	as[]="acdefghiklmnpqrstvwy";
char	**name,translist[1000][20],emilist[1000][20],OF_NAME[1000],**locates;
int	**seq,*prpart,**MH,trdb[1000],NS,NT,NP,drb,*n,nmax,nas,nma,**LOCATE,locnum;
int	FUT=1,NFH=1000,START=PSEUDO,PMODE,OPP,OPC,OPL;
int  	HMMTOP_ITERATION=1;
double	***P,***T,**I,**E,**H,PM,PMM,TM=10000;
STATE   state[1000];
FILE	*LOG=NULL,*LOCF=NULL;
int first_fasta=1;

int main(int argc, char **argv)
{
int i,j=0;
char *sf,*ifi,*pi,*is,*ofi,*lf;
char *locf;

  fprintf(stderr,"\n\n");
  fprintf(stderr,"HMMTOP %d.%d, Copyright (C) 2000, Gabor E. Tusnady\n\n",
  	HMMTOP_VERSION_NUMBER,HMMTOP_MINOR_NUMBER);
  srand48(time(NULL)%97);
  OPP=OPC=OPL=0;
  nmax=nma=nas=0;
  sf=ifi=ofi=pi=is=lf=locf=NULL;
  locnum=0;
  for (i=1; i<argc; i++) if (argv[i][0]=='-')
  { if (!strcmp(argv[i],"--help")) hmmtop_usage(argv[0]);
    if (!strncmp(argv[i],"--sequence_format=",18)) sf=&argv[i][18];
    if (!strncmp(argv[i],"--input_file=",13)) ifi=&argv[i][13];
    if (!strncmp(argv[i],"--output_file=",14)) ofi=&argv[i][13];
    if (!strncmp(argv[i],"--log_file=",11)) lf=&argv[i][11];
    if (!strncmp(argv[i],"--process_inputfile=",19)) pi=&argv[i][19];
    if (!strncmp(argv[i],"--pseudo_size=",14)) TM=atof(&argv[i][14]);
    if (!strncmp(argv[i],"--iteration_start=",18)) is=&argv[i][18];
    if (!strncmp(argv[i],"--iteration_number=",19)) FUT=atoi(&argv[i][19]);
    if (!strncmp(argv[i],"--print_probabilities",21)) OPP=1;
    if (!strncmp(argv[i],"--print_pseudocount",19)) OPC=1;
    if (!strncmp(argv[i],"--print_longprediction",22)) OPL=1;
    if (!strncmp(argv[i],"--short_help",12)) hmmtop_shortusage(argv[0]);
    if (!strncmp(argv[i],"--locate=",9)) hmmtop_addlocate(&argv[i][9]);
    if (!strncmp(argv[i],"--loc_file=",11)) locf=&argv[i][11];
    if (!strcmp(argv[i],"--noiteration")) HMMTOP_ITERATION=0;
    if (!strcmp(argv[i],"-sh")) hmmtop_shortusage(argv[0]);
    if (!strcmp(argv[i],"-h")) hmmtop_usage(argv[0]);
    if (!strncmp(argv[i],"-sf=",4)) sf=&argv[i][4];
    if (!strncmp(argv[i],"-if=",4)) ifi=&argv[i][4];
    if (!strncmp(argv[i],"-of=",4)) ofi=&argv[i][4];
    if (!strncmp(argv[i],"-lf=",4)) lf=&argv[i][4];
    if (!strncmp(argv[i],"-pi=",4)) pi=&argv[i][4];
    if (!strncmp(argv[i],"-ps=",4)) TM=atof(&argv[i][4]);
    if (!strncmp(argv[i],"-is=",4)) is=&argv[i][4];
    if (!strncmp(argv[i],"-in=",4)) FUT=atoi(&argv[i][4]);
    if (!strncmp(argv[i],"-pp",3)) OPP=1;
    if (!strncmp(argv[i],"-pc",3)) OPC=1;
    if (!strncmp(argv[i],"-pl",3)) OPL=1;
    if (!strncmp(argv[i],"-loc=",5)) hmmtop_addlocate(&argv[i][5]);
    if (!strncmp(argv[i],"-locf=",6)) locf=&argv[i][6];
    if (!strcmp(argv[i],"-noit")) HMMTOP_ITERATION=0;
    j++;
  }
  if (j==0) hmmtop_shortusage(argv[0]);

  if (ifi==NULL)
  { printf("No input sequence file is given, see usage\n"); hmmtop_usage(argv[0]);}
  
  if (lf!=NULL)
  { if (!strcmp(lf,"--")) LOG=stdout; else
    if ((LOG=fopen(lf,"w"))==NULL) perror(lf), hmmtop_error("Bad file name\n");
    setlinebuf(LOG);
    fprintf(LOG,"Log file of HMMTOP %d.%d\n",HMMTOP_VERSION_NUMBER,HMMTOP_MINOR_NUMBER);
    fprintf(LOG,"Parameters were the following:\n");
    for (i=0; i<argc; i++) fprintf(LOG,"%d \"%s\"\n",i,argv[i]);
    fprintf(LOG,"End of parameter list\n");
  }
  
  if (locf!=NULL) {
    if ((LOCF=fopen(locf,"r"))==NULL) perror(locf), hmmtop_error("Bad file name");
  }
  
  hmmtop_main(ifi,ofi,sf,pi,is);
  
  if (LOG!=NULL) 
  { fprintf(LOG,"Program terminated normally\n");
    fclose(LOG);
  }

  return 0;
}

void hmmtop_usage(char *pn)
{
  printf("Usage: \"%s options\", where options are the following:\n",pn);
  printf("\n    -if=filename, --input_file=filename\n");
  printf("          path and name of the input sequence file.  If  filename\n");
  printf("          is -- program reads from the standard input.\n");
  printf("\n    -of=filename, --output_file=filename\n");
  printf("          path and name of the output sequence file. If  filename\n");
  printf("          is -- program writes to the standard output.\n");
  printf("\n    -lf=filename, --log_file=filename\n");
  printf("          path and name of the log file.  Use  only for debugging\n");
  printf("          purpose.\n");
  printf("\n    -sf=XXX, --sequence_format=XXX\n"); 
  printf("          format of sequence(s).  XXX may be FAS for fasta format\n");
  printf("          (default), PIR for NBRF/PIR format or  SWP  for  SWISS-\n");
  printf("          PROT format.\n");
  printf("\n    -pi=xxx, --proces_inputfile=xxx\n");
  printf("          treating sequences in input file as single or homologue\n");
  printf("          sequences. xxx may be spred or mpred.  In case of spred\n");
  printf("          for each sequence in input_file will be made prediction\n");
  printf("          (default). In case of mpred only for the first sequence\n");
  printf("          in the input file will  be made  prediction,  remaining\n");
  printf("          sequences will be treated as helpers.\n");
  printf("\n    -ps=xxx, --pseudo_size=xxx\n");
  printf("          size of the pseudo count vector.  xxx may be from 0 (no\n");
  printf("          pseudo count vector used) to 10000 (default=10000).\n");
  printf("\n    -loc=bpos-epos-spart, --locate=bpos-epos-spart\n");
  printf("          localize a portion of the sequence in a given structur-\n");
  printf("          al part.  bpos and  epos are the begin and end position\n");
  printf("          in the sequence, respectively;  spart is the structural\n");
  printf("          part and may be i,I,H,o and O.  epos  may be the letter\n");
  printf("          'E', which means the C terminal end of the sequence\n");
  printf("\n    -locf=filename, --loc_file=filename\n");
  printf("          name of the  'locate file',  which contains the locates\n");
  printf("          for each protein  if multiple  sequences with  pi=spred\n");
  printf("          are given.  Each line contains one or more locates. The\n");
  printf("          format of the locates are as given in the -loc= switch.\n");
  printf("\n    -noit, --noiteration\n");
  printf("          do not optimize the parameters on the given sequence(s)\n");
  printf("          just calculate the best path if the model is given.  It\n");
  printf("          is useful for rapid but not so accurate prediction.\n");
  printf("\n    -is=xxx, --iteration_start=xxx\n");
  printf("          starting point of the iteration(s).  xxx  may be pseudo\n");
  printf("          or random.  In case of pseudo iteration starts from the\n");
  printf("          pseudo  count  vector  (default).  In  case  of  random\n");
  printf("          iteration starts from random values.\n");
  printf("\n    -in=xxx, --iteration_number=xxx\n");
  printf("          maximum allowed iterations (if is=random).\n");
  printf("\n    -pp, --print_probabilities\n");
  printf("          print the optimized probabilies.\n");
  printf("\n    -pc, --print_pseudocount\n");
  printf("          print the pseudocount vector used.\n");
  printf("\n    -pl, --print_longprediction\n");
  printf("          print prediction in a long format,  ie. localization of\n");
  printf("          each amino acid in the input sequence.\n");
  printf("\n    -h, --help\n");
  printf("          this help message.\n");
  printf("\n    -sh, --short_help\n");
  printf("          print a short help message.\n");
  printf("\n\nReferences:\n");
  printf("          G.E. Tusnady and I. Simon (1998)\n");
  printf("          Principles Governing Amino Acid Composition of Integral\n");
  printf("          Membrane Proteins:  Applications to topology prediction\n");
  printf("          J. Mol. Biol. 283, 489-506\n");
  printf("          http://www.enzim.hu/hmmtop\n\n");
  exit(0);
}

void hmmtop_shortusage(char *pn)
{
  printf("Usage: \"%s options\", where options are the following:\n",pn);
  printf("-if=filename    name of input file, -- stands stdin\n");
  printf("-of=filename    name of output file, -- stands stdout\n");
  printf("-lf=filename    name of the logfile, use only for debugging purpose\n");
  printf("-sf=XXX         format of sequence(s).  XXX=PIR or FAS or SWP\n");
  printf("-pi=xxx         process input xxx=spred or mpred\n");
  printf("-ps=xxx         size of pseudo count vector xxx=0-10000\n");
  printf("-loc=b-e-sp     localize the sequence between b and e in sp\n");
  printf("-locf=filename  name of the loc_file containing the locates for each sequence\n");
  printf("-noit           do not iterate for optimization of the model\n");
  printf("-is=xxx         starting point of the iteration(s). xxx=pseudo or random\n");
  printf("-in=xxx         number of iterations (only if is=random)\n");
  printf("-pp             print the optimized probabilies\n");
  printf("-pc             print the pseudocount vector used\n");
  printf("-pl             print prediction in a long format\n");
  printf("-h              print a long help message\n");
  printf("-sh             print this (short) help message\n");
  printf("References:     J. Mol. Biol. (1998), 283, 489-506\n");
  printf("                http://www.enzim.hu/hmmtop\n");
  exit(0);
}

void hmmtop_main(char *ifn, char *ofn, char *sf, char *pm, char *is)
{
int i;
char *mysf;
char *mypm;
char *myis;
FILE *IF=NULL,*OF=NULL;

  if (ifn==NULL) hmmtop_error("No input file");
  else
  { if (!strcmp(ifn,"--")) IF=stdin; else 
    { if ((IF=fopen(ifn,"r"))==NULL) perror(ifn), hmmtop_error("Bad input file name");}
  }
  
  if (ofn==NULL) 
  {	OF=stdout; strcpy(OF_NAME,"stdout");}
  else
  { if (!strcmp(ofn,"--")) 
    { OF=stdout; strcpy(OF_NAME,"stdout");}
    else
    { if ((OF=fopen(ofn,"w"))==NULL) perror(ofn), hmmtop_error("Bad output filename");
      setlinebuf(OF);
      strcpy(OF_NAME,ofn);
    }
  }
  
  if (sf==NULL) 
  { hmmtop_info("No input file format, using fasta format");
    mysf=calloc(4,sizeof(char)); strcpy(mysf,"FAS");
  } else mysf=sf;
  if (pm==NULL) 
  { hmmtop_info("No process mode, using spred");
    mypm=calloc(10,sizeof(char)); strcpy(mypm,"spred");
  } else mypm=pm;
  if (is==NULL) 
  { hmmtop_info("No iteration start point. Start from pseudo vector");
    myis=calloc(10,sizeof(char)); strcpy(myis,"pseudo");
  } else myis=is;
  

  
  
  if (!strcasecmp("FAS",mysf) && !strcasecmp("PIR",mysf) && !strcasecmp("SWP",mysf))
  hmmtop_error("Unknown input file format");
  
  if (!strcasecmp("mpred",mypm) && !strcasecmp("spred",mypm))
  hmmtop_error("Unknown process mode");
  
  
  hmmtop_readarchitecture();
  hmmtop_getpsv();
  
  
  if (!strcasecmp(myis,"pseudo")) { FUT=1; START=PSEUDO;} else START=RANDOM;
  
  if (!strcasecmp(mypm,"mpred"))
  { PMODE=MPRED;
    drb=0;
    while (hmmtop_readseq(IF,sf));
    fclose(IF);
    hmmtop_optim();
    hmmtop_bestpath(0);
    hmmtop_output(OF);
  }
  else if (!strcasecmp(mypm,"spred"))
  { PMODE=SPRED;
    drb=0; hmmtop_newsequence();
    while (hmmtop_readseq(IF,sf))
    if (n[0]>20)
    { hmmtop_optim();
      hmmtop_bestpath(0);
      hmmtop_output(OF);
      drb=0; n[0]=0;
    }
    fclose(IF);
  }
  fclose(OF);  
}

void hmmtop_error(char *s)
{ if (LOG!=NULL) fprintf(LOG,"HMMTOP error: %s\n",s);
  fprintf(stderr,"HMMTOP error: %s\n",s); exit(-1);
exit(-1);
}

void hmmtop_info(char *s)
{ if (LOG!=NULL) fprintf(LOG,"HMMTOP info: %s\n",s);
  fprintf(stderr,"HMMTOP info: %s\n",s);
}

void hmmtop_setmem()
{
int i,j;

  if (LOG!=NULL) fprintf(LOG,"Entering setmem %d %d\n",nma,nmax);
  if (nmax>nma) 
  { for (i=0; i<NS; i++)
    { E[i]=realloc(E[i],nmax*sizeof(double));
      H[i]=realloc(H[i],nmax*sizeof(double));
      MH[i]=realloc(MH[i],nmax*sizeof(int));
      LOCATE[i]=realloc(LOCATE[i],nmax*sizeof(int));
    }
    prpart=realloc(prpart,nmax*sizeof(int));
    nma=nmax;
  }
  if (LOG!=NULL) fprintf(LOG,"Setmem OK %5d %5d\n",nma,nmax);
}

void hmmtop_newsequence()
{
int *ns,i;
int **sseq;
char **sname;
 
  if (LOG!=NULL) fprintf(LOG,"Newseq: %d\n",drb);
  if ((ns=calloc(drb+1,sizeof(int)))==NULL) hmmtop_error("Not enough memory");
  if ((sseq=calloc(drb+1,sizeof(short*)))==NULL) hmmtop_error("Not enough memory");
  if ((sname=calloc(drb+1,sizeof(char*)))==NULL) hmmtop_error("Not enough memory");
  for (i=0; i<drb; i++) ns[i]=n[i],sseq[i]=seq[i],sname[i]=name[i];
  if (n!=NULL) free(n);
  if (seq!=NULL) free(seq);
  if (name!=NULL) free(name);
  n=ns; seq=sseq; name=sname;
  n[drb]=0; nas=0;
}

void hmmtop_newaminoacid(char aa)
{
char a;

  a='`';
  if (strchr(as,aa)!=NULL) a=(char)(strchr(as,aa)-as);
  else if (strchr(AS,aa)!=NULL) a=(char)(strchr(AS,aa)-AS);
  if (a!='`')
  { if (n[drb]==nas)
    { nas+=50;
      if ((seq[drb]=realloc(seq[drb],nas*sizeof(int)))==NULL) hmmtop_error("Not enough memory");
    }
    seq[drb][n[drb]++]=a;
  }
}

void hmmtop_getname(FILE *IF)
{
int i;
char c,s[100];
    for (c=' ',i=0; (c!='\n' && !feof(IF) && i<100); i++)
    { c=getc(IF); s[i]=c;}
    s[i-1]='\0';
    while (c!='\n' && !feof(IF)) c=getc(IF);
    if ((name[drb]=realloc(name[drb],i*sizeof(char)))==NULL) 
      hmmtop_error("Not enough memory");
    strcpy(name[drb],s);
    if (LOG!=NULL) fprintf(LOG,"New name: %d %d %s %p\n",drb,i,name[drb],name);

}

int hmmtop_readfasta(FILE *IF)
{
char c;
int i;

	if (first_fasta) {
		for (c=' '; (c!='>' && !feof(IF)); c=getc(IF));
    first_fasta=0;
	}
  if (!feof(IF))
  { if (PMODE==MPRED) hmmtop_newsequence(); 
    else {drb=0; n[drb]=0;}
    hmmtop_getname(IF);
    for (c=' '; (c!='>' && !feof(IF)); )
    { c=getc(IF); hmmtop_newaminoacid(c);}
    if (n[drb]>nmax) nmax=n[drb];
    if (n[drb]>0) 
    { if (LOG!=NULL) 
      { fprintf(LOG,"Sequence in fasta format is the following:\n");
        fprintf(LOG,">%s\n",name[drb]);
	for (i=0; i<n[drb]; i++) 
	{ fprintf(LOG,"%c",AS[seq[drb][i]]);
	  if ((i+1)%10==0) fprintf(LOG," ");
	  if ((i+1)%50==0) fprintf(LOG,"\n");
	}
	fprintf(LOG,"\n");
      }
      drb++;
      return 1;
    }
    else return 0;
  }
  else return 0;
}

int hmmtop_readpir(FILE *IF)
{
char c;
int i;

  for (c=' '; (c!='>' && !feof(IF)); c=getc(IF));
  if (!feof(IF))
  { if ((c=getc(IF))=='P')
    { if ((c=getc(IF))=='1')
      { if ((c=getc(IF))==';')
	{ if (PMODE==MPRED) hmmtop_newsequence(); 
	  else {drb=0; n[drb]=0;}
	  hmmtop_getname(IF);
	  for (c=' '; (c!='\n' && !feof(IF)); c=getc(IF));
	  while (c!='*' && !feof(IF))
	  { c=getc(IF); hmmtop_newaminoacid(c);}
	  if (n[drb]>nmax) nmax=n[drb];
	  if (n[drb]>0) 
	  { if (LOG!=NULL) 
	    { fprintf(LOG,"Sequence in pir format is the following:\n");
              fprintf(LOG,">P1;%s\n\n",name[drb]);
	      for (i=0; i<n[drb]; i++) 
	      { fprintf(LOG,"%c",AS[seq[drb][i]]);
		if ((i+1)%10==0) fprintf(LOG," ");
		if ((i+1)%50==0) fprintf(LOG,"\n");
	      }
	      fprintf(LOG,"*\n");
	    }
	    drb++;
	    return 1;
	  }
	  else return 0;
	} else return 0;
      } else return 0;
    } else return 0;
  }
  return 0;
}

int hmmtop_readswissprot(FILE *IF)
{
char s[100],c;
int i;
  for (strcpy(s,""); (!feof(IF) && strncmp(s,"ID",2)); fgets(s,100,IF));
  if (!feof(IF))
  { if (PMODE==MPRED) hmmtop_newsequence();
    else {drb=0; n[drb]=0;}
    if ((name[drb]=realloc(name[drb],12*sizeof(char)))==NULL) 
      hmmtop_error("Not enough memory");
    for (i=0; (s[i+5]!=' ' && i<11); i++) name[drb][i]=s[i+5]; name[drb][i]='\0';
    for (strcpy(s,""); (!feof(IF) && strncmp(s,"SQ",2)); fgets(s,100,IF));
    c=' ';
    while (c!='/' && !feof(IF))
    { c=getc(IF); hmmtop_newaminoacid(c);}
    if (n[drb]>nmax) nmax=n[drb];
    if (n[drb]>0) 
    { if (LOG!=NULL) 
      { fprintf(LOG,"Sequence in Swissprot format is the following:\n");
        fprintf(LOG,"ID   %s\n",name[drb]);
	fprintf(LOG,"SQ\n");
	for (i=0; i<n[drb]; i++) 
	{ fprintf(LOG,"%c",AS[seq[drb][i]]);
	  if ((i+1)%10==0) fprintf(LOG," ");
	  if ((i+1)%50==0) fprintf(LOG,"\n");
	}
	fprintf(LOG,"\n//\n");
      }
      drb++;
      return 1;
    }
  }
  return 0;
}

int hmmtop_readseq(FILE *IF, char *sf)
{ if (!feof(IF))
  { if (LOCF!=NULL) hmmtop_readlocates();
    if (sf==NULL) return hmmtop_readfasta(IF);
    else
    { if (!strcasecmp(sf,"FAS")) return hmmtop_readfasta(IF);
      if (!strcasecmp(sf,"PIR")) return hmmtop_readpir(IF);
      if (!strcasecmp(sf,"SWP")) return hmmtop_readswissprot(IF);
    }
  } 
  return 0;
}

void hmmtop_getpsv()
{
int i,j;
FILE *f;

	if (LOG!=NULL) fprintf(LOG,"Entering Pseudocount reading\n");
	if (getenv("HMMTOP_PSV")==NULL) {
		hmmtop_info("HMMTOP_PSV is not set, trying in the current directory...");
#ifdef sgi
		putenv("HMMTOP_PSV=hmmtop.psv");
#else
		setenv("HMMTOP_PSV","hmmtop.psv",0);
#endif
	}
    	if ((f=fopen(getenv("HMMTOP_PSV"),"r"))==NULL) 
		hmmtop_error("psv file not found");
		
	for (i=0; i<NS; i++) fscanf(f,"%lf",&I[HCUR][i]); 
	fscanf(f,"\n");
	for (i=0; i<NP; i++) {
		for (j=0; j<20; j++) {
			fscanf(f,"%lf",&P[HCUR][i][j]);
		}
		fscanf(f,"\n");
	}
	for (i=0; i<NT; i++) {
		for (j=0; j<trdb[i]; j++) {
			fscanf(f,"%lf",&T[HCUR][i][j]);
		}
		fscanf(f,"\n");
	}
	fclose(f);
	hmmtop_norm();
	hmmtop_copy(HCUR,HINI);
	for (i=0; i<NP; i++) for (j=0; j<20; j++) P[HINI][i][j]*=TM;
}

double hmmtop_log(double q) {
double p; if (q>1e-20) p=log(q); else p=-1000; return(p);}

void hmmtop_onenorm(double *x, int h) {
int i;
double q;
  for (i=0; i<h; i++) if (x[i]<1e-20) x[i]=0;
  for (i=0,q=0; i<h; i++) q+=x[i];
  if (q>0) for (i=0; i<h; i++) x[i]/=q;
}

void hmmtop_onelog(double *x, double *y, int h) {
int i; for (i=0; i<h; i++) y[i]=hmmtop_log(x[i]);}

void hmmtop_norm() {
int i;
  hmmtop_onenorm(I[HCUR],NS); hmmtop_onelog(I[HCUR],I[HLOG],NS);
  for (i=0; i<NP; i++) hmmtop_onenorm(P[HCUR][i],20), hmmtop_onelog(P[HCUR][i],P[HLOG][i],20);
  for (i=0; i<NT; i++) hmmtop_onenorm(T[HCUR][i],trdb[i]), hmmtop_onelog(T[HCUR][i],T[HLOG][i],trdb[i]);
  if (LOG!=NULL) fprintf(LOG,"norm OK\n");
}
void hmmtop_copy(int from, int to) {
int i,j;
  for (i=0; i<NS; i++) I[to][i]=I[from][i];
  for (i=0; i<NP; i++) for (j=0; j<20; j++) P[to][i][j]=P[from][i][j];
  for (i=0; i<NT; i++) for (j=0; j<trdb[i]; j++) T[to][i][j]=T[from][i][j];
  if (LOG!=NULL) fprintf(LOG,"copy %d %d OK\n",from,to);
}
void hmmtop_zero(int what) {
int i,j;
  for (i=0; i<NS; i++) I[what][i]=0;
  for (i=0; i<NP; i++) for (j=0; j<20; j++) P[what][i][j]=0;
  for (i=0; i<NT; i++) for (j=0; j<trdb[i]; j++) T[what][i][j]=0;
  if (LOG!=NULL) fprintf(LOG,"zero %d OK\n",what);
}
double hmmtop_forward(int f) {
int i,j,k;
double q,max,S[1000];

  if (LOG!=NULL) fprintf(LOG,"Entering forward %d\n",f);
  for (j=0; j<NS; j++) E[j][0]=I[HLOG][j]+P[HLOG][state[j].e][seq[f][0]]-1e5*LOCATE[j][0];
  for (i=1; i<n[f]; i++) for (j=0; j<NS; j++) {
    for (k=0, max=-1e30; k<state[j].nb; k++) {
      S[k]=E[state[j].bt[k]][i-1]+T[HLOG][state[state[j].bt[k]].t][state[j].bm[k]];
      if (max<S[k]) max=S[k];
    }
    for (k=0,q=0; k<state[j].nb; k++) q+=exp(S[k]-max);
    E[j][i]=log(q)+max+P[HLOG][state[j].e][seq[f][i]]-1e5*LOCATE[j][i];
  }
  for (i=0; i<NS; i++) E[i][n[f]-1]+=state[i].end;
  for (i=0,max=-1e30; i<NS; i++) 
    if (max<E[i][n[f]-1]) max=E[i][n[f]-1]; 
  for (i=0,q=0; i<NS; i++) q+=exp(E[i][n[f]-1]-max);
  q=log(q)+max;
  if (LOG!=NULL) fprintf(LOG,"forward OK %.2f\n",q);
  return(q);
}

void hmmtop_backward(int f) {
int i,j,k,to;
double q,max,S[1000];

  if (LOG!=NULL) fprintf(LOG,"Entering backward\n");
  for (i=0; i<NS; i++) H[i][n[f]-1]=state[i].end-1e5*LOCATE[i][n[f]-1];
  for (i=n[f]-2; i>=0; i--) for (j=0; j<NS; j++) {
    for (k=0, max=-1e30; k<state[j].nf; k++) {
      to=state[j].f[k];
      S[k]=H[to][i+1]+T[HLOG][state[j].t][k]+P[HLOG][state[to].e][seq[f][i+1]];
      if (max<S[k]) max=S[k];
    }
    for (k=0,q=0; k<state[j].nf; k++) q+=exp(S[k]-max);
    H[j][i]=log(q)+max-1e5*LOCATE[j][i];
  }
  if (LOG!=NULL) fprintf(LOG,"backward OK\n");
}

double hmmtop_core() {
int   i,j,k,fsz,to;
double   sum,NN,LOP;

  hmmtop_zero(HNEW);
  NN=0; sum=0;
  for (fsz=0; fsz<drb; fsz++) {  
    LOP=hmmtop_forward(fsz); 
    hmmtop_backward(fsz); 
    sum+=LOP; NN+=n[fsz];
    for (i=0; i<NS; i++) I[HNEW][i]+=exp(E[i][0]+H[i][0]-LOP);
    for (i=0; i<n[fsz]; i++) for (j=0; j<NS; j++)
      P[HNEW][state[j].e][seq[fsz][i]]+=exp(E[j][i]+H[j][i]-LOP);
    for (i=0; i<n[fsz]-1; i++) for (j=0; j<NS; j++)
    for (k=0; k<state[j].nf; k++) {
      to=state[j].f[k];
      T[HNEW][state[j].t][k]+=
      exp(E[j][i]+T[HLOG][state[j].t][k]+
          P[HLOG][state[to].e][seq[fsz][i+1]]+H[to][i+1]-LOP);
    }
  }
  for (i=0; i<NP; i++) for (j=0; j<20; j++) P[HNEW][i][j]+=(P[HINI][i][j]);
  for (i=0; i<NP; i++) for (j=0; j<20; j++)
  if (P[HLOG][i][j]>-100) {
    sum+=P[HLOG][i][j]*(P[HINI][i][j]);
    NN+=(P[HINI][i][j]);
  }
  hmmtop_copy(HNEW,HCUR);
  hmmtop_norm();
  return exp(-sum/NN);
}

void hmmtop_optim() {
double p,rp,elt;
int i,j,k,sz=0,jo=0,r;
  hmmtop_setmem();
  hmmtop_setlocate();
  for (k=0,PM=1000; k<FUT; k++) {
    p=999; jo=0; sz=r=0; 
    hmmtop_copy(HINI,HCUR);
    if (START==RANDOM) {
      hmmtop_zero(HINI);
      for (i=0; i<NP; i++) for (j=0; j<20; j++) P[HCUR][i][j]=drand48();
      for (i=0; i<NT; i++) for (j=0; j<trdb[i]; j++) T[HCUR][i][j]=drand48();
      for (i=0; i<NS; i++) I[HCUR][i]=(state[i].ini==0?drand48():0);
      if (LOG!=NULL) fprintf(LOG,"random initialization\n");
    }
    hmmtop_norm(); 
    if (HMMTOP_ITERATION)
    while (jo==0) {
     	rp=p; sz++;
	    p=hmmtop_core(); elt=rp-p;
	    if (elt<LIMIT) jo=1;
	    if (elt<-1e-10) jo=0;
	    if (sz>IFUT) jo=1;
	if (LOG!=NULL) fprintf(LOG,"Core: %5d %20.8f%20.8f\n",sz,p,elt);
	if (elt<-1e-6) hmmtop_error("Serious error, delta probability negative");
    }
    if (FUT>1) {
      if (p<PM) {
        PM=p; hmmtop_copy(HCUR,HSAV);
      }
    }
    else PM=p;
  }
  if (FUT>1) hmmtop_copy(HSAV,HCUR);
}

void hmmtop_bestpath(int f) {
int i,j,k,h;
double q,max;

  if (LOG!=NULL) fprintf(LOG,"Entering bestpath\n");
  for (i=0; i<NS; i++) {
    E[i][0]=I[HLOG][i]+P[HLOG][state[i].e][seq[f][0]]-1e10*LOCATE[i][0];
    MH[i][0]=0;
  }
  for (i=1; i<n[f]; i++) for (j=0; j<NS; j++) {
    for (k=0, max=-1e30; k<state[j].nb; k++) {
      q=E[state[j].bt[k]][i-1]+T[HLOG][state[state[j].bt[k]].t][state[j].bm[k]];
      if (max<q) max=q,h=state[j].bt[k];
    }
    E[j][i]=max+P[HLOG][state[j].e][seq[f][i]]-1e10*LOCATE[j][i];
    MH[j][i]=h;
  }
  for (i=0,max=-1e30; i<NS; i++) {
    q=E[i][n[f]-1]+state[i].end;
    if (max<q) max=q,h=i; 
  }
  q=n[f];
  for (i=0; i<NP; i++) for (j=0; j<20; j++)
  if (P[HLOG][i][j]>-100)
  {	max+=P[HLOG][i][j]*(P[HINI][i][j]);
	  q+=(P[HINI][i][j]);
  }
  PMM=exp(-max/q); 
  for (i=n[f]-1; i>=0; i--) {
    prpart[i]=h; h=MH[h][i];
  }
  if (LOG!=NULL) fprintf(LOG,"Best path ok\n");
}

void hmmtop_output(FILE *OF) {
int i,j,ke,ve,td;

	if (LOG!=NULL)
	{	fprintf(LOG,">HP: %d %s\n",n[0],name[0]);
		fprintf(LOG,"Output file: %s\n",OF_NAME);
	}
	for (i=0,td=0; i<n[0]-1; i++)
	 	if (state[prpart[i]].let!='H' && state[prpart[i+1]].let=='H') td++;
	fprintf(OF,">HP: %d %s ",n[0],name[0]);
	if (state[prpart[0]].let=='o'||state[prpart[0]].let=='O') 
		fprintf(OF,"OUT "); else fprintf(OF," IN ");
	fprintf(OF,"%3d ",td);
	for (i=0; i<n[0]-1; i++)
	{	if (state[prpart[i]].let!='H' && state[prpart[i+1]].let=='H') fprintf(OF,"%4d ",i+2);
		if (state[prpart[i]].let=='H' && state[prpart[i+1]].let!='H') fprintf(OF,"%4d ",i+1);
	}
	fprintf(OF,"\n");
	
	if (OPL)
	{	fprintf(OF,"The best model:\n\n");
		ve=0;
        	while (ve<n[0])
        	{       ke=ve; ve+=50; if (ve>n[0]) ve=n[0];
                	fprintf(OF,"     seq  ");
                	for (j=ke; j<ve; j++) 
                	{       fprintf(OF,"%c",AS[seq[0][j]]);
                        	if ((int)((j+1)/10)*10==j+1) fprintf(OF," ");
                	}
                	fprintf(OF,"%5d\n",ve);
                	fprintf(OF,"     pred ");
                	for (j=ke; j<ve; j++) 
                	{       fprintf(OF,"%c",state[prpart[j]].let);
                        	if ((int)((j+1)/10)*10==j+1) fprintf(OF," ");
                	}
                	fprintf(OF,"\n\n");
        	}
	}
	if (OPP)
        {	fprintf(OF,"\nTotal entropy of the model: %8.4f\n",PM);
        	fprintf(OF,"Entropy of the best path: %8.4f\n\n",PMM);
        	for (i=0; i<NS; i++) fprintf(OF,"%8.4f",100.0*I[HCUR][i]);
		fprintf(OF,"\n");
		for (i=0; i<NP; i++) {
			for (j=0; j<20; j++)
   				fprintf(OF,"%8.4f",P[HCUR][i][j]*100);
        	    fprintf(OF,"\n");
        	}
		for (i=0; i<NT; i++) {
			for (j=0; j<trdb[i]; j++)
   				fprintf(OF,"%8.4f",T[HCUR][i][j]*100);
        	    fprintf(OF,"\n");
        	}
	}
	if (OPC)
        {	fprintf(OF,"\nSize of pseudocount vector:%5.0f\n",TM);
        	fprintf(OF,"Pseudocount vector used:\n");
        	for (i=0; i<NP; i++) fprintf(OF,"%8s",emilist[i]); fprintf(OF,"\n");
        	for (j=0; j<20; j++)
        	{   fprintf(OF,"%3c",AS[j]);
        	    for (i=0; i<NP; i++) fprintf(OF,"%8.2f",P[HINI][i][j]);
        	    fprintf(OF,"\n");
        	}
	}
}

void hmmtop_addlocate(char *s)
{
	if ((locates=realloc(locates,(locnum+1)*sizeof(char *)))==NULL)
		hmmtop_error("No enough memory for locates");
	if ((locates[locnum]=calloc(strlen(s),sizeof(char)))==NULL)
		hmmtop_error("No enough memory for locatesn");
	strcpy(locates[locnum],s);
	if (LOG!=NULL) fprintf(LOG,"Locate has been added: %d %s\n",locnum,locates[locnum]);
	locnum++;
}

void hmmtop_setlocate() {
int i,j,k,d,bpos,epos;
char s[100],spart;
	
	if (LOG!=NULL) fprintf(LOG,"Entering setlocate\n");
	for (i=0; i<NS; i++) for (j=0; j<n[0]; j++) LOCATE[i][j]=0;
	for (i=0; i<locnum; i++)
	{	if (strchr(locates[i],'-')!=NULL)
		{	j=(int)(strchr(locates[i],'-')-locates[i]);}
		else
		{	hmmtop_error("--locates=... syntaxis are incorrect (1)");}
		if (strchr(&locates[i][j+1],'-')!=NULL)
		{	k=(int)(strchr(&locates[i][j+1],'-')-&locates[i][j+1]);}
		else
		{	hmmtop_error("--locates=... syntaxis are incorrect (2)");}
		strncpy(s,locates[i],j); s[j]='\0';
		if (!strcmp(s,"E")) bpos=n[0]-1; else bpos=atoi(s)-1;
		strncpy(s,&locates[i][j+1],k); s[k]='\0';
		if (!strcmp(s,"E")) epos=n[0]; else epos=atoi(s);
		if (j+k+2>=strlen(locates[i]))
		{	hmmtop_error("--locates=... syntaxis are incorrect (3)");}
		spart=locates[i][j+k+2];
		d=0;
		for (j=0; j<NS; j++) if (state[j].let!=spart) {
		  for (k=bpos; k<epos; k++) LOCATE[j][k]=1;
		  d++;
		}
		if (d==NS) hmmtop_error("No such letter for locate");
	}
	if (LOG!=NULL) fprintf(LOG,"Setlocate OK\n");
}

void hmmtop_readarchitecture() {
FILE *f;
int numline,i,j,b,ll,k,l,trd,ns,num1,num2,db;
char line[200],s[200],ss[100];
char trlist[100][20];

  if (getenv("HMMTOP_ARCH")==NULL) {
    hmmtop_info("HMMTOP_ARCH is not set, trying in the current directory ...");
#ifdef sgi
    putenv("HMMTOP_ARCH=hmmtop.arch");
#else
    setenv("HMMTOP_ARCH","hmmtop.arch",0);
#endif
  }
  if ((f=fopen(getenv("HMMTOP_ARCH"),"r"))==NULL)
    hmmtop_error("Architecture file is not found");
    
  NS=0;
  numline=0;
  NT=0;
  NP=0;
  while (!feof(f)) {
    for (strcpy(line,""); (!feof(f)&&strstr(line,"state")==NULL); 
       numline++) fgets(line,200,f);
    if (strchr(line,'{')==NULL&&!feof(f)) {
      sprintf(s,"Syntax error in architecture file (1), line %d\n%s",numline,line);
      hmmtop_error(s);
    }
    if (!feof(f)) {
      for (i=0, j=(int)(strchr(line,'{')-line)+1; line[j]!='\n'; j++)
        if (line[j]!='\n'&&line[j]!=' ') s[i++]=line[j];
      s[i]='\0';
      state[NS].code=calloc(i+1,sizeof(char));
      strcpy(state[NS].code,s);
      state[NS].ini=state[NS].end=-1000;
    }
    for (;(!feof(f)&&strchr(line,'}')==NULL);
      fgets(line,200,f), numline++) 
    if (!feof(f)&&strchr(line,'{')==NULL&&line[0]!='#') {
      line[strlen(line)-1]='\0';
      ll=strlen(line);
      if (strchr(line,'=')!=NULL) {
        b=(int)(strchr(line,'=')-line)+1;
        for (j=b,i=0; (line[j]!=','&&j<ll); j++)
          if (line[j]!=','&&line[j]!=' ') s[i++]=line[j]; 
        s[i]='\0';
        if (strstr(line,"transition")!=NULL) {
          for (i=0; (i<NT&&strcmp(translist[i],s)); i++);
          if (i==NT) strcpy(translist[NT++],s);
          state[NS].t=i;
        }
        else if (strstr(line,"emission")!=NULL) {
          for (i=0; (i<NP&&strcmp(emilist[i],s)); i++);
          if (i==NP) strcpy(emilist[NP++],s);
          state[NS].e=i;
        }
        if (strstr(line,"name")!=NULL) {
          state[NS].name=calloc(strlen(line)-b+3,sizeof(char));
          strcpy(state[NS].name,&line[b]);
        }
        if (strstr(line,"letter")!=NULL) {
          if (strchr(line,'"')==NULL) {
            sprintf(s,"Syntax error: no \" in line %d\n%s",numline,line);
            hmmtop_error(s);
          }
          state[NS].let=line[(int)(strchr(line,'"')-line)+1];
        }
        if (strstr(line,"ends")!=NULL) {
          if (strstr(&line[b],"start")!=NULL) state[NS].ini=0;
          if (strstr(&line[b],"end")!=NULL) state[NS].end=0;
        }
      }
    }
    if (!feof(f)) NS++;
  }
  rewind(f);
  if (LOG!=NULL) fprintf(LOG,"Second reading\n");
  ns=numline=0;
  while (!feof(f)) {
    for (strcpy(line,""); (!feof(f)&&strstr(line,"state")==NULL); 
       numline++) fgets(line,200,f);
    if (!feof(f)) {
      for (i=0, j=(int)(strchr(line,'{')-line)+1; line[j]!='\n'; j++)
        if (line[j]!='\n'&&line[j]!=' ') s[i++]=line[j];
      s[i]='\0';
      if (strcmp(s,state[ns].code)) { 
        sprintf(s,"second read failed:%d\n%s",numline,line);
        hmmtop_error(s);
      }
    }
    for (;(!feof(f)&&strchr(line,'}')==NULL);
      fgets(line,200,f), numline++) 
    if (!feof(f)&&strchr(line,'{')==NULL&&line[0]!='#') {
      line[strlen(line)-1]='\0';
      ll=strlen(line);
      if (strchr(line,'=')!=NULL) {
        b=(int)(strchr(line,'=')-line)+1;
        if (strstr(line,"transition")!=NULL) {
          for (i=b,trd=j=0; i<=ll; i++) if (line[i]!=' ') {
            if (line[i]!=','&&i!=ll) s[j++]=line[i];
            if (line[i]==','||i==ll) {
              s[j]='\0';
              if (strchr(s,'-')!=NULL) {
                if (strchr(s,'.')==NULL) 
                  hmmtop_error("Syntax error in architecture file (state.num-num)");
                k=(int)(strchr(s,'.')-s)+1;
                l=(int)(strchr(s,'-')-s)+1;
                strncpy(ss,&s[k],l-k-1); ss[l-k-1]='\0';
                num1=atoi(ss);
                num2=atoi(&s[l]);
                strncpy(ss,s,k); ss[k]='\0';
                if (LOG!=NULL) fprintf(LOG,"Transition list found: %s %d %d\n",s,num1,num2);
                for (j=num1; j<=num2; j++, trd++) sprintf(trlist[trd],"%s%d",ss,j);
              }
              else strcpy(trlist[trd++],s);
              j=0;
            }
          }
          if (trdb[state[ns].t]==0) trdb[state[ns].t]=trd-1;
          if (trdb[state[ns].t]!=trd-1) hmmtop_error("Tanslist number wrong");
          state[ns].nf=trd-1;
          state[ns].f=calloc(trd,sizeof(int));
          if (LOG!=NULL) fprintf(LOG,"Transition list of state %d %s (.e=%c .t=%s): ",ns,state[ns].code,state[ns].let,translist[state[ns].t]);
          for (i=1; i<trd; i++) {
            for (j=0; (j<NS&&strcmp(state[j].code,trlist[i])); j++);
            if (j==NS) hmmtop_error("Transition state not found");
            state[ns].f[i-1]=j; 
            if (LOG!=NULL) fprintf(LOG," %d",j);
          }
          if (LOG!=NULL) fprintf(LOG,"\n");
        }
      }
    }
    if (!feof(f)) ns++;
  }

  fclose(f);
  if (LOG!=NULL) fprintf(LOG,"Architecture have been successfully read!\n");
  for (i=0; i<NS; i++) {
    for (j=0,db=0; j<NS; j++) for (k=0; k<state[j].nf; k++)
    if (state[j].f[k]==i) db++;
    if ((state[i].bt=calloc(db,sizeof(int)))==NULL) hmmtop_error("Can't alloc state[i].bt");
    if ((state[i].bm=calloc(db,sizeof(int)))==NULL) hmmtop_error("Can't alloc state[i].bm");
    for (j=0,db=0; j<NS; j++) for (k=0; k<state[j].nf; k++)
    if (state[j].f[k]==i) {
      state[i].bt[db]=j;
      state[i].bm[db]=k;
      db++;
    }
    state[i].nb=db;
  }

  if (LOG!=NULL) fprintf(LOG,"Allocating memory for the given architecture\n");
  if ((I=calloc(ATA,sizeof(double*)))==NULL) hmmtop_error("Can't alloc I");
  for (i=0; i<ATA; i++)
  if ((I[i]=calloc(NS,sizeof(double)))==NULL) hmmtop_error("Can't alloc I[i]");
  if ((P=calloc(ATA,sizeof(double**)))==NULL) hmmtop_error("Can't alloc P");
  for (i=0; i<ATA; i++) {
    if ((P[i]=calloc(NP,sizeof(double*)))==NULL) hmmtop_error("Can't alloc P[i]");
    for (j=0; j<NP; j++)
      if ((P[i][j]=calloc(20,sizeof(double)))==NULL) hmmtop_error("Can't alloc P[i][j]");
  }
  if ((T=calloc(ATA,sizeof(double**)))==NULL) hmmtop_error("Can't alloc T");
  for (i=0; i<ATA; i++) {
    if ((T[i]=calloc(NT,sizeof(double*)))==NULL) hmmtop_error("Can't alloc T[i]");
    for (j=0; j<NT; j++)
      if ((T[i][j]=calloc(trdb[j],sizeof(double)))==NULL) hmmtop_error("Can't alloc T[i][j]");
  }
  E=calloc(NS,sizeof(double*));
  H=calloc(NS,sizeof(double*));
  MH=calloc(NS,sizeof(int*));
  LOCATE=calloc(NS,sizeof(int*));
  if (LOG!=NULL) fprintf(LOG,"Memory Allocation for the given architecture is done\n");
}

void hmmtop_readlocates() {
  char s[100],c;
  locnum=0;
  for (c=' ';(c!='\n'&&!feof(LOCF));c=getc(LOCF)) {
  	fscanf(LOCF,"%s",s);
	hmmtop_addlocate(s);
  }
  if (LOG!=NULL) fprintf(LOG,"Line from loc_file has been read\n");
}
