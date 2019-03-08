#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include <mpi.h>

#define PI 3.14159265
#define tiny 1e-16
#define fmin 0.1
#define fmax 100.0
#define send_data_tag 2001
#define finish_data_tag 2002
#define pcyrtokms 0.9778

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


char* concat(const char *s1, const char *s2)
  {
      char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the null-terminator
      //in real code you would check for errors in malloc here
      strcpy(result, s1);
      strcat(result, s2);
      return result;
  }

double uniform(double minimo, double maximo){
  return minimo+(maximo-minimo)*rand()/(1.0*RAND_MAX);
}

double walker(char *name, char *snapshot,double n1_min , double n1_max, double pc_min, double pc_max,
  double dv_min, double dv_max,
  double pmin, double pmax, int Nbins, double dx,float *sigma, float *delta, float *error, int steps, int id){
        printf("This is the %03d\t sampler \n ", id);


    FILE *output;
    char outname[sizeof "_third100.txt"];
    sprintf(outname, "_third%03d.txt", id);
    char* outpar = concat("Files_model_2/", name);
    outpar       = concat(outpar,outname);
    output       = fopen(outpar,"w");

    FILE *arx;
    char* parameters;
    char filename[sizeof "_sampler100.txt"];

    float pc,pc_b;
  	double n1  , fc  , dv  , prob  ;
    double n1_b, fc_b, dv_b, prob_b, accept_rate;
	double alpha, rnd, numerator, denominator;
  	int i,counter,counter_b, ntrial, nacc;
    int i1,i2;
    int l1,l2;
    int u1,u2;
    int number_file;
  	int N=60;
    double aux;
    double* tab_1 = malloc(Nbins*sizeof(double*));
    double* tab_2 = malloc(Nbins*sizeof(double*));
    double* tab_3 = malloc(Nbins*sizeof(double*));
    double* tab_4 = malloc(Nbins*sizeof(double*));
    double* temp  = malloc(Nbins*sizeof(double*));


    double mu1,mu2,mu3;
    double sd1,sd2,sd3;
    double n11,pc1;
    double n12,pc2;
    double n13,pc3;
    double n14,pc4;

    float v1,v2,v3,v4,vtot;

    char line[700];
    int iline;

	/***********************
	 * First random numbers
	************************/

    n1=uniform(n1_min,n1_max);
    pc=uniform(pc_min,pc_max);
    dv=uniform(dv_min,dv_max);
    counter=0;
    ntrial=0;
    nacc=0;

    mu1=0.0*(n1_min+n1_max);
    mu2=0.0*(pc_min+pc_max);
    mu3=0.0*(dv_min+dv_max);
    sd1=0.00*(n1_max-n1_min);
    sd2=0.00*(pc_max-pc_min);
    sd3=0.00*(dv_max-dv_min);

    counter_b=-100;

    //printf("Start %03d\t th sampler \n ", id);
    while(counter<steps){
        ntrial++;
        if(counter<200){
            if(counter%3==0)n1=uniform(n1_min,n1_max);
            if(counter%3==1)pc=uniform(pc_min,pc_max);
            if(counter%3==2)dv=uniform(dv_min,dv_max);}

        else{
            do{n1=n1_b+0.6*sqrt(-2.0*sd1*log(uniform(0,1)))*cos(6.2831853*uniform(0,1));}while((n1<n1_min)||(n1>n1_max));
            do{pc=pc_b+0.6*sqrt(-2.0*sd2*log(uniform(0,1)))*cos(6.28318530718*uniform(0,1));}while((pc<pc_min)||(pc>pc_max));
            do{dv=dv_b+0.6*sqrt(-2.0*sd3*log(uniform(0,1)))*cos(6.2831853*uniform(0,1));}while((dv<dv_min)||(dv>dv_max));
            }

        aux= 1.0*(N-1)*1.0*(n1-n1_min)/(n1_max-n1_min);

        l1 = (int) aux;
        u1 = l1+1;

        number_file = l1/3;

        if(l1%3==0) i1 = 0;
        if(l1%3==1) i1 = N;
        if(l1%3==2) i1 = 2*N;

        aux = (N-1)*1.0*(log(1.0*pc/pc_min)/log(1.0*pc_max/pc_min)) +1;

        l2 = (int) aux;
        i2 = l2;

        sprintf(filename, "_sampler%03d.txt", number_file);
        parameters = concat("Files_model_2/", snapshot);
        parameters       = concat(parameters,filename);
        arx = fopen(parameters,"r");

        iline = 0;
        while (fgets(line, sizeof(line), arx)) {
            iline++;
            //printf("ints %d\t %d\n",iline,i1+i2+i3);

            //// PRIMER NODO ////

            if(iline == i1+i2){
            sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                    &n11,&pc1,&tab_1[0],&tab_1[1],&tab_1[2],&tab_1[3],&tab_1[4],&tab_1[5],&tab_1[6],&tab_1[7],&tab_1[8],&tab_1[9],
                    &tab_1[10],&tab_1[11],&tab_1[12],&tab_1[13],&tab_1[14],&tab_1[15],&tab_1[16],&tab_1[17],&tab_1[18],&tab_1[19],
                    &tab_1[20],&tab_1[21],&tab_1[22],&tab_1[23],&tab_1[24],&tab_1[25],&tab_1[26],&tab_1[27],&tab_1[28],&tab_1[29],
                    &tab_1[30],&tab_1[31],&tab_1[32],&tab_1[33],&tab_1[34],&tab_1[35],&tab_1[36],&tab_1[37],&tab_1[38],&tab_1[39],
                    &tab_1[40],&tab_1[41],&tab_1[42],&tab_1[43],&tab_1[44],&tab_1[45],&tab_1[46],&tab_1[47],&tab_1[48],&tab_1[49],
                    &tab_1[50]);}

            //// SEGUNDO NODO ////

            if(iline == i1+i2+1){
            sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                    &n12,&pc2,&tab_2[0],&tab_2[1],&tab_2[2],&tab_2[3],&tab_2[4],&tab_2[5],&tab_2[6],&tab_2[7],&tab_2[8],&tab_2[9],
                    &tab_2[10],&tab_2[11],&tab_2[12],&tab_2[13],&tab_2[14],&tab_2[15],&tab_2[16],&tab_2[17],&tab_2[18],&tab_2[19],
                    &tab_2[20],&tab_2[21],&tab_2[22],&tab_2[23],&tab_2[24],&tab_2[25],&tab_2[26],&tab_2[27],&tab_2[28],&tab_2[29],
                    &tab_2[30],&tab_2[31],&tab_2[32],&tab_2[33],&tab_2[34],&tab_2[35],&tab_2[36],&tab_2[37],&tab_2[38],&tab_2[39],
                    &tab_2[40],&tab_2[41],&tab_2[42],&tab_2[43],&tab_2[44],&tab_2[45],&tab_2[46],&tab_2[47],&tab_2[48],&tab_2[49],
                    &tab_2[50]);}



            if(iline == i1+i2+1)break;

            }

        fclose(arx);

        number_file = (l1+1)/3;

        if((l1+1)%3==0) i1 = 0;
        if((l1+1)%3==1) i1 = N;
        if((l1+1)%3==2) i1 = 2*N;

        sprintf(filename, "_sampler%03d.txt", number_file);
        parameters = concat("Files_model_2/", snapshot);
        parameters = concat(parameters,filename);
        arx = fopen(parameters,"r");
        iline = 0;

        /*if(id==0){
            printf("filename %s\n", parameters);
            printf("indices %d\t %d\t",i1,i2);}
        */

        while (fgets(line, sizeof(line), arx)) {
            iline++;
            //printf("ints %d\t %d\n",iline,i1+i2+i3);

            //// PRIMER NODO ////

            if(iline == i1+i2){
            sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                    &n13,&pc3,&tab_3[0],&tab_3[1],&tab_3[2],&tab_3[3],&tab_3[4],&tab_3[5],&tab_3[6],&tab_3[7],&tab_3[8],&tab_3[9],
                    &tab_3[10],&tab_3[11],&tab_3[12],&tab_3[13],&tab_3[14],&tab_3[15],&tab_3[16],&tab_3[17],&tab_3[18],&tab_3[19],
                    &tab_3[20],&tab_3[21],&tab_3[22],&tab_3[23],&tab_3[24],&tab_3[25],&tab_3[26],&tab_3[27],&tab_3[28],&tab_3[29],
                    &tab_3[30],&tab_3[31],&tab_3[32],&tab_3[33],&tab_3[34],&tab_3[35],&tab_3[36],&tab_3[37],&tab_3[38],&tab_3[39],
                    &tab_3[40],&tab_3[41],&tab_3[42],&tab_3[43],&tab_3[44],&tab_3[45],&tab_3[46],&tab_3[47],&tab_3[48],&tab_3[49],
                    &tab_3[50]);}

            //// SEGUNDO NODO ////

            if(iline == i1+i2+1){
            sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                    &n14,&pc4,&tab_4[0],&tab_4[1],&tab_4[2],&tab_4[3],&tab_4[4],&tab_4[5],&tab_4[6],&tab_4[7],&tab_4[8],&tab_4[9],
                    &tab_4[10],&tab_4[11],&tab_4[12],&tab_4[13],&tab_4[14],&tab_4[15],&tab_4[16],&tab_4[17],&tab_4[18],&tab_4[19],
                    &tab_4[20],&tab_4[21],&tab_4[22],&tab_4[23],&tab_4[24],&tab_4[25],&tab_4[26],&tab_4[27],&tab_4[28],&tab_4[29],
                    &tab_4[30],&tab_4[31],&tab_4[32],&tab_4[33],&tab_4[34],&tab_4[35],&tab_4[36],&tab_4[37],&tab_4[38],&tab_4[39],
                    &tab_4[40],&tab_4[41],&tab_4[42],&tab_4[43],&tab_4[44],&tab_4[45],&tab_4[46],&tab_4[47],&tab_4[48],&tab_4[49],
                    &tab_4[50]);}


            if(iline == i1+i2+1)break;

            }

        fclose(arx);
        v1=1.0/fabs((n1-n11)*(pc-pc1));
        v2=1.0/fabs((n1-n12)*(pc-pc2));
        v3=1.0/fabs((n1-n13)*(pc-pc3));
        v4=1.0/fabs((n1-n14)*(pc-pc4));
        vtot=v1+v2+v3+v4;
        v1/=vtot;
        v2/=vtot;
        v3/=vtot;
        v4/=vtot;

        //printf("n1 %lf\t %lf\t %lf\t %lf\t %lf\n",n1,n11,n12,n13,n14);
        //printf("pc %lf\t %lf\t %lf\t %lf\t %lf\n",pc,pc1,pc2,pc3,pc4);

        for(i=0;i<Nbins;i++)temp[i]=v1*tab_1[i]+v2*tab_2[i]+v3*tab_3[i]+v4*tab_4[i];

        numerator=0.0;
        denominator=0.0;

        for(i=0;i<Nbins;i++){
            numerator+=sigma[i]*sigma[i]/(error[i]*error[i]);
            denominator+=temp[i]*sigma[i]/(error[i]*error[i]);}
        fc=numerator/denominator;
        fc+=dv;

        prob=0.0;

        for(i=0;i<Nbins;i++){
            prob=prob-0.5*(pow((sigma[i]-fc*temp[i])/error[i],2)+2*log(error[i]));}

        if(counter==0)prob_b=prob-1.0e5;
        alpha=prob-prob_b;

        if (alpha>0.0){
            prob_b=prob;
            n1_b=n1;
            pc_b=pc;
            fc_b=fc;
            dv_b=dv;
            if(counter_b>=0){
            sd1 = (1.0*counter_b/(1.0*counter_b+1))*(sd1 +(n1_b-mu1)*(n1_b-mu1)/(1.0*counter_b+1));
            sd2 = (1.0*counter_b/(1.0*counter_b+1))*(sd2 +(pc_b-mu2)*(pc_b-mu2)/(1.0*counter_b+1));
            sd3 = (1.0*counter_b/(1.0*counter_b+1))*(sd3 +(dv_b-mu3)*(dv_b-mu3)/(1.0*counter_b+1));

            mu1 = (1.0*counter_b*mu1+n1_b)/(1.0*counter_b+1);
            mu2 = (1.0*counter_b*mu2+pc_b)/(1.0*counter_b+1);
            mu3 = (1.0*counter_b*mu3+dv_b)/(1.0*counter_b+1);}
            counter_b++;
            counter++;
            nacc++;

            if(counter>1000)fprintf(output, "%f\t %f\t %10.8f\t %f\t %f\n", n1, pc, pcyrtokms*fc, accept_rate);}
        else{
            rnd=log(uniform(0.0,1.0));
            if(rnd<alpha){
                prob_b=prob;
                n1_b=n1;
                pc_b=pc;
                fc_b=fc;
                dv_b=dv;

                if(counter_b>=0){
                sd1 = (1.0*counter_b/(1.0*counter_b+1))*(sd1 +(n1_b-mu1)*(n1_b-mu1)/(1.0*counter_b+1));
                sd2 = (1.0*counter_b/(1.0*counter_b+1))*(sd2 +(pc_b-mu2)*(pc_b-mu2)/(1.0*counter_b+1));
                sd3 = (1.0*counter_b/(1.0*counter_b+1))*(sd3 +(dv_b-mu3)*(dv_b-mu3)/(1.0*counter_b+1));

                mu1 = (1.0*counter_b*mu1+n1_b)/(1.0*counter_b+1);
                mu2 = (1.0*counter_b*mu2+pc_b)/(1.0*counter_b+1);
                mu3 = (1.0*counter_b*mu3+dv_b)/(1.0*counter_b+1);}

                counter++;
                counter_b++;
                nacc++;
            if(counter>1000)fprintf(output, "%f\t %f\t %10.8f\t %f\t %f\n", n1, pc, pcyrtokms*fc, accept_rate);

            }}

        accept_rate=1.0*nacc/(1.0*ntrial);
        //if(counter>200){
        //if(id==0)printf("acceptance rate %f\n",accept_rate);
        //if(id==0)printf("%f\t %f\t %.15f\t %f\t %f\t %d\n",sd1,sd2,sd3,sd4, pc_b,counter_b);}



    }


	//return 0.0;

}


main(int argc, char* argv[]) {
MPI_Status status;
int my_id, an_id, root_process, ierr, i, num_procs,rec_id, fin;
root_process=0;
fin=0;
int N, Nres, Nbins, steps;
double L, out;
double pmin,pmax,dx,n1_min,n1_max,pc_min,pc_max,dv_min,dv_max;
char *name, *snapshot;
//clock_t start, end;
//double cpu_time_used;


name=argv[1];
snapshot=argv[2];
sscanf(argv[3],"%lf",&L);
sscanf(argv[4],"%d",&N);
sscanf(argv[5],"%lf",&n1_min);
sscanf(argv[6],"%lf",&n1_max);
sscanf(argv[7],"%lf",&pc_min);
sscanf(argv[8],"%lf",&pc_max);
sscanf(argv[9],"%lf",&dv_min);
sscanf(argv[10],"%lf",&dv_max);

Nbins=(argc-11)/3;
//printf("Nbins %d", Nbins);
steps=10000;
float sigma[Nbins];
float delta[Nbins];
float error[Nbins];

for(i = 0; i < Nbins; i++){
 sscanf(argv[11+i],"%f",&delta[i]);
 sscanf(argv[11+Nbins+i],"%f",&sigma[i]);
 sscanf(argv[11+2*Nbins+i],"%f",&error[i]);
}

pmin=4.0/L;
dx=L/(1.0*N);
pmax=1.0/(4.0*dx);

ierr = MPI_Init(&argc, &argv);

ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);    // Get id of the current processor
ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);// Get number of processes
srand(time(NULL)+ my_id);

if(my_id==root_process){

	// Send respective id
	for(an_id=1 ; an_id<num_procs ; an_id++){
    		ierr = MPI_Send(&an_id, 1, MPI_INT, an_id, send_data_tag, MPI_COMM_WORLD);
	}

	out=walker(name,snapshot,n1_min , n1_max, pc_min, pc_max, dv_min, dv_max,
	  	pmin, pmax, Nbins, dx, sigma, delta, error, steps, 0);
	for(an_id=1 ; an_id<num_procs ; an_id++)
	ierr = MPI_Recv(&fin, 1, MPI_INT, MPI_ANY_SOURCE, finish_data_tag, MPI_COMM_WORLD, &status);
}

else{
    	ierr = MPI_Recv(&rec_id, 1, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD, &status);
	out=walker(name,snapshot,n1_min , n1_max, pc_min, pc_max, dv_min, dv_max,
                pmin, pmax, Nbins, dx, sigma, delta, error, steps, rec_id);

	ierr = MPI_Send(&fin, 1, MPI_INT, 0, finish_data_tag, MPI_COMM_WORLD);
	}
ierr=MPI_Finalize();
return 0;
}
