#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define imax 10000000
#define nThreads 2

void prodMVet(double *x, double *r, int TAM, int NNZERO, int *columns, int *rows, double *entries) { //versao paralelizada
	int i, j, id;
	double *r_private[nThreads];
	omp_set_num_threads(nThreads);

	for(i=0;i<TAM;i++){
		r[i]=0;
	}
	for(i = 0;i < nThreads;i++){
		r_private[i] = (double *)malloc(TAM*sizeof(double));
	}

	for(i = 0;i < nThreads;i++){
		for(j=0; j<TAM; j++) {
		    r_private[i][j] = 0;
		}
	}

	#pragma omp parallel private(j,i,id)
	{
		id = omp_get_thread_num();
		for (i = id; i < TAM; i+=nThreads) {
			for (j = columns[i]; j < columns[i+1]; j++) {
			    r_private[id][rows[j]] += entries[j] * x[i];
			    if (rows[j] != i) {
				r_private[id][i] += entries[j] * x[rows[j]];
			    }
			}
		}
	}


	for(i = 0;i < nThreads;i++){
		for(j=0; j<TAM; j++) {
		    r[j] += r_private[i][j];
		}
		free(r_private[i]);
	}
}

void prodMVet_seq(double *x, double *r, int TAM, int NNZERO, int *columns, int *rows, double *entries) { //versao sequencial
    int i, j;

    for(i=0;i<TAM;i++){
        r[i]=0;
    }

    for (i = 0; i < TAM; i++) {
        for (j = columns[i]; j < columns[i+1]; j++) {
            if (rows[j] == i) {
                r[rows[j]] += entries[j] * x[i];
            } else {
                r[rows[j]] += entries[j] * x[i];
                r[i] += entries[j] * x[rows[j]];
            }
        }
    }
}


void prodMVet_sem(double *x, double *r, int TAM, int NNZERO, int *columns, int *rows, double *entries) { //versao com semaforos
	int i, j, tid;

	for(i=0;i<TAM;i++){
		r[i]=0;
	}

	omp_lock_t simple_lock[TAM];
	for(i=0; i<TAM; i++){
		omp_init_lock(&simple_lock[i]);
	}
	omp_set_num_threads(2);

	#pragma omp parallel for private (i,j)
	for (i = 0; i < TAM; i++) {
		for (j = columns[i]; j < columns[i+1]; j++) {
		    omp_set_lock(&simple_lock[rows[j]]);
		    r[rows[j]] += entries[j] * x[i];
		    omp_unset_lock(&simple_lock[rows[j]]);	
		    if (rows[j] != i) {
 		        omp_set_lock(&simple_lock[i]);
			r[i] += entries[j] * x[rows[j]];

        		omp_unset_lock(&simple_lock[i]);			
		    }
		}
	}

    	for(i=0; i<TAM; i++){
		omp_destroy_lock(&simple_lock[i]);
    	}
}

void subVet(double a[], double b[], double r[], int TAM)
{
    int i;
    for (i=0;i<TAM;i++) r[i]=a[i]-b[i];
}

double prodEsc(double a[], double b[], int TAM)
{
    double soma=0;
    int i;
    for (i=0;i<TAM;i++)
        soma+=a[i]*b[i];
    return soma;
}

int gc(double *b, double *x, int TAM, int NNZERO, int *rows, int *columns, double *entries)
{
    double erro=0.00001;
    int i;
    for (i=0;i<TAM;i++) x[i]=0;
    int it=0; // iteracao corrente
    double *aux=(double *)malloc(TAM*sizeof(double));
    double *r=(double *)malloc(TAM*sizeof(double));
    subVet(b,aux,r,TAM);
    double *d=(double *)malloc(TAM*sizeof(double));
    for (i=0;i<TAM;i++) d[i]=r[i];
    double sigma_novo=prodEsc(r,r,TAM);
    double sigma0=sigma_novo;
    double *q=(double *)malloc(TAM*sizeof(double));
    while (it<imax && sigma_novo>erro*erro*sigma0)
    {
         //if(it % 500 == 0)
           //  printf("it=%d sigma=%lf\n",it,sigma_novo);
        prodMVet(d,q,TAM, NNZERO, rows, columns, entries);
        double Alpha=sigma_novo/prodEsc(d,q,TAM);
        for (i=0;i<TAM;i++) x[i]+=Alpha*d[i];
        if (it%50==0)
        {
            prodMVet(x,aux,TAM,NNZERO,rows,columns,entries);
            for (i=0;i<TAM;i++) r[i]=b[i]-aux[i];
        }
        else
            for (i=0;i<TAM;i++) r[i]-=Alpha*q[i];
        double sigma_velho = sigma_novo;
        sigma_novo = prodEsc(r,r,TAM);
        double beta = sigma_novo / sigma_velho;
        for (i=0;i<TAM;i++) d[i]=r[i]+beta*d[i];
        it=it+1;
    }
    return it;
}

void mostra(double **A, int TAM)
{
    int i,j;
    for (i=0;i<TAM;i++)
    {
        for (j=0;j<TAM;j++)
            printf("%10lf ",A[i][j]);
        printf("\n");
    }
}

int gauss(double **A, double *b, double *x, int n)
{
    int i,j,k;
    for(j=0;j<n-1; j++)
    {
        for(i=j+1; i<n; i++)
        {
            if (A[j][j]==0) return 0;
            double c=A[i][j]/A[j][j];
            for(k=0; k<n; k++)
                A[i][k]=A[i][k]-c*A[j][k];
            b[i]=b[i]-c*b[j];
        }
    }

    if (A[n-1][n-1]==0) return 0;
    x[n-1]=b[n-1]/A[n-1][n-1];
    for(i=n-2; i>=0; i--)
    {
        double sum=0;
        for(j=i+1; j<n; j++)
            sum=sum+A[i][j]*x[j];
        x[i]=(b[i]-sum)/A[i][i];
    }
    return 1;
}

int main() {
    FILE *f;
    int i, j;

    FILE *arqnomes = fopen("arquivos.txt", "rt");
    if (arqnomes == NULL) {
        printf("Erro na abertura do arquivo de nomes\n");
        return 0;
    }
    char nomearq[20];
    fscanf(arqnomes, "%s", nomearq);
    while (strcmp(nomearq, "fim")) {
        FILE *arqin = fopen(nomearq, "rt");
        char linha[100];

        printf("Analisando arquivo %s\n", nomearq);

        fgets(linha, 100, arqin);

        int TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD;
        char MXTYPE[3];
        int NROW, NCOL, NNZERO, NELTVL;
        char PTRFMT[10], INDFMT[10], VALFMT[10], RHSFMT[10];

        fscanf(arqin, "%d %d %d %d %d", &TOTCRD, &PTRCRD, &INDCRD, &VALCRD, &RHSCRD);
        fscanf(arqin, "%s %d %d %d %d", MXTYPE, &NROW, &NCOL, &NNZERO, &NELTVL);
        fscanf(arqin, "%s %s %s", PTRFMT, INDFMT, VALFMT);

        int n_linhas1, n_elementos1;
        sscanf(PTRFMT, "(%dI%d)", &n_linhas1, &n_elementos1);

        int n_linhas2, n_elementos2;
        sscanf(INDFMT, "(%dI%d)", &n_linhas2, &n_elementos2);

        int n_elementos_linhas, n_elementos_coluna, tam_elemento;
        sscanf(VALFMT, "(%dE%d.%d)", &n_elementos_linhas, &n_elementos_coluna, &tam_elemento);

        int *columns = (int *) calloc(NROW+1, sizeof(int));
        int *rows = (int *) calloc(NNZERO, sizeof(int));
        double *entries = (double *) calloc(NNZERO, sizeof(double));

        char line[100];
        char dest[6];
        int k = 0;
        fgets(line, 100, arqin);

        char st[100], staux[n_elementos1 + 1];
        char *p = line;
        int int_aux;

        for (i = 0; i < NROW+1; i++) {
            if (i % n_linhas1 == 0) {
                fgets(line, 100, arqin);
                p = line;
            }
            strncpy(staux, p, n_elementos1);
            staux[n_elementos1] = '\0';
            sscanf(staux, "%d", &int_aux);
            columns[k] = int_aux;
            columns[k]--;
            k++;
            p += n_elementos1;
        }
        k = 0;

        for (i = 0; i < NNZERO; i++) {
            if (i % n_linhas2 == 0) {
                fgets(line, 100, arqin);
                p = line;
            }
            strncpy(staux, p, n_elementos2);
            staux[n_elementos2] = '\0';
            sscanf(staux, "%d", &int_aux);
            rows[k] = int_aux;
            rows[k]--;
            k++;
            p += n_elementos2;
        }
        k = 0;

        char staux2[n_elementos_coluna];
        double aux2 = 0;
        for (i = 0; i < NNZERO; i++) {
            if (i % n_elementos_linhas == 0) {
                fgets(line, 100, arqin);
                p = line;
            }

            strncpy(staux2, p, n_elementos_coluna);
            staux[n_elementos_coluna] = '\0';
            sscanf(staux2, "%lf", &aux2);
            entries[i] = aux2;
            p += n_elementos_coluna;
        }

        double *xgc = (double *) calloc(NROW, sizeof(double));
        double *xgauss = (double *) calloc(NROW, sizeof(double));
        double *b = (double *) calloc(NROW, sizeof(double));
        for (i = 0; i < NROW; i++) b[i] = 1.0;

        double *aux = (double *) calloc(NROW, sizeof(double));

        float ini, fim;
        ini = omp_get_wtime();
        int itera = gc(b, xgc, NROW, NNZERO, columns, rows, entries);
        subVet(aux, b, aux, NROW);
        double erro = prodEsc(aux, aux, NROW);
        fim = omp_get_wtime();
	printf("Tempo: %f\n",fim-ini);
        if (itera == imax)
            printf("%s nao convergiu em %d iteracoes\n", nomearq, imax);
        else
            printf("%s convergiu em %d iteracoes\n\n", nomearq, itera);

        fclose(arqin);
        fscanf(arqnomes, "%s\n", nomearq);
    }
}
