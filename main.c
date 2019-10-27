#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define imax 10000000

void prodMVet(double **M, double x[], double r[], int TAM)
{
    int i,j;
    for (i=0;i<TAM;i++)
    {
        double soma=0;
        for (j=0;j<TAM;j++)
            soma+=M[i][j]*x[j];
        r[i]=soma;
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

int gc(double **A, double *b, double *x, int TAM)
{
    double erro=0.00001;
    int i;
    for (i=0;i<TAM;i++) x[i]=0;
    int it=0; // iteracao corrente
    double *aux=(double *)malloc(TAM*sizeof(double));
    prodMVet(A,x,aux,TAM);
    double *r=(double *)malloc(TAM*sizeof(double)); // vetor do res�duo
    subVet(b,aux,r,TAM);
    double *d=(double *)malloc(TAM*sizeof(double));
    for (i=0;i<TAM;i++) d[i]=r[i];
    double sigma_novo=prodEsc(r,r,TAM);
    double sigma0=sigma_novo;
    while (it<imax && sigma_novo>erro*erro*sigma0)
    {
        if(it % 500 == 0)
            printf("it=%d sigma=%lf\n",it,sigma_novo);
        double *q=(double *)malloc(TAM*sizeof(double));
        prodMVet(A,d,q,TAM);
        double Alpha=sigma_novo/prodEsc(d,q,TAM);
        for (i=0;i<TAM;i++) x[i]+=Alpha*d[i];
        if (it%50==0)
        {
            prodMVet(A,x,aux,TAM);
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
//    printf("x[%d]=%lf\n",n-1,x[n-1]);
    for(i=n-2; i>=0; i--)
    {
        double sum=0;
        for(j=i+1; j<n; j++)
            sum=sum+A[i][j]*x[j];
        x[i]=(b[i]-sum)/A[i][i];
    }
    return 1;
}

int main(){
    FILE *f;
    int i, j;

    FILE *arqnomes=fopen("arquivos.txt","rt");
    if (arqnomes==NULL){printf("Erro na abertura do arquivo de nomes\n");return 0;}
    char nomearq[20];
    fscanf(arqnomes,"%s",nomearq);
    while (strcmp(nomearq,"fim")) {
//        getchar();
        FILE *arqin = fopen(nomearq, "rt");
        char linha[100];

        printf("Analisando arquivo %s\n",nomearq);

        fgets(linha, 100, arqin);
//        printf("%s\n", linha);

        int TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD;
        char MXTYPE[3];
        int NROW, NCOL, NNZERO, NELTVL;
        char PTRFMT[10], INDFMT[10], VALFMT[10], RHSFMT[10];

        fscanf(arqin, "%d %d %d %d %d", &TOTCRD, &PTRCRD, &INDCRD, &VALCRD, &RHSCRD);
        fscanf(arqin, "%s %d %d %d %d", MXTYPE, &NROW, &NCOL, &NNZERO, &NELTVL);
        fscanf(arqin, "%s %s %s", PTRFMT, INDFMT, VALFMT);

        //
        int n_linhas1, n_elementos1;
        sscanf(PTRFMT, "(%dI%d)", &n_linhas1, &n_elementos1);

        int n_linhas2, n_elementos2;
        sscanf(INDFMT, "(%dI%d)", &n_linhas2, &n_elementos2);

        int n_elementos_linhas, n_elementos_coluna, tam_elemento;
        sscanf(VALFMT, "(%dE%d.%d)", &n_elementos_linhas, &n_elementos_coluna, &tam_elemento);

        int *rows = (int *)calloc(NROW,sizeof(int));
        int *columns = (int *)calloc(NNZERO,sizeof(int));
        double *entries = (double *)calloc(NNZERO,sizeof(double));

//        printf("LINHAS\n");
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
            rows[k] = int_aux;
            rows[k]--;
            k++;
            p += n_elementos1;
        }
        k = 0;
        //printf("COLUNAS\n");
//        for(int i =0;i<NROW+1;i++){
//            printf("%d ",rows[i]);
//            fflush(stdout);
//        }

        for (i = 0; i < NNZERO; i++) {
            if (i % n_linhas2 == 0) {
                fgets(line, 100, arqin);
                p = line;
            }
            strncpy(staux, p, n_elementos2);
            staux[n_elementos2] = '\0';
            sscanf(staux, "%d", &int_aux);
            columns[k] = int_aux;
            columns[k]--;
            k++;
            p += n_elementos2;
        }
        k = 0;
        //        printf("Columns\n");
//        for(int i = 0;i<NNZERO;i++){
//            printf("%d ",columns[i]);
//            fflush(stdout);
//        }

        char staux2[n_elementos_coluna];
        double aux2;
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

        double **matriz=(double **)malloc(NROW*sizeof(double *));

        for(int i = 0;i < NROW;i++){
            matriz[i] = (double *)calloc(NCOL, sizeof(double));
        }
        int a=0;

        int contEntries = 0;
        for (i = 0; i < NROW + 1; i++) {
            int limite = 0;
            if (i != NROW) {
                limite = rows[i + 1];
            }
            int finalRow = i + 1;
            for (j = rows[i]; j < limite; j++) {
                matriz[i][columns[j]] = entries[contEntries];
                matriz[columns[j]][i] = entries[contEntries]; //como a matriz é simétrica, coloca tanto no [i][j] quanto no [j][i]
                contEntries++;
            }
        }


//        printf("MATRIZ \n");
//        for(int i=0;i<NROW;i++){
//            for(int j=0;j<NCOL;j++){
//                printf("%lf ",matriz[i][j]);
//                fflush(stdout);
//            }
//            printf("\n");
//        }

        double *xgc=(double *)calloc(NROW,sizeof(double));
        double *xgauss=(double *)calloc(NROW,sizeof(double));
        double *b=(double *)calloc(NROW,sizeof(double));
        for (i=0;i<NROW;i++) b[i]=1.0;

        double *aux=(double *)calloc(NROW,sizeof(double));

//        for(int i = 0;i < TOTCRD;i++){
//            for(int j = 0;j<n_elementos_linhas;j++){
//                printf("%ld ",matriz[i][j]);
//            }
//            printf("\n");
//        }

        int itera=gc(matriz,b,xgc,NROW);
        prodMVet(matriz,xgc,aux,NROW);
        subVet(aux,b,aux,NROW);
        double erro=prodEsc(aux,aux,NROW);
//        printf("Erro do gc: %lf\n",erro);
        if (itera==imax)
            printf("%s nao convergiu em %d iteracoes\n",nomearq,imax);
        else
            printf("%s convergiu em %d iteracoes\n\n",nomearq,itera);
//        getchar();
        if (!gauss(matriz,b,xgauss,NROW)){printf("Deu pau no Gauss!!!\n");return 0;}
//        for (i=0;i<NROW;i++)
//            printf("xgc[%d]=%lf xgauss[%d]=%lf\n",i,xgc[i],i,xgauss[i]);

        fclose(arqin);
        fscanf(arqnomes,"%s\n",nomearq);
    }

}
