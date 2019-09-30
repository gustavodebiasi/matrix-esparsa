#include <stdio.h>
#include <stdlib.h>

int main(){
    FILE *f;
    int i, j;
    //vetor de solução começa com um chute inicial
    //haverá também o vetor de termos independentes, verificar com cadinho como será preenchido

    f = fopen("bcsstk01.rsa", "rt");

    char linha[100];
    fgets(linha,100,f);
    printf("%s\n",linha);

    int TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD;
    char MXTYPE[3];
    int NROW, NCOL, NNZERO, NELTVL;
    char PTRFMT[10], INDFMT[10], VALFMT[10], RHSFMT[10];

    fscanf(f, "%d %d %d %d %d", &TOTCRD, &PTRCRD, &INDCRD, &VALCRD, &RHSCRD);
    fscanf(f, "%s %d %d %d %d", MXTYPE, &NROW, &NCOL, &NNZERO, &NELTVL);
    fscanf(f, "%s %s %s", PTRFMT, INDFMT, VALFMT);

    //alocação dos vetores
    int *rows = (int)malloc((NROW+1)*sizeof(int));
    int *columns = (int)malloc(NNZERO*sizeof(int));
    long long float *entries = (long long float)malloc(NNZERO*sizeof(long long float));

    printf("LINHAS\n");
    for(i = 0;i < NROW + 1;i++){
        fscanf(f,"%d",&rows[i]);
        printf("%d ",rows[i]);
    }
    printf("\n");

    printf("COLUNAS\n");
    for(i = 0;i < NNZERO;i++){
        fscanf(f,"%d",&columns[i]);
        printf("%d ",columns[i]);
    }
    printf("\n");

    printf("ENTRADAS\n");
    for(i = 0;i < NNZERO + 1;i++){
        fscanf(f,"%llf",&entries[i]);
        printf("%llf ",entries[i]);
    }
    printf("\n");
}
