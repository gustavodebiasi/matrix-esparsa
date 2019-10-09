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
    printf("%d \n", PTRCRD);
    printf("%d \n", INDCRD);
    printf("%d \n", VALCRD);
    int rows[NROW];
    int columns[NNZERO];
    double entries[NNZERO];

    printf("LINHAS\n");
    for(i = 0;i < NROW + 1;i++){
        fscanf(f,"%d",&rows[i]);
        rows[i]--;
        printf("%d ",rows[i]);
    }
    printf("\n");

    printf("COLUNAS\n");
    for(i = 0;i < NNZERO;i++){
        fscanf(f,"%d",&columns[i]);
        columns[i]--;
        printf("%d ",columns[i]);
    }
    printf("\n");

    printf("ENTRADAS\n");
    for(i = 0;i < NNZERO;i++){
        fscanf(f,"%lf",&entries[i]);
        printf("%lf ",entries[i]);
    }
    printf("\n");

    double matriz[NROW][NCOL];
    for (i = 0; i < NROW; i++) {
        for (j = 0; j < NCOL; j++) {
            matriz[i][j] = 0;
        }
    }

    int contEntries = 0;
    for (i = 0; i < NROW + 1; i++) {
        int limite = 0;
        if (i != NROW) {
            limite = rows[i+1];
        }
        int finalRow = i+1;
        for (j = rows[i]; j < limite; j++) {
            matriz[i][columns[j]] = entries[contEntries];
            matriz[columns[j]][i] = entries[contEntries]; //como a matriz é simétrica, coloca tanto no [i][j] quanto no [j][i]
            contEntries++;
        }
    }
    printf("MATRIZ \n");
    for (i = 0; i<NROW; i++) {
        for (j = 0; j < NCOL; j++) {
            if(matriz[i][j] != 0)
                //printf("%lf ", matriz[i][j]);
                printf("x ");
            else
                printf("  ");
        }
        printf("\n");
    }
}
