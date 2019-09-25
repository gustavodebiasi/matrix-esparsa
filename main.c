#include <stdio.h>
#include <stdlib.h>

int main(){
  FILE *f;

  //o tamanho dos 3 arrays
  //ler do arquivo e alimentar os arrays

  f = fopen("bcsstk01.rsa", "rt");

  char linha[100];
  fgets(linha,100,f);
  printf("%s\n",linha);

  int TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD;
  char MXTYPE[3];
  int NROW, NCOL, NNZERO, NELTVL;
  //while(f != EOF){
    fscanf(f, "%d %d %d %d %d", &TOTCRD, &PTRCRD, &INDCRD, &VALCRD, &RHSCRD);
    fscanf(f, "%s %d %d %d %d", MXTYPE, &NROW, &NCOL, &NNZERO, &NELTVL);

    printf("%d %d %d %d %d\n",TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD);
    printf("%s %d %d %d %d\n", MXTYPE, NROW, NCOL, NNZERO, NELTVL);
  //}

}
