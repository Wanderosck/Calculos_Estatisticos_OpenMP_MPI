#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>
#include <time.h>

void bubble_sort (double *matriz, int tam) {
    int i, j; 
    double aux;
    for (i = tam - 1; i > 0; i--) {
        for (j = 0; j < i; j++) {
            if (matriz[j] > matriz[j+1]) {
                aux=matriz[j];
                matriz[j]=matriz[j+1];
                matriz[j+1]=aux;
            }
        }
    }
}

void ordena_colunas(double *matriz,int lin, int col) {
    int i;
    for (i = 0; i < col; i++) { 
        //manda o endereco do primeiro elemento da coluna, limites inf e sup e a largura da matriz
        bubble_sort(&matriz[i*lin], lin);
    }
} 

void calcula_media(double *matriz, double *vet,int lin, int col){
    int i,j;
    double soma;
    for(i=0;i<col;i++){
        soma=0;
        for(j=0;j<lin;j++){
            soma+=matriz[i*lin+j];
        }
        vet[i]=soma/lin; 
    }   
}

void calcula_media_harmonica(double *matriz, double *vet,int lin, int col){
    int i,j;
    double soma;
    for(i=0;i<col;i++){
        soma=0;
        for(j=0;j<lin;j++){
            soma+=(1/(matriz[i*lin+j]));
        }
        vet[i]=lin/soma; 
    }   
}

void calcula_mediana(double *matriz, double *vet, int lin, int col) {  
     int i=0;
     if((lin-1)%2){ //Quantidade par de elementos
         for (i = 0; i < col; i++) {
            vet[i]=(matriz[(((i*lin)+(i+1)*lin)-1)/2]+matriz[((((i*lin)+(i+1)*lin)-1)/2)+1]);
            vet[i]*=0.5;
            
         }
     }
     else{ //Quantidade ímpar de elementos
         for (i = 0; i < col; i++) {
            vet[i]=matriz[(((i*lin)+(i+1)*lin)-1)/2];
         }
     }
    
 } 

//Adaptado de https://www.clubedohardware.com.br/forums/topic/1291570-programa-em-c-que-calcula-moda-media-e-mediana/
double moda_aux(double *matriz,int lin){
    int i, j; 
    double *cont;
    cont=(double*)malloc(lin*sizeof(double));
	float conta=0, moda=0;
	
	for(i=0;i<lin;i++){
        for(j=i+1;j<lin;j++){
        	
			if(matriz[i]==matriz[j]){
				cont[i]++;
					if(cont[i]>conta){
						conta=cont[i];
						moda=matriz[i];
					}
			}

        }
        cont[i]=0;
    }
    free(cont);
    if(conta == 0){
    	return -1;
	}
	else{
		return moda;
	}
}


void calcula_moda(double *matriz,double *moda,int lin, int col){
    int i;
    for(i=0;i<col;i++){
        moda[i]=moda_aux(matriz+(i*lin),lin);
    }
}

void calcula_variancia(double *matriz, double *media,double *variancia, int lin, int col){
    int i,j;
    double soma;
    for(i=0;i<col;i++){
        soma=0;
        for(j=0;j<lin;j++){
            soma+=pow((matriz[i*lin+j]-media[i]),2);
        }
        variancia[i]=soma/(lin-1); 
    }
}

void calcula_desvio_padrao(double *variancia,double *dp, int col){
    int i;
    for(i=0;i<col;i++){
        dp[i]=sqrt(variancia[i]);
    }  
}

void calcula_coeficiente_variacao(double *media,double *dp,double *cv, int col){
    int i;
    for(i=0;i<col;i++){
        cv[i]=dp[i]/media[i];
    }  
}

int main(int argc,char **argv){
    clock_t t1, t2, inicio, fim, a, b;
    double *matriz_local, *moda_local;
    double *matriz, *mediana,*media,*media_har,*moda,*variancia,*dp,*cv; // Define a matriz (forma linear), vetores de medidas estatísticas
    int rank, processos, i, j, k, colPorProcesso, lin, col;
    FILE *f, *saida;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processos);
    MPI_Status status;
    if (rank == 0) {
        //Leitura do arquivo
        f = fopen("./input","r");
        if(!f) {printf("Erro");}

        fscanf(f, "%d ", &lin); // Lê a quantidade de linhas da matriz
        fscanf(f, "%d\n", &col); // Lê a quantidade de colunas da matriz

        colPorProcesso = col/processos;

        // Alocações
        matriz=(double *)malloc(lin*col * sizeof(double)); // Aloca a matriz
        media=(double *)malloc(col * sizeof(double)); // Aloca o vetor de media
        media_har=(double *)malloc(col * sizeof(double)); // Aloca o vetor de media harmônica
        mediana=(double *)malloc(col * sizeof(double)); // Aloca o vetor de mediana
        moda=(double *)malloc(col * sizeof(double)); // Aloca o vetor de moda
        variancia=(double *)malloc(col * sizeof(double)); // Aloca o vetor de variância
        cv=(double *)malloc(col * sizeof(double)); // Aloca o vetor de coeficiente de variação
        dp=(double *)malloc(col * sizeof(double)); // Aloca o vetor de desvio padrão
        matriz_local=(double *)malloc(lin*colPorProcesso * sizeof(double)); // Aloca a matriz local

        for(i=0;i<lin;i++){
            for(j=0;j<col;j++){
                fscanf(f, "%lf ",&(matriz[j*lin+i])); // Lê os dados transpostos em uma matriz de entrada
            }
        };

        fclose(f);
        inicio = clock();
        calcula_media(matriz,media,lin,col);
        calcula_media_harmonica(matriz,media_har,lin,col);
        calcula_variancia(matriz,media,variancia,lin,col);
        calcula_desvio_padrao(variancia,dp,col);
        calcula_coeficiente_variacao(media,dp,cv,col);
        for(i=1; i<processos; i++){
            MPI_Send(&colPorProcesso, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&lin, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&(matriz[i*colPorProcesso*lin]), colPorProcesso*lin, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        for(i=0; i<lin*colPorProcesso; i++){
            matriz_local[i] = matriz[i]; 
        }
    }
    else {
        a = clock();
        MPI_Recv(&colPorProcesso, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&lin, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        matriz_local=(double *)malloc(lin*colPorProcesso * sizeof(double)); // Aloca a matriz local 
        MPI_Recv(matriz_local, colPorProcesso*lin, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        b = clock();
        printf("Processo %d: %f s\n", rank, (float)(b-a)/CLOCKS_PER_SEC);
    }
    moda_local=(double *)malloc(colPorProcesso * sizeof(double)); // Aloca o vetor de moda local
    ordena_colunas(matriz_local,lin,colPorProcesso);
    calcula_moda(matriz_local,moda_local,lin,colPorProcesso);

    if (rank != 0) {
        MPI_Send(&(moda_local[0]), colPorProcesso, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&(matriz_local[0]), colPorProcesso*lin, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
        for(i=0;i<colPorProcesso;i++) {
            moda[i] = moda_local[i];
        }
        for(i=0;i<colPorProcesso*lin;i++) {
            matriz[i] = matriz_local[i];
        }
        for(i=1;i<processos;i++) {
            MPI_Recv(moda_local, colPorProcesso, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(matriz_local, colPorProcesso*lin, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
            k = 0;
            for(j=colPorProcesso*i;j<colPorProcesso*(i+1);j++) {
                moda[j] = moda_local[k];
                k++;
            }
            k = 0;
            for(j=colPorProcesso*lin*i;j<colPorProcesso*lin*(i+1);j++) {
                matriz[j] = matriz_local[k];
                k++;
            }
        }
    }
    if (rank == 0) {
        calcula_mediana(matriz,mediana,lin,col);
        fim = clock();
        printf("Tempo final do processo 0: %f s\n", (float)(fim-inicio)/CLOCKS_PER_SEC);
        char nomeArquivo[11] = "saidaP.txt";
        saida = fopen(nomeArquivo, "w");
        for(i=0;i<col;i++)
            fprintf(saida,"%.1lf ",media[i]); // Imprime as médias aritméticas de cada coluna
        fprintf(saida,"\n");
        for(i=0;i<col;i++)
            fprintf(saida,"%.1lf ",media_har[i]); // Imprime as médias harmônicas de cada coluna
        fprintf(saida,"\n");
        for(i=0;i<col;i++)
            fprintf(saida,"%.1lf ",mediana[i]); // Imprime as medianas de cada coluna
        fprintf(saida,"\n");
        for(i=0;i<col;i++)
            fprintf(saida,"%.1lf ",moda[i]); // Imprime as modas de cada coluna
        fprintf(saida,"\n");
        for(i=0;i<col;i++)
            fprintf(saida,"%.1lf ",variancia[i]); // Imprime as variâncias amostrais de cada coluna
        fprintf(saida,"\n");
        for(i=0;i<col;i++)
            fprintf(saida,"%.1lf ",dp[i]); // Imprime os desvios padrão de cada coluna
        fprintf(saida,"\n");
        for(i=0;i<col;i++)
            fprintf(saida,"%.1lf ",cv[i]); // Imprime os coeficientes de variação de cada coluna
        fclose(saida);

        // Desaloca memória
        free(matriz);
        free(media);
        free(media_har);
        free(mediana);
        free(moda);
        free(variancia);
        free(dp);
        free(cv);
    }
    MPI_Finalize();
}