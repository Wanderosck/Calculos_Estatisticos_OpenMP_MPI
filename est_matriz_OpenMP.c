#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>

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
    #pragma omp parallel shared(col,lin,matriz) private(i) 
	{
        #pragma omp for schedule(static)
            for (i = 0; i < col; i++) {    
                //manda o endereco do primeiro elemento da coluna, limites inf e sup e a largura da matriz
                bubble_sort(&matriz[i*lin], lin);
            }
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
    #pragma omp parallel shared(col,lin,matriz,moda) private(i)
	{
        #pragma omp for schedule(static)
        for(i=0;i<col;i++){
            moda[i]=moda_aux(matriz+(i*lin),lin);
        }
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
    clock_t inicio, fim, t1,t2;
    int lin,col,i,j; // Define as variáveis de índices e dimensões
    double *matriz,*mediana,*media,*media_har,*moda,*variancia,*dp,*cv; // Define a matriz (forma linear), vetores de medidas estatísticas
    FILE *f, *saida;

    //Leitura do arquivo
    f = fopen("./input","r");

    if(!f) {printf("Erro");}

    fscanf(f, "%d ", &lin); // Lê a quantidade de linhas da matriz
    fscanf(f, "%d\n", &col); // Lê a quantidade de colunas da matriz

    // Alocações
    matriz=(double *)malloc(lin*col * sizeof(double)); // Aloca a matriz
    media=(double *)malloc(col * sizeof(double)); // Aloca o vetor de media
    media_har=(double *)malloc(col * sizeof(double)); // Aloca o vetor de media harmônica
    mediana=(double *)malloc(col * sizeof(double)); // Aloca o vetor de mediana
    moda=(double *)malloc(col * sizeof(double)); // Aloca o vetor de moda
    variancia=(double *)malloc(col * sizeof(double)); // Aloca o vetor de variância
    cv=(double *)malloc(col * sizeof(double)); // Aloca o vetor de coeficiente de variação
    dp=(double *)malloc(col * sizeof(double)); // Aloca o vetor de desvio padrão

    for(i=0;i<lin;i++){
        for(j=0;j<col;j++){
            fscanf(f, "%lf ",&(matriz[j*lin+i])); // Lê os dados transpostos em uma matriz de entrada
        }
    };

    fclose(f);
    inicio = clock();
    calcula_media(matriz,media,lin,col);
    calcula_media_harmonica(matriz,media_har,lin,col);

    t1 = clock();
    ordena_colunas(matriz,lin,col);
    t2 = clock();
    printf("Ordena colunas: %f s\n", (float)(t2-t1)/CLOCKS_PER_SEC);

    calcula_mediana(matriz,mediana,lin,col);

    t1 = clock();
    calcula_moda(matriz,moda,lin,col);
    t2 = clock();
    printf("Calcula Moda: %f s\n", (float)(t2-t1)/CLOCKS_PER_SEC);

    calcula_variancia(matriz,media,variancia,lin,col);
    calcula_desvio_padrao(variancia,dp,col);
    calcula_coeficiente_variacao(media,dp,cv,col);
    fim = clock();
    printf("Tempo final: %f s\n", (float)(fim-inicio)/CLOCKS_PER_SEC);

    char nomeArquivo[11] = "saidaS.txt";
    #ifdef _OPENMP
    nomeArquivo[5] = 'P';
    #endif

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