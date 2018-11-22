#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <err.h>

#include <fcntl.h> //ham 
#include <unistd.h> //ham close

#include "string.h"
#include <math.h>

# ifndef PI
# define PI	3.14159265358979323846264338327950288
# endif

typedef float real;
typedef struct complex_t {
		real re; 
		real im;
} complex;
//Doc header file Wav

typedef struct {
    char     chunkID[4];
    uint32_t chunkSize;
    char     format[4];
    char     subchunk1ID[4];
    uint32_t subchunk1Size;
    uint16_t audioFormat;
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bps;
    char     datachunkID[4];
    uint32_t datachunkSize;
}WavHeader;
WavHeader *header;

////////////////////////////////////////////////////////////////////////
complex conv_from_polar(double r, double radians);
complex add(complex left, complex right);
complex multiply(complex left, complex right);
complex* DFT_naive(complex* x, int N);
complex* FFT_CooleyTukey(complex* input, int N, int N1, int N2);
void fft(complex *v, int n, complex *tmp);
complex *wavread(char *fileName, char *fileToSave, int N);
////////////////////////////////////////////////////////////////////////

complex conv_from_polar(double r, double radians) {
    complex result;
    result.re = r * cos(radians);
    result.im = r * sin(radians);
    return result;
}

complex add(complex left, complex right) {
    complex result;
    result.re = left.re + right.re;
    result.im = left.im + right.im;
    return result;
}

complex multiply(complex left, complex right) {
    complex result;
    result.re = left.re*right.re - left.im*right.im;
    result.im = left.re*right.im + left.im*right.re;
    return result;
}

complex* DFT_naive(complex* x, int N) {
    complex* X = (complex*) malloc(sizeof(struct complex_t) * N);
    int k, n;
    for(k = 0; k < N; k++) {
        X[k].re = 0.0;
        X[k].im = 0.0;
        for(n = 0; n < N; n++) {
            X[k] = add(X[k], multiply(x[n], conv_from_polar(1, -2*PI*n*k/N)));
        }
    }
    return X;
}

complex* FFT_CooleyTukey(complex* input, int N, int N1, int N2) {
    int k1, k2;

    /* Cap phat cot ma tran */
    complex** columns = (complex**) malloc(sizeof(struct complex_t*) * N1);
    for(k1 = 0; k1 < N1; k1++) {
        columns[k1] = (complex*) malloc(sizeof(struct complex_t) * N2);
    }
    
    /* Cap phat hang ma tran */
    complex** rows = (complex**) malloc(sizeof(struct complex_t*) * N2);
    for(k2 = 0; k2 < N2; k2++) {
        rows[k2] = (complex*) malloc(sizeof(struct complex_t) * N1);
    }
    
    /* Dinh dang lai dau vao N cot */
    for (k1 = 0; k1 < N1; k1++) {
        for(k2 = 0; k2 < N2; k2++) {
            columns[k1][k2] = input[N1*k2 + k1];
        }
    }

    /* Tinh toan N1 DFTs cua do dai N2 su dung phuong phap thong thuong */
    for (k1 = 0; k1 < N1; k1++) {
        columns[k1] = DFT_naive(columns[k1], N2);
    }
    
    /* Nhan 2 nhan to nghich ( e^(-2*pi*j/N * k1*k2)) va chuyen vi */
    for(k1 = 0; k1 < N1; k1++) {
        for (k2 = 0; k2 < N2; k2++) {
            rows[k2][k1] = multiply(conv_from_polar(1, -2.0*PI*k1*k2/N), columns[k1][k2]);
        }
    }
    
    /* Tinh toan N2 DFTs cua do dai N1 su dung phuong phap thong thuong */
    for (k2 = 0; k2 < N2; k2++) {
        rows[k2] = DFT_naive(rows[k2], N1);
    }
    
    /* Flatten into single output */
    complex* output = (complex*) malloc(sizeof(struct complex_t) * N);
    for(k1 = 0; k1 < N1; k1++) {
        for (k2 = 0; k2 < N2; k2++) {
            output[N2*k1 + k2] = rows[k2][k1];
        }
    }

    /* Xoa cap phat bo nho */
    for(k1 = 0; k1 < N1; k1++) {
        free(columns[k1]);
    }
    for(k2 = 0; k2 < N2; k2++) {
        free(rows[k2]);
    }
    free(columns);
    free(rows);
    return output;
}

complex *wavread(char *fileName, char *fileToSave, int N)
{
    int fd;
    
    FILE *fin = fopen(fileName, "rb");
    
    if (!fileName)
        errx(1, "Filename not specified");
    if ((fd = open(fileName, O_RDONLY)) < 1)
        errx(1, "Error opening file");
    if (!header)
        header = (WavHeader*)malloc(sizeof(WavHeader));
    if (read(fd, header, sizeof(WavHeader)) < sizeof(WavHeader))
        errx(1, "File broken: header");
    if (strncmp(header->chunkID, "RIFF", 4) ||
        strncmp(header->format, "WAVE", 4))
        errx(1, "Not a wav file");
    if (header->audioFormat != 1)
        errx(1, "Only PCM encoding supported");
    
    int chunkRealSize = header->chunkSize - 36;
    // printf("%d \n",chunkRealSize);
    close(fd);
    
    //Number of samples
    int sample_size = (header->bps)/ 8;
    int samples_count = chunkRealSize * 8 / (header->bps);
    //printf("Samples count = %i\n", samples_count);
    
    short int *value = (short int*)malloc((samples_count-22)* sizeof(short int));
    short int *value_test = (short int*)malloc(samples_count* sizeof(short int));
        //Reading data
    int i = 0;
    for (i = 0; i < samples_count; i++) {
        value_test[i] = 0;
        fread(&value_test[i], sample_size, 1, fin);
    }
    
    samples_count = samples_count - 22;
    
    for (i = 0; i < samples_count; i++) {
        value[i] = value_test[i + 22];
    }
    
    /*//Write data into the file
    FILE *fout = fopen(fileToSave, "w");
    i = 0;
    for (i = 0; i < samples_count; i++) {
        fprintf(fout, "%d\n", value[i]);
    }*/
    
    fclose(fin);
   // fclose(fout);
   
    i = 0;
    complex  *result1;
    
    complex *v = (complex*)malloc(N* sizeof(complex));
    for (i = 0; i < N; i++) {
        v[i].re = value[i];
        v[i].im = 0;
    }
    
    result1 = FFT_CooleyTukey(v, N, N, 1);
    return result1;
    free(v);
}

float DOAS_Beam(complex *channel1, complex *channel2, int num_mic, int sample_fft){
    complex **X = (complex**)malloc(sample_fft* sizeof(complex*));
    int i = 0, j = 0, k = 0;
    for (i = 0; i < sample_fft; i++) {
		X[i]=(complex*)malloc(sizeof(complex)*num_mic);
    }//X[sample_fft][num_mic];
    i = 0;
    for (i = 0; i < sample_fft; i++) {
    	X[i][0].re = channel1[i].re;
        X[i][0].im = channel1[i].im;
        X[i][1].re = channel2[i].re;
        X[i][1].im = channel2[i].im;
    }
    
    i = 0;
    complex **Y = (complex**)malloc(num_mic* sizeof(complex*));
    for (i = 0; i < num_mic; i++) {
		Y[i]=(complex*)malloc(sizeof(complex)*sample_fft);
    }//complex Y[num_mic][sample_fft];
    for (i = 0; i < sample_fft; i++)
    for (j = 0; j < num_mic; j++) {
    	Y[j][i].re = X[i][j].re;
        Y[j][i].im = -X[i][j].im;
    }
    //Tinh ma tran hiep phuong sai
    complex **R = (complex**)malloc(num_mic* sizeof(complex*));
    for (i = 0; i < num_mic; i++) {
		R[i]=(complex*)malloc(sizeof(complex)*num_mic);
    } //complex R[num_mic][num_mic];
    i = 0; j = 0;
    for (i = 0; i < num_mic; i++)
    for (j = 0; j < num_mic ; j++) {
            R[i][j].re = 0;
            R[i][j].im = 0;
	    k = 0;
            for (k = 0; k < sample_fft; k++) {
                R[i][j].re += (Y[i][k].re*X[k][j].re - Y[i][k].im*X[k][j].im)/sample_fft;
                R[i][j].im += (Y[i][k].re*X[k][j].im + Y[i][k].im*X[k][j].re)/sample_fft;
                }
    }
    
    //Tinh pho Beamforming
    int daipho = 180/0.01;
    float *angles = (float*)malloc(daipho* sizeof(float*)); //float angles[daipho];
    float *angles_dgree = (float*)malloc(daipho* sizeof(float*)); //float angles_dgree[daipho];
    i = 0;
    for (i = 0; i <= daipho; i++) {
    	angles_dgree[i] = -90 + i*0.01;
    	angles[i] = angles_dgree[i]*PI/180;
    }
    
    double **theta_al = (double**)malloc(daipho* sizeof(double*));
    for (i = 0; i < daipho; i++) {
		theta_al[i]=(double*)malloc(sizeof(double)*num_mic);
    }//double theta_al[daipho][num_mic];
    float d = 0.04;
    i = 0;
    j = 0;
    for (i = 0; i < daipho; i++)
    for (j = 0; j < num_mic; j++) {
    	theta_al[i][j] = 2*PI*d*j*sin(angles[i]);
    }

    complex **al = (complex**)malloc(daipho* sizeof(complex*));
    for (i = 0; i < daipho; i++) {
		al[i]=(complex*)malloc(sizeof(complex)*num_mic);
    }//complex al[daipho][num_mic];
    i = 0;
    j = 0;
    for (i = 0; i< daipho; i++)
	for (j = 0; j < num_mic; j++) {	
		al[i][j].re = cos(theta_al[i][j]);
		al[i][j].im = -sin(theta_al[i][j]);
	}
	
	//Ma tran dao cua vector trong luong
    complex **al_re = (complex**)malloc(num_mic* sizeof(complex*));
    for (i = 0; i < num_mic; i++) {
		al_re[i]=(complex*)malloc(sizeof(complex)*daipho);
    }//complex al_re[num_mic][daipho];
    i = 0;
    j = 0;
    for (i = 0; i < daipho; i++)
    for (j = 0; j < num_mic; j++) {
    	al_re[j][i].re = al[i][j].re;
        al_re[j][i].im = -al[i][j].im;
    }
	
    //Pho Beamforming
    complex **Classical_1 = (complex**)malloc(num_mic* sizeof(complex*));
    for (i = 0; i < num_mic; i++) {
		Classical_1[i]=(complex*)malloc(sizeof(complex)*daipho);
    }//complex Classical_1[num_mic][daipho];
    i = 0;
    j = 0;
    for (i = 0; i < num_mic; i++)
    for (j = 0; j < daipho ; j++) {
            Classical_1[i][j].re = 0;
            Classical_1[i][j].im = 0;
	    k = 0;
            for (k = 0; k < num_mic; k++) {
                Classical_1[i][j].re += R[i][k].re*al_re[k][j].re - R[i][k].im*al_re[k][j].im;
                Classical_1[i][j].im += R[i][k].re*al_re[k][j].im + R[i][k].im*al_re[k][j].re;
            }
    }
    
    complex *Classical = (complex*)malloc(daipho* sizeof(complex*)); //complex Classical[daipho];
    i = 0;
    for (i = 0; i< daipho; i++) {	
      
	    Classical[i].re = 0;
            Classical[i].im = 0;
	    k = 0;
            for (k = 0; k < num_mic; k++){
                Classical[i].re += al[i][k].re*Classical_1[k][i].re - al[i][k].im*Classical_1[k][i].im;
                Classical[i].im += al[i][k].re*Classical_1[k][i].im + al[i][k].im*Classical_1[k][i].re;
            }
    }
    
    float max_beamforming = 0;
    int i_max = 0;
    i = 0;
    for (i = 0; i < daipho; i++) {
        if (Classical[i].re >= max_beamforming) { 
		max_beamforming = Classical[i].re;
		i_max = i;
	}
    }
    printf("Bien do Beamforming lon nhat %f \nVi tri lon nhat %d \nGoc den tin hieu %f \n", max_beamforming, i_max, angles_dgree[i_max]);
    return angles_dgree[i_max];
}

int main(int argc, char **argv){
	
    FILE *fd = fopen("datatest.txt", "w");
     
    int K = 800, N = 2;
    complex *result1, *result2;
    result1 = wavread("a.wav", "open1.txt", K*2);
    result2 = wavread("b.wav", "open2.txt", K*2);
	
    float angle;
    angle = DOAS_Beam(result1, result2 , N, K);

    fprintf(fd, "%d      %f\n", K*2, angle + 20); 
   
    fclose(fd);
    return 0;
}
