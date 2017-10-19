#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include "svpng/svpng.inc"
#include <complex.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include <mpi.h>
#include <omp.h>

#define FILE_NAME "record.txt"

//escape tolerance when doing the convergence calculation
const double tolerance = 1e-7;
//list of polynomial terms with their constants in front, change it to {0,4,3,0,0,1} (0+4z+3z^2+0*z^3+0*z^4+1*z^5) per PIAZZA
const double polynomial[]={0,4,3,0,0,1};
// Number of terms in the polynomial
const unsigned int N = sizeof(polynomial)/sizeof(*polynomial);

// Get the result (also the numerator when doing the iteration) based on the polynomial given
double complex get_numerator(double complex z) {
    double complex numerator = polynomial[0]+polynomial[1]*z+polynomial[2]*z*z+polynomial[3]*z*z*z;
    for (unsigned int i = 4; i < N; ++i) {
        if (polynomial[i]==0)
            continue;
        numerator+=polynomial[i]*cpow(z,i);
    }
    return numerator;
}

//Assume N>=3, apply z-=f(z)/f'(z)
void poly_calc(double complex *z) {
    double complex numerator = polynomial[0]+polynomial[1]*(*z)+polynomial[2]*(*z)*(*z)+polynomial[3]*(*z)*(*z)*(*z);
    double complex denominator = polynomial[1]+2*polynomial[2]*(*z)+3*polynomial[3]*(*z)*(*z);
    for (unsigned int i = 4; i < N; ++i) {
        if (polynomial[i]==0) continue;
        numerator+=polynomial[i]*cpow(*z,i);
        denominator+=i*polynomial[i]*cpow(*z,i-1);
    }
    *z -= numerator/denominator;
}

void get_root(unsigned int precision, unsigned int xPrecision, unsigned int yPrecision, mpc_t root) {

    mpfr_t root8799, ui81, div13, div23, p15,m82_3;
    mpc_t c81_r8799,first_root1, first_root2,base,base_root;
    mpfr_init2(root8799,precision);
    mpfr_init2(div13,precision);
    mpfr_init2(m82_3,precision);
    mpfr_init2(div23,precision);
    mpfr_init2(p15,precision);
    mpfr_init2(ui81,precision);
    mpc_init3(c81_r8799,xPrecision,yPrecision);
    mpc_init3(base,xPrecision,yPrecision);
    mpc_init3(base_root,xPrecision,yPrecision);
    mpc_init3(first_root1,xPrecision,yPrecision);
    mpc_init3(first_root2,xPrecision,yPrecision);

    mpfr_set_ui(div13,3,MPFR_RNDZ);
    mpfr_ui_div(div13,1,div13,MPFR_RNDZ);
    mpfr_ui_div(div23,2,div13,MPFR_RNDZ);
    mpfr_set_ui(ui81,81,MPFR_RNDZ);
    mpfr_sqrt_ui(root8799,8799,MPFR_RNDZ);
    mpc_set_fr_fr(c81_r8799,ui81,root8799,MPC_RNDZZ);
    mpc_mul_ui(first_root1,c81_r8799,2,MPC_RNDZZ);
    mpc_pow_fr(first_root1,first_root1,div13,MPC_RNDZZ);
    mpfr_set_ui(p15,15*15,MPFR_RNDZ);
    mpfr_cbrt(p15,p15,MPFR_RNDZ);
    mpc_div_fr(first_root1,first_root1,p15,MPC_RNDZZ);

    mpfr_set_ui(m82_3,2048,MPFR_RNDZ);
    mpfr_cbrt(m82_3,m82_3,MPFR_RNDZ);
    mpc_mul_ui(first_root2,c81_r8799,15,MPC_RNDZZ);
    mpc_pow_fr(first_root2,first_root2,div13,MPC_RNDZZ);
    mpc_fr_div(first_root2,m82_3,first_root2,MPC_RNDZZ);

    mpc_add(base,first_root1,first_root2,MPC_RNDZZ);
    mpc_sqrt(base_root,base,MPC_RNDZZ);
    mpc_mul_ui(root,base_root,5,MPC_RNDZZ);
    mpc_ui_div(root,12,root,MPC_RNDZZ);
    mpc_sub(root,root,base,MPC_RNDZZ);
    mpc_sqrt(root,root,MPC_RNDZZ);
    mpc_sub(root,root,base_root,MPC_RNDZZ);
    mpc_div_ui(root,root,2,MPC_RNDZZ);

    mpfr_clear(root8799);
    mpfr_clear(ui81);
    mpfr_clear(div13);
    mpfr_clear(div23);
    mpfr_clear(p15);
    mpfr_clear(m82_3);
    mpc_clear(c81_r8799);
    mpc_clear(first_root1);
    mpc_clear(first_root2);
    mpc_clear(base);
    mpc_clear(base_root);
}

//Assume N>=3, apply z-=f(z)/f'(z), using MPC (based on mpfr which is based on gmp)
void poly_calc_mpc(mpfr_t mul, mpc_t zc, mpc_t term, mpc_t numerator, mpc_t denominator) {
    mpc_set_ui(numerator,0,MPC_RNDZZ);
    mpc_set_ui(denominator,0,MPC_RNDZZ);

    // the follow for loop gets the numerator and denominator of the polynomial
    for (unsigned int j = 0; j < N; ++j) {
        if (polynomial[j]==0) continue;
        mpfr_set_d(mul,polynomial[j],MPFR_RNDZ);
        mpc_pow_ui(term,zc,j,MPC_RNDZZ);
        mpc_mul_fr(term,term,mul,MPC_RNDZZ);
        mpc_add(numerator,numerator,term,MPC_RNDZZ);
        if (j==0) continue;
        mpc_pow_ui(term,zc,j-1,MPC_RNDZZ);
        mpc_mul_fr(term,term,mul,MPC_RNDZZ);
        mpc_mul_ui(term,term,j,MPC_RNDZZ);
        mpc_add(denominator,denominator,term,MPC_RNDZZ);
    }
    mpc_div(term,numerator,denominator,MPC_RNDZZ);
    mpc_sub(zc,zc,term,MPC_RNDZZ);
    mpc_set_ui(numerator,0,MPC_RNDZZ);

    // the follow for loop gets the numerator of the polynomial
    for (unsigned int j = 0; j < N; ++j) {
        if (polynomial[j]==0) continue;
        mpfr_set_d(mul,polynomial[j],MPFR_RNDZ);
        mpc_pow_ui(term,zc,j,MPC_RNDZZ);
        mpc_mul_fr(term,term,mul,MPC_RNDZZ);
        mpc_add(numerator,numerator,term,MPC_RNDZZ);
    }

    mpc_abs(mul,numerator,MPFR_RNDZ);
}


// Apply the iterations, 1024 iterations are used for double) and 
// linearly increase when MPC is used to prevent from black screen
long unsigned int newton(double complex *z, mpc_t zc) {
    if (z!=NULL) {
        for (unsigned int i = 0; i < 102400; ++i) {
            poly_calc(z);
            if (cabs(get_numerator(*z)) < tolerance)
                return i;
        }
    } else {
        long int pr,pi;
        mpfr_t mul;
        mpc_get_prec2(&pr,&pi,zc);
        mpfr_init2(mul,pr);
        mpc_t numerator,denominator,term;
        mpc_init3(numerator,pr,pi);
        mpc_init3(denominator,pr,pi);
        mpc_init3(term,pr,pi);
        for (long unsigned int i = 0; i < 128*pr/N; ++i) {
            poly_calc_mpc(mul,zc,term,numerator,denominator);
            if (mpfr_cmp_d(mul,tolerance)<0) {
                mpfr_clear(mul);
                mpc_clear(numerator);
                mpc_clear(denominator);
                mpc_clear(term);
                return i;
            }
        }
        mpfr_clear(mul);
        mpc_clear(numerator);
        mpc_clear(denominator);
        mpc_clear(term);
    }
    return 0;
}


//Read file from the spec file called spec, does a bit of error handling
bool read_spec(FILE *spec, unsigned int *width, unsigned int *height, double *leftx, double *lefty, double *rightx,
               double *righty, double *zoom, unsigned int *frames) {
    char *line;
    size_t len = 0;
    ssize_t read;

    for (unsigned int i = 0; i < 5; ++i) {
        if ((read = getline(&line, &len, spec)) == -1)
            return false;

        switch (i) {
            case 0:
                sscanf(line, "%u %u", width, height);
                break;
            case 1:
                sscanf(line, "%lf %lf", leftx, lefty);
                break;
            case 2:
                sscanf(line, "%lf %lf", rightx, righty);
                break;
            case 3:
                sscanf(line, "%lf", zoom);
                break;
            case 4:
                sscanf(line, "%u", frames);
                break;
        }
    }
    return true;
}


// Generate random colours, for fun
uint8_t** random_colours(void) {
    //to seed rand() which is used to generate a set of random colours
    time_t t;
    srand((unsigned) time(&t));
    uint8_t** colours = malloc(N*sizeof(*colours));
    for (unsigned int i = 0; i < N; ++i) {
        colours[i] = malloc(3*sizeof(**colours));
    }

    //Generate a set of random colours based on the number of roots
    for (unsigned int i = 0; i < N; ++i) {
        switch (i) {
            case 0:
                colours[i][0] = 244;
                colours[i][1] = 173;
                colours[i][2] = 66;
                break;
            case 1:
                colours[i][0] = 75;
                colours[i][1] = 0;
                colours[i][2] = 130;
                break;
            case 2:
                colours[i][0] = 255;
                colours[i][1] = 255;
                colours[i][2] = 0;
                break;
            case 3:
                colours[i][0] = 0;
                colours[i][1] = 191;
                colours[i][2] = 255;
                break;
            default:
                colours[i][0] = rand() % 256;
                colours[i][1] = rand() % 256;
                colours[i][2] = rand() % 256;
        }
    }
    return colours;
}

int main(int argc, char const *argv[]) {

    // MPI initialisation
    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (argc!=2)
        puts("Usage: ./fract SPEC_X");

    FILE* spec = fopen(argv[1], "r");
    if (spec == NULL) {
        puts("save the spec to a file, also check permission");
        return 0;
    }

    unsigned int width, height, frames;
    double leftx, lefty, rightx, righty, zoom;

    if (!read_spec(spec, &width, &height, &leftx, &lefty, &rightx, &righty, &zoom, &frames)) {
        puts("Error reading spec");
        return 0;
    }

    fclose(spec);

    FILE* file_ptr = fopen(FILE_NAME, "a");
    fprintf(file_ptr,"start,%lu\n", time(NULL));
    fclose(file_ptr);
    uint8_t** colours = random_colours();
    const double xDistance = rightx - leftx;
    const double yDistance = lefty - righty;

    mpfr_t zoom_64, current_zoom;
    mpfr_init2(zoom_64,64);
    mpfr_set_d (zoom_64, zoom, MPFR_RNDZ);
    mpfr_init2(current_zoom,64);


    for (unsigned int frame = 0; frame < frames; ++frame) {
        // Split each frame to each corresponding node
        if(frame%world_size==world_rank) {
            uint8_t *data = malloc(sizeof(*data) * height * width * 3);
            //long unsigned int *iteration_counts = malloc(sizeof(*iteration_counts) * height * width);
            char name[20];
            //double complex current_zoom_centre=0+0*I;
            // Increase the precision with the the current area being rendered
            const long unsigned int precision = (10+log10(1/zoom)*8/xDistance*frame*(N-1));

            const double current_zoom_double = pow(zoom,frame);
            const long unsigned int xPrecision = (precision+log(width)/log(2)+1);
            const long unsigned int yPrecision = (precision+log(height)/log(2)+1);

            // zoom in to where 4 + 6 z + 5 z^4 == 0
            // -0.703844398150689 + 0.2629949912237024 I
            
            mpc_t root;
            mpfr_t xMid_mpfr, yMid_mpfr;
            mpc_init3(root,xPrecision,yPrecision);
            mpfr_init2(xMid_mpfr,xPrecision);
            mpfr_init2(yMid_mpfr,yPrecision);
            get_root(precision,xPrecision,yPrecision,root);
            mpc_real(xMid_mpfr,root,MPFR_RNDZ);
            mpc_imag(yMid_mpfr,root,MPFR_RNDZ);
            mpc_clear(root);


            // Use standard double precision when precision needed is <= 64
            if (precision<=64) {
                double xMid = mpfr_get_d(xMid_mpfr,MPFR_RNDZ);
                double yMid = mpfr_get_d(yMid_mpfr,MPFR_RNDZ);

                #pragma omp parallel
                {
                    #pragma omp for
                    for (unsigned int i = 0; i < height; ++i) {
                        unsigned int temp = 3*i*width;
                        double yTemp = yMid + yDistance*(0.5- (double)i/height)*current_zoom_double;
                        for (unsigned int j = 0; j < width; ++j) {
                            double xTemp = xMid + xDistance*((double)j/width - 0.5)*current_zoom_double;
                            double complex z = xTemp + yTemp * I;
                            long unsigned int n = newton(&z,NULL);
                            if (n) {
                            	data[temp++] = colours[n%N][0];
	                            data[temp++] = colours[n%N][1];
	                            data[temp++] = colours[n%N][2];
                            } else {
                            	data[temp++] = 0;
                            	data[temp++] = 0;
                            	data[temp++] = 0;
                            }
                        }
                    }
                }
            } else{
                // Use mpc with proper precision calculated above
                mpfr_set_prec(current_zoom,precision);
                mpfr_pow_ui(current_zoom,zoom_64,frame,MPFR_RNDZ);
                printf("Precision for x %lu, y %lu, zoom value ", xPrecision,yPrecision);
                mpfr_out_str(stdout,10,0,current_zoom,MPFR_RNDZ);
                puts("\nUsing MPC");
                fflush(stdout);
                mpfr_t half;

                mpfr_init2(half,precision);
                mpfr_set_ui(half,2,MPFR_RNDZ);
                mpfr_ui_div(half,1,half,MPFR_RNDZ);


                #pragma omp parallel
                {
                    #pragma omp for
                    for (unsigned int i = 0; i < height; ++i) {
                        unsigned int temp = 3*i*width;
                        mpfr_t xTemp, yTemp;
                        mpc_t z;
                        mpfr_init2(xTemp,xPrecision);
                        mpfr_init2(yTemp,yPrecision);
                        mpc_init3(z,xPrecision,yPrecision);

                        mpfr_set_ui(yTemp,height,MPFR_RNDZ);
                        mpfr_set_ui(xTemp,width,MPFR_RNDZ);
                        mpfr_ui_div(yTemp,i,yTemp,MPFR_RNDZ);
                        mpfr_sub(yTemp,half,yTemp,MPFR_RNDZ);
                        mpfr_mul_d(yTemp,yTemp,yDistance, MPFR_RNDZ);

                        mpfr_mul(yTemp,current_zoom,yTemp,MPFR_RNDZ);
                        mpfr_add(yTemp,yTemp,yMid_mpfr,MPFR_RNDZ);

                        for (unsigned int j = 0; j < width; ++j) {
                            /*mpfr_ui_div(xTemp,j,xTemp,MPFR_RNDZ);
                            mpfr_sub(xTemp,xTemp,half,MPFR_RNDZ);
                            mpfr_mul_d(xTemp,xTemp,xDistance,MPFR_RNDZ);
                            mpfr_mul(xTemp,xTemp,current_zoom,d)*/
                            mpfr_mul_d(xTemp,current_zoom, xDistance*((double)j/width - 0.5),MPFR_RNDZ);
                            mpfr_add(xTemp,xTemp,xMid_mpfr,MPFR_RNDZ);
                            mpc_set_fr_fr(z,xTemp,yTemp,MPC_RNDZZ);
                            long unsigned int n = newton(NULL,z);

                            if (n) {
                            	data[temp++] = colours[n%N][0];
	                            data[temp++] = colours[n%N][1];
	                            data[temp++] = colours[n%N][2];
                            } else {
                            	data[temp++] = 0;
                            	data[temp++] = 0;
                            	data[temp++] = 0;
                            }
                        }
                        mpfr_clear(xTemp);
                        mpfr_clear(yTemp);
                        mpc_clear(z);
                    }
                }
                mpfr_clear(half);
            }

            
            mpfr_clear(xMid_mpfr);
            mpfr_clear(yMid_mpfr);
            //Add padding to the file name
            sprintf(name, "%06d.png", frame);
            FILE *file = fopen(name, "wb");

            //svpng save the data as a whole to a png file, faster than CImg
            svpng(file, width, height, data, 0);
            fclose(file);
            FILE* file_ptr = fopen(FILE_NAME, "a");
            fprintf(file_ptr,"%s,%lu\n", name, time(NULL));
            fclose(file_ptr);
            free(data);
            //free(iteration_counts);
        }
    }

    // Free all the ram and cache
    mpfr_clear(zoom_64);
    mpfr_clear(current_zoom);
    mpfr_free_cache();
    for (int i = 0; i < N; ++i) {
        free(colours[i]);
    }
    free(colours);
    MPI_Finalize();
}
