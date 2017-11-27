//Sikender Shahid - 0981476
//programming assignment1

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

#define TOLERANCE 0.001 

double * poly_array = (void*) 0; 

double horners_algorithm(double *array, int N, double x_value){
	if(N ==0){
	 	printf("Error: horners_algorithm requires a polynomial\n");
		exit(EXIT_FAILURE); 
	}

	poly_array = (double*) malloc(sizeof(double)* (N));
	poly_array[N-1] = array[0];  
	for(int i = 1 ; i < N ; i++){
		poly_array[N-1-i] = (poly_array[N-i] * (double) x_value )+ array[i]; 
	}

	return poly_array[0]; 
}

double horners_algorithm_diff(double *array, int N, double x_value){
	if(poly_array == NULL){
		horners_algorithm(array, N, x_value);  
	}

	double sum = 0.00; 
	for(int i = N-1 ; i > 0; i--){
		sum = sum + ( poly_array[i] * pow(x_value, (double) i-1) ); 
	}
	return sum; 
}

int N_MAX(double a, double b){
	return (log10(b-a)-log10(TOLERANCE * .01))/(log10(2));
}

int bisection_method_ERROR_TEST(double a, double b){
	return (b - a)/2.0 <= TOLERANCE;
}

double bisection_method(double * function, double coef_num,double a_bound, double b_bound){
	int N_max = N_MAX(a_bound, b_bound);  
	int N = 0;
	double c = 0; 
	while (N < N_max){
		c = (a_bound + b_bound)/2.00;
		if(horners_algorithm(function, coef_num, c) == 0 || \
			bisection_method_ERROR_TEST(a_bound, b_bound)){
			printf("B. Bisection Method number of iteration: %d\n", N);  
			return c; 
		}
		N++;
		if(horners_algorithm(function, coef_num, a_bound)* \
			horners_algorithm(function, coef_num,c) > 0)
			a_bound = c; 
		else
			b_bound = c; 
	} 
	printf("B. Bisection Method: error exited N_MAX\n");
	return 0; 
}

int newtons_method_ERROR_TEST(double x, double x_init){
	return fabs(x - x_init) < TOLERANCE; 
}

double newtons_method(double *function, double coef_num, double a_bound, double b_bound){
	int N_max = N_MAX(a_bound, b_bound); 
	double x_init = (a_bound + b_bound)/2.0;
	double x_next, y, y_prime;

	int N = 0; 
	
	while (N < N_max){
		y = horners_algorithm(function, coef_num, x_init);
		y_prime = horners_algorithm_diff(function, coef_num, x_init);
		x_next = x_init - (y/y_prime); 
		if(newtons_method_ERROR_TEST(x_next, x_init)){
			printf("D. Newton's Method number of iteration: %d\n", N);
			return x_next; 
		}
		x_init = x_next; 
		N++; 
	} 
	printf("D. Newton's Method: error exited N_MAX\n");
	return 0; 
}

int main(int argc, char *argv[]){ 
	if(argc == 1){
		printf("enter more arguments\n"); 
		exit(EXIT_FAILURE); 
	}

	double a_bound = atof(argv[argc-2]);
	double b_bound = atof(argv[argc-1]);
	int coef_num = argc -3;
	double * coef_array = (double*) malloc(sizeof(double)* (coef_num));
	for(int i = 0 ; i< coef_num; i++){
		coef_array[i] = atof(argv[1+i]);
	}

	//test print 1
	printf("A. Polynomial: ");
	for(int i = 0 ; i < coef_num ; i++){
		printf("%fx^(%d) ", coef_array[i], coef_num-(i+1)); 
	}
	printf("[%f,%f] \n", a_bound, b_bound);

	// //test print 2
	// printf(" R = %f \n", horners_algorithm(coef_array, coef_num, 3)); 
	// printf(" differential q1(x) = %f \n", horners_algorithm_diff(coef_array, coef_num, 3));

	//test print 3
	printf("C. Bisection Method root: %f\n", bisection_method(coef_array, coef_num, a_bound, b_bound));
	printf("E. Newton's Method root:  %f\n", newtons_method(coef_array, coef_num, a_bound, b_bound)); 

	
	return 0; 
}