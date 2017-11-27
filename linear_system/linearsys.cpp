//Sikender Shahid
#include <iostream>
#include <fstream>
#include <string.h> 
#include <cmath>
using namespace std;

class MatrixLUDecomposition{
public:
	MatrixLUDecomposition(int aRowNum, int aColNum){
		rowNum = aRowNum; 
		colNum = aColNum;
		l_vector = (int*)malloc((rowNum)*sizeof(int));
		for(int i = 0; i < rowNum ; i++){
			l_vector[i] = i; 
		}
	};
	~MatrixLUDecomposition(){};

	// use: locating the maximum value in a col to set as pivot
	// input: matrix and col number to find max value in
	// output: no return, sets a maintained row vector e.g. P Matrix 
	void partial_pivot_locator(double** matrix, int col){
		double max = abs(matrix[col][col]); 
		int max_position =0; 
		int temp_position = 0; 
		int i; 
		for(i =col+1; i<rowNum; i++){ 
			if(max < abs(matrix[i][col])){
				max = abs(matrix[i][col]);
				max_position = i; 	 
			}
		}
		if(max_position){
			temp_position = l_vector[col];
			l_vector[col] = l_vector[max_position];
			l_vector[max_position] = temp_position;
		}
	};

	// use: performing gaussian elemination to obtain L and U factorization
	// input: main matrix, l factoring matrix, pivot control (for partial pivoting)
	// output: return is the U factorized matrix
	double** gaussian_elemination(double** A_matrix, double** L_matrix, bool pivot){
		//creating U matrix 
		double ** U_matrix = (double **)malloc((rowNum)*sizeof(double*)); 
		for(int i = 0; i < colNum; i ++)
			U_matrix[i] = (double*)malloc((colNum)*sizeof(double));
		
		//creating L matrix reinitizlization
		for(int i = 0 ; i < rowNum ; i++){
			for(int j = 0; j < colNum ; j++){
				if(i == j)
					L_matrix[i][j] = 1; 
				else
					L_matrix[i][j] = 0; 
			}
		}

		// gaussian elemination row = row - mP 
		double multiplier = 0;
		for(int col = 0; col < colNum-1 ; col++){
			if(pivot)
				partial_pivot_locator(A_matrix, col);
			for(int row = col+1; row < rowNum; row++){
				multiplier = A_matrix[l_vector[row]][col]/A_matrix[l_vector[col]][col];
				L_matrix[row][col] = multiplier; 
				for(int colU = 0; colU < colNum; colU++)
					A_matrix[l_vector[row]][colU] = A_matrix[l_vector[row]][colU] - multiplier*(A_matrix[l_vector[col]][colU]);
			} 
		}
		
		// setting values for the U matix
		double value; 
		for(int i = 0 ; i < rowNum; i++){
			for(int j= 0; j < colNum; j++){
				value = A_matrix[l_vector[i]][j]; 
				U_matrix[i][j] = value; 
			}
		}
		return U_matrix; 
	};

	// use: algorithm for calculating roots through back subsitution
	// input: matrix and vector 
	// output: root vector which contains the solution
	double* back_subsitution(double ** matrix, double * vector){
		double * vectorR = (double*)malloc((rowNum)*sizeof(double)); 
		vectorR[rowNum-1] = vector[rowNum-1]/matrix[rowNum-1][colNum-1];
		double sum = 0; 
		for(int i = rowNum-2; i >=0; i--){
			sum = vector[i];
			for(int j = i+1; j < rowNum; j++){
				sum = sum - (matrix[i][j] * vectorR[j]);
			}
			vectorR[i] = sum/matrix[i][i]; 
		}
		return vectorR; 
	};

	// use: algorithm for calculating roots through forward subsitution
	// input: matrix and vector 
	// output: root vector which contains the solution
	double* forward_subsitution(double ** matrix, double * vector){
		double * vectorR = (double*)malloc((rowNum)*sizeof(double)); 
		vectorR[0] = vector[0]/matrix[0][0]; 
		double sum = 0; 
		for(int i = 1; i < rowNum; i++){
			sum = vector[i]; 
			for(int j = i-1; j >= 0; j--){
				sum = sum - (matrix[i][j] * vectorR[j]);
			}
			vectorR[i] = sum/matrix[i][i];  
		}
		return vectorR; 
	};

	// use: printing matrix to consol, for debugging or showcasing results
	void print_matrix(double ** matrix, string type){
		cout << "status - printing "<<type <<" matrix to consol" << endl;
		for(int i=0; i< rowNum; i++){
			for(int j = 0; j<colNum; j++)
				cout << matrix[i][j] <<"  ";
			cout << endl; 
		}
		cout << "\nstatus - printing matrix to console complete"<< endl;
	};

	// use: printing vector to consol, for debugging or showcasing results
	void print_vector(double * vector, string type){
		cout << "status - printing "<<type <<" vector to console" << endl; 
		for(int i = 0; i< rowNum; i++)
			cout << type<<i<<": "<<vector[i]<<"\t"; 
		cout << "\nstatus - printing vector to console complete" <<endl;
	};

	// use: to get permutation matrix from the maintained row vector
	// input: no input
	// output: permutation matrix
	double** get_p_matrix(){
		double ** matrix = (double **)malloc((rowNum)*sizeof(double*)); 
		for(int i = 0; i < colNum; i ++)
			matrix[i] = (double*)malloc((colNum)*sizeof(double)); 
		for(int i = 0; i < rowNum ; i++){
			for(int j = 0; j < colNum; j++){
				matrix[i][j] = 0;
				if(j==l_vector[i])
					matrix[i][j] = 1; 
			} 
		}
		return matrix; 
	}

private:
	int rowNum; 
	int colNum; 
	int * l_vector;
};

class MatrixParser{
public:
	MatrixParser(){
		lineNum =0;
	};

	// use: store filename and extract size info of the file
	void read_from_file(string afilename){
		cout <<"status - reading from file: " << afilename <<endl;
		filename = afilename; 
		ifstream infile(filename);
		string file; 
		while(getline(infile, file)){lineNum++;} 
		infile.close();
	};

	// use: parse file to obtain the matrix and b_vector from Ax=b
	void create_matrix(void){
		cout <<"status - creating matrix and vector from file: "<<filename<< endl;
		matrix = (double**)malloc(lineNum * sizeof(double *));
		vector = (double*)malloc(lineNum* sizeof(double)); 
		for(int i =0; i < lineNum; i++)
		 	matrix[i] = (double *)malloc((lineNum) * sizeof(double));

		ifstream infile(filename);
		double number; 

		for(int row = 0; row < lineNum; row++){
			for(int col = 0; col < lineNum+1; col++){
				infile >> number;
				if(col == lineNum)
					vector[row] = number; 
				else
					matrix[row][col] = number;
			}
		}
		infile.close(); 
	};

	// use: print matrix for verfication 
	void print_matrix(void){
		cout << "status - printing matrix to console from file: "<< filename << endl; 
		for(int row = 0; row < lineNum; row++){
			for(int col = 0; col < lineNum+1; col++){
				if(col==lineNum)
					cout << "| " << vector[row] << endl;
				else
					cout << matrix[row][col] << " ";

			}
		}
	}; 

	// use: print matrix to designated file
	void print_matrix_to_file(string afilename, string type, double** matrix){
		cout << "status - writing " << type <<" matrix to file: "<<afilename << endl; 
		ofstream outfile(afilename, ofstream::app); 
		for(int row = 0 ; row < lineNum ; row++){
			for(int col = 0; col < lineNum; col++){
				outfile <<matrix[row][col]<< " "; 
			}
			outfile <<"\n"; 
		} 
		outfile.close();
		cout << "status - writing matrix to file complete" << endl;  
	};

	int get_lineNum(){
		return lineNum; 
	};
	double ** get_matrix(){
		return matrix; 
	};
	double * get_vector(){
		return vector; 
	};
	~MatrixParser(){
		// free matrix;
		cout << "Program end ... freeing memory" << endl;
		for(int i = 0; i < lineNum; i++){
			free(matrix[i]); 
		}
		free(matrix);
		free(vector);  
	}; 
private:
	int lineNum;
	string filename;  
	double** matrix; 
	double* vector; 
};

int main(int argv , char** argc){
	if(argv!=3){
		cout << "Error : Specify Input and Output file \a" << endl; 
		exit(1); 
	}
	// extract data from file
	MatrixParser parse; 
	parse.read_from_file(argc[1]);
	parse.create_matrix();
	parse.print_matrix(); 
	 
	// obtain a matrix from parse, create l matrix, and perform operation to obtain u matrix
	double** a_matrix = parse.get_matrix(); 
	double ** l_matrix = (double **)malloc((parse.get_lineNum())*sizeof(double*));
	for(int i = 0; i < parse.get_lineNum(); i ++)
		l_matrix[i] = (double*)malloc((parse.get_lineNum())*sizeof(double));

	MatrixLUDecomposition LU(parse.get_lineNum(), parse.get_lineNum());
	// 3 param is the PIVOT control 
	double ** u_matrix = LU.gaussian_elemination(a_matrix, l_matrix, 1); 

	//write L, U, P to file ** needs work 
	parse.print_matrix_to_file(argc[2], "L matrix", l_matrix); 
	parse.print_matrix_to_file(argc[2], "U matrix", u_matrix); 
	parse.print_matrix_to_file(argc[2], "P matrix", LU.get_p_matrix()); 
		 
	// solve for root, using LU decomposition, forward_sub returns vector that is used by back_sub
	double* x_vector = LU.back_subsitution(u_matrix, LU.forward_subsitution(l_matrix, parse.get_vector())); 

	//printing solution to console
	LU.print_vector(x_vector, "LUP-X");

	

	MatrixParser parse1; 
	parse1.read_from_file(argc[1]);
	parse1.create_matrix();
	parse1.print_matrix(); 
	 
	// // obtain a matrix from parse, create l matrix, and perform operation to obtain u matrix
	a_matrix = parse1.get_matrix(); 
	MatrixLUDecomposition LU1(parse1.get_lineNum(), parse1.get_lineNum());
	// // 3 param is the PIVOT control 
	u_matrix = LU1.gaussian_elemination(a_matrix, l_matrix, 0); 
	// // solve for root, using LU decomposition, forward_sub returns vector that is used by back_sub
	x_vector = LU1.back_subsitution(u_matrix, LU1.forward_subsitution(l_matrix, parse1.get_vector())); 

	// //LU.print_vector(z_vector); 
	LU.print_vector(x_vector, "LU-X");

	return 0; 
}