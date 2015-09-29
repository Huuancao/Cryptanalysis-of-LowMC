#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <bitset>
#include <algorithm>

using namespace std;
//////////////////
//     MAIN     //
//////////////////
const unsigned blocksize = 16;
typedef std::bitset<blocksize> block;

block MultiplyWithGF2Matrix
        (const std::vector<block> matrix, const block message) {
    block temp = 0;
    for (unsigned i = 0; i < blocksize; ++i) {
        temp[i] = (message & matrix[i]).count() % 2;
    }
    return temp;
}

std::vector<block> invert_Matrix (const std::vector<block> matrix) {
    std::vector<block> mat; //Copy of the matrix 
    for (auto u : matrix) {
        mat.push_back(u);
    }
    std::vector<block> invmat(blocksize, 0); //To hold the inverted matrix
    for (unsigned i = 0; i < blocksize; ++i) {
        invmat[i][i] = 1;
    }

    unsigned size = mat[0].size();
    //Transform to upper triangular matrix
    unsigned row = 0;
    for (unsigned col = 0; col < size; ++col) {
        if ( !mat[row][col] ) {
            unsigned r = row+1;
            while (r < mat.size() && !mat[r][col]) {
                ++r;
            }
            if (r >= mat.size()) {
                continue;
            } else {
                auto temp = mat[row];
                mat[row] = mat[r];
                mat[r] = temp;
                temp = invmat[row];
                invmat[row] = invmat[r];
                invmat[r] = temp;
            }
        }
        for (unsigned i = row+1; i < mat.size(); ++i) {
            if ( mat[i][col] ) {
                mat[i] ^= mat[row];
                invmat[i] ^= invmat[row];
            }
        }
        ++row;
    }

    //Transform to identity matrix
    for (unsigned col = size; col > 0; --col) {
        for (unsigned r = 0; r < col-1; ++r) {
            if (mat[r][col-1]) {
                mat[r] ^= mat[col-1];
                invmat[r] ^= invmat[col-1];
            }
        }
    }

    return invmat;
}

void printMatrix(std::vector<block> mat){
	cout << "Printing Matrix !" << endl;
	for(int i = 0; i<mat.size(); ++i){
			cout << mat[i];
			cout << endl;
	}
	
}

std::vector<block> init_Matrix(){
	std::vector<block> mat;
	block tempMatrixLine;
    tempMatrixLine.reset();
    tempMatrixLine[9]=1;
    cout << "tempMatrixLine: " << tempMatrixLine << endl; 
    for(int i=0; i<blocksize; ++i){
    	
        mat.push_back(tempMatrixLine);
        tempMatrixLine <<= 1;
    }
        
    
    return mat;

}

int main () {
	block message;
	message.reset();
	message[1]=message[4]=message[6]=message[8]=message[11]=1;
	cout << message << endl;
	block temp;


	std::vector<block> identMatrix;

	identMatrix = init_Matrix();
	temp = MultiplyWithGF2Matrix(identMatrix, message);
	cout << "Result mult matrix : " << temp << endl;

	printMatrix(identMatrix);






	

    return 0;
}
