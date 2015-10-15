#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

const unsigned numofboxes = 3; // Number of Sboxes
const unsigned boxsize = 3; //Number of bits in Sbox
const unsigned blocksize = 16; // Block size in bits
const unsigned keysize = 6; // Key size in bits
const unsigned rounds = 6; // Number of rounds
const unsigned partialRounds = 4; // Number of rounds to compute high order constants
const unsigned tail = 7; // Number of bits in tail
const unsigned dimension = 12; //Dimension of vector space
const unsigned maxpermut= 4080; //Bound binary:111111110000
const unsigned numSubspaces=495; //Combination C(12,8)
const unsigned nummonomials=283;
const unsigned numPartialCiphertexts=4096;
const std::vector<unsigned> Sbox = {0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02}; // Sboxes
const std::vector<unsigned> invSbox = {0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04}; // Invers Sboxes



const string plainPath = "../LowMC/plaintexts.txt";
const string cipherPath = "../LowMC/ciphertexts.txt";
const string partialCipherPath = "../LowMC/partialCiphertexts.txt";
const string freeCoefPath= "a0.txt";
const string monomialsPath = "monomials.txt";

const unsigned identitysize = blocksize - 3*numofboxes;

typedef std::bitset<blocksize> block; // Store messages and states
typedef std::bitset<keysize> keyblock;
typedef std::bitset<dimension> vecspace;
typedef std::bitset<nummonomials> monomatrix;
typedef std::bitset<numSubspaces> freeCoef;




//////////////////
//   FUNCTIONS  //
//////////////////
/*
Computes GCD.
*/
unsigned long long
gcd(unsigned long long x, unsigned long long y)
{
    while (y != 0)
    {
        unsigned long long t = x % y;
        x = y;
        y = t;
    }
    return x;
}
/*
Compute number of combinations.
*/
unsigned long long
choose(unsigned long long n, unsigned long long k){

    unsigned long long r(1);
    for (unsigned long long d=1; d <= k; ++d, --n)
    {
        unsigned long long g = gcd(r, d);
        r /= g;
        unsigned long long t = n / (d / g);
        r *= t;
    }
    return r;
}
/*
Compute Rank of a Matrix.
*/
unsigned rank_of_Matrix (const std::vector<block> matrix) {
    std::vector<block> mat; //Copy of the matrix 
    for (auto u : matrix) {
        mat.push_back(u);
    }
    unsigned size = mat[0].size();
    //Transform to upper triangular matrix
    unsigned row = 0;
    for (unsigned col = 1; col <= size; ++col) {
        if ( !mat[row][size-col] ) {
            unsigned r = row;
            while (r < mat.size() && !mat[r][size-col]) {
                ++r;
            }
            if (r >= mat.size()) {
                continue;
            } else {
                auto temp = mat[row];
                mat[row] = mat[r];
                mat[r] = temp;
            }
        }
        for (unsigned i = row+1; i < mat.size(); ++i) {
            if ( mat[i][size-col] ) mat[i] ^= mat[row];
        }
        ++row;
        if (row == size) break;
    }
    return row;
}
/*
Compute
*/
void preprocessingFreeCoef(freeCoef& a0, const vector<block>& partialCiphertexts, 
                        const std::vector<block>& base, const vector<vecspace>& subspaces, 
                        const unsigned int targetBit){
    for(int i=0; i< subspaces.size(); ++i){
        vector<block> tempSubspace;
        for(int j=0; j < subspaces[i].size(); ++j){
            if (subspaces[i][j]){
                tempSubspace.push_back(base[j]);
            }
        }
        for(int k=0; k<partialCiphertexts.size(); ++k){
            tempSubspace.push_back(partialCiphertexts[k]);
            if(rank_of_Matrix(tempSubspace)==8){
                a0[i]=a0[i]^partialCiphertexts[k][targetBit];
            }
            tempSubspace.pop_back();
        }
    }
}
/*
Generate next lexicographically permutation.
*/
unsigned int nextPermut(unsigned int currentPermut){
    unsigned int v(currentPermut); // current permutation of bits 
    unsigned int w; // next permutation of bits

    unsigned int t = v | (v - 1); // t gets v's least significant 0 bits set to 1
    // Next set to 1 the most significant bit to change, 
    // set to 0 the least significant ones, and add the necessary 1 bits.
    w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));  
    vecspace v1(v);
    vecspace v2(w);
    //cout << "current: " << v1 << " next: " << v2 << endl;
    return w;
}
/*
Generate C(12,8) subspaces.
*/
void setSubspaces(vector<vecspace>& subspaces){
    unsigned int current(255);
    do {
            subspaces.push_back(current);
            current = nextPermut(current);
    }while(current <= maxpermut);
}
/*
Generate vector space 12x16.
*/
void setVectorSpace(std::vector<block>& base){
    block tempVector(0);
    tempVector[0]=1;
    for(int i=0; i<dimension; ++i){
        base.push_back(tempVector);
        tempVector = tempVector << 1;
    }
}
/*
Print vector of vectors
*/
void printVectorVectors(vector<vector<double>>& vector){
    for(int i=0; i< vector.size(); ++i){
        cout << "Entry n" <<i << ": ";
        for(int j=0; j< vector[i].size(); ++j){
            cout << vector[i][j] << " ";
        }
        cout << endl;
    }
}
/*
Print Bitset freeCoef.
*/
void printFreeCoef(freeCoef& vector){
    for(int i=0; i<vector.size();++i){
        cout<<vector[i];
    }
    cout << endl;
}

/*
Print vector of integers.
*/
void printVector(vector<double>& vector){
    for(int i=0; i<vector.size();++i){
        cout<<vector[i];
    }
    cout << endl;
}

/*
Print vector of sequences of blocks.
*/
void printSequencesBlocks(const std::vector<block>& sequences){
    for(int i=0; i< sequences.size(); ++i){
        cout << "Entry n" <<i << ": " << sequences[i] << endl;
    }
}

/*
Print vector of sequences of vecspaces.
*/
void printSequencesVecspaces(const std::vector<vecspace>& sequences){
    for(int i=0; i< sequences.size(); ++i){
        cout << "Entry n" <<i << ": " << sequences[i] << endl;
    }
}

/*
Print vector of sequences of vecspaces.
*/
void printSequencesMonoMatrices(const std::vector<monomatrix>& sequences){
    for(int i=0; i< sequences.size(); ++i){
        cout << "Entry n" <<i << ": " << sequences[i] << endl;
    }
}

/*
Read file and set inputs in vector of blocks.
*/
void initInputs(vector<block>& textfile, string filePath){
    ifstream myFile(filePath.c_str());
    if(myFile)
    {
        string bitLine;
        
        while (getline(myFile, bitLine))
        {
            block b(bitLine);
            textfile.push_back(b);
        }
    }
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    myFile.close();
}

/*
Read file and set inputs in vector of blocks.
*/
void initInputs(freeCoef& inputs, string filePath){
    ifstream myFile(filePath.c_str());
    if(myFile)
    {
        string bitLine;
        
        while (getline(myFile, bitLine))
        {
            freeCoef b(bitLine);
            inputs=b;
        }
    }
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    myFile.close();
}

/*
Write preprocessed free coefs in file.
*/
void writeFreeCoef(const freeCoef& a0){
    ofstream myFile;
    //myFile.open("ciphertexts.txt");
    myFile.open(freeCoefPath.c_str());
    for(int i=0; i< a0.size(); ++i){
        myFile << a0[i];
    }
    myFile.close();
}

/*
bool isInSubspace(vector<vecspace>& subspaces, block v){
    bool isIn= true;
    for (int i = 0; i < subspaces.size(); ++i){
        for (int j = 0; j < blocksize; ++j){
            if (subspaces[i][j]*v[j]){
                isIn=false;
                return isIn;
                break;
            }
        }
    }
    return isIn;
}*/
/*
Bitset multiplication functions
bool fullAdder(bool bit1, bool bit2, bool& carry){
    bool sum = (bit1^bit2)^carry;
    carry = ((bit1 && bit2) || (bit1 && carry) || (bit2 && carry));
    return sum;
}

void bitsetAdd(block& x, const block& y){
    bool carry = false;
    for(int i=0; i < x.size(); ++i){
        x[i]= fullAdder(x[i], y[i], carry);
    }
}
void bitsetMultiply(block& result, const block& x, const block& y){
    block temp = x;
    result.reset();
    if(temp.count() < y.count()){
        for (int i=0; i < result.size(); ++i){
            if (x[i]) {
                bitsetAdd(result , y << i);
            }
        }
    }
    else{
        for (int i=0; i < result.size(); i++){
            if (y[i]) {
                bitsetAdd(result, temp << i);
            }
        }
    }
}*/
/*
Functions to multiply bitsets.
*/
void bitsetMultiply(block& result, const block& x, const block& y){
    result = x|y;
}
/*
Generate monomials.
*/
void generateMonomials(vector<block>& monomials){
    //u1 generation
    block firstMonomial(1);
    for(int i=0; i<blocksize; ++i){
        monomials.push_back(firstMonomial<<i);
    }
    //u2 generation
    for(int j=0; j<numofboxes; ++j){
        int permut(3);
        for(int k=0; k<boxsize; ++k){
            block tempMono(permut);
            tempMono <<= (numofboxes-1-j)*boxsize;
            monomials.push_back(tempMono);
            permut=nextPermut(permut);
        }
    }

    //u3 generation
    unsigned int sizeMonomials = monomials.size();
    for(int l=0;l<sizeMonomials; ++l){
        for (int m=l+1; m<sizeMonomials; ++m){
            bool alreadyIn(false);
            block resultMult(0);
            bitsetMultiply(resultMult, monomials[l], monomials[m]);
            for(int n=0; n<monomials.size(); ++n){
                if(resultMult==monomials[n]){
                    alreadyIn=true;
                }
            }
            if(resultMult!=0 && !alreadyIn){
                monomials.push_back(resultMult);
            }
        }
    }
}

/*
Write blocks in file.
*/
void writeVectorsBlocks(const vector<block>& vectorBlocks, const string fileName){
    ofstream myFile;
    //myFile.open("ciphertexts.txt");
    myFile.open(fileName.c_str());
    for(int i=0; i< vectorBlocks.size(); ++i){
        myFile << vectorBlocks[i] << endl;
    }
    myFile.close();
}

/*
Generate Matrix A, Prod c_i^u_i
*/
void generateMatrixA(vector<block>& monomials, vector<block>& ciphertexts, vector<monomatrix>& matrixA){
    for (int i = 0; i < ciphertexts.size(); ++i){
        matrixA.push_back(0);
        for (int j = 0; j < monomials.size(); ++j){
            for (int k = 0; k < blocksize; ++k){
                //cout << j << " " << k << endl;
                if (!ciphertexts[i][k] && monomials[j][k]){
                    matrixA[i][j]=0;
                    //A[i].push_back(0);
                    break;
                }else{    
                    matrixA[i][j]=1;
                    //A[i].push_back(1);
                }
            }
        }
    }
}
/*
Generate Matrix E by testing if ciphertext J belong to subspace i and then adding the corresponding jth row of A in the ith-row of E
*/
void generateMatrixE(const vector<monomatrix>& A, const vector<block>& ciphertexts, 
                    const vector<vecspace>& subspaces,const std::vector<block>& base, 
                    vector<monomatrix>& E){
    for (int i = 0; i < subspaces.size(); ++i){
        vector<block> tempSubspace;
        E.push_back(0);
        for(int k=0; k < subspaces[i].size(); ++k){
            if (subspaces[i][k]){
                tempSubspace.push_back(base[k]);
            }
        }
        for (int j = 0; j < ciphertexts.size(); ++j){
            tempSubspace.push_back(ciphertexts[j]);
            if(rank_of_Matrix(tempSubspace)==8){
                E[i]=E[i]^A[j];
            }
        tempSubspace.pop_back();
        }        
    }      
}

void setUpEquation(vector<monomatrix>& E, vector<vector<double>>& linearSystem, const freeCoef& a0){
    for (int i = 0; i < E.size(); ++i){
        for(int j =0; j < E[i].size(); ++j){
            if(E[i][j]){
                linearSystem[i].push_back(1);
            }
            else linearSystem[i].push_back(0);
        }
        if(a0[i]){
            linearSystem[i].push_back(1);
        }
        else linearSystem[i].push_back(0);
    }
}

void gauss(vector< vector<double> >& A) {
    int n = A.size();

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }
/*
    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;

*/
}

void swapRow(vector<block>& E,int i, int j){
    block temp= E[i];
    E[i]=E[j];
    E[j]=temp;
}

void solveEquation(vector<block>& E){
    for (int i = 0; i < E.size(); ++i){
        int k = 0;

        while(!E[k][i]) k++;
        swapRow(E,k,i);
        for (int j = i+1; j < E.size(); ++j){
            if (E[j][i]){
                E[j]=E[i]^E[j];
            }  
        }
    }
    for (int i = E.size(); i >= 0; --i){
        int k=0;
        while(!E[i][k] && k<E.size()) k++; 
        if (k< E.size()){
            for (int j = i-1;  j>=0; --j){
                if (E[j][k]){
                    E[j]=E[i]^E[j]; 
                }
            }
        }       
    }
}



/*
Generate all components of a given integer such that n&i == i
*/
vector<unsigned> components(int n){
    vector<unsigned> list;
    for(int i=0; i < n+1; i++){
        if ((n&i) == i){
            list.push_back(i);
        }
    }
    return list;
}
/*
Transform Sbox to ANF
*/
vector<vector<unsigned>> sboxToANF(vector<unsigned> Box){
    unsigned num = Box.size();
    unsigned dim = log2(num);
    if(pow(2,dim) != num){
        cout << "Invalid size of sbox!" << endl;
    }
    vector<vector<unsigned>> ls;
    for(int i=0; i < dim; ++i){
        vector<unsigned> dc;
        for(int j=0; j< num; ++j){
            unsigned sign(0);
            vector<unsigned> compo = components(j);
            for(int k=0; k<compo.size(); ++k){
                sign=sign^((Box[compo[k]]>>i) & 0x1);
            }
            dc.push_back(sign);
        }
        ls.push_back(dc);
    }
    return ls;
}
/*
Print ANF
*/
void printANF(string mode){
    vector<unsigned> Box;
    if(mode == "reverse"){
        Box=invSbox;
    }
    else{
        Box=Sbox;
    }
    cout << "Printing our the ANF equations." << endl;
    cout << "x0 = c, x1 = b, x2 = a" << endl;
    vector<vector<unsigned>> anf = sboxToANF(Box);
    for(int i=0; i < anf.size(); ++i){
        cout << "y" << i << " = ";
        bool firstTime(true);
        vector<unsigned> ls = anf[i];
        for(int j=0; j<ls.size(); ++j){
            if(ls[j]){
                if(firstTime){
                    firstTime=false;
                }
                else{
                    cout << " + ";
                }
                if(j==0){
                    cout << "1";
                }
                else{
                    unsigned mask = 0x1;
                    while(mask <= j){
                        if((mask&j) > 0){
                            cout << "x" << int(log2(mask));
                        }
                        mask = mask << 1;
                    }
                }
            }
        }
        if(firstTime){
            cout << "0";
        }
        cout << endl;
    }
}


//////////////////
//     MAIN     //
//////////////////




int main(int argc, const char * argv[]) {
    vector<block> plaintexts;
    vector<block> ciphertexts;
    vector<block> partialCiphertexts;

    vector<block> base;
    vector<vecspace> subspaces;

    vector<block> monomials;
    //vector<int> a0(numPartialCiphertexts ,0);
    freeCoef a0(0);

    vector<monomatrix> matrixA;
    vector<monomatrix> matrixE;
    vector<vector <double>> linearSystem(numSubspaces);
    vector<double> linearSystemSolution;
    

    unsigned int targetBit(9);
    initInputs(plaintexts, plainPath);
    initInputs(ciphertexts, cipherPath);
    initInputs(partialCiphertexts, partialCipherPath);
    initInputs(a0, freeCoefPath);
    initInputs(monomials, monomialsPath);
    setVectorSpace(base);
    setSubspaces(subspaces);


    //generateMonomials(monomials);


    //printANF("");

    //printSequencesBlocks(monomials);
    //writeVectorsBlocks(monomials, monomialsPath);


    generateMatrixA(monomials, ciphertexts, matrixA);
    //printSequencesMonoMatrices(matrixA);
    generateMatrixE(matrixA,ciphertexts,subspaces, base, matrixE);
    //printSequencesMonoMatrices(matrixE);

    setUpEquation(matrixE, linearSystem, a0);
    gauss(linearSystem);

    //linearSystemSolution=gauss(linearSystem);
    printVectorVectors(linearSystem);
    
    //solveEquation(A);
    //printSequencesBlocks(A);


    //preprocessingFreeCoef(a0, partialCiphertexts, base, subspaces, targetBit);
    //writeFreeCoef(a0);
    
    //cout << a0 << endl;
    //printVector(a0);
    //printSequencesVecspaces(subspaces);
    //printSequencesBlocks(base);
    //printSequencesBlocks(plaintexts);
    //printSequencesBlocks(ciphertexts);
    //printSequencesBlocks(partialCiphertexts);

    return 0;
}