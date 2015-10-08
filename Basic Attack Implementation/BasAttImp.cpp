#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

const unsigned numofboxes = 3; // Number of Sboxes
const unsigned blocksize = 16; // Block size in bits
const unsigned keysize = 6; // Key size in bits
const unsigned rounds = 6; // Number of rounds
const unsigned partialRounds = 4; // Number of rounds to compute high order constants
const unsigned tail = 7; // Number of bits in tail
const unsigned dimension = 12; //Dimension of vector space
const unsigned maxpermut= 4080; //Bound binary:111111110000
const unsigned numSubspaces= 495; //Combination C(12,8)
const unsigned numPartialCiphertexts=4096;
const std::vector<unsigned> Sbox = {0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02}; // Sboxes



const string plainPath = "../LowMC/plaintexts.txt";
const string cipherPath = "../LowMC/ciphertexts.txt";
const string partialCipherPath = "../LowMC/partialCiphertexts.txt";
const string freeCoefPath= "a0.txt";
const string monomialsPath = "monomials.txt";

const unsigned identitysize = blocksize - 3*numofboxes;

typedef std::bitset<blocksize> block; // Store messages and states
typedef std::bitset<keysize> keyblock;
typedef std::bitset<dimension> vecspace;
typedef std::bitset<numPartialCiphertexts> freeCoef;




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
void preprocessingFreeCoef(vector<int>& a0, const vector<block>& partialCiphertexts, 
                        const std::vector<block>& base, const vector<vecspace>& subspaces, 
                        const unsigned int targetBit){
    for(int i=0; i< partialCiphertexts.size(); ++i){
        for(int j=0; j< subspaces.size(); ++j){
            vector<block> tempSubspace;
            for(int k=0; k < subspaces[j].size(); ++k){
                if (subspaces[j][k]){
                    tempSubspace.push_back(base[k]);
                }
            }
            tempSubspace.push_back(partialCiphertexts[i]);
            if(rank_of_Matrix(tempSubspace)==8){
                a0[i]=(a0[i]+partialCiphertexts[i][blocksize-7])%2;
            }
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
Print vector of integers.
*/
void printVector(vector<int>& vector){
    for(int i=0; i<vector.size();++i){
        cout<<vector[i];
    }
    cout << endl;
}

/*
Print vector of sequences of blocks.
*/
void printSequencesBlocks(std::vector<block>& sequences){
    for(int i=0; i< sequences.size(); ++i){
        cout << "Entry n" <<i << ": " << sequences[i] << endl;
    }
}

/*
Print vector of sequences of vecspaces.
*/
void printSequencesVecspaces(std::vector<vecspace>& sequences){
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
Write preprocessed free coef in file.
*/
void writeFreeCoef(const vector<int>& a0){
    ofstream myFile;
    //myFile.open("ciphertexts.txt");
    myFile.open(freeCoefPath.c_str());
    for(int i=0; i< a0.size(); ++i){
        myFile << a0[i];
    }
    myFile.close();
}


bool isInSubspace(vector<vecspace>& subspaces, block v){
    bool isIn= true;
    for (int i = 0; i < subspaces.size(); ++i)
    {
        for (int j = 0; j < blocksize; ++j)
        {
            if (subspaces[i][j]*v[j]==1)
            {
                isIn=false;
                return isIn;
                break;
            }
        }
    }
    return isIn;
}
/*
Functions to multiply bitsets.
*/
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
}
/*
Generate monomials.
*/
void generateMonomials(vector<block>& monomials){
    block firstMonomial(1);
    for(int i=0; i<blocksize; ++i){
        monomials.push_back(firstMonomial<<i);
    }
    for(int j=0; j<numofboxes; ++j){
        for(int k=0; k<numofboxes; ++k){
            block temp(0);
            temp ^= Sbox[ ((monomials[3*j+k] >> 3*j)
                      & block(0x7)).to_ulong()];
            temp = temp << 3*j;
            if(temp.count() == 2){
                monomials.push_back(temp);
            }
        }
    }
    /*
    unsigned int current(3);
    unsigned int max(49152);
    std::vector<block> combinationsMono;
    do{
        combinationsMono.push_back(current);
        nextPermut(current);
    }while(current<=max);
    for(int l=0; l< combinationsMono.size(); ++l){
        vector<block> tempMonomial;
        for(int m=0; m<blocksize; ++m){
            if(combinationsMono[l][m]){
                tempMonomial.push_back(monomials[m]);            
            }
        }
        block resultMult(0);
        bitsetMultiply(resultMult, tempMonomial[0], tempMonomial[1]);
        if(resultMult!=0){
            monomials.push_back(resultMult);
        }
    }
    */
    unsigned int sizeMonomials = monomials.size();
    for(int l=0;l<sizeMonomials; ++l){
        for (int m=l+1; m<sizeMonomials; ++m){
            block resultMult(0);
            bitsetMultiply(resultMult, monomials[l], monomials[m]);
            if(resultMult!=0)
                monomials.push_back(resultMult);
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


void generateMatrixA(vector<block>& monomials, vector<block>& ciphertexts, vector<block>& A){
    for (int i = 0; i < ciphertexts.size(); ++i)
    {
        for (int j = 0; j < monomials.size(); ++j)
        {
            for (int k = 0; k < blocksize; ++k)
            {
                if (ciphertexts[i][k]==0 && monomials[j][k]==1)
                {
                    A[i][j]=0;
                    //A[i].push_back(0);
                    break;
                }else{
                    A[i][j]=1;
                    //A[i].push_back(1);
                }
            }
        }
    }
}

void generateMatrixE(vector<block>& A, const vector<vecspace>& subspaces,const std::vector<block>& base, vector<block>& E){
    for (int i = 0; i < subspaces.size(); ++i)
    {
        vector<block> tempSubspace;
        for(int k=0; k < subspaces[i].size(); ++k){
                if (subspaces[i][k]){
                    tempSubspace.push_back(base[k]);
                }
            }
        for (int j = 0; j < A.size(); ++j)
        {
            tempSubspace.push_back(A[j]);
            if(rank_of_Matrix(tempSubspace)==8){
                E[i]=E[i]^A[j];
            }
            tempSubspace.pop_back();
        }
            
        }
        
}

void setUpEquation(vector<block>& E, const vector<int>& a0){
    for (int i = 0; i < E.size(); ++i)
    {
        E[i].set(E.size(),a0[i]);
    }
}
void swapRow(vector<block>& E,int i, int j){
    block temp= E[i];
    E[i]=E[j];
    E[j]=temp;
}

void solveEquation(vector<block>& E){
    for (int i = 0; i < E.size(); ++i)
    {
        int k = 0;
        while(E[k][i] != 1) k++;
        swapRow(E,k,i);
        for (int j = i+1; j < E.size(); ++j)
        {
            if (E[j][i]==1)
            {
                E[j]=E[i]^E[j];
            }
            
        }
    }
    for (int i = E.size()-1; i >= 0; --i)
    {
        for (int j = i-1;  j>=0; --j)
        {
            if (E[j][i]==1)
            {
                E[j]=E[i]^E[j];
            }
        }
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
    freeCoef a0;
    unsigned int targetBit(9);
    initInputs(plaintexts, plainPath);
    initInputs(ciphertexts, cipherPath);
    initInputs(partialCiphertexts, partialCipherPath);
    initInputs(a0, freeCoefPath);
    initInputs(monomials, monomialsPath);
    setVectorSpace(base);
    setSubspaces(subspaces);

    generateMonomials(monomials);


    printSequencesBlocks(monomials);
    //writeVectorsBlocks(monomials, monomialsPath);



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