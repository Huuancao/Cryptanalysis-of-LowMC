#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <fstream>
#include <cmath>
#include <set>


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
const unsigned relationLength = 22;
const unsigned identitysize = blocksize - 3*numofboxes;
const std::vector<unsigned> Sbox = {0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02}; // Sboxes
const std::vector<unsigned> invSbox = {0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04}; // Invers Sboxes



const string plainPath = "../LowMC/plaintexts.txt";
const string cipherPath = "../LowMC/ciphertexts.txt";
const string partialCipherPath = "../LowMC/partialCiphertexts.txt";
const string linMatPath = "../LowMC/linmatrices.txt";
const string keyMatPath = "../LowMC/keymatrices.txt";
const string roundConstPath = "../LowMC/roundconstants.txt";
const string freeCoefPath= "a0.txt";
const string monomialsPath = "monomials.txt";
const string pythonPath1 ="python1.txt";
const string pythonPath2 ="python2.txt";
const string invLinMatPath ="invlinmatrices.txt";
const string peelOffCipherPath ="peeledOffCiphertexts.txt";

typedef std::bitset<blocksize> block; // Store messages and states
typedef std::bitset<keysize> keyblock;
typedef std::bitset<dimension> vecspace;
typedef std::bitset<nummonomials> monomatrix;
typedef std::bitset<numSubspaces> freeCoef;
typedef std::bitset<relationLength> relationRepresentation;



//////////////////
//     CLASS    //
//////////////////

class relationRepComp{
   public:
    bool operator()( relationRepresentation const & r1, relationRepresentation const & r2 ) const{
         return r1.to_ulong() < r2.to_ulong();
    }
};

typedef set<relationRepresentation, relationRepComp>  relationSetType;

class blockComp{
   public:
    bool operator()( block const & r1, block const & r2 ) const{
         return r1.to_ulong() < r2.to_ulong();
    }
};

typedef set<block, blockComp>  blockSetType;

//////////////////
//  USELESS FUN //
//////////////////
/*
Set up the linear equations system by appending a0.
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

/*
Gaussian elimination.

void gauss(vector< vector<double> >& A) {
    int n = A.size();
    int m = A[0].size();

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
        for (int k=i; k<m;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = A[k][i];
            for (int j=i; j<m; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] = fmod(A[k][j] + A[i][j], 2);
                }
            }
        }


    }
    for (int i = A.size()-1; i >= 0; --i){
        int k=0;
        while(A[i][k] == 0  && k< A[i].size()) k++; 
        if (k < A[i].size()){
            for (int j = i-1;  j>=0; --j){
                if (A[j][k]==1){
                    for (int a = k; a < A[i].size() ; ++a)
                    {
                        A[j][a]=fmod(A[i][a]+A[j][a],2);
                    }
                     
                }
            }
        }       
    }
    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
    
}
*/
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


//////////////////
//   FUNCTIONS  //
//////////////////
/*
Computes GCD.
*/
unsigned long long
gcd(unsigned long long x, unsigned long long y){
    while (y != 0){
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
    for (unsigned long long d=1; d <= k; ++d, --n){
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
Compute preprocessed free coefs.
*/
void preprocessingFreeCoef(vector<freeCoef>& a0, const vector<block>& partialCiphertexts, 
                        const vector<block>& plaintexts,
                        const std::vector<block>& base, const vector<vecspace>& subspaces){
    for(int i=0; i< subspaces.size(); ++i){
        vector<block> tempSubspace;
        for(int j=0; j < subspaces[i].size(); ++j){
            if (subspaces[i][j]){
                tempSubspace.push_back(base[j]);
            }
        }
        for(int k=0; k<plaintexts.size(); ++k){
            tempSubspace.push_back(plaintexts[k]);
            if(rank_of_Matrix(tempSubspace)==8){
                for(int l=0; l< a0.size(); ++l){
                    a0[l][i]=a0[l][i]^partialCiphertexts[k][l];
                }
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
            cout << vector[i][j];
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
Print set of blocks.
*/
void printSequencesBlocks(const blockSetType& blockSet){
    for(blockSetType::iterator i=blockSet.begin(); i!=blockSet.end();++i){
        cout << "Entry n" << distance(blockSet.begin(),i) << ": " << *i << endl;
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
Print vector of vectors
*/
void printVectorVectorsBlock(vector<vector<block>>& vector){
    for(int i=0; i< vector.size(); ++i){
        cout << "Entry n" <<i << ": "<<endl;
        for(int j=0; j< vector[i].size(); ++j){
            cout << vector[i][j] << endl;
        }
        cout << endl;
    }
}
void printVectorVectorsKeyBlock(vector<vector<keyblock>>& vector){
    for(int i=0; i< vector.size(); ++i){
        cout << "Entry n" <<i << ": " << endl;
        for(int j=0; j< vector[i].size(); ++j){
            cout << vector[i][j] << endl;
        }
        cout << endl;
    }
}
/*
Print Relations map.
*/
void printRelationMap(vector<vector<vector<relationRepresentation>>>& relationMap){
    for(int i=0; i<relationMap.size(); ++i){
        cout << "Round " << i << endl;
        for(int j=0; j< blocksize; ++j){
            cout << "Bit "<< j << endl;
            for(int k=0; k<relationMap[0][j].size(); ++k){ 
                cout << "Relation Element " << k << ": " << relationMap[0][j][k] << endl;
            }
        }
    }
}
/*
Read file and set inputs in vector of vector of blocks vector<vector<block> linearMatrices.
*/
void initInputsLinearMatrices(vector<vector<block>>& linearMatrices, string filePath){
    ifstream myFile(filePath.c_str());
    if(myFile){
        block temp(0);
        vector<block> tempVector;
        tempVector.push_back(temp);
        string bitLine;
        int increment(0);
        for (int i=0; i < rounds; i++){
            linearMatrices.push_back(tempVector);
        }
        while (getline(myFile, bitLine)){   
            if(bitLine.empty()){
                linearMatrices[increment].erase(linearMatrices[increment].begin());
                ++increment;
            }else{
                block b(bitLine);
                linearMatrices[increment].push_back(b);
            }
        }         
    }
    else{
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    myFile.close();
}
/*
Read file and set inputs in vector of vector of blocks vector<vector<block> linearMatrices.
*/
void initInputsKeyMatrices(vector<vector<keyblock>>& keyMatrices, string filePath){
    ifstream myFile(filePath.c_str());
    if(myFile){
        keyblock temp(0);
        vector<keyblock> tempVector;
        tempVector.push_back(temp);
        string bitLine;
        int increment(0);
        for (int i=0; i < 7; i++){
            keyMatrices.push_back(tempVector);
        }
        while (getline(myFile, bitLine)){
            if(bitLine.empty()){
                keyMatrices[increment].erase(keyMatrices[increment].begin());
                ++increment;
            }else{
                keyblock b(bitLine);
                keyMatrices[increment].push_back(b);
            }
        }      
    }
    else{
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    myFile.close();
}
/*
Read file and set inputs in vector of blocks.
*/
void initInputs(vector<block>& textfile, string filePath){
    ifstream myFile(filePath.c_str());
    if(myFile){
        string bitLine;
        while (getline(myFile, bitLine))
        {
            block b(bitLine);
            textfile.push_back(b);
        }
    }
    else{
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    myFile.close();
}
/*
Read file and set inputs in vector of blocks.
*/
void initInputs(vector<freeCoef>& inputs, string filePath){
    ifstream myFile(filePath.c_str());
    if(myFile){
        string bitLine;
        int increment(0);
        
        while (getline(myFile, bitLine)){
            freeCoef b(bitLine);
            inputs[increment]=b;
            ++increment;
        }
    }
    else{
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    myFile.close();
}
/*
Read file and set inputs in set of blocks.
*/
void initInputs(blockSetType& monomials, string filePath){
    ifstream myFile(filePath.c_str());
    if(myFile){
        string bitLine;
        while (getline(myFile, bitLine))
        {
            block b(bitLine);
            monomials.insert(b);
        }
    }
    else{
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    myFile.close();
}

/*
Write preprocessed free coefs in file.
*/
void writeFreeCoef(const vector<freeCoef>& a0){
    ofstream myFile;
    remove(freeCoefPath.c_str());
    myFile.open(freeCoefPath.c_str());
    for(int i=0; i< a0.size(); ++i){
        for(int j=0; j< a0[i].size(); ++j){
            myFile << a0[i][j];
        }
        myFile << endl;
    }
    myFile.close();
}
/*
Write Inputs for Sage.
*/
writePython(vector<monomatrix>& matrixE, vector<freeCoef>& a0){
    ofstream myFile;
    myFile.open(pythonPath1.c_str());
    myFile << "[";
    for(int i=0; i< matrixE.size(); ++i){
        myFile << "[";
        for(int j=0; j< matrixE[0].size(); ++j){
            if (j==matrixE[0].size()){
                    myFile << matrixE[i][j];
            }else{
                 myFile << matrixE[i][j]<< " ";
            }
        }
        if(i != matrixE.size()) 
            {
                myFile << "],";
        }else{
            myFile << "]";
        }
    }
    myFile << "]";
    myFile.close();

    myFile.open(pythonPath2.c_str());
    myFile << "[";
    for(int i=0; i< a0.size(); ++i){
        myFile << "[";
        for(int j=0; j< a0[0].size(); ++j){
            if (j==a0[0].size()){
                    myFile << a0[i][j];
            }else{
                 myFile << a0[i][j]<< " ";
            }
        }
        if(i != a0.size()-1) 
            {
                myFile << "],";
        }else{
            myFile << "]";
        }
    }
    myFile << "]";
    myFile.close();
}
/*
Functions to multiply bitsets.
*/
void bitsetMultiply(block& result, const block& x, const block& y){
    result = x|y;
}
/*
Generate monomials.
*/
void generateMonomials(blockSetType& monomials){
    //u1 generation
    block firstMonomial(1);
    for(int i=0; i<blocksize; ++i){
        monomials.insert(firstMonomial<<i);
    }
    //u2 generation
    for(int j=0; j<numofboxes; ++j){
        int permut(3);
        for(int k=0; k<boxsize; ++k){
            block tempMono(permut);
            tempMono <<= (numofboxes-1-j)*boxsize;
            monomials.insert(tempMono);
            permut=nextPermut(permut);
        }
    }
    //u3 generation
    blockSetType temp;
    for(blockSetType::iterator k=monomials.begin(); k!=monomials.end(); ++k){
        for(blockSetType::iterator l=monomials.begin(); l!=monomials.end(); ++l){
            block resultMult(0);
            bitsetMultiply(resultMult, *k, *l);
            if(monomials.find(resultMult)==monomials.end() && temp.find(resultMult)==temp.end()){
                temp.insert(resultMult);
            }
        }
    }
    for(auto element : temp){
        monomials.insert(element);
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
Write blocks in file.
*/
void writeBlockSet(const blockSetType& blockSet, const string fileName){
    ofstream myFile;
    //myFile.open("ciphertexts.txt");
    myFile.open(fileName.c_str());
    for(auto block : blockSet){
        myFile << block << endl;
    }
    myFile.close();
}
/*
Write Matrices in file.
*/
void writeMatrices(std::vector<std::vector<block>>& matrix, string fileName){
    std::ofstream myFile;
    myFile.open(fileName.c_str());

    for(int i=0; i<matrix.size();++i){
        for(int j=0; j<matrix[i].size(); ++j){
            myFile << matrix[i][j] << std::endl;
        }
        myFile << std::endl;
    }
    myFile.close();
}
/*
Generate Matrix A, Prod c_i^u_i.
*/
void generateMatrixA(vector<block>& monomials, vector<block>& peeledOffCiphertexts, vector<monomatrix>& matrixA){
    for (int i = 0; i < peeledOffCiphertexts.size(); ++i){
        for (int j = 0; j < monomials.size(); ++j){
            for (int k = 0; k < blocksize; ++k){
                if (!peeledOffCiphertexts[i][k] && monomials[j][k]){
                    matrixA[i][j]=0;
                    break;
                }else{    
                    matrixA[i][j]=1;
                }
            }
        }
    }
}
/*
Generate Matrix E by testing if ciphertext J belong to subspace i and then adding the corresponding jth row of A in the ith-row of E
*/
void generateMatrixE(const vector<monomatrix>& A,const vector<block>& plaintexts,  const vector<block>& ciphertexts, 
                    const vector<vecspace>& subspaces,const std::vector<block>& base, 
                    vector<monomatrix>& E){
    for (int i = 0; i < subspaces.size(); ++i){
        vector<block> tempSubspace;
        for(int k=0; k < subspaces[i].size(); ++k){
            if (subspaces[i][k]){
                tempSubspace.push_back(base[k]);
            }
        }
        for (int j = 0; j < ciphertexts.size(); ++j){
            tempSubspace.push_back(plaintexts[j]);
            if(rank_of_Matrix(tempSubspace)==8){
                E[i]=E[i]^A[j];
            }
        tempSubspace.pop_back();
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
/*
Invert Matrix.
*/
vector<block> invertMatrix(const vector<block>& matrix) {
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
/*
Initialize inverted matrices of linear matrices.
*/
void initInvMatrices(vector<vector<block>>& linearMatrices, vector<vector<block>>& invLinearMatrices){
    for(int i=0; i < linearMatrices.size(); ++i){
        invLinearMatrices.push_back(invertMatrix(linearMatrices[i]));
    }
}
/*
Multiply with matrix in GF(2).
*/
block MultiplyWithGF2Matrix(const std::vector<block> matrix, const block message) {
    block temp = 0;
    for (unsigned i = 0; i < blocksize; ++i) {
        temp[i] = (message & matrix[i]).count() % 2;
    }
    return temp;
}
/*
Peel off last layer by adding the last round constants and multiplying with inverse of last linear matrix.
*/
void peelingOffCiphertexts(const vector<block>& ciphertexts, const block& roundConstant, 
                        const vector<block>& invLinearMatrix, vector<block>& peeledOffCiphertexts){
    for(int i=0; i< ciphertexts.size(); ++i){
        block temp(0);
        temp=ciphertexts[i]^roundConstant;
        peeledOffCiphertexts.push_back(MultiplyWithGF2Matrix(invLinearMatrix, temp));
    }
}
/*
Linear layer function.
*/
void lineaLayerMixing(vector<vector<vector<relationRepresentation>>>& relationMap,
                      const vector<block>& linearMatrix){
    vector<vector<relationRepresentation>> tempRelationVectorVectors;
    vector<relationRepresentation> tempRelationVectors;
    for(int i=0; i<blocksize; ++i){
        for(int k=0; k<blocksize; ++k){
            if(linearMatrix[i][k]){
                for(int l=0; l<relationMap[0][i].size(); ++l){
                        relationRepresentation tempOutMono(relationMap[0][k][l]<<6);
                        relationRepresentation tempOutKey(relationMap[0][k][l]>>16);
                        tempOutMono >>=6;
                        tempOutKey <<=16;
                        bool alreadyIn(false);
                    if(k==i && linearMatrix[i][k]){
                    break;
                    }
                    else {
                        for(int m=0; m<tempRelationVectors.size(); ++m){
                            relationRepresentation tempInMono(tempRelationVectors[m]<<6);
                            tempInMono >>=6;
                            if(tempInMono.to_ullong() == tempOutMono.to_ullong()){
                                alreadyIn=true;
                                tempRelationVectors[m]=tempRelationVectors[m]^tempOutKey;
                            }  
                        }
                        if(!alreadyIn){
                            tempRelationVectors.push_back(relationMap[0][k][l]);
                        }
                    }
                }
            }
        }
        tempRelationVectorVectors.push_back(tempRelationVectors);
        tempRelationVectors.clear();
    }
    relationMap.clear();
    relationMap.push_back(tempRelationVectorVectors);
}
/*
Adding key function.
*/
void keyRoundAdd(vector<vector<relationRepresentation>>& tempRelation, const vector<keyblock>& keyMatrix){
    for(int i=0; i< blocksize; ++i){
        relationRepresentation tempKey(keyMatrix[i].to_ullong());
        tempKey<<=blocksize;
        for(int j=0; j< tempRelation[i].size(); ++j){
           tempRelation[i][j] = tempRelation[i][j]^tempKey;
        }
    }
}
/*
Init iniput first round & key whitening.
*/
void initRelationWhitening(vector<vector<vector<relationRepresentation>>>& relationMap,
                        const vector<vector<keyblock>>& keyMatrices){
    vector<vector<relationRepresentation>> tempRelationVectorVectors;  
    vector<relationRepresentation> tempRelationVectors;
    relationRepresentation tempRelation(1);
    for(int i=0; i<blocksize; ++i){
        tempRelationVectors.clear();
        tempRelationVectors.push_back(tempRelation);
        tempRelation <<=1;
        tempRelationVectorVectors.push_back(tempRelationVectors);
    }
    keyRoundAdd(tempRelationVectorVectors, keyMatrices[0]);
    relationMap.push_back(tempRelationVectorVectors);
}
void insertRemastered(vector<relationRepresentation>& InsertionResult, vector<relationRepresentation>& toInsert){
    for (int i = 0; i < toInsert.size(); ++i){
        InsertionResult.push_back(toInsert[i]);
    }
}
/*
SBoxes function for a vector of bitset of size relationLength 
*/
void SBoxRelation(vector<relationRepresentation>& a, 
    vector<relationRepresentation>& b, 
    vector<relationRepresentation>& c){
    vector<relationRepresentation> tempSBoxA;
    vector<relationRepresentation> tempSBoxB;
    vector<relationRepresentation> tempSBoxC;
    
    for(int i = 0; i < a.size(); ++i){
        for (int j = 0; j < b.size(); ++j){
            for (int k = 0; k < c.size(); ++k){
                tempSBoxA.push_back(b[j]|c[k]);
                tempSBoxB.push_back(a[i]|c[k]);
                tempSBoxC.push_back(a[i]|b[j]);
            }
        }
    }
    insertRemastered(b,a);
    insertRemastered(c,b);
    insertRemastered(a,tempSBoxA);
    insertRemastered(b,tempSBoxB);
    insertRemastered(c,tempSBoxC);   
}
/*
Relation mapping creation.
*/
void relationMapping(vector<vector<vector<relationRepresentation>>>& relationMap, 
                        const vector<vector<block>>& linearMatrices,
                        const vector<vector<keyblock>>& keyMatrices){
    initRelationWhitening(relationMap, keyMatrices);
    lineaLayerMixing(relationMap, linearMatrices[0]);

}

//////////////////
//     MAIN     //
//////////////////




int main(int argc, const char * argv[]) {
    vector<block> plaintexts;
    vector<block> ciphertexts;
    vector<block> partialCiphertexts;
    vector<block> peeledOffCiphertexts;

    vector<block> base;
    vector<vecspace> subspaces;

    //vector<block> monomials;
    blockSetType monomials;
    blockSetType monomialsv1;
    vector<freeCoef> a0(blocksize, 0);

    vector<monomatrix> matrixA(numPartialCiphertexts,0);
    vector<monomatrix> matrixE(numSubspaces, 0);

    vector<vector<block>> linearMatrices;
    vector<vector<block>> invLinearMatrices;
    vector<vector<keyblock>> keyMatrices;
    vector<block> roundConstants;

    //vector<vector<vector<relationRepresentation>>> relationMap;
    relationSetType relationMap;
 /*   relationRepresentation temp(22);
    relationRepresentation temp2(1);
    relationRepresentation temp21(1942);
    relationRepresentation temp1(0);

    relationMap.insert(temp);
    relationMap.insert(temp2);
    relationMap.insert(temp21);
    relationMap.insert(temp1);
    for(relationSetType::iterator i = relationMap.begin();i!=relationMap.end(); ++i){
        cout << *i << endl;
    }
*/

    initInputs(plaintexts, plainPath);
    initInputs(ciphertexts, cipherPath);
    initInputs(partialCiphertexts, partialCipherPath);
    initInputs(a0, freeCoefPath);
    initInputs(monomials, monomialsPath);
    initInputs(peeledOffCiphertexts, peelOffCipherPath);
    setVectorSpace(base);
    setSubspaces(subspaces);
    initInputsLinearMatrices(linearMatrices, linMatPath);
    initInputsKeyMatrices(keyMatrices, keyMatPath);
    initInputs(roundConstants, roundConstPath);
    initInputsLinearMatrices(invLinearMatrices, invLinMatPath);
    

    //relationMapping(relationMap, linearMatrices, keyMatrices);
    //printRelationMap(relationMap);
    



    //peelingOffCiphertexts(ciphertexts, roundConstants[5], invLinearMatrices[5], peeledOffCiphertexts);
    //printSequencesBlocks(peeledOffCiphertexts);
    //writeVectorsBlocks(peeledOffCiphertexts, peelOffCipherPath);
    
    //initInvMatrices(linearMatrices, invLinearMatrices);
    //printVectorVectorsBlock(invLinearMatrices);
    //writeMatrices(invLinearMatrices, invLinMatPath);

    //printVectorVectorsBlock(linearMatrices);
    //printVectorVectorsKeyBlock(keyMatrices);
    //printVectorVectorsBlock(invLinearMatrices);
    //preprocessingFreeCoef(a0, partialCiphertexts, plaintexts, base, subspaces);
    //writeFreeCoef(a0);
    //generateMonomials(monomials);
    printSequencesBlocks(monomials);
    //cout << "Previous monomials equal to new monomials set? " << (monomials == monomialsv1) << endl;
    //writeBlockSet(monomials, monomialsPath);

    //printANF("");

    //generateMatrixA(monomials, ciphertexts, matrixA);
    //generateMatrixE(matrixA, plaintexts, ciphertexts,subspaces, base, matrixE);

    //printSequencesMonoMatrices(matrixA);
    //printSequencesMonoMatrices(matrixE);

    //writePython(matrixE, a0);
    
    //printSequencesVecspaces(subspaces);
    //printSequencesBlocks(base);
    //printSequencesBlocks(plaintexts);
    //printSequencesBlocks(ciphertexts);
    //printSequencesBlocks(partialCiphertexts);

    return 0;
}
