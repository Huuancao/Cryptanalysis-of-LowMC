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
const unsigned blocksize = 21; // Block size in bits
const unsigned keysize = 6; // Key size in bits
const unsigned rounds = 5; // Number of rounds
const unsigned partialRounds = 2; // Number of rounds to compute high order constants
const unsigned tail = 12; // Number of bits in tail
const unsigned dimension = 12; //Dimension of vector space
const unsigned subDimension = 4; //Dimension of vector subspace
const unsigned firstpermut = 15; // 000000001111
const unsigned maxpermut = 3840; //Bound binary:111100000000
const unsigned numSubspaces = 495; //Combination C(12,4)
const unsigned nummonomials = 423;
const unsigned numPartialCiphertexts = 4096;
const unsigned relationLength = 27;
const unsigned identitysize = blocksize - 3*numofboxes;
const unsigned targetBit = 9;
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
const string pythonPath1 = "python1.txt";
const string pythonPath2 = "python2.txt";
const string invLinMatPath = "invlinmatrices.txt";
const string peelOffCipherPath = "peeledOffCiphertexts.txt";
const string peeledOffPartialCiphertextsPath = "peeledOffPartialCiphertexts.txt";
const string relationRepresentationPath = "relationRepresentation.txt";
const string relationRepresentationTargetPath ="relationRepresentationTarget.txt";
const string keysMonomialsPath = "keysMonomials.txt";

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
//  PROTOTYPES  //
//////////////////

vector<block> invertMatrix(const vector<block>& matrix);

//////////////////
//  UNUSED FUN  //
//////////////////

/*
Merge keys with same monomial.
*/
/*
void mergeMonomialsKeys(relationSetType& relationMapBit){
    relationSetType relationMapMonoKeys;
    relationSetType::iterator iter1=relationMapBit.begin();
    while(iter1 != relationMapBit.end()){
        relationRepresentation currentMonomial(*iter1);
        currentMonomial>>=keysize;
        relationRepresentation lowerBound(currentMonomial.to_ulong());
        lowerBound<<=keysize;
        relationRepresentation upperBound(currentMonomial.to_ulong()+1);
        upperBound<<=keysize;
        upperBound=upperBound.to_ulong()-1;
        //cout << *iter1 << " " << lowerBound << " " << upperBound << endl;
        relationSetType::iterator it1=relationMapBit.lower_bound(lowerBound);
        relationSetType::iterator it2=relationMapBit.upper_bound(upperBound);
        relationRepresentation tempRelaRep(lowerBound);

        //cout << "Outer: " << tempRelaRep << endl;
        while(it1!=it2){
        //    cout << "Loop" << endl;
            relationRepresentation tempMonoKeys(63); // All keybits set 111111
            tempMonoKeys&=*it1;
        //    cout << tempMonoKeys<< endl;
            tempRelaRep^=tempMonoKeys;
        ///    cout << tempRelaRep << endl;;
            ++it1;
            ++iter1;
        }
        relationMapMonoKeys.insert(tempRelaRep);
    }
    relationMapBit.clear();
    insertRemastered(relationMapBit, relationMapMonoKeys);
}
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
*/
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
block Substitution (const block message) {
    block temp = 0;

    //std::cout<<"temp Mask : "<<temp<<std::endl;
    //Get the identity part of the message
    temp ^= (message >> 3*numofboxes);

    //std::cout<<"Shifted message : "<<(message >> 3*numofboxes)<<std::endl;
    //std::cout<<"identity Part : "<<temp<<std::endl;
    //Get the rest through the Sboxes
    for (unsigned i = 1; i <= numofboxes; ++i) {
        temp <<= 3;
        //std::cout<<"Sbox Part " << i << ": "<< (message >> 3*(numofboxes-i)) <<std::endl;
        //std::cout<<"Sbox Part " << i << ": "<< ((message >> 3*(numofboxes-i))& block(0x7)).to_ulong() <<std::endl;
        temp ^= Sbox[ ((message >> 3*(numofboxes-i))
                      & block(0x7)).to_ulong()];
    }
    return temp;
}

void testSubstitution(const int val){
    block test(val);
    test <<= 6;
    cout << "Test: " << endl << test << endl;
    block substitutionTest = Substitution(test);
    cout << "Result of permutation: " << endl << substitutionTest << endl;
}
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
void preprocessingFreeCoef(vector<freeCoef>& a0, const vector<block>& peeledOffPartialCiphertexts, 
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
            if(rank_of_Matrix(tempSubspace)==subDimension){
                for(int l=0; l< a0.size(); ++l){
                    a0[l][i]=a0[l][i]^peeledOffPartialCiphertexts[k][l];
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
    unsigned int current(firstpermut);
    do {
            subspaces.push_back(current);
            current = nextPermut(current);
    }while(current <= maxpermut);
}
/*
Generate vector space 12x16.
*/
void setVectorSpace(vector<block>& base){
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
void printSequencesBlocks(const vector<block>& sequences){
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
        for (int i=0; i < rounds+1; i++){
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
Read file and set inputs in vector of blocks.
*/
void initInputs(relationSetType& textfile, string filePath){
    ifstream myFile(filePath.c_str());
    if(myFile){
        string bitLine;
        while (getline(myFile, bitLine))
        {
            relationRepresentation b(bitLine);
            textfile.insert(b);
        }
    }
    else{
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
    myFile.close();
}
/*
Write Relations map.
*/
void writeRelationMap(vector<relationSetType>& relationMap){
    ofstream myFile;
    remove(relationRepresentationPath.c_str());
    myFile.open(relationRepresentationPath.c_str());

    for(int i=0; i<relationMap.size(); ++i){
        myFile << "Bit "<< i << endl;
        int indexJ(0);
        for(relationSetType::iterator j = relationMap[i].begin(); j!=relationMap[i].end(); ++j, ++indexJ){
            myFile << "Relation Element " << indexJ << ": " << *j << endl;
        }
    }
    myFile.close();
}
/*
Write Relations map.
*/
void writeRelationMapTarget(relationSetType& relationMap){
    ofstream myFile;
    remove(relationRepresentationTargetPath.c_str());
    myFile.open(relationRepresentationTargetPath.c_str());
    int indexJ(0);
    for(relationSetType::iterator j = relationMap.begin(); j!=relationMap.end(); ++j, ++indexJ){
        myFile <<  *j << endl;
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
void writePython(vector<monomatrix>& matrixE, vector<freeCoef>& a0){
    ofstream myFile;
    myFile.open(pythonPath1.c_str());
    myFile << "[";
    for(int i=0; i< matrixE.size(); ++i){
        myFile << "[";
        for(int j=0; j< matrixE[0].size(); ++j){
            if (j==matrixE[0].size()-1){
                    myFile << matrixE[i][j];
            }else{
                 myFile << matrixE[i][j]<< " ";
            }
        }
        if(i != matrixE.size()-1) 
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
        for(int j=0; j< a0[targetBit].size(); ++j){
            if (j==a0[targetBit].size()-1){
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
Write Inputs for Sage.
*/
void writePython(vector<keyblock>& keys){
    ofstream myFile;
    myFile.open(keysMonomialsPath.c_str());
    myFile << "[";
    for(int i=0; i< keys.size(); ++i){
        myFile << "[";
        for(int j=0; j<keysize; ++j){
            if (j==keys[0].size()-1){
                    myFile << keys[i][j];
            }else{
                 myFile << keys[i][j]<< " ";
            }
        }
        if(i != keys.size()-1) 
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
Functions to multiply bitsets.
*/
void relationRepresentationMultiply(relationRepresentation& result, 
                                    const relationRepresentation& x, 
                                    const relationRepresentation& y){
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
void generateMatrixA(blockSetType& monomials, vector<block>& peeledOffCiphertexts, vector<monomatrix>& matrixA){
    for (int i = 0; i < peeledOffCiphertexts.size(); ++i){
        for (blockSetType::iterator j = monomials.begin(); j != monomials.end(); ++j){
            int indexJ=distance(monomials.begin(), j);
            block currentMonomial(*j);
            for (int k = 0; k < blocksize; ++k){
                if (!peeledOffCiphertexts[i][k] && currentMonomial[k]){
                    matrixA[i][indexJ]=0;
                    break;
                }else{    
                    matrixA[i][indexJ]=1;
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
            if(rank_of_Matrix(tempSubspace)==subDimension){
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
Generate inverted matrices of linear matrices.
*/
void generateInvMatrices(vector<vector<block>>& linearMatrices, vector<vector<block>>& invLinearMatrices){
    for(int i=0; i < linearMatrices.size(); ++i){
        invLinearMatrices.push_back(invertMatrix(linearMatrices[i]));
    }
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
vector<block> invertMatrix(const vector<block>& matrix){
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
void setInsert(relationSetType& set, relationRepresentation element){
    if(set.find(element)!=set.end()){
        set.erase(element);
    }
    else{
        set.insert(element);
    }
}
/*
Linear layer function.
*/
void linearLayerMixing(vector<relationSetType>& relationMap,
                      const vector<block>& invLinearMatrices){
    vector<relationSetType> tempRelationMap;
    relationSetType tempRelationSet;
    tempRelationMap.clear();
    tempRelationSet.clear();
    for(int h=0; h<blocksize; ++h){
        tempRelationMap.push_back(tempRelationSet);
    }
    for(int i=0; i<blocksize; ++i){
        /*if(round == 2 && i >= 9){
            break;
        }
        else{*/
            for(int j=0; j<blocksize; ++j){
                if(invLinearMatrices[i][j]){
                    for(auto element : relationMap[j]){
                        setInsert(tempRelationMap[i], element);
                    }
                }
            }
        }/*
    }*/
    relationMap.clear();
    for(int k=0; k<blocksize; ++k){
        relationMap.push_back(tempRelationMap[k]);
    }
}
/*
Adding key function.
*/
void keyRoundAdd(vector<relationSetType>& relationMap, const vector<keyblock>& keyMatrix){
    vector<relationSetType> tempRelation;
    for(int i=0; i< blocksize; ++i){
        relationRepresentation tempKey(keyMatrix[i].to_ullong());
        relationSetType tempRelationSet;
        for(auto element : relationMap[i]){
            relationRepresentation temp(element);
            temp=temp^tempKey;
            tempRelationSet.insert(temp);
        }
        tempRelation.push_back(tempRelationSet);
    }
    relationMap.clear();
    for(int k=0; k<blocksize; ++k){
        relationMap.push_back(tempRelation[k]);
    }
}
/*
Init iniput first round & key whitening.
*/
void initRelationWhitening(vector<relationSetType>& relationMap,
                        const vector<vector<keyblock>>& keyMatrices,
                        string mode){
    relationMap.clear();
    relationRepresentation tempRelation(1);
    tempRelation <<=keysize;
    relationSetType bitSet;
    for(int i=0; i<blocksize; ++i){
        bitSet.clear();
        bitSet.insert(tempRelation);
        tempRelation <<=1;
        relationMap.push_back(bitSet);
    }
    if(mode == ""){
        keyRoundAdd(relationMap, keyMatrices[0]);
    }
}/*
void insertRemastered(vector<relationRepresentation>& InsertionResult, vector<relationRepresentation>& toInsert){
    for (int i = 0; i < toInsert.size(); ++i){
        InsertionResult.push_back(toInsert[i]);
    }
}*/
void insertRemastered(relationSetType& InsertionResult, const relationSetType& toInsert){
    for (auto element : toInsert){
        setInsert(InsertionResult, element);
    }
}
/*
SBoxes function for a vector of bitset of size relationLength 
*/
void SBoxRelation(vector<relationSetType>& relationMap, string mode){
    vector<relationSetType> tempRelationMap;
    tempRelationMap.clear();
    for(int h=0; h<blocksize; ++h){
        tempRelationMap.push_back(relationMap[h]);
    }
    for(int i=numofboxes-1; i>=0; --i){      
        int x0(3*i);
        int x1(x0+1);
        int x2(x1+1);
        for(auto elementX0 : relationMap[x0]){
            for(auto elementX1 : relationMap[x1]){
                relationRepresentation tempResult(0);
                relationRepresentationMultiply(tempResult, elementX0, elementX1);
                setInsert(tempRelationMap[x2], tempResult);
            }
            for(auto elementX2 : relationMap[x2]){
                relationRepresentation tempResult(0);
                relationRepresentationMultiply(tempResult, elementX0, elementX2);
                setInsert(tempRelationMap[x1], tempResult);
            }
        }
        for(auto elementX1 : relationMap[x1]){
            for(auto elementX2 : relationMap[x2]){
                relationRepresentation tempResult(0);
                relationRepresentationMultiply(tempResult, elementX1, elementX2);
                setInsert(tempRelationMap[x0], tempResult);
            }
        }
        insertRemastered(tempRelationMap[x0], relationMap[x1]);
        insertRemastered(tempRelationMap[x0], relationMap[x2]);
        if(mode=="reverse"){
            insertRemastered(tempRelationMap[x2], relationMap[x1]);
        }
        else{
            insertRemastered(tempRelationMap[x1], relationMap[x2]);
        }
    }
    relationMap.clear();
    for(int k=0; k<blocksize; ++k){
        relationMap.push_back(tempRelationMap[k]);
    }
}
/*
Relation mapping creation.
*/
void relationMapping(vector<relationSetType>& relationMap,
                    const vector<vector<block>>& invLinearMatrices,
                    const vector<vector<keyblock>>& keyMatrices){
    initRelationWhitening(relationMap, keyMatrices, "reverse");
    /*
    for(int i=0; i<3; ++i){
    //    cout << i << endl;
        SBoxRelation(relationMap, "");
        /*cout << "Sbox"<< endl;
        for(int k =0; k < blocksize; ++k){
            cout << "Bit " << k << " : " <<  relationMap[k].size() << endl;
        }
        linearLayerMixing(relationMap, linearMatrices[i], i);

        /*cout << "Linear layer"<< endl;
        for(int l =0; l < blocksize; ++l){
            cout << "Bit " << l << " : " <<  relationMap[l].size() << endl;
        }
        keyRoundAdd(relationMap, keyMatrices[i+1]);
        /*for(int m =0; m < blocksize; ++m){
            cout << "Bit " << m << " : " <<  relationMap[m].size() << endl;
        }
    }*/
    for(int j=rounds-1; j>rounds-2; --j){
        keyRoundAdd(relationMap, keyMatrices[j+1]);
        /*cout << "Key Round"<< endl;
        for(int m =0; m < blocksize; ++m){
            //cout << "Bit " << m << " : " <<  relationMap[m].size() << endl;
            cout << "Bit " << m << " : " << endl;
            for (auto element1 : relationMap[m]){
                cout << element1 << endl;
            }
        }*/
        linearLayerMixing(relationMap, invLinearMatrices[j]);
        /*cout << "Linear layer"<< endl;
        for(int n =0; n < blocksize; ++n){
            //cout << "Bit " << n << " : " <<  relationMap[n].size() << endl;
            cout << "Bit " << n << " : " << endl;
            for (auto element1 : relationMap[n]){
                cout << element1 << endl;
            }
        }*/
        //SBoxRelation(relationMap, "reverse");
        /*cout << "Sbox"<< endl;
        for(int o =0; o < blocksize; ++o){
            //cout << "Bit " << o << " : " <<  relationMap[o].size() << endl;
            cout << "Bit " << o << " : " << endl;
            for (auto element1 : relationMap[o]){
                cout << element1 << endl;
            }
        }*/
    }
}
/*
Extract Key information according to monomials precomputed.
*/
void extractMonomialsKeys(const relationSetType& relationMapTarget, 
                            relationSetType& relationMapMonoKeys, 
                            const blockSetType& monomials){
    for(blockSetType::iterator iter1=monomials.begin(); iter1 != monomials.end(); ++iter1){
        block currentMonomial(*iter1);
        relationRepresentation lowerBound(currentMonomial.to_ulong());
        lowerBound<<=keysize;
        relationRepresentation upperBound(currentMonomial.to_ulong()+1);
        upperBound<<=keysize;
        upperBound=upperBound.to_ulong()-1;
        //cout << *iter1 << " " << lowerBound << " " << upperBound << endl;
        relationSetType::iterator it1=relationMapTarget.lower_bound(lowerBound);
        relationSetType::iterator it2=relationMapTarget.upper_bound(upperBound);
        relationRepresentation tempRelaRep(lowerBound);

        //cout << "Outer: " << tempRelaRep << endl;
        for(it1; it1!=it2; ++it1){
        //    cout << "Loop" << endl;
            relationRepresentation tempMonoKeys(63); // All keybits set 111111
            tempMonoKeys&=*it1;
        //    cout << tempMonoKeys<< endl;
            tempRelaRep^=tempMonoKeys;
        ///    cout << tempRelaRep << endl;;
        }
        relationMapMonoKeys.insert(tempRelaRep);
    }
}
/*
Setup linear equation system Keys per monomial and alpha_u per monomial.
*/
void setUpLinearEquationKeyAlphas(vector<keyblock>& keysMonomials, const relationSetType& relationMapMonoKeys){
    for(auto element : relationMapMonoKeys){
        relationRepresentation temp(element);
        temp <<= blocksize;
        temp >>= blocksize;
        keysMonomials.push_back(temp.to_ulong());
    }
    writePython(keysMonomials);
}
//////////////////
//     MAIN     //
//////////////////




int main(void) {
    vector<block> plaintexts;
    vector<block> ciphertexts;
    vector<block> partialCiphertexts;
    vector<block> peeledOffCiphertexts;
    vector<block> peeledOffPartialCiphertexts;

    vector<block> base;
    vector<vecspace> subspaces;

    blockSetType monomials;
    vector<freeCoef> a0(blocksize, 0);

    vector<monomatrix> matrixA(numPartialCiphertexts,0);
    vector<monomatrix> matrixE(numSubspaces, 0);

    vector<vector<block>> linearMatrices;
    vector<vector<block>> invLinearMatrices;
    vector<vector<keyblock>> keyMatrices;
    vector<block> roundConstants;
    
    vector<relationSetType> relationMap;
    relationSetType relationMapTarget;
    relationSetType relationMapMonoKeys;
    vector<keyblock> keysMonomials;

    //Pre-generating variables Functions
    //generateMonomials(monomials);
    setVectorSpace(base);
    setSubspaces(subspaces);

    //Initialize variables functions
    initInputs(plaintexts, plainPath);
    initInputs(ciphertexts, cipherPath);
    initInputs(partialCiphertexts, partialCipherPath);
    //initInputs(a0, freeCoefPath);
    initInputs(monomials, monomialsPath);
    //initInputs(peeledOffCiphertexts, peelOffCipherPath);
    //initInputs(peeledOffPartialCiphertexts, peeledOffPartialCiphertextsPath);
    //initInputs(relationMapTarget, relationRepresentationTargetPath);
    initInputsLinearMatrices(linearMatrices, linMatPath);
    initInputsKeyMatrices(keyMatrices, keyMatPath);
    initInputs(roundConstants, roundConstPath);
    initInputsLinearMatrices(invLinearMatrices, invLinMatPath);


    //Post-generating elements functions
    //generateInvMatrices(linearMatrices, invLinearMatrices);
    //peelingOffCiphertexts(ciphertexts, roundConstants[rounds-1], invLinearMatrices[rounds-1], peeledOffCiphertexts);
    //peelingOffCiphertexts(partialCiphertexts, roundConstants[rounds-3], invLinearMatrices[rounds-3], peeledOffPartialCiphertexts);
    //preprocessingFreeCoef(a0, peeledOffPartialCiphertexts, plaintexts, base, subspaces);
    //generateMatrixA(monomials, ciphertexts, matrixA);
    //generateMatrixE(matrixA, plaintexts, ciphertexts,subspaces, base, matrixE);
    relationMapping(relationMap, invLinearMatrices, keyMatrices);


    //Operational functions
    extractMonomialsKeys(relationMap[targetBit], relationMapMonoKeys, monomials);
    setUpLinearEquationKeyAlphas(keysMonomials, relationMapMonoKeys);
    
    //Printing Functions
    //printANF("");
    //printVectorVectorsBlock(linearMatrices);
    //printVectorVectorsKeyBlock(keyMatrices);
    //printVectorVectorsBlock(invLinearMatrices);
    //printSequencesBlocks(peeledOffCiphertexts);
    //printSequencesBlocks(peeledOffPartialCiphertexts);
    //printSequencesBlocks(monomials);
    //printSequencesMonoMatrices(matrixA);
    //printSequencesMonoMatrices(matrixE);
    //printSequencesVecspaces(subspaces);
    //printSequencesBlocks(base);
    //printSequencesBlocks(plaintexts);
    //printSequencesBlocks(ciphertexts);
    //printSequencesBlocks(partialCiphertexts);

    //Writing Functions
    //writeBlockSet(monomials, monomialsPath);
    //writeVectorsBlocks(peeledOffPartialCiphertexts, peeledOffPartialCiphertextsPath);
    //writeVectorsBlocks(peeledOffCiphertexts, peelOffCipherPath);
    //writeMatrices(invLinearMatrices, invLinMatPath);
    //writeFreeCoef(a0);
    //writePython(matrixE, a0);
    writeRelationMap(relationMap);
    writeRelationMapTarget(relationMap[targetBit]);

    

    //Testing functions

    /*keyblock tempKey(27);
    cout << "Key: " << tempKey << endl;
    relationSetType::iterator iter2=reverseRelationMapMonoKeys.begin();
    for(relationSetType::iterator iter1=relationMapMonoKeys.begin(); iter1 != relationMapMonoKeys.end(); ++iter1, ++iter2){
        relationRepresentation tempMono(0);
        relationRepresentation reverseTempMono(0);
        relationRepresentation MonoKey(*iter1);
        relationRepresentation ReverseMonoKey(*iter2);
        for(int j=0; j<keysize; ++j){
            if(MonoKey[j]){
                tempMono[0] = tempMono[0]^tempKey[j];
            }
            else if(ReverseMonoKey[j]){

                reverseTempMono[0] = reverseTempMono[0]^tempKey[j];
            }
        }
        cout << "Monokey: " << MonoKey << " ReverseMonoKey" << ReverseMonoKey << "Xor result: " << (tempMono[0] == reverseTempMono[0]) << endl;
    }

    //cout << "Previous monomials equal to new monomials set? " << (monomials == monomialsv1) << endl;
    
    //testSubstitution(3);

    /*relationSetType relationMap1;
    relationSetType relationMap2;
    relationRepresentation temp(22);
    relationRepresentation temp2(24);
    relationRepresentation temp21(1942);
    relationRepresentation temp1(0);
    relationRepresentation temp3(26);

    relationMap1.insert(temp);
    relationMap2.insert(temp);
    relationMap1.insert(temp2);

    cout << (relationMap1==relationMap2) << endl;
    relationMap1.insert(temp2);
    relationMap1.insert(temp21);
    relationMap1.insert(temp1);
    relationMap1.insert(temp3);
    for(relationSetType::iterator i = relationMap1.begin();i!=relationMap.end(); ++i){
        cout << *i << endl;
    }
    /*
    relationSetType::iterator it1=relationMap1.lower_bound(22);
    relationSetType::iterator it2=relationMap1.upper_bound(30);
    for(it1; it1!=it2; ++it1){
        cout << *it1 << endl;
    }*/

        /*cout << "Reverse Relation" <<endl;
    for(int i=0; i< relationMap.size(); ++i){
        for(auto element : relationMap[targetBit]){
            cout << element << endl;
        }
    }
    cout << "Reverse Relation" <<endl;
    for(int i=0; i< reverseRelationMap.size(); ++i){
        for(auto element : reverseRelationMap[targetBit]){
            cout << element << endl;
        }
    }
    /*
    for(auto element : reverseRelationMap[targetBit]){
        cout << element << endl;
    }

    /*cout << relationMapMonoKeys.size() << endl;

    cout << reverseRelationMapMonoKeys.size() << endl;

    cout << (relationMapMonoKeys == reverseRelationMapMonoKeys) << endl;
    /*for(auto element : relationMapMonoKeys){
        //cout << relationMapMonoKeys.size() << endl;
        cout << element << endl;
    }*/

    return 0;
}
