#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <fstream>

using namespace std;

const unsigned numofboxes = 3;    // Number of Sboxes
const unsigned blocksize = 16;   // Block size in bits
const unsigned keysize = 6; // Key size in bits
const unsigned rounds = 6; // Number of rounds
const unsigned tail = 7; // Number of bits in tail
const unsigned dimension = 12; //Dimension of vector space

const string plainPath = "../LowMC/plaintexts.txt";
const string cipherPath = "../LowMC/ciphertexts.txt";

const unsigned identitysize = blocksize - 3*numofboxes;

typedef std::bitset<blocksize> block; // Store messages and states
typedef std::bitset<keysize> keyblock;
typedef std::bitset<dimension> vecspace;




//////////////////
//   FUNCTIONS  //
//////////////////
/*
Generate combination C(12,8) subspaces
*/
void setSubspaces(std::vector<vecspace>& subspaces){

}
/*
Generate vector space 12x16
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
Print vector of sequences of bitsets.
*/
void printSequences(std::vector<block>& sequences){
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



//////////////////
//     MAIN     //
//////////////////




int main(int argc, const char * argv[]) {
    vector<block> plaintexts;
    vector<block> ciphertexts;
    vector<block> base;
    vector<vecspace> subspaces;
    //initInputs(plaintexts, plainPath);
    //initInputs(ciphertexts, cipherPath);
    setVectorSpace(base);

    
    //printSequences(base);
    //printSequences(plaintexts);
    //printSequences(ciphertexts);

    return 0;
}