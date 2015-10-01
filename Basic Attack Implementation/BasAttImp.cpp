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

const string plainPath = "../LowMC/plaintexts.txt";
const string cipherPath = "../LowMC/ciphertexts.txt";

const unsigned identitysize = blocksize - 3*numofboxes;

typedef std::bitset<blocksize> block; // Store messages and states
typedef std::bitset<keysize> keyblock;




//////////////////
//   FUNCTIONS  //
//////////////////
void printSequences(std::vector<block>& sequences){
    for(int i=0; i< sequences.size(); ++i){
        cout << "Entry n" <<i << ": " << sequences[i] << endl;
    }
}

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
    initInputs(plaintexts, plainPath);
    initInputs(ciphertexts, cipherPath);
    //printSequences(plaintexts);
    //printSequences(ciphertexts);

    return 0;
}