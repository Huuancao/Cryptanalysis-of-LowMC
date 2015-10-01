#include "LowMC.h"
#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <fstream>

#define MAX 128
using namespace std;
//////////////////
//     MAIN     //
//////////////////
unsigned int reverseBits(unsigned int num){;
    unsigned int reverseNum = 0;
    for (int i = 0; i < blocksize; i++){
            if((num & (1 << i)))
               reverseNum |= 1 << ((blocksize - 1) - i);  
    }
    return reverseNum;
}

void generatePlaintexts(std::vector<block>& plaintexts, string mode){

    for(int i=0; i < MAX; ++i){
        if(mode == "reverse"){      
            plaintexts.push_back(reverseBits(i));
        }
        else{
            block input(i);
            plaintexts.push_back(input << blocksize - tail);
        }
        
        
    }
    /*for(int j=0; j< plaintexts.size(); ++j){
        std::cout << "Input n" << j <<": " << plaintexts[j]<< std::endl;
    }*/
}

void printSequences(std::vector<block>& sequences){
    for(int i=0; i< sequences.size(); ++i){
        cout << "Entry n" <<i << ": " << sequences[i] << endl;
    }
}
void generateCiphertexts(const vector<block> plaintexts, vector<block>& ciphertexts, LowMC cipher){
    for(int i=0; i< plaintexts.size(); ++i){
        block a = cipher.encrypt(plaintexts[i]);
        ciphertexts.push_back(a);
        /*
        std::cout << "Plaintext:" << std::endl;
        std::cout << plaintexts[i] << std::endl;
        std::cout << "Ciphertext:" << std::endl;
        std::cout << a << std::endl;
        a = cipher.decrypt( a );
        std::cout << "Encryption followed by decryption of plaintext:" << std::endl;
        std::cout << a << std::endl;
        */
    }
}

void writePlaintexts(const vector<block>& plaintexts){
    ofstream myFile;
    myFile.open("plaintexts.txt");
    for(int i=0; i< plaintexts.size(); ++i){
        myFile << plaintexts[i] << endl;
    }
    myFile.close();
}

void writeCiphertexts(const vector<block>& ciphertexts){
    ofstream myFile;
    myFile.open("ciphertexts.txt");
    for(int i=0; i< ciphertexts.size(); ++i){
        myFile << ciphertexts[i] << endl;
    }
    myFile.close();
}

int main () {
    LowMC cipher(0);
    std::vector<block> plaintexts;
    std::vector<block> ciphertexts;
    string mode("reverse");
    generatePlaintexts(plaintexts, mode);
    writePlaintexts(plaintexts);

    //printSequences(plaintexts);
    generateCiphertexts(plaintexts, ciphertexts, cipher);
    writeCiphertexts(ciphertexts);



    /*block a = 0xABCD;

    std::cout << "Plaintext:" << std::endl;
    std::cout << a << std::endl;
    a = cipher.encrypt( a );
    std::cout << "Ciphertext:" << std::endl;
    std::cout << a << std::endl;
    a = cipher.decrypt( a );
    std::cout << "Encryption followed by decryption of plaintext:" << std::endl;
    std::cout << a << std::endl;
    */
    return 0;
}
