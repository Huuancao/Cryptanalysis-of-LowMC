#include "LowMC.h"
#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <fstream>
#include <cmath>


using namespace std;
//////////////////
//     MAIN     //
//////////////////

/*
Reverse the bits order in identity part of inputs.
*/
unsigned int reverseBits(unsigned int num){;
    unsigned int reverseNum = 0;
    for (int i = 0; i < blocksize; i++){
            if((num & (1 << i)))
               reverseNum |= 1 << ((blocksize - 1) - i);  
    }
    return reverseNum;
}

/*
No brain generator of bits sequences in reverse order. Shift left 9 to feed identity part.
*/
void generatePlaintexts(std::vector<block>& plaintexts,const string mode, const int maxPlaintexts){

    for(int i=0; i < maxPlaintexts; ++i){
        if(mode == "reverse"){      
            plaintexts.push_back(reverseBits(i));
        }
        else{
            block input(i);
            plaintexts.push_back(input<<=9);
        }
        
        
    }
    /*for(int j=0; j< plaintexts.size(); ++j){
        std::cout << "Input n" << j <<": " << plaintexts[j]<< std::endl;
    }*/
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
Generate ciphertexts.
*/
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


int main () {
    LowMC cipher(0); //Set key to 101101
    std::vector<block> plaintexts;
    std::vector<block> ciphertexts;
    int maxPlaintexts(pow(2,12));
    string mode("no"); //Type "reverse" to reverse inputs, else type anything different
    generatePlaintexts(plaintexts, mode, maxPlaintexts);
    writeVectorsBlocks(plaintexts, "plaintexts.txt");

    //printSequences(plaintexts);
    generateCiphertexts(plaintexts, ciphertexts, cipher);
    printSequences(ciphertexts);
    


    writeVectorsBlocks(ciphertexts, "ciphertexts.txt");
    //writeVectorsBlocks(ciphertexts, "partialCiphertexts.txt");

    /*
    block a = 0xABCD;

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
