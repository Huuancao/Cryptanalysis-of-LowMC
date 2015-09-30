#include "LowMC.h"
#include <iostream>
#include <vector>
#include <bitset>

#define MAX 512
//////////////////
//     MAIN     //
//////////////////

void generateInputs(std::vector<block>& inputs){
    block input(0);
    for(int i=0; i < MAX; ++i){
        input=i;
        inputs.push_back(input << blocksize - tail);

    }
    /*for(int j=0; j< inputs.size(); ++j){
        std::cout << "Input n" << j <<": " << inputs[j]<< std::endl;
    }*/
}

int main () {
    LowMC cipher(0);
    std::vector<block> Inputs;
    generateInputs(Inputs);

    block a = 0xABCD;

    std::cout << "Plaintext:" << std::endl;
    std::cout << a << std::endl;
    a = cipher.encrypt( a );
    std::cout << "Ciphertext:" << std::endl;
    std::cout << a << std::endl;
    a = cipher.decrypt( a );
    std::cout << "Encryption followed by decryption of plaintext:" << std::endl;
    std::cout << a << std::endl;

    return 0;
}
