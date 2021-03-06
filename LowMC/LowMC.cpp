#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>

#include "LowMC.h"


/////////////////////////////
//     LowMC functions     //
/////////////////////////////

block LowMC::encrypt (const block message) {
    block c = message ^ roundkeys[0];
    for (unsigned r = 1; r <= rounds; ++r) {
        c =  Substitution(c);
        c =  MultiplyWithGF2Matrix(LinMatrices[r-1], c);
        c ^= roundconstants[r-1];
        c ^= roundkeys[r];
    }
    return c;
}


block LowMC::decrypt (const block message) {
    block c = message;
    for (unsigned r = rounds; r > 0; --r) {
        c ^= roundkeys[r];
        c ^= roundconstants[r-1];
        c =  MultiplyWithGF2Matrix(invLinMatrices[r-1], c);
        c =  invSubstitution(c);
    }
    c ^= roundkeys[0];
    return c;
}


void LowMC::set_key (keyblock k) {
    key = k;

    std::cout << "Key: " << key << std::endl;
    keyschedule();
}

//////////////////////////////
// Custom private functions //
//////////////////////////////

/*
Set last 7 lines to identity matrix' lines in order to "eat" a layer in degree.
*/
void setLastPartMatrix( std::vector<block>& mat){
    block tempMatrixLine;
    tempMatrixLine.reset();
    tempMatrixLine[9]=1;
    for(int i=0; i<12; ++i){
        mat.push_back(tempMatrixLine);       
        tempMatrixLine <<= 1;
    }
}

/*
Write data functions.
*/
void writeMatrices(std::vector<std::vector<block>> matrix, std::string fileName){
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
void writeMatrices(std::vector<std::vector<keyblock>> matrix, std::string fileName){
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
void writeConstants(std::vector<block>& roundconstants, std::string fileName){
    std::ofstream myFile;
    myFile.open(fileName.c_str());

    for(int i=0; i<roundconstants.size();++i){
        myFile << roundconstants[i] << std::endl;
    }
    myFile.close();
}
/*
Write Roundkeys.
*/
void writeRoundKeys(std::vector<block> roundkeys){
    std::ofstream myFile;
    myFile.open("roundkeys.txt");

    for(int i=0; i<roundkeys.size();++i){
            myFile << roundkeys[i] << std::endl;
    }
    myFile.close();
}

/////////////////////////////
// LowMC private functions //
/////////////////////////////

block LowMC::Substitution (const block message) {
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

block LowMC::invSubstitution (const block message) {
    block temp = 0;
    //Get the identity part of the message
    temp ^= (message >> 3*numofboxes);
    //Get the rest through the invSboxes
    for (unsigned i = 1; i <= numofboxes; ++i) {
        temp <<= 3;
        temp ^= invSbox[ ((message >> 3*(numofboxes-i))
                         & block(0x7)).to_ulong()];
    }
    return temp;
}



block LowMC::MultiplyWithGF2Matrix
        (const std::vector<block> matrix, const block message) {
    block temp = 0;
    for (unsigned i = 0; i < blocksize; ++i) {
        temp[i] = (message & matrix[i]).count() % 2;
    }
    return temp;
}


block LowMC::MultiplyWithGF2Matrix_Key
        (const std::vector<keyblock> matrix, const keyblock k) {
    block temp = 0;
    for (unsigned i = 0; i < blocksize; ++i) {
        temp[i] = (k & matrix[i]).count() % 2;
    }
    return temp;
}

void LowMC::keyschedule () {
    roundkeys.clear();
    for (unsigned r = 0; r <= rounds; ++r) {
        roundkeys.push_back( MultiplyWithGF2Matrix_Key (KeyMatrices[r], key) );
    }
    writeRoundKeys(roundkeys);
    return;
}



void LowMC::instantiate_LowMC () {
    // Create LinMatrices and invLinMatrices
    LinMatrices.clear();
    invLinMatrices.clear();
    for (unsigned r = 0; r < rounds; ++r) {
        // Create matrix
        std::vector<block> mat;
        // Fill matrix with random bits
        do {
            mat.clear();
            for (unsigned i = 0; i < blocksize; ++i) {
                if(r==2 && i == 9){
                    setLastPartMatrix(mat);
                    break;
                }
                mat.push_back( getrandblock () );
            }
        // Repeat if matrix is not invertible
        } while ( rank_of_Matrix(mat) != blocksize );
        LinMatrices.push_back(mat);
        invLinMatrices.push_back(invert_Matrix (LinMatrices.back()));
    }

    writeMatrices(LinMatrices, "linmatrices.txt");

    //printMatrix(LinMatrices);

    // Create roundconstants
    roundconstants.clear();
    for (unsigned r = 0; r < rounds; ++r) {
        roundconstants.push_back( getrandblock () );
    }

    writeConstants(roundconstants, "roundconstants.txt");

    // Create KeyMatrices
    KeyMatrices.clear();
    for (unsigned r = 0; r <= rounds; ++r) {
        // Create matrix
        std::vector<keyblock> mat;
        // Fill matrix with random bits
        do {
            mat.clear();
            for (unsigned i = 0; i < blocksize; ++i) {
                mat.push_back( getrandkeyblock () );
            }
        // Repeat if matrix is not of maximal rank
        } while ( rank_of_Matrix_Key(mat) < std::min(blocksize, keysize) );
        KeyMatrices.push_back(mat);
    }
    writeMatrices(KeyMatrices, "keymatrices.txt");
    //printMatrix(KeyMatrices);
    
    return;
}


/////////////////////////////
// Binary matrix functions //
/////////////////////////////


unsigned LowMC::rank_of_Matrix (const std::vector<block> matrix) {
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


unsigned LowMC::rank_of_Matrix_Key (const std::vector<keyblock> matrix) {
    std::vector<keyblock> mat; //Copy of the matrix 
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


std::vector<block> LowMC::invert_Matrix (const std::vector<block> matrix) {
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

///////////////////////
// Pseudorandom bits //
///////////////////////


block LowMC::getrandblock () {
    block tmp = 0;
    for (unsigned i = 0; i < blocksize; ++i) tmp[i] = getrandbit ();
    return tmp;
}

keyblock LowMC::getrandkeyblock () {
    keyblock tmp = 0;
    for (unsigned i = 0; i < keysize; ++i) tmp[i] = getrandbit ();
    return tmp;
}


// Uses the Grain LSFR as self-shrinking generator to create pseudorandom bits
// Is initialized with the all 1s state
// The first 160 bits are thrown away
bool LowMC::getrandbit () {
    static std::bitset<80> state; //Keeps the 80 bit LSFR state
    bool tmp = 0;
    //If state has not been initialized yet
    if (state.none ()) {
        state.set (); //Initialize with all bits set
        //Throw the first 160 bits away
        for (unsigned i = 0; i < 160; ++i) {
            //Update the state
            tmp =  state[0] ^ state[13] ^ state[23]
                       ^ state[38] ^ state[51] ^ state[62];
            state >>= 1;
            state[79] = tmp;
        }
    }
    //choice records whether the first bit is 1 or 0.
    //The second bit is produced if the first bit is 1.
    bool choice = false;
    do {
        //Update the state
        tmp =  state[0] ^ state[13] ^ state[23]
                   ^ state[38] ^ state[51] ^ state[62];
        state >>= 1;
        state[79] = tmp;
        choice = tmp;
        tmp =  state[0] ^ state[13] ^ state[23]
                   ^ state[38] ^ state[51] ^ state[62];
        state >>= 1;
        state[79] = tmp;
    } while (! choice);
    return tmp;
}



