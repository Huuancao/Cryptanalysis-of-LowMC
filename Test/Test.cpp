#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <bitset>

const unsigned maxbits = 512;
const unsigned tail = 7;
const unsigned blocksize = 16;

typedef std::bitset<blocksize> block;



using namespace std;

/*void increase(vector<int> &current, int values) {
    current.back() += 1;
    long index = current.size() - 1;
    while (current[index] > values && index >= 0) {
        current[index] = 0;
        index -= 1;
        if (index >= 0) {
            current[index] += 1;
        }
    }
}
 
vector<block> allSequences(int values, int depth) {
    vector<block> result;
    vector<int> current(depth,0);
    vector<int> limit(depth,values);
    while (current != limit) {
        result.push_back(current);
        increase(current,values);
    }

    return result;

}*/
unsigned int reverseBits(unsigned int num){;
    unsigned int reverseNum = 0;
    int i;
    for (i = 0; i < blocksize; i++){
            if((num & (1 << i)))
               reverseNum |= 1 << ((blocksize - 1) - i);  
    }
    return reverseNum;
}
void generateInputs(std::vector<block>& inputs){
    block input(0);
    for(int i=0; i < maxbits; ++i){
        input=reverseBits(i);
        inputs.push_back(input);
    }
    for(int j=0; j< inputs.size(); ++j){
        std::cout << "Input n" << j <<": " << inputs[j]<< std::endl;
    }
}




int main(int argc, const char * argv[]) {
    std::vector<block> Inputs;
    generateInputs(Inputs);




}