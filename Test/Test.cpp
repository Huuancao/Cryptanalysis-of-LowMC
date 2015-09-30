#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <bitset>




using namespace std;
/*
void increase(vector<int> &current, int values) {
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
        int tempSum(0);
        for(int i=0; i< depth; ++i){
            tempSum+=current[i];
        }
        result.push_back(current);
        increase(current,values);
    }

    return result;
}
*/
int main(int argc, const char * argv[]) {


}