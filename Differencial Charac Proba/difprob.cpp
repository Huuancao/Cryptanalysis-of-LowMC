#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#define NUMSBOXES 3
#define DATACOMPLEXITY 4
#define BLOCKSIZECONST 65536
#define TAIL 128
//////////////////
//     MAIN     //
//////////////////

using namespace std;

unsigned long long
gcd(unsigned long long x, unsigned long long y)
{
    while (y != 0)
    {
        unsigned long long t = x % y;
        x = y;
        y = t;
    }
    return x;
}

unsigned long long
choose(unsigned long long n, unsigned long long k){

    unsigned long long r(1);
    for (unsigned long long d=1; d <= k; ++d, --n)
    {
        unsigned long long g = gcd(r, d);
        r /= g;
        unsigned long long t = n / (d / g);
        r *= t;
    }
    return r;
}

double V(int a){
    return choose(NUMSBOXES,a)*pow(7,a)*TAIL;
}


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
 
vector<vector<int> > allSequences(int values, int depth) {
    vector<vector<int>> result;
    vector<int> current(depth,0);
    vector<int> limit(depth,values);
    while (current != limit) {
        int tempSum(0);
        for(int i=0; i< depth; ++i){
            tempSum+=current[i];
        }
        if(tempSum <= DATACOMPLEXITY){
            result.push_back(current);
        }
        increase(current,values);
    }

    return result;
}

void printSequences(vector<vector<int>> sequences){
    for (int i = 0; i < sequences.size(); ++i) {
        for (int j = 0; j < sequences[i].size(); ++j) {
            cout << sequences[i][j] << " ";
        }
        cout << endl;
    }
}
 
double computeDifferentialCharacProba(vector<vector<int>> sequences, int round){
    double tempProba(0);

    for (int i = 0; i < sequences.size(); ++i) {
        double tempProd(1);
        double tempSum(0);
        for (int j = 0; j < sequences[i].size(); ++j) {
            tempSum+=sequences[i][j];
            tempProd*=V(sequences[i][j]);
            /*if(j== sequences[i].size()-1){
                cout << "Sum alphas = " << tempSum << endl;
                cout << "Produit V = " << tempProd << endl;
            }*/
        }
        tempProba+=tempProd*pow(4,tempSum)/pow(BLOCKSIZECONST-1,round-1);
        //cout << "Diffenrential Characteristic Proba = " << tempProba << endl;
    }
    return tempProba;
}
 
int main(int argc, const char * argv[]) {
    string answer("No");
    double threshold = pow(2, -DATACOMPLEXITY);


    int round(6);
    double difProba(0);
    vector<vector<int> >sequences = allSequences(NUMSBOXES, round);
    //printSequences(sequences);
    difProba=computeDifferentialCharacProba(sequences, round);
    if(difProba<threshold){
        answer = "Yes";
    }

    cout << "Diffenrential Characteristic of Probability for " << round << " rounds = " << difProba << endl;
    cout << "Is " << round << " rounds enough? " << answer << endl;

   
    return 0;
}
