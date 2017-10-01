#include <interface.hpp>
#include <iostream>

#define VECTOR_LENGTH 100

int main(int argc, char** argv){
    nornir::Parameters p;
    std::vector<uint> v;
    v.resize(VECTOR_LENGTH);
    nornir::parallel_for(0, VECTOR_LENGTH, 1, 2, 4, &p, [&v](long long int idx, uint id){
        v[idx] = idx;
    });

    for(uint x : v){
        std::cout << x << std::endl;
    }
    return 0;
}