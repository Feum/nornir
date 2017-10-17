#include <interface.hpp>
#include <iostream>

int main(int argc, char** argv){
    int nworkers, startloop, endloop, step, chunksize;

    if (argc != 6) {
        std::cerr << "use: "
                  << argv[0]
                  << " nworkers startloop endloop step chunksize\n";
        return -1;
    }
    nworkers = atoi(argv[1]);
    startloop = atoi(argv[2]);
    endloop = atoi(argv[3]);
    step = atoi(argv[4]);
    chunksize = atoi(argv[5]);

    nornir::Parameters p("parameters.xml");
    std::vector<uint> v;
    v.resize(endloop - startloop);
    nornir::parallel_for(startloop, endloop, step, chunksize, nworkers, &p, [&v](long long int idx, uint id){
        v[idx] = idx;
    });

    for(uint x : v){
        std::cout << x << std::endl;
    }
    return 0;
}
