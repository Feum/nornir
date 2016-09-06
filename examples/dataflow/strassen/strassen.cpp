/*
 * Dataflow implementaion of Strassen algorithm for matrix multiplication.
 *
 *  Created on: Sep 6, 2016
 *  Author: Daniele De Sensi (d.desensi.software@gmail.com)
 *          Massimo Torquati (torquati@di.unipi.it)
 */

/*
 * graph.cpp
 *
 * Created on: 30/04/2016
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of nornir.
 *
 *  nornir is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  nornir is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with nornir.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

#include <iostream>
#include <stdlib.h>
#include "../../../src/dataflow/interpreter.hpp"

using namespace nornir::dataflow;

const double THRESHOLD = 0.001;

long CheckResults(long m, long n, const double *C, const double *C1)
{
    for (long i=0; i<m; i++)
        for (long j=0; j<n; j++) {
            long idx = i*n+j;
            if (fabs(C[idx] - C1[idx]) > THRESHOLD) {
                printf("ERROR %ld,%ld %f != %f\n", i, j, C[idx], C1[idx]);
                return 1;
            }
        }
    printf("OK.\n");
    return 0;
}

void random_init (long M, long N, long P, double *A, double *B) {
    for (long i = 0; i < M; i++) {
    for (long j = 0; j < P; j++) {
        A[i*P+j] = 5.0 - ((double)(rand()%100) / 10.0);
    }
    }

    for (long i = 0; i < P; i++) {
    for (long j = 0; j < N; j++) {
        B[i*N+j] = 5.0 - ((double)(rand()%100) / 10.0);
    }
    }
}

void printarray(const double *A, long m, long n, long N) {
    for (long i=0; i<m; i++) {
        for (long j=0; j<n; j++)
            printf("%f\t", A[i*N+j]);
        printf("\n");
    }
}

inline void seqMatMult(long m, long n, long p,
                       const double* A, const long AN,
                       const double* B, const long BN,
                       double* C, const long CN)  {
    for (long i = 0; i < m; i++)
        for (long j = 0; j < n; j++) {
            C[i*CN+j] = 0.0;
            for (long k = 0; k < p; k++)
                C[i*CN+j] += A[i*AN+k]*B[k*BN+j];
        }
}


// m by n with row stride XN for X YN for Y and CN for C
inline void mmsum(const double *X, long XN, const double *Y, long YN,
                   double *C, long CN,  long m, long n) {

    for(long i = 0; i < m; i++){
        for (long j=0; j<n; j++){
            C[i*CN+j] = X[i*XN+j] + Y[i*YN + j];
        }
    };
}

// m by n with row stride XN for X YN for Y and CN for C
inline void mmsub (const double *X, long XN,  const double *Y, long YN,
                   double *C, long CN,  long m, long n) {
    for(long i = 0; i < m; i++){
        for (long j=0; j<n; j++){
            C[i*CN+j] = X[i*XN+j] - Y[i*YN + j];
        }
    };
}

typedef struct{
    const double* A11;
    const double* A12;
    const double* A21;
    const double* A22;

    const double* B11;
    const double* B12;
    const double* B21;
    const double* B22;

    double* C;

    double* P1;
    double* P2;
    double* P3;
    double* P4;
    double* P5;
    double* P6;
    double* P7;
}Task;

class MatrixStreamRate: public nornir::dataflow::InputStreamRate{
private:
    long _m, _n, _p, _AN, _BN, _matrices;
public:
    inline MatrixStreamRate(std::string fileName, long m, long n, long p,
                            long AN, long BN,
                            long matrices = 1000):InputStreamRate(fileName),
                            _m(m), _n(n), _p(p), _AN(AN), _BN(BN), _matrices(matrices){
        ;
    }

protected:
    std::vector<void*> loadObjects(){
        std::vector<void*>  objects;
        for(size_t i = 0; i < (size_t) _matrices; i++){
            Task* t = new Task();
            memset(t, 0, sizeof(Task));
            const double *A = (double*)malloc(_m*_p*sizeof(double));
            const double *B = (double*)malloc(_p*_n*sizeof(double));
            assert(A); assert(B);

            random_init(_m, _n, _p, const_cast<double*>(A), const_cast<double*>(B));
            double *C = (double*)malloc(_m*_n*sizeof(double));

            t->A11 = &A[0];
            t->A12 = &A[_p/2];
            t->A21 = &A[_m/2*_AN];
            t->A22 = &A[_m/2*_AN+_p/2];

            t->B11 = &B[0];
            t->B12 = &B[_n/2];
            t->B21 = &B[_p/2*_BN];
            t->B22 = &B[_p/2*_BN+_n/2];

            t->C = C;
        }
        return objects;
    }
};

class DemoOutputStream: public nornir::dataflow::OutputStream{
private:
    long _m, _n, _p;
    bool _check;
public:
    DemoOutputStream(long m, long n, long p, bool check = false):_m(m), _n(n), _p(p), _check(check){;}
    void put(void* a){
        Task* t = (Task*) a;
        if(_check){
            double *C2 = (double*)malloc(_m*_n*sizeof(double));
            seqMatMult(_m, _n, _p, t->A11, _p, t->B11, _n, C2, _n);
            CheckResults(_m, _n, t->C, C2);
            free(C2);
        }
        delete t->A11;
        delete t->B11;
        delete t->C;
        delete t->P1;
        delete t->P2;
        delete t->P3;
        delete t->P4;
        delete t->P5;
        delete t->P6;
        delete t->P7;
        delete t;
    }
};

Computable *IN, *P1Ins, *P2Ins, *P3Ins, *P4Ins, *P5Ins, *P6Ins, *P7Ins, *C11Ins, *C22Ins, *C12Ins, *C21Ins, *OUT;

class INCode: public Computable{
public:
    void compute(Data* d){
        Task* t = (Task*) d->getInput();

        d->setOutput((void*) t, P1Ins);
        d->setOutput((void*) t, P2Ins);
        d->setOutput((void*) t, P3Ins);
        d->setOutput((void*) t, P4Ins);
        d->setOutput((void*) t, P5Ins);
        d->setOutput((void*) t, P6Ins);
        d->setOutput((void*) t, P7Ins);
    }
};

class P1Code: public Computable{
private:
    long _m2, _n2, _p2;
    long _AN, _BN, _CN;
public:
    P1Code(long m2, long n2, long p2,
           long AN, long BN, long CN):
               _m2(m2), _n2(n2), _p2(p2),
               _AN(AN), _BN(BN), _CN(CN){;}

    void compute(Data* d){
        Task* t = (Task*) d->getInput();

        double *sumA= (double*)malloc(_m2*_p2*sizeof(double));
        double *sumB= (double*)malloc(_p2*_n2*sizeof(double));
        double *P1  = (double*)malloc(_m2*_n2*sizeof(double));

        mmsum(t->A11, _AN, t->A22, _AN, sumA, _p2, _m2, _p2); // S1
        mmsum(t->B11, _BN, t->B22, _BN, sumB, _n2, _p2, _n2); // S2
        seqMatMult(_m2, _n2, _p2, sumA, _p2, sumB, _n2, P1, _n2);  // P1

        free(sumA);free(sumB);

        t->P1 = P1;

        d->setOutput((void*) t, C11Ins);
        d->setOutput((void*) t, C22Ins);
    }
};

class P2Code: public Computable{
private:
    long _m2, _n2, _p2;
    long _AN, _BN, _CN;
public:
    P2Code(long m2, long n2, long p2,
           long AN, long BN, long CN):
               _m2(m2), _n2(n2), _p2(p2),
               _AN(AN), _BN(BN), _CN(CN){;}

    void compute(Data* d){
        Task* t = (Task*) d->getInput();

        double *P2 = (double*)malloc(_m2*_n2*sizeof(double));
        double *sumA = (double*)malloc(_m2*_p2*sizeof(double));

        mmsum(t->A21, _AN, t->A22, _AN, sumA, _p2, _m2, _p2); // S3
        seqMatMult(_m2, _n2, _p2, sumA, _p2, t->B11, _BN, P2, _CN); // P2

        free(sumA);

        t->P2 = P2;

        d->setOutput((void*) t, C22Ins);
        d->setOutput((void*) t, C21Ins);
    }
};


class P3Code: public Computable{
private:
    long _m2, _n2, _p2;
    long _AN, _BN, _CN;
public:
    P3Code(long m2, long n2, long p2,
           long AN, long BN, long CN):
               _m2(m2), _n2(n2), _p2(p2),
               _AN(AN), _BN(BN), _CN(CN){;}

    void compute(Data* d){
        Task* t = (Task*) d->getInput();

        double *P3 = (double*)malloc(_m2*_n2*sizeof(double));
        double *sumB= (double*)malloc(_p2*_n2*sizeof(double));

        mmsub(t->B12, _BN, t->B22, _BN, sumB, _n2, _p2, _n2); // S4
        seqMatMult(_m2, _n2, _p2, t->A11, _AN, sumB, _n2, P3, _CN); // P3

        free(sumB);

        t->P3 = P3;

        d->setOutput((void*) t, C22Ins);
        d->setOutput((void*) t, C12Ins);
    }
};

class P4Code: public Computable{
private:
    long _m2, _n2, _p2;
    long _AN, _BN, _CN;
public:
    P4Code(long m2, long n2, long p2,
           long AN, long BN, long CN):
               _m2(m2), _n2(n2), _p2(p2),
               _AN(AN), _BN(BN), _CN(CN){;}

    void compute(Data* d){
        Task* t = (Task*) d->getInput();

        double *P4 = (double*)malloc(_m2*_n2*sizeof(double));
        double *sumB= (double*)malloc(_p2*_n2*sizeof(double));

        mmsub(t->B21, _BN, t->B11, _BN, sumB, _n2, _p2, _n2); // S5
        seqMatMult(_m2, _n2, _p2, t->A22, _AN, sumB, _n2, P4, _n2); // P4

        free(sumB);

        t->P4 = P4;

        d->setOutput((void*) t, C11Ins);
        d->setOutput((void*) t, C21Ins);
    }
};

class P5Code: public Computable{
private:
    long _m2, _n2, _p2;
    long _AN, _BN, _CN;
public:
    P5Code(long m2, long n2, long p2,
           long AN, long BN, long CN):
               _m2(m2), _n2(n2), _p2(p2),
               _AN(AN), _BN(BN), _CN(CN){;}

    void compute(Data* d){
        Task* t = (Task*) d->getInput();

        double *P5 = (double*)malloc(_m2*_n2*sizeof(double));
        double *sumA= (double*)malloc(_m2*_p2*sizeof(double));

        mmsum(t->A11, _AN, t->A12, _AN, sumA, _p2, _m2, _p2); // S6
        seqMatMult(_m2, _n2, _p2, sumA, _p2, t->B22, _BN, P5, _n2); // P5

        free(sumA);

        t->P5 = P5;

        d->setOutput((void*) t, C11Ins);
        d->setOutput((void*) t, C12Ins);
    }
};

class P6Code: public Computable{
private:
    long _m2, _n2, _p2;
    long _AN, _BN, _CN;
public:
    P6Code(long m2, long n2, long p2,
           long AN, long BN, long CN):
               _m2(m2), _n2(n2), _p2(p2),
               _AN(AN), _BN(BN), _CN(CN){;}

    void compute(Data* d){
        Task* t = (Task*) d->getInput();

        double *P6 = (double*)malloc(_m2*_n2*sizeof(double));
        double *sumA= (double*)malloc(_m2*_p2*sizeof(double));
        double *sumB= (double*)malloc(_p2*_n2*sizeof(double));

        mmsub(t->A21, _AN, t->A11, _AN, sumA, _p2, _m2, _p2); // S7
        mmsum(t->B11, _BN, t->B12, _BN, sumB, _n2, _p2, _n2); // S8
        seqMatMult(_m2, _n2, _p2, sumA, _p2, sumB, _n2, P6, _CN);// P6

        free(sumA);free(sumB);

        t->P6 = P6;

        d->setOutput((void*) t, C22Ins);
    }
};

class P7Code: public Computable{
private:
    long _m2, _n2, _p2;
    long _AN, _BN, _CN;
public:
    P7Code(long m2, long n2, long p2,
           long AN, long BN, long CN):
               _m2(m2), _n2(n2), _p2(p2),
               _AN(AN), _BN(BN), _CN(CN){;}

    void compute(Data* d){
        Task* t = (Task*) d->getInput();

        double *P7 = (double*)malloc(_m2*_n2*sizeof(double));
        double *sumA= (double*)malloc(_m2*_p2*sizeof(double));
        double *sumB= (double*)malloc(_p2*_n2*sizeof(double));

        mmsub(t->A12, _AN, t->A22, _AN, sumA, _p2, _m2, _p2); // S9
        mmsum(t->B21, _BN, t->B22, _BN, sumB, _n2, _p2, _n2); // S10
        seqMatMult(_m2, _n2, _p2, sumA, _p2, sumB, _n2, P7, _CN);// P7

        free(sumA);free(sumB);

        t->P7 = P7;

        d->setOutput((void*) t, C11Ins);
    }
};

class C11Code: public Computable{
private:
    long _m2, _n2, _CN;
public:
    C11Code(long m2, long n2, long CN):_m2(m2), _n2(n2), _CN(CN){;}

    void compute(Data* d){
        Task* P1t = (Task*) d->getInput(P1Ins);
        Task* P4t = (Task*) d->getInput(P4Ins);
        Task* P5t = (Task*) d->getInput(P5Ins);
        Task* P7t = (Task*) d->getInput(P7Ins);
        double *C11 = &(P1t->C[0]);
        for(long i = 0; i < _m2; i++){
            for(long j = 0; j< _n2; ++j) {
                C11[i*_CN + j] = P1t->P1[i*_n2 + j] + P4t->P4[i*_n2 + j] - P5t->P5[i*_n2 + j] + P7t->P7[i*_n2 + j];
            }
        }
        P1t->P4 = P4t->P4;
        P1t->P5 = P5t->P5;
        P1t->P7 = P7t->P7;
        d->setOutput((void*) P1t, OUT);
    }
};

class C12Code: public Computable{
private:
    long _m2, _n2, _CN;
public:
    C12Code(long m2, long n2, long CN):_m2(m2), _n2(n2), _CN(CN){;}

    void compute(Data* d){
        Task* P3t = (Task*) d->getInput(P3Ins);
        Task* P5t = (Task*) d->getInput(P5Ins);
        double *C12 = &(P3t->C[_n2]);
        for(long i = 0; i < _m2; i++){
            for(long j = 0; j< _n2; ++j) {
                C12[i*_CN + j] = P3t->P3[i*_n2 + j] + P5t->P5[i*_n2 + j];
            }
        }
        P3t->P5 = P5t->P5;
        d->setOutput((void*) P3t, OUT);
    }
};

class C21Code: public Computable{
private:
    long _m2, _n2, _CN;
public:
    C21Code(long m2, long n2, long CN):_m2(m2), _n2(n2), _CN(CN){;}

    void compute(Data* d){
        Task* P2t = (Task*) d->getInput(P2Ins);
        Task* P4t = (Task*) d->getInput(P4Ins);
        double *C21 = &(P2t->C[_m2*_CN]);
        for(long i = 0; i < _m2; i++){
            for(long j = 0; j< _n2; ++j) {
                C21[i*_CN + j] = P2t->P2[i*_n2 + j] + P4t->P4[i*_n2 + j];
            }
        }
        P2t->P4 = P4t->P4;
        d->setOutput((void*) P2t, OUT);
    }
};

class C22Code: public Computable{
private:
    long _m2, _n2, _CN;
public:
    C22Code(long m2, long n2, long CN):_m2(m2), _n2(n2), _CN(CN){;}

    void compute(Data* d){
        Task* P1t = (Task*) d->getInput(P1Ins);
        Task* P2t = (Task*) d->getInput(P2Ins);
        Task* P3t = (Task*) d->getInput(P3Ins);
        Task* P6t = (Task*) d->getInput(P6Ins);
        double *C22 = &(P1t->C[_m2*_CN+_n2]);
        for(long i = 0; i < _m2; i++){
            for(long j = 0; j< _n2; ++j) {
                C22[i*_CN + j] = P1t->P1[i*_n2 + j] - P2t->P2[i*_n2 + j] + P3t->P3[i*_n2 + j] + P6t->P6[i*_n2 + j];
            }
        }
        P1t->P2 = P2t->P2;
        P1t->P3 = P3t->P3;
        P1t->P6 = P6t->P6;
        d->setOutput((void*) P1t, OUT);
    }
};

class OUTCode: public Computable{
public:
    void compute(Data* d){
        Task* C11t = (Task*) d->getInput(C11Ins);
        Task* C22t = (Task*) d->getInput(C22Ins);

        // P1 already in C11t
        C11t->P2 = C22t->P2;
        C11t->P3 = C22t->P3;
        // P4 already in C11t
        // P5 already in C11t
        C11t->P6 = C22t->P6;
        // P7 already in C11t
        d->setOutput((void*) C11t);
    }
};

int main(int argc, char** argv){
    if(argc < 5){
        std::cerr << "Usage: " << argv[0] << " ratefile m n p" << std::endl;
        return -1;
    }

    long m=atol(argv[2]),n=atol(argv[3]),p=atol(argv[4]);
    long AN=p,BN=n,CN=n;
    long m2 = m/2, n2 = n/2, p2 = p/2;

    /* Create streams. */
    MatrixStreamRate inp(argv[1], m, n, p, AN, BN);
    DemoOutputStream out(m, n, p, true);

    /* Create instructions. */
    IN = new INCode();
    P1Ins = new P1Code(m2, n2, p2, AN, BN, CN);
    P2Ins = new P2Code(m2, n2, p2, AN, BN, CN);
    P3Ins = new P3Code(m2, n2, p2, AN, BN, CN);
    P4Ins = new P4Code(m2, n2, p2, AN, BN, CN);
    P5Ins = new P5Code(m2, n2, p2, AN, BN, CN);
    P6Ins = new P6Code(m2, n2, p2, AN, BN, CN);
    P7Ins = new P7Code(m2, n2, p2, AN, BN, CN);
    C11Ins = new C11Code(m2, n2, CN);
    C22Ins = new C22Code(m2, n2, CN);
    C12Ins = new C12Code(m2, n2, CN);
    C21Ins = new C21Code(m2, n2, CN);
    OUT = new OUTCode();

    /* Link instructions. */
    Mdfg graph;
    graph.link(IN, P1Ins); graph.link(IN, P2Ins); graph.link(IN, P3Ins);
    graph.link(IN, P4Ins); graph.link(IN, P5Ins); graph.link(IN, P6Ins); graph.link(IN, P7Ins);
    graph.link(P1Ins, C11Ins); graph.link(P1Ins, C22Ins);
    graph.link(P2Ins, C22Ins); graph.link(P2Ins, C21Ins);
    graph.link(P3Ins, C22Ins); graph.link(P3Ins, C12Ins);
    graph.link(P4Ins, C11Ins); graph.link(P4Ins, C21Ins);
    graph.link(P5Ins, C11Ins); graph.link(P5Ins, C12Ins);
    graph.link(P6Ins, C22Ins);
    graph.link(P7Ins, C11Ins);
    graph.link(C11Ins, OUT); graph.link(C22Ins, OUT); graph.link(C12Ins, OUT); graph.link(C21Ins, OUT);

    /* Create interpreter. */
    nornir::Parameters par("parameters.xml");
    nornir::dataflow::Interpreter inter(&par, &graph, &inp, &out);
    inter.start();
    inter.wait();
}






