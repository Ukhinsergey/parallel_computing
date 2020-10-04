#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>


void Generate(int Nx, int Ny, int k1, int k2, int &N, std::vector<int> &ia, std::vector<int> &ja){
    N = (Nx + 1) * (Ny + 1);
    int numsinrow = 0;
    ia.push_back(0);
    for(int i = 0 ; i < N; ++i) {
        int r = i / (Nx + 1);
        int c = i % (Nx + 1);
        //up
        if (r != 0) {
            ja.push_back(i - (Nx + 1));
            numsinrow++;
        }
        //upright
        if (r != 0 && c != Nx && ( ((r - 1) * Nx + c) % (k1 + k2) >= k1)) {
            ja.push_back(i - Nx);
            numsinrow++;
        }
        //left
        if (c != 0) {
            ja.push_back(i - 1);
            numsinrow++;
        }
        //middle
        ja.push_back(i);
        numsinrow++;
        //right
        if (c != Nx) {
            ja.push_back(i + 1);
            numsinrow++;
        }
        //downleft
        if (c != 0 && r != Ny && ((r * Nx + c - 1) % (k1 + k2) >= k1) ) {
            ja.push_back(i + Nx);
            numsinrow++;          
        }
        //down
        if(r != Ny) {
            ja.push_back(i + Nx + 1);
            numsinrow++;
        }
        ia.push_back(numsinrow);
    }

}


int xysize(int Nx, int Ny,int k1,int k2) {
    return 4 * 3 + (Nx - 1) * 8 + (Ny - 1) * 8 + (Nx - 1) * (Ny - 1) * 5 + Nx * Ny / (k1 + k2) * 2 * k2 + std::max(0, Nx * Ny % (k1 + k2) - k1) * 2;
}

void Generate1(int Nx, int Ny, int k1, int k2, int &N, std::vector<int> &ia, std::vector<int> &ja){
    N = (Nx + 1) * (Ny + 1);
    ia.resize(N + 1);
    int jasize = xysize(Nx, Ny, k1, k2);
    ia[N] = jasize;
    ja.resize(jasize);
    for(int i = 0 ; i < N; ++i) {
        int r = i / (Nx + 1);
        int c = i % (Nx + 1);
        int pos = 0;
        if (r != 0) {
            pos += 3 * 2 + (r - 1) * 4 * 2 + (Nx - 1) * 4 + (r - 1) * (Nx- 1) * 5;
            if (c != 0 && r != Ny) {
                pos += 4 + (c - 1) * 5;
            } else if (c != 0) {
                pos += 3 + (c - 1) * 4;
            }
        } else {
            if (c != 0) {
                pos += 3 + (c - 1) * 4;
            } else {
                pos = 0;
            }
        }
        //^ pos without diag

        //now upright 
        if (r != 0) { 
            int temp = (i - (Nx + 1) - (r - 1)) / (k1 + k2) * k2 + std::max(0, (i - (Nx + 1) - (r - 1)) % (k1 + k2) - k1);
            pos += temp;
        }

        //now downleft
        if (r != Ny) {
            int c0 = 1;
            if (c == 0 ) {
                c0 = 0;
            }
            int temp = (i - r - c0) / (k1 + k2) * k2 + std::max(0, (i - r - c0) % (k1 + k2) -k1);
            pos += temp;

        } else {
            int temp = (i - c - r) / (k1 + k2) * k2 + std::max(0, (i - r - c) % (k1 + k2) -k1);
            pos += temp;
        }


        ia[i] = pos;
        //up
        if (r != 0) {
            ja[pos++] = (i - (Nx + 1));
        }
        //upright
        if (r != 0 && c != Nx && ( ((r - 1) * Nx + c) % (k1 + k2) >= k1)) {
            ja[pos++] = (i - Nx);
        }
        //left
        if (c != 0) {
            ja[pos++] = (i - 1);
        }
        //middle
        ja[pos++] = (i);
        //right
        if (c != Nx) {
            ja[pos++] = (i + 1);
        }
        //downleft
        if (c != 0 && r != Ny && ((r * Nx + c - 1) % (k1 + k2) >= k1) ) {
            ja[pos++] =(i + Nx);          
        }
        //down
        if(r != Ny) {
            ja[pos++] = (i + Nx + 1);
        }
    }
}


void Fill(int N, std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &A, std::vector<double> &b) {
    A.resize(ia[N]);
    b.resize(N);
    int curpos = 0;
    for(int i = 0 ; i < N; ++i) {
        double sum = 0;
        int pos = 0;
        for(int j = ia[i]; j < ia[i + 1]; ++j) {
            if (ja[j] == i) {
                pos = curpos;
                
            } else {
                A[curpos] = std::cos(i * ja[j] + i + ja[j]);
                sum += std::abs(A[curpos]);
                
            }
            curpos++;
        }
        sum *=1.234;
        A[pos] = sum;
        b[i] = std::sin(i);
    }
}

int main(int argc, char **argv) {
    int Nx, Ny;
    int k1, k2;
    bool debug = false;
    if (argc != 5 && argc != 6) {
        std::cout << "input: Nx, Ny, k1, k2, debug(optional)" << std::endl;
        return 1;
    } else {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        k1 = atoi(argv[3]);
        k2 = atoi(argv[4]);
        if (argc == 6) {
            debug = atoi(argv[5]);
        }
    }
    std::vector<int> ia;
    std::vector<int> ja;
    int N;
    Generate1(Nx, Ny, k1, k2, N, ia, ja);
    std::vector<double> A;
    std::vector<double> b;
    Fill(N, ia, ja, A, b);
    if (debug) {
        std::ofstream fout("output.txt");
        fout << "N = " << N << std::endl;
        fout << "JA : ";
        for(const auto &j:ja) {
            fout << j << ' ';
        }
        fout << std::endl;
        fout << "IA : ";
        for(const auto &i: ia) {
            fout << i << ' ';
        }
        fout << std::endl;
        fout << "A : ";
        for(const auto &i: A) {
            fout << i << ' ';
        }
        fout << std::endl;
        fout << "b : ";
        for(const auto &i: b) {
            fout << i << ' ';
        }
        fout << std::endl;
    }
    return 0;
}