#include <fstream>
#include <iostream>
#include <vector>


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
    ia.resize(N);
    int jasize = xysize(Nx, Ny, k1, k2);
    // std::cout << jasize;
    // return;
    ja.resize(jasize);
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
    }
    return 0;
}