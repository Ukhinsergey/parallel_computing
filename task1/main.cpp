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

        //upright 
        if (r != 0) { 
            int temp = (i - (Nx + 1) - (r - 1)) / (k1 + k2) * k2 + std::max(0, (i - (Nx + 1) - (r - 1)) % (k1 + k2) - k1);
            pos += temp;
        }

        //downleft
        if (r != Ny) {
            int c0 = 1;
            if (c == 0 ) {
                c0 = 0;
            }
            int temp = (i - r - c0) / (k1 + k2) * k2 + std::max(0, (i - r - c0) % (k1 + k2) - k1);
            pos += temp;

        } else {
            int temp = (i - c - r) / (k1 + k2) * k2 + std::max(0, (i - r - c) % (k1 + k2) - k1);
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

void Fill(int N, std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &A, std::vector<double> &b, std::vector<double> &M) {
    A.resize(ia[N]);
    b.resize(N);
    M.resize(N);
    for(int i = 0 ; i < N; ++i) {
        double sum = 0;
        int pos = 0;
        int curpos = ia[i];
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
        M[i] = 1.0 / sum;
        b[i] = std::sin(i);
    }
}

double dot(const std::vector<double> &x, const std::vector<double> &y) {
    double res = 0;
    int n = x.size();
    int n2 = y.size();
    if(n != n2) {
        std::cout << "mismatch vector sizes for dot " << n << ' ' << n2 << std::endl;
        exit(0);
    }
    for(int i = 0 ; i < n; ++i) {
        res += x[i] * y[i];
    }
    return res;
}

void axpby(const std::vector<double> &x, const std::vector<double> &y, double a, double b, std::vector<double> &res) {
    int n = x.size();
    int n2 = y.size();
    if(n != n2) {
        std::cout << "mismatch vector sizes for axpby " << n << ' ' << n2 << std::endl;
        exit(0);
    }
    res.resize(n);
    for(int i = 0 ; i < n; ++i) {
        res[i] = a * x[i] + b * y[i];
    }
}

void SpMv(const std::vector<int> &ia, const std::vector<int> &ja, const std::vector<double> &a, const  std::vector<double> &b, std::vector<double> &res) {
    int n = ia.size();
    res.resize(n - 1);
    for(int i = 0; i < n - 1; ++i){
        double sum = 0.0;
        const int jb = ia[i];
        const int je = ia[i+1];
        for(int j=jb; j<je; ++j) sum += a[j]*b[ja[j]];
        res[i] = sum;
    }
}


void Solve(int N, const std::vector<int> &ia, const std::vector<int> &ja, const std::vector<double> &A, const std::vector<double> &b, const std::vector<double> &M, double tol, std::vector<double> &x, int &k, double &res) {
    std::vector<double> r(N);
    std::vector<double> tmp(N);
    std::vector<double> z(N);
    std::vector<double> p(N);
    std::vector<double> q(N);
    int maxiter = 10;
    double beta;
    double rhonow;
    double rhoprev;
    double alpha;
    SpMv(ia, ja, A, x, tmp);
    for(int i = 0; i < N; ++i) {
        r[i] = b[i] - tmp[i];
    }
    bool convergence = false;
    k = 1;
    while(!convergence) {
        for(int i = 0 ; i < N; ++i) {
            z[i] = M[i] * r[i];
        }
        rhonow = dot(r,z);
        std::cout << std::sqrt(dot(r,r)) << ' ' <<rhonow << std::endl;
        if (k == 1) {
            p = z;
        } else {
            beta = rhonow / rhoprev;
            for(int i = 0 ; i < N; ++i) {
                p[i] = z[i] + beta * p[i];
            }
        }
        SpMv(ia, ja, A, p, q);
        alpha = rhonow / dot(p, q);
        for(int i = 0; i < N; ++i) {
            x[i] = x[i] + alpha * p[i];
            r[i] = r[i] - alpha * q[i];
        }
        if (rhonow < tol || k >= maxiter) {
            convergence = true;
        } else {
            k++;
        }
        rhoprev = rhonow;
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
    std::vector<double> M;
    Fill(N, ia, ja, A, b, M);
    double tol = 0.0001;
    std::vector<double> x(N);
    for(int i = 0 ; i < N; ++i) {
        x[i] = 0.0;
    }
    int k;
    double res;
    Solve(N, ia, ja, A, b, M, tol, x, k, res);
    // std::vector<double> x {1, 3.4, 2.5};
    // std::vector<double> y {2.1, 2.3, 4};
    // std::cout << dot(x,y) << std::endl;
    // std::vector<double> res;
    // axpby(x,y,1,2.5,res);
    // for(auto i: res) {
    //     std::cout << i << ' ';
    // }
    for(auto i: x) {
        std::cout << i << ' ';
    }
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