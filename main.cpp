#include <iostream>

#include "munkres_algorithm.hpp"
#include "matrix.hpp"






int main(int /*argc*/, const char* /*argv*/[]) {
  //  const Matrix<int> m = {{250, 400, 350}, {400, 600, 350}, {200, 400, 250}};
 //  const Matrix<int> s = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};

 const Matrix<int> s = {{0,0,0,0,1,0},{0,0,1,0,0,0},{0,0,0,0,0,1},{0,0,0,1,0,0},{0,1,0,0,0,0},{1,0,0,0,0,0}};

  const Matrix<int> m = {{68,35,37,10,47,31},{10,70,26,52,58,74},{71,59,86,65,84,40},{65,87,53,4,69,77},
                         {33,28,31,68,67,38},{5,34,72,93,95,18}};


    const auto res = run_munkres_algorithm(m);

    if(s.nrows() != res.nrows() || s.ncols() != res.ncols()) {
        std::cerr << "Cannot compare matrices" << std::endl;
        return 1;
    }

    for(size_t i = 0; i < s.nrows(); i++) {
        for(size_t j = 0; j < s.ncols(); j++) {
            if(s(i,j) != res(i,j)){
                std::cerr << "Found invalid result" << std::endl;
                return 1;
            }
        }
    }



    return 0;
}
