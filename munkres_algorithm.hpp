#ifndef MUNKRES_ALGORITHM_HPP
#define MUNKRES_ALGORITHM_HPP

#include "matrix.hpp"


#include <algorithm>
#include <iostream>
#include <list>
#include <set>
#include <stack>
#include <stdexcept>
#include <vector>
#include <iterator>






class Munkres {

public:
    //user enters square matrix
    explicit Munkres(Matrix<int>& x)
                :C(x),
                n_row(C.nrows()),
                n_col(C.ncols()),
                cover_row(n_row, false),
                cover_col(n_col, false),
                star_matrix(n_row, std::vector<bool>(n_col, false)),
                prime_matrix(n_row, std::vector<bool>(n_col, false)){};


    size_t STEP_0();

    size_t STEP_1();

    size_t STEP_2();

    size_t STEP_3();

    size_t STEP_4();


    bool find_zero(size_t&, size_t&);

    bool col_starred (size_t&);

    bool row_starred (size_t&);


    std::pair<size_t, size_t> get_star_zero (size_t&);

    bool uncovered_zero ();

    int smallest_element();

    void check_validity ();

    Matrix<int> output() {

        /* To get the fnal assignment matrix */

        Matrix<int> assignment(n_row, n_col, 0);

        for (size_t row = 0; row < n_row; ++row) {
            for (size_t col = 0; col < n_col; ++col) {

                if (star_matrix[row][col]) {
                    assignment(row, col) = 1;
                }
            }
        }

        return assignment;
    };


private:

    Matrix<int> C;
    size_t n_row;
    size_t n_col;
    std::vector<bool>cover_row;
    std::vector<bool>cover_col;
    std::vector<std::vector<bool>> star_matrix;
    std::vector<std::vector<bool>>  prime_matrix;



    int emin;
    size_t row_z0;
    size_t col_z0;



};


Matrix<int> run_munkres_algorithm(Matrix<int> c);

#endif //MUNKRES_ALGORITHM_HPP
