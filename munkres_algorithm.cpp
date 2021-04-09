#include "munkres_algorithm.hpp"



void Munkres::check_validity() {

    /* FUNCTION to check if given matrix is valid*/

    if(n_row != n_col){

        std::cerr << "Invalid Matrix. \n";
        exit(0);

    }


    for (size_t row = 0; row < n_row; row++) {

        for (size_t col = 0; col < n_col; col++) {

            if(C(row,col) < 0){

                std::cerr << "Invalid Matrix. \n";
                exit(0);
            }

        }
    }
}


bool Munkres::row_starred (size_t& row){

    /*FUNCTION to check if row has starred zero*/

    size_t size = star_matrix.size();

    for (size_t col = 0; col < size; col++) {

        if(star_matrix[row][col]) { return true; };

    }

    return false;
}



bool Munkres::col_starred (size_t& col){

    /*FUNCTION to check if column has starred zero*/

    size_t size = star_matrix.size();

    for (size_t row = 0; row < size; row++) {

        if(star_matrix[row][col]) { return true; };

    }


    return false;
}


bool Munkres:: find_zero (size_t&  row, size_t& col){

    /* FUNCTION to find if given cell holds zero*/

    if(C(row,col) == 0){ return true; }

    return false;

}



bool Munkres:: uncovered_zero (){

    /* FUNCTION to find if cost matrix holds uncovered zero*/


    for (size_t row = 0; row < n_row; row++) {

        for (size_t col = 0; col < n_col; col++) {

            if ((!cover_col[col] && !cover_row[row]) && find_zero(row, col)) {

                return true;
            }
        }
    }

    return false;

}



std::pair<size_t, size_t> Munkres::get_star_zero (size_t& col){

    /*FUNCTION to find a starred zero in a column*/

    for(size_t row =0; row < n_row; row++) {

        if (star_matrix[row][col] == 1) { return std::make_pair(row,col); }

    }

    return std::pair<size_t, size_t>();
}




int Munkres::smallest_element () {

    /* FUNCTION returns the smallest uncovered element*/

    int min = INT16_MAX;
    for (size_t row = 0; row < n_row; row++) {

        if (!cover_row[row]) {

            for (size_t col = 0; col < n_col; col++) {

                if (!cover_col[col]) {

                    if (C(row, col) < min) { min = C(row, col); }
                }
            }
        }
    }

    return min;
}



size_t Munkres::STEP_0() {


    /*
     *  For each row r in C, subtract its smallest element from every element in r
     *For each column c in C, subtract its smallest element from every element in c
     */


    for (size_t row = 0; row < n_row; row++) {

        int row_min = INT16_MAX;

        for (size_t col = 0; col < n_col; col++) {

            if (C(row, col) < row_min) { row_min = C(row, col); }
        }


        for (size_t col = 0; col < n_col; col++) {

            C(row, col) = C(row, col) - row_min;
        }
    }


    for (size_t col = 0; col < n_col; col++) {

        int col_min = INT16_MAX;

        for (size_t row = 0; row < n_row; row++) {

            if (C(row, col) < col_min) { col_min = C(row, col); }
        }


        for (size_t row = 0; row < n_row; row++) {

            C(row, col) = C(row, col) - col_min;
        }
    }


    /* For all zeros z_i in C, mark z_i with a star if there is no starred zero in its row or column*/

    for (size_t row = 0; row < n_row; row++) {

        for (size_t col = 0; col < n_col; col++) {

            if (!row_starred(row)  && !col_starred(col) && find_zero(row,col)) {

                star_matrix[row][col] = true;
            }
        }
    }

    return 1;
}


/******************************************/


size_t Munkres::STEP_1() {


    /*
     * for Each column containing a starred zero, cover this column
     * if n columns are covered then GOTO DONE
     */


    for (size_t row = 0; row < n_row; row++) {

        for (size_t col = 0; col < n_col; col++) {

            if (col_starred(col)) {

                cover_col[col] = true;
            }
        }
    }




    const bool covered = std::all_of(cover_col.begin(), cover_col.end(), [](bool v) {return v; });




    if(covered){

        return 5;
    }

    return 2;
}


/*******************************************************/


size_t Munkres::STEP_2() {

    /*
     * find uncovered zero
     * if found:
     * Prime it
     * check for Z* in the row
     * if: No Z* -> step 3
     * else: cover this row, uncover the col -> step 2
     * else:
     * save the smallest uncovered element, e_min -> step 4
    */



    if(uncovered_zero()) {


        for (size_t row = 0; row < n_row; row++) {

            for (size_t col = 0; col < n_col; col++) {

                if ((!cover_col[col] && !cover_row[row]) && find_zero(row, col)) { //find uncovered zero

                    prime_matrix[row][col] = true; //prime it


                    if (!row_starred(row)) { //if no starred zero in row

                        row_z0 = row;
                        col_z0 = col;


                        return 3;

                    } else { //if starred zero

                        cover_row[row] = true; //cover row

                        for (size_t column = 0; column < n_col; column++) {

                            if (star_matrix[row][column]) { //uncover col containing star zero

                                cover_col[column] = false;

                                return 2;
                            }
                        }
                    }
                }
            }
        }

    } else { //save the smallest uncovered element, e_min

        emin = smallest_element();
    }


    return 4;
}


/**************************************************/

size_t Munkres::STEP_3() {


    /*
 * Construct a series S of alternating primed and starred zeros as follows:
 * Insert Z 0 into S
 * while In the column of Z 0 exists a starred zero Z 1 do
 * Insert Z 1 into S
 * Replace Z 0 with the primed zero in the row of Z 1 . Insert Z 0 into S
 * end while
 * */


    std::pair<size_t ,size_t > Z0; //primed zero
    std::pair<size_t ,size_t > Z1; // starred zero
    Z0 = std::make_pair(row_z0,col_z0);
    std::vector<std::pair<size_t, size_t>> S;
    S.emplace_back(Z0);


    while(col_starred(Z0.second)){ //if column of Z0 has a starred zero Z1

            Z1 = get_star_zero(Z0.second);
            S.emplace_back(Z1); // insert Z1 to S

            for(size_t col =0; col < n_col;  col++){

                if(prime_matrix[Z1.first][col]){ //Replace Z0 with the primed zero in the row of Z1

                    Z0 = std::make_pair(Z1.first,col);
                    S.emplace_back(Z0); //Insert Z0 into S

                    break;
                }
            }
    }


/*
 * Unstar each starred zero in S and replace all primes with stars.
 * Erase all other primes and uncover every line in C GOTO STEP 1
*/

for(size_t index =0; index < S.size(); index=index+2){


    size_t z0_row = S[index].first;
    size_t z0_col = S[index].second;


    star_matrix[z0_row][z0_col] = true; //replace all primes with stars

    if(index != S.size() -1) {

        size_t z1_row = S[index + 1].first;
        size_t z1_col = S[index + 1].second;

        star_matrix[z1_row][z1_col] = false; //Unstar each starred zero in S
    }
}


for(size_t i =0; i < n_row; i++){

    for(size_t j =0; j < n_row; j++){

        prime_matrix[i][j] = false; //Erase all other primes
    }
}



//uncover every line in C
for(auto row:cover_row){

    row = false;
}


for(auto col:cover_col){

    col= false;
}


return 1;

}


/***************************************************/


size_t Munkres::STEP_4() {

    /*Add e min to every element in covered rows
     * subtract it from every element in uncovered columns.
     * GOTO STEP 2
     */



    for (size_t row = 0; row < n_row; row++) {

        if (cover_row[row]) {

            for (size_t col = 0; col < n_col; col++) {

                C(row, col) = C(row,col) + emin;
            }
        }
    }


    for (size_t col = 0; col < n_row; col++) {

        if (!cover_col[col]) {

            for (size_t row = 0; row < n_col; row++) {

                C(row, col) = C(row,col) - emin;
            }
        }
    }


    return 2;

}


/***************************************/







Matrix<int> run_munkres_algorithm(Matrix<int> C) {



    Munkres munkres_matrix(C);

    munkres_matrix.check_validity();

    size_t step = munkres_matrix.STEP_0();

    bool DONE = false;

    while (!DONE) {

        switch (step) {
            case 1:
                step = munkres_matrix.STEP_1();
                break;
            case 2:
                step = munkres_matrix.STEP_2();
                break;
            case 3:
                step = munkres_matrix.STEP_3();
                break;
            case 4:
                step = munkres_matrix.STEP_4();
                break;
            default:
                DONE = true;
                break;
        }
    }


    return munkres_matrix.output();
}





