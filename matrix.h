#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


#define LENGTH_PRINT 10

typedef float matrix_content;


struct matrix_s {
    int nb_of_lines;
    int nb_of_columns;
    matrix_content** content;
};


/**
 * @brief Check if a matrix is correctly initialized.
 *
 * @param matrix The matrix to check.
 * @return An integer: 0 if the matrix is not initialized, else 1.
 */
int is_matrix_init(struct matrix_s matrix);


/**
 * @brief Check if a matrix is a square matrix.
 *
 * @param matrix The matrix to test.
 * @return An integer: 0 if the matrix is not square, else 1.
 */
int is_matrix_square(struct matrix_s matrix);


/**
 * @brief Create a matrix.
 *
 * @param nb_of_lines The number of lines of the matrix.
 * @param nb_of_col The number of columns of the matrix.
 * @result A new matrix.
*/
struct matrix_s create_matrix(int nb_of_lines, int nb_of_col);


/**
 * @brief Delete a matrix.
 *
 * @param matrix The address of the matrix to delete.
*/
void delete_matrix(struct matrix_s* matrix);


/**
 * @brief Copy an existing matrix.
 *
 * @param matrix The matrix to copy.
 * @return A new matrix that is the copy of the original matrix.
 */
struct matrix_s copy_matrix(struct matrix_s matrix);


/**
 * @brief A helper function to ease the use of "matrix_term_to_term_opp" and
 * "matrix_all_terms_opp". It makes a simple addition.
 *
 * @param a The first operand.
 * @param b The second operand.
 */
matrix_content simple_addition(matrix_content a, matrix_content b);


/**
 * @brief A helper function to ease the use of "matrix_term_to_term_opp" and
 * "matrix_all_terms_opp". It makes a simple substraction.
 *
 * @param a The first operand.
 * @param b The second operand.
 */
matrix_content simple_substraction(matrix_content a, matrix_content b);


/**
 * @brief A helper function to ease the use of "matrix_term_to_term_opp" and
 * "matrix_all_terms_opp". It makes a simple multiplication.
 *
 * @param a The first operand.
 * @param b The second operand.
 */
matrix_content simple_product(matrix_content a, matrix_content b);


/**
 * @brief A helper function to ease the use of "matrix_term_to_term_opp" and
 * "matrix_all_terms_opp". It makes a simple division.
 *
 * @param a The first operand.
 * @param b The second operand.
 */
matrix_content simple_division(matrix_content a, matrix_content b);


/**
 * @brief Calculate a matrix addition.
 *
 * @param matrix_a The first matrix.
 * @param matrix_b The second matrix.
 * @param operation The function to apply.
 * @result A new matrix, the result of the operation.
*/
struct matrix_s matrix_term_to_term_opp(
    struct matrix_s matrix_a,
    struct matrix_s matrix_b,
    matrix_content (*operation)(matrix_content, matrix_content)
);


/**
 * @brief Calculate a matrix addition.
 *
 * @param matrix The matrix on which the operation will be performed.
 * @param nb The number with which the operations shall be performed
 * @param operation The opperation to apply. The first argument is nb and the
 * second is the matrix content.
 * @return A new matrix with the result.
*/
struct matrix_s matrix_all_terms_opp(
    struct matrix_s matrix,
    matrix_content nb,
    matrix_content (*operation)(matrix_content, matrix_content)
);


/**
 * @brief Calculate the product of 2 matrix.
 *
 * Note that for a product the number of columns of the first matrix shall be
 * equal to the number of line of the second matrix.
 *
 * @param matrix_a The first matrix.
 * @param matrix_b The second matrix.
 * @result A new matrix (two dimension table), the product of 2 matrix.
*/
struct matrix_s matrix_prod(struct matrix_s matrix_a, struct matrix_s matrix_b);


/**
* @brief Print a matrix.
*
* @param matrix The matrix to print.
*/
void print_matrix(struct matrix_s matrix);


/**
 * @brief Calculate a matrix trace.
 *
 * The matrix MUST have the same number of line and column.
 * If not, 0 will be returned.
 *
 * @param matrix The matrix of which the trace will be computed.
 * @result The trace of the matrix.
*/
matrix_content matrix_trace(struct matrix_s matrix);


/**
 * @brief Calculate a matrix transpose.
 *
 * @param matrix The matrix to transpose.
 * @result A new matrix that is the transposed of the given one.
*/
struct matrix_s matrix_transpose(struct matrix_s matrix);


/**
 * @brief Calculate a matrix determinant.
 *
 * The matrix MUST have the same number of lines and columns. Else 0 will be
 * returned.
 *
 * @param matrix The matrix of which the determinant will be computed.
 * @result The determinant of the matrix.
*/
matrix_content matrix_det(struct matrix_s matrix);


/**
 * @brief Extract a matrix from an other.
 *
 * @param matrix The original matrix.
 * @param line_begin The first line to extract.
 * @param col_begin The first column to extract.
 * @param line_end The last line to extract.
 * @param col_end The last column to extract.
 * @result A new matrix that corresponds to the extracted matrix.
*/
struct matrix_s matrix_extract(
    struct matrix_s matrix,
    int line_begin,
    int col_begin,
    int line_end,
    int col_end
);


/**
 * @brief Remove a line from a matrix.
 *
 * @param matrix The matrix to remove the line .
 * @param line_to_remove The line number to remove (starting from 0).
 * @return A new matrix that correspond to the given matrix without the
 * given line.
 */
struct matrix_s remove_line_from_matrix(
    struct matrix_s matrix,
    int line_to_remove
);


/**
 * @brief Remove a column from a matrix.
 *
 * @param matrix The matrix to remove the column from.
 * @param column_to_remove The column number to remove (starting from 0).
 * @return A new matrix that correspond to the given matrix without the
 * given column.
 */
struct matrix_s remove_column_from_matrix(
    struct matrix_s matrix,
    int column_to_remove
);


/**
 * @brief Remove a line and a column from a matrix.
 *
 * @param matrix The matrix to remove the line and column from.
 * @param line_to_remove The line number to remove (starting from 0).
 * @param column_to_remove The column number to remove (starting from 0).
 * @return A new matrix that correspond to the given matrix without the
 * given line and given column.
 */
struct matrix_s remove_line_and_column_from_matrix(
    struct matrix_s matrix,
    int line_to_remove,
    int column_to_remove
);


/**
 * @brief Calculate a matrix co matrix.
 *
 * @param matrix The number of lines of the matrix.
 * @param nb_col: the number of columns of the matrix.
 * @result A new matrix that correspond to theco matrix of the given matrix.
*/
struct matrix_s co_matrix(struct matrix_s matrix);


/**
 * @brief Calculate a matrix invert.
 *
 * @param matrix The matrix to invert.
 * @result A new metrix that correspond to the invert of the given matrix.
*/
struct matrix_s matrix_invert(struct matrix_s matrix);


/**
 * @brief Calculate the left pseudo invert.
 *
 * @param matrix The matrix to invert.
 * @result A new matrix that corresponds to the pseudo invert left of the
 * matrix.
*/
struct matrix_s left_pseudo_inv(struct matrix_s matrix);


/**
 * @brief Calculate the left pseudo invert.
 *
 * @param matrix The matrix to invert.
 * @result A new matrix that corresponds to the pseudo invert right of the
 * matrix.
*/
struct matrix_s right_pseudo_inv(struct matrix_s matrix);


#endif //MATRIX_H
