#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "matrix.h"
#include "utils.h"


int matrix_is_init(struct matrix_s matrix) {
    return matrix.content == NULL ? 0 : 1;
}


int matrix_is_square(struct matrix_s matrix) {
    return matrix.row_nb == matrix.col_nb ? 1 : 0;
}


matrix_content mypow(matrix_content a, int exp) {
    if (exp == 0)
        return 1;

    matrix_content res = 1;
    for (int i = 0; i < exp; i++)
        res*=a;
    return res;
}


struct matrix_s matrix_create(unsigned int row_nb, unsigned int nb_col) {
    struct matrix_s matrix = {
        .row_nb = row_nb,
        .col_nb = nb_col,
        .content = NULL
    };

    matrix_content** mat = malloc(row_nb * sizeof(matrix_content*));

    if (mat == NULL)
        return matrix;

    for (unsigned int i = 0; i < row_nb; i++) {
        mat[i] = calloc(nb_col, sizeof(matrix_content));
        if (mat[i] == NULL) {
            for (unsigned int j = 0; j < i; j++)
                free(mat[j]);
            free(mat);
            return matrix;
        }
    }

    matrix.content = mat;
    return matrix;
}


void matrix_init(struct matrix_s* matrix, matrix_content* content) {
    if (!matrix_is_init(*matrix))
        return;

    for (unsigned int i = 0; i < matrix->row_nb; i++)
        for (unsigned int j = 0; j < matrix->col_nb; j++)
            matrix->content[i][j] = content[i * matrix->row_nb + j];
}



void matrix_delete(struct matrix_s* matrix) {
    if (!matrix_is_init(*matrix))
        return;

    for (unsigned int i = 0; i < matrix->row_nb; i++)
        free(matrix->content[i]);

    free(matrix->content);
    matrix->content = NULL;
}


struct matrix_s matrix_copy(struct matrix_s matrix) {
    struct matrix_s res =
        matrix_create(matrix.row_nb, matrix.col_nb);

    if (!matrix_is_init(matrix) || !matrix_is_init(res))
        return res;

    for (unsigned int i = 0; i < res.row_nb; i++)
        for (unsigned int j = 0; j < res.col_nb; j++)
            res.content[i][j] = matrix.content[i][j];

    return res;
}


struct matrix_s matrix_term_to_term_opp(
    struct matrix_s matrix_a,
    struct matrix_s matrix_b,
    matrix_content (*operation)(matrix_content, matrix_content)
) {
    unsigned int row_nb = matrix_a.row_nb;
    unsigned int nb_of_cols = matrix_a.col_nb;

    struct matrix_s res = {
        .row_nb = matrix_a.row_nb,
        .col_nb = matrix_a.col_nb,
        .content = NULL
    };

    if (
        matrix_a.row_nb != matrix_b.row_nb ||
        matrix_a.col_nb != matrix_b.col_nb
    )
        return res;

    res = matrix_create(row_nb, nb_of_cols);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int i = 0; i < row_nb; i++)
        for (unsigned int j = 0; j< nb_of_cols; j++)
            res.content[i][j] =
                operation(matrix_a.content[i][j], matrix_b.content[i][j]);

    return res;
}


struct matrix_s matrix_all_terms_opp(
    struct matrix_s matrix,
    matrix_content nb,
    matrix_content (*operation)(matrix_content, matrix_content)
) {
    struct matrix_s res = {
        .row_nb = matrix.row_nb,
        .col_nb = matrix.col_nb,
        .content = NULL
    };

    res = matrix_create(matrix.row_nb, matrix.col_nb);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int i = 0; i < matrix.row_nb; i++)
        for (unsigned int j = 0; j < matrix.col_nb; j++)
            res.content[i][j] += operation(nb, matrix.content[i][j]);

    return res;
}


struct matrix_s matrix_product(
    struct matrix_s matrix_a,
    struct matrix_s matrix_b
) {
    struct matrix_s res = {
        .row_nb = matrix_a.row_nb,
        .col_nb = matrix_b.col_nb,
        .content = NULL
    };

    if (matrix_a.col_nb != matrix_b.row_nb)
        return res;

    res = matrix_create(res.row_nb, res.col_nb);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int i = 0; i < matrix_a.row_nb; i++)
        for (unsigned int j = 0; j < matrix_b.col_nb; j++) {
            matrix_content cij = 0;
            for (unsigned int k = 0; k < matrix_a.col_nb; k++)
                cij += matrix_a.content[i][k] * matrix_b.content[k][j];

            res.content[i][j] = cij;
        }

    return res;
}


void matrix_print(struct matrix_s matrix) {
    for (unsigned int i = 0; i < matrix.row_nb; i++) {
        printf("| ");
        for (unsigned int j = 0; j < matrix.col_nb; j++)
            printf("\t%.2f", matrix.content[i][j]);
        printf("\t|\n");
    }
}


matrix_content matrix_trace(struct matrix_s matrix) {
    if (!matrix_is_square(matrix))
        return 0;

    matrix_content trace = 0;
    for (unsigned int i = 0; i < matrix.row_nb; i++)
        trace += matrix.content[i][i];

    return trace;
}


struct matrix_s matrix_transpose(struct matrix_s matrix) {
    struct matrix_s res = {
        .row_nb = matrix.col_nb,
        .col_nb = matrix.row_nb,
        .content = NULL
    };

    res = matrix_create(res.row_nb, res.col_nb);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int i = 0; i < matrix.row_nb; i++)
        for (unsigned int j = 0; j < matrix.col_nb; j++)
            res.content[j][i] = matrix.content[i][j];

    return res;
}


matrix_content matrix_det(struct matrix_s matrix) {
    if (!matrix_is_square(matrix))
        return 0;

    if (matrix.row_nb == 1)
        return matrix.content[0][0];

    if (matrix.row_nb == 2)
        return (
            matrix.content[0][0] * matrix.content[1][1] -
            matrix.content[1][0]* (matrix.content[0][1])
        );

    int det = 0;
    for (unsigned int i = 0; i < matrix.row_nb; i++) {
        struct matrix_s under_mat =
            matrix_create(matrix.row_nb - 1, matrix.col_nb - 1);

        if (!matrix_is_init(under_mat))
            return 0;

        for (unsigned int k = 0; k < matrix.row_nb - 1; k++) {
            for (unsigned int l = 0; l < i; l++)
                under_mat.content[k][l] = matrix.content[k+1][l];
            for (unsigned int l = i; l < matrix.col_nb-1; l++)
                under_mat.content[k][l] = matrix.content[k+1][l+1];
        }

        det += mypow(-1,i) * matrix.content[0][i] * matrix_det(under_mat);

        matrix_delete(&under_mat);
    }

    return det;
}


struct matrix_s matrix_extract(
    struct matrix_s matrix,
    unsigned int row_begin,
    unsigned int col_begin,
    unsigned int row_end,
    unsigned int col_end
) {
    struct matrix_s res = {
        .row_nb = 0,
        .col_nb = 0,
        .content = NULL
    };

    if (row_begin > matrix.row_nb || row_end > matrix.row_nb)
        return res;

    if (col_begin > matrix.col_nb || col_end > matrix.col_nb)
        return res;

    if (row_end < row_begin || col_end < col_begin)
        return res;

    res = matrix_create(row_end - row_begin + 1, col_end - col_begin + 1);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int i = row_begin; i <= row_end; i++)
        for (unsigned int j = col_begin; j <= col_end; j++)
            res.content[i][j] = matrix.content[i][j];

    return res;
}


struct matrix_s matrix_remove_line(
    struct matrix_s matrix,
    unsigned int row_to_remove
) {
    struct matrix_s res = {
        .row_nb = matrix.row_nb - 1,
        .col_nb = matrix.col_nb,
        .content = NULL
    };

    if (row_to_remove > matrix.row_nb)
        return res;

    res = matrix_create(matrix.row_nb - 1, matrix.col_nb);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int i = 0; i < matrix.row_nb; i++) {
        if (i == row_to_remove)
            continue;

        for (unsigned int j = 0; j < matrix.col_nb; j++)
            if (i > row_to_remove)
                res.content[i-1][j] = matrix.content[i][j];
            else
                res.content[i][j] = matrix.content[i][j];
    }

    return res;
}


struct matrix_s matrix_remove_column(
    struct matrix_s matrix,
    unsigned int column_to_remove
) {
    struct matrix_s res = {
        .row_nb = matrix.row_nb,
        .col_nb = matrix.col_nb - 1,
        .content = NULL
    };

    if (column_to_remove > matrix.col_nb)
        return res;

    res = matrix_create(matrix.row_nb, matrix.col_nb - 1);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int j = 0; j < matrix.col_nb; j++) {
        if (j == column_to_remove)
            continue;
        for (unsigned int i = 0; i < matrix.row_nb; i++)
            if (j > column_to_remove)
                res.content[i][j-1] = matrix.content[i][j];
            else
                res.content[i][j] = matrix.content[i][j];
    }

    return res;
}


struct matrix_s matrix_remove_line_and_column(
    struct matrix_s matrix,
    unsigned int row_to_remove,
    unsigned int column_to_remove
) {
    // Do not use other functions for performance reasons

    struct matrix_s res = {
        .row_nb = matrix.row_nb - 1,
        .col_nb = matrix.col_nb - 1,
        .content = NULL
    };

    if (
        column_to_remove > matrix.col_nb ||
        row_to_remove > matrix.row_nb
    )
        return res;

    res = matrix_create(matrix.row_nb - 1, matrix.col_nb - 1);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int i = 0; i < matrix.row_nb; i++) {
        if (i == row_to_remove)
            continue;
        unsigned int line_nb_in_new = i > row_to_remove ? i-1 : i;

        for (unsigned int j = 0; j < matrix.col_nb; j++) {
            if (j == column_to_remove)
                continue;

            if (j > column_to_remove)
                res.content[line_nb_in_new][j-1] = matrix.content[i][j];
            else
                res.content[line_nb_in_new][j] = matrix.content[i][j];
        }
    }

    return res;
}


struct matrix_s matrix_co(struct matrix_s matrix) {
    struct matrix_s res = {
        .row_nb = matrix.row_nb,
        .col_nb = matrix.col_nb,
        .content = NULL
    };

    if (!matrix_is_square(matrix))
        return res;

    res = matrix_create(matrix.row_nb, matrix.col_nb);

    if (!matrix_is_init(res))
        return res;

    for (unsigned int i = 0; i < matrix.row_nb; i++)
        for (unsigned int j = 0; j < matrix.col_nb; j++) {
            struct matrix_s under_mat =
                matrix_remove_line_and_column(matrix, i, j);

            if (!matrix_is_init(under_mat)) {
                matrix_delete(&res);
                return res;
            }

            res.content[i][j] = mypow(-1,i+j) * matrix_det(under_mat);
            matrix_delete(&under_mat);
        }

    return res;
}


matrix_content simple_addition(matrix_content a, matrix_content b) {
    return a + b;
}


matrix_content simple_substraction(matrix_content a, matrix_content b) {
    return a - b;
}


matrix_content simple_product(matrix_content a, matrix_content b) {
    return a * b;
}


matrix_content simple_division(matrix_content a, matrix_content b) {
    return a / b;
}


struct matrix_s matrix_invert(struct matrix_s matrix) {
    struct matrix_s pre_res = {
        .row_nb = matrix.col_nb,
        .col_nb = matrix.row_nb,
        .content = NULL
    };

    if (!matrix_is_square(matrix))
        return pre_res;

    struct matrix_s mat_inv_trans = matrix_co(matrix);

    if (!matrix_is_init(mat_inv_trans))
        return pre_res;

    float det = matrix_det(matrix);

    if ( det == 0 )
        return pre_res;

    matrix_content coef = 1 / det;

    pre_res = matrix_transpose(mat_inv_trans);
    matrix_delete(&mat_inv_trans);

    if (!matrix_is_init(pre_res))
        return pre_res;

    struct matrix_s res = matrix_all_terms_opp(pre_res, coef, simple_product);
    matrix_delete(&pre_res);

    return res;
}


struct matrix_s matrix_left_pseudo_inv(struct matrix_s matrix) {
    struct matrix_s res = {
        .row_nb = matrix.col_nb,
        .col_nb = matrix.row_nb,
        .content = NULL
    };

    if (matrix.col_nb > matrix.row_nb)
        return res;

    struct matrix_s mat_trans = matrix_transpose(matrix);
    if (!matrix_is_init(mat_trans))
        return res;

    struct matrix_s term_1 = matrix_product(mat_trans, matrix);
    if (!matrix_is_init(term_1)) {
        matrix_delete(&mat_trans);
        return res;
    }

    struct matrix_s term_1_inv = matrix_invert(term_1);
    matrix_delete(&term_1);

    if (!matrix_is_init(term_1_inv)) {
        matrix_delete(&mat_trans);
        return res;
    }

    res = matrix_product(term_1_inv, mat_trans);
    matrix_delete(&mat_trans);
    matrix_delete(&term_1_inv);

    return res;
}


struct matrix_s matrix_right_pseudo_inv(struct matrix_s matrix) {
    struct matrix_s res = {
        .row_nb = matrix.col_nb,
        .col_nb = matrix.row_nb,
        .content = NULL
    };

    if (matrix.col_nb < matrix.row_nb)
        return res;

    struct matrix_s mat_trans = matrix_transpose(matrix);
    if (!matrix_is_init(mat_trans))
        return res;

    struct matrix_s term_1 = matrix_product(matrix, mat_trans);
    if (!matrix_is_init(term_1)) {
        matrix_delete(&mat_trans);
        return res;
    }

    struct matrix_s term_1_inv = matrix_invert(term_1);
    matrix_delete(&term_1);
    if (!matrix_is_init(term_1_inv)) {
        matrix_delete(&mat_trans);
        return res;
    }

    res = matrix_product(mat_trans, term_1_inv);

    matrix_delete(&mat_trans);
    matrix_delete(&term_1_inv);

    return res;
}
