#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "matrix.h"
#include "utils.h"


int is_matrix_init(struct matrix_s matrix) {
    return matrix.content == NULL ? 0 : 1;
}


int is_matrix_square(struct matrix_s matrix) {
    return matrix.nb_of_lines == matrix.nb_of_columns ? 1 : 0;
}


matrix_content mypow(matrix_content a, int exp) {
    if (exp == 0)
        return 1;

    matrix_content res = 1;
    for (int i = 0; i < exp; i++)
        res*=a;
    return res;
}


struct matrix_s create_matrix(int nb_lines, int nb_col) {
    struct matrix_s matrix = {
        .nb_of_lines = nb_lines,
        .nb_of_columns = nb_col,
        .content = NULL
    };

    matrix_content** mat = malloc(nb_lines*sizeof(matrix_content*));

    if (mat == NULL)
        return matrix;

    for (int i = 0; i < nb_lines; i++) {
        mat[i] = calloc(nb_col, sizeof(matrix_content));
        if (mat[i] == NULL) {
            for (int j = 0; j < i; j++)
                free(mat[j]);
            free(mat);
            return matrix;
        }
    }

    matrix.content = mat;
    return matrix;
}


void delete_matrix(struct matrix_s* matrix) {
    if (!is_matrix_init(*matrix))
        return;

    for (int i = 0; i < matrix->nb_of_lines; i++)
        free(matrix->content[i]);

    free(matrix->content);
    matrix->content = NULL;
}


struct matrix_s copy_matrix(struct matrix_s matrix) {
    struct matrix_s res =
        create_matrix(matrix.nb_of_lines, matrix.nb_of_columns);

    if (!is_matrix_init(matrix) || !is_matrix_init(res))
        return res;

    for (int i = 0; i < res.nb_of_lines; i++)
        for (int j = 0; j < res.nb_of_columns; j++)
            res.content[i][j] = matrix.content[i][j];

    return res;
}


struct matrix_s matrix_term_to_term_opp(
    struct matrix_s matrix_a,
    struct matrix_s matrix_b,
    matrix_content (*operation)(matrix_content, matrix_content)
) {
    int nb_of_lines = matrix_a.nb_of_lines;
    int nb_of_cols = matrix_a.nb_of_columns;

    struct matrix_s res = {
        .nb_of_lines = matrix_a.nb_of_lines,
        .nb_of_columns = matrix_a.nb_of_columns,
        .content = NULL
    };

    if (
        matrix_a.nb_of_lines != matrix_b.nb_of_lines ||
        matrix_a.nb_of_columns != matrix_b.nb_of_columns
    )
        return res;

    res = create_matrix(nb_of_lines, nb_of_cols);

    if (!is_matrix_init(res))
        return res;

    for (int i = 0; i < nb_of_lines; i++)
        for (int j = 0; j< nb_of_cols; j++)
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
        .nb_of_lines = matrix.nb_of_lines,
        .nb_of_columns = matrix.nb_of_columns,
        .content = NULL
    };

    res = create_matrix(matrix.nb_of_lines, matrix.nb_of_columns);

    if (!is_matrix_init(res))
        return res;

    for (int i = 0; i < matrix.nb_of_lines; i++)
        for (int j = 0; j < matrix.nb_of_columns; j++)
            res.content[i][j] += operation(nb, matrix.content[i][j]);

    return res;
}


struct matrix_s matrix_prod(
    struct matrix_s matrix_a,
    struct matrix_s matrix_b
) {
    struct matrix_s res = {
        .nb_of_lines = matrix_a.nb_of_lines,
        .nb_of_columns = matrix_b.nb_of_columns,
        .content = NULL
    };

    if (matrix_a.nb_of_columns != matrix_b.nb_of_lines)
        return res;

    res = create_matrix(res.nb_of_lines, res.nb_of_columns);

    if (!is_matrix_init(res))
        return res;

    for (int i = 0; i < matrix_a.nb_of_lines; i++)
        for (int j = 0; j < matrix_b.nb_of_columns; j++) {
            matrix_content cij = 0;
            for (int k = 0; k < matrix_a.nb_of_columns; k++)
                cij += matrix_a.content[i][k] * matrix_b.content[k][j];

            res.content[i][j] = cij;
        }

    return res;
}


void print_matrix(struct matrix_s matrix) {
    for (int i = 0; i < matrix.nb_of_lines; i++) {
        printf("| ");
        for (int j = 0; j < matrix.nb_of_columns; j++)
            printf("\t%.2f", matrix.content[i][j]);
        printf("\t|\n");
    }
}


matrix_content matrix_trace(struct matrix_s matrix) {
    if (!is_matrix_square(matrix))
        return 0;

    matrix_content trace = 0;
    for (int i = 0; i < matrix.nb_of_lines; i++)
        trace += matrix.content[i][i];

    return trace;
}


struct matrix_s matrix_transpose(struct matrix_s matrix) {
    struct matrix_s res = {
        .nb_of_lines = matrix.nb_of_columns,
        .nb_of_columns = matrix.nb_of_lines,
        .content = NULL
    };

    res = create_matrix(res.nb_of_lines, res.nb_of_columns);

    if (!is_matrix_init(res))
        return res;

    for (int i = 0; i < matrix.nb_of_lines; i++)
        for (int j = 0; j < matrix.nb_of_columns; j++)
            res.content[j][i] = matrix.content[i][j];

    return res;
}


matrix_content matrix_det(struct matrix_s matrix) {
    if (!is_matrix_square(matrix))
        return 0;

    if (matrix.nb_of_lines == 1)
        return matrix.content[0][0];

    if (matrix.nb_of_lines == 2)
        return (
            matrix.content[0][0] * matrix.content[1][1] -
            matrix.content[1][0]* (matrix.content[0][1])
        );

    int det = 0;
    for (int i = 0; i < matrix.nb_of_lines; i++) {
        struct matrix_s under_mat =
            create_matrix(matrix.nb_of_lines - 1, matrix.nb_of_columns - 1);

        if (!is_matrix_init(under_mat))
            return 0;

        for (int k = 0; k < matrix.nb_of_lines - 1; k++) {
            for (int l = 0; l < i; l++)
                under_mat.content[k][l] = matrix.content[k+1][l];
            for (int l = i; l < matrix.nb_of_columns-1; l++)
                under_mat.content[k][l] = matrix.content[k+1][l+1];
        }

        det += mypow(-1,i) * matrix.content[0][i] * matrix_det(under_mat);

        delete_matrix(&under_mat);
    }

    return det;
}


struct matrix_s matrix_extract(
    struct matrix_s matrix,
    int line_begin,
    int col_begin,
    int line_end,
    int col_end
) {
    struct matrix_s res = {
        .nb_of_lines = 0,
        .nb_of_columns = 0,
        .content = NULL
    };

    if (line_begin < 0 || line_end < 0 || col_begin < 0 || col_end < 0)
        return res;

    if (line_begin > matrix.nb_of_lines || line_end > matrix.nb_of_lines)
        return res;

    if (col_begin > matrix.nb_of_columns || col_end > matrix.nb_of_columns)
        return res;

    if (line_end < line_begin || col_end < col_begin)
        return res;

    res = create_matrix(line_end - line_begin + 1, col_end - col_begin + 1);

    if (!is_matrix_init(res))
        return res;

    for (int i = line_begin; i <= line_end; i++)
        for (int j = col_begin; j <= col_end; j++)
            res.content[i][j] = matrix.content[i][j];

    return res;
}


struct matrix_s remove_line_from_matrix(
    struct matrix_s matrix,
    int line_to_remove
) {
    struct matrix_s res = {
        .nb_of_lines = matrix.nb_of_lines - 1,
        .nb_of_columns = matrix.nb_of_columns,
        .content = NULL
    };

    if (line_to_remove < 0 || line_to_remove > matrix.nb_of_lines)
        return res;

    res = create_matrix(matrix.nb_of_lines - 1, matrix.nb_of_columns);

    if (!is_matrix_init(res))
        return res;

    for (int i = 0; i < matrix.nb_of_lines; i++) {
        if (i == line_to_remove)
            continue;

        for (int j = 0; j < matrix.nb_of_columns; j++) {
            if (i > line_to_remove)
                res.content[i-1][j] = matrix.content[i][j];
            else
                res.content[i][j] = matrix.content[i][j];
        }
    }

    return res;
}


struct matrix_s remove_column_from_matrix(
    struct matrix_s matrix,
    int column_to_remove
) {
    struct matrix_s res = {
        .nb_of_lines = matrix.nb_of_lines,
        .nb_of_columns = matrix.nb_of_columns - 1,
        .content = NULL
    };

    if (column_to_remove < 0 || column_to_remove > matrix.nb_of_columns)
        return res;

    res = create_matrix(matrix.nb_of_lines, matrix.nb_of_columns - 1);

    if (!is_matrix_init(res))
        return res;

    for (int j = 0; j < matrix.nb_of_columns; j++) {
        if (j == column_to_remove)
            continue;
        for (int i = 0; i < matrix.nb_of_lines; i++)
            if (j > column_to_remove)
                res.content[i][j-1] = matrix.content[i][j];
            else
                res.content[i][j] = matrix.content[i][j];
    }

    return res;
}


struct matrix_s remove_line_and_column_from_matrix(
    struct matrix_s matrix,
    int line_to_remove,
    int column_to_remove
) {
    // Do not use other functions for performance reasons

    struct matrix_s res = {
        .nb_of_lines = matrix.nb_of_lines - 1,
        .nb_of_columns = matrix.nb_of_columns - 1,
        .content = NULL
    };

    if (
        column_to_remove < 0 || line_to_remove < 0 ||
        column_to_remove > matrix.nb_of_columns ||
        line_to_remove > matrix.nb_of_lines
    )
        return res;

    res = create_matrix(matrix.nb_of_lines - 1, matrix.nb_of_columns - 1);

    if (!is_matrix_init(res))
        return res;

    for (int i = 0; i < matrix.nb_of_lines; i++) {
        if (i == line_to_remove)
            continue;
        int line_nb_in_new = i > line_to_remove ? i-1 : i;

        for (int j = 0; j < matrix.nb_of_columns; j++) {
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


struct matrix_s co_matrix(struct matrix_s matrix) {
    struct matrix_s res = {
        .nb_of_lines = matrix.nb_of_lines,
        .nb_of_columns = matrix.nb_of_columns,
        .content = NULL
    };

    if (!is_matrix_square(matrix))
        return res;

    res = create_matrix(matrix.nb_of_lines, matrix.nb_of_columns);

    if (!is_matrix_init(res))
        return res;

    for (int i = 0; i < matrix.nb_of_lines; i++)
        for (int j = 0; j < matrix.nb_of_columns; j++) {
            struct matrix_s under_mat =
                remove_line_and_column_from_matrix(matrix, i, j);

            if (!is_matrix_init(under_mat)) {
                delete_matrix(&res);
                return res;
            }

            res.content[i][j] = mypow(-1,i+j) * matrix_det(under_mat);
            delete_matrix(&under_mat);
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
        .nb_of_lines = matrix.nb_of_columns,
        .nb_of_columns = matrix.nb_of_lines,
        .content = NULL
    };

    if (!is_matrix_square(matrix))
        return pre_res;

    struct matrix_s mat_inv_trans = co_matrix(matrix);

    if (!is_matrix_init(mat_inv_trans))
        return pre_res;

    float det = matrix_det(matrix);

    if ( det == 0 )
        return pre_res;

    matrix_content coef = 1 / det;

    pre_res = matrix_transpose(mat_inv_trans);
    delete_matrix(&mat_inv_trans);

    if (!is_matrix_init(pre_res))
        return pre_res;

    struct matrix_s res = matrix_all_terms_opp(pre_res, coef, simple_product);
    delete_matrix(&pre_res);

    return res;
}


struct matrix_s left_pseudo_inv(struct matrix_s matrix) {
    struct matrix_s res = {
        .nb_of_lines = matrix.nb_of_columns,
        .nb_of_columns = matrix.nb_of_lines,
        .content = NULL
    };

    if (matrix.nb_of_columns > matrix.nb_of_lines)
        return res;

    struct matrix_s mat_trans = matrix_transpose(matrix);
    if (!is_matrix_init(mat_trans))
        return res;

    struct matrix_s term_1 = matrix_prod(mat_trans, matrix);
    if (!is_matrix_init(term_1)) {
        delete_matrix(&mat_trans);
        return res;
    }

    struct matrix_s term_1_inv = matrix_invert(term_1);
    delete_matrix(&term_1);

    if (!is_matrix_init(term_1_inv)) {
        delete_matrix(&mat_trans);
        return res;
    }

    res = matrix_prod(term_1_inv, mat_trans);
    delete_matrix(&mat_trans);
    delete_matrix(&term_1_inv);

    return res;
}


struct matrix_s right_pseudo_inv(struct matrix_s matrix) {
    struct matrix_s res = {
        .nb_of_lines = matrix.nb_of_columns,
        .nb_of_columns = matrix.nb_of_lines,
        .content = NULL
    };

    if (matrix.nb_of_columns < matrix.nb_of_lines)
        return res;

    struct matrix_s mat_trans = matrix_transpose(matrix);
    if (!is_matrix_init(mat_trans))
        return res;

    struct matrix_s term_1 = matrix_prod(matrix, mat_trans);
    if (!is_matrix_init(term_1)) {
        delete_matrix(&mat_trans);
        return res;
    }

    struct matrix_s term_1_inv = matrix_invert(term_1);
    delete_matrix(&term_1);
    if (!is_matrix_init(term_1_inv)) {
        delete_matrix(&mat_trans);
        return res;
    }

    res = matrix_prod(mat_trans, term_1_inv);

    delete_matrix(&mat_trans);
    delete_matrix(&term_1_inv);

    return res;
}
