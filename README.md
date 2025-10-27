# Matrix library

A simple implementation of matrixes in C.

This library contains simple operations, matrix product, transpose, traces,
etc...

If you need to change the type of the matrix content, simply change the
```c
typedef float matrix_content;
```
to
```c
typedef <your-type> matrix_content;
```

## Example

A very simple example that illustrate how the library shall be used:

```c
struct matrix_s my_mat = create_matrix(5, 5);
    matrix_content init[25] = {
            1. ,2. ,3. , 4. ,5.,
            6. ,7. ,8. ,9. ,10.,
            11.,12.,13.,14.,15.,
            16.,17.,18.,19.,20.,
            21.,22.,23.,24.,25.,
        };

    init_matrix(&my_mat, (matrix_content*) init);
    struct matrix_s my_mat2 = create_matrix(5, 5);
    init_matrix(&my_mat2, (matrix_content*) init);


    print_matrix(my_mat);

    struct matrix_s prod_res = matrix_prod(my_mat, my_mat2);
    print_matrix(prod_res);

    // Do not forget to delete the matrix to free memory.
    delete_matrix(&prod_res);
    delete_matrix(&my_mat2);
    delete_matrix(&my_mat);
```
