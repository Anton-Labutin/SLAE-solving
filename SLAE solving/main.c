//
//  SLAE solving
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

enum { INPUT_FILE = 1, INPUT_FORMULAS = 2 };
typedef double MX_TYPE;

int n = 25, m = 10;
long long int addG = 0, multG = 0, divG = 0;
long long int iter_cnt = 0, min_iter_cnt = -1, max_iter_cnt = 20000;
double fast_omega = 0.0;
long long int addR = 0, multR = 0, divR = 0;
long long int addRfast = 0, multRfast = 0, divRfast = 0;


// выделяем память под матрицу размера row_cnt на col_cnt
MX_TYPE**
init(int row_cnt, int col_cnt)
{
    MX_TYPE* *matrix = (MX_TYPE**) malloc(row_cnt * sizeof(MX_TYPE*));

    for (int i = 0; i < row_cnt; ++i) {
        matrix[i] = (MX_TYPE*) malloc(col_cnt * sizeof(MX_TYPE));
    }

    return matrix;
}


// освобождаем память под матрицу
void
free_matrix_space(MX_TYPE* *matrix, int N)
{
    for (int i = 0; i < N; ++i) {
        free(matrix[i]);
    }

    if (N > 1) {
        free(matrix);
    }
}

//печать матрицы
void
print_matrix(MX_TYPE* *matrix, int row_cnt, int col_cnt)
{
    for (int i = 0; i < row_cnt; ++i) {
        for (int j = 0; j < col_cnt; ++j) {
            printf("%.10g ", matrix[i][j]);
        }

        printf("\n");
    }
}

// копировать числовые значения матрицы old_matrix в матрицу new_matrix
void
copy_matrix_values(MX_TYPE* *old_matrix, int row_cnt, int col_cnt, MX_TYPE* *new_matrix)
{
    for(int i = 0; i < row_cnt; ++i) {
        for (int j = 0; j < col_cnt; ++j) {
            new_matrix[i][j] = old_matrix[i][j];
        }
    }
}

// сделать матрицу единичной
void
init_identity_matrix(MX_TYPE* *matrix, int mxsize)
{
    for (int i = 0; i < mxsize; ++i) {
        for (int j = 0; j < mxsize; ++j) {
            if (j != i) {
                matrix[i][j] = 0;
            } else {
                matrix[i][j] = 1;
            }
        }
    }
}

// ищем первый после строки с номером cur_row ненулевой элемент и возвращаем номер строки
MX_TYPE
get_leading_matrix_el_idx(MX_TYPE* *matrix, int mxsize, int cur_row)
{
    int index = cur_row;

    for (int k = cur_row + 1; k < mxsize; ++k) {
        if (matrix[k][cur_row] != 0) {
            index = k;
            break;

        }
    }

    return index;
}

// меняем строки с номерами row1 и row2 местами
void
swap_matrix_rows(MX_TYPE* *matrix, int col_cnt, int row1, int row2)
{
    for (int k = 0; k < col_cnt; ++k) {
        MX_TYPE temp = matrix[row1][k];
        matrix[row1][k] = matrix[row2][k];
        matrix[row2][k] = temp;
    }
}

// ищем ведущий элемент и меняем строки местами
int
Gauss_transform_matrix(MX_TYPE* *matrix, int mxsize, int cur_row)
{
        // ищем первый ненулевой элемент
    int lead_idx = get_leading_matrix_el_idx(matrix, mxsize, cur_row);

        //меняем строки с номерами cur_row и lead_idx местами
    swap_matrix_rows(matrix, mxsize, cur_row, lead_idx);

    return lead_idx;
}

//нахождение обратной матрицы методом Гаусса
void
inverse_matrix(MX_TYPE* *matr, int mxsize, MX_TYPE* *invmatrix)
{
    MX_TYPE* *matrix = init(mxsize, mxsize);
    copy_matrix_values(matr, mxsize, mxsize, matrix);

    // сделаем обратную матрицу единичной
    init_identity_matrix(invmatrix, mxsize);

    // прямой ход методом Гаусса
    for (int i = 0; i < mxsize - 1; ++i) {
        for (int j = i + 1; j < mxsize; ++j) {
                //если ведущий элемент равен 0
            if (matrix[i][i] == 0) {
                // ищем ненулевой элемент
                int lead_idx = get_leading_matrix_el_idx(matrix, mxsize, 0);

                //меняем строки местами
                swap_matrix_rows(matrix, mxsize, i, lead_idx);
                swap_matrix_rows(invmatrix, mxsize, i, lead_idx);
            }

            //вычисление коэффициента
            double coef = matrix[j][i] / matrix[i][i];

            //вычитание из j-той строчки i-той строки, умноженной на коэффициент
            for (int k = 0; k < mxsize; ++k) {
                matrix[j][k] -= coef * matrix[i][k];
                invmatrix[j][k] -= coef * invmatrix[i][k];
            }
        }
    }

    //обратный ход вычисления элементов обратной матрицы
    for (int i = mxsize - 1; i >= 0; --i) {
        for (int j = i; j > 0; --j) {
            double coef = matrix[j - 1][i] / matrix[i][i];

            for (int k = mxsize - 1; k >= 0; --k) {
                invmatrix[j - 1][k] -= invmatrix[i][k] * coef;
            }
        }
    }

    for (int i = 0; i < mxsize; i++) {
        for (int j = 0; j < mxsize; j++) {
            invmatrix[i][j] /= matrix[i][i];
        }
    }
}

//прямой ход метода Гаусса
double
Gauss_forward_elimination(MX_TYPE* *matrix, MX_TYPE mxsize, MX_TYPE *f)
{
    double matrix_det = 1;
    int sgn = 1; // sign of the determinant

    for (int i = 0; i < mxsize - 1; ++i) { // current row
        for (int j = i + 1; j < mxsize; ++j) { // next rows
                                               //если ведущий элемент равен 0
            if (matrix[i][i] == 0) {
                int leading_el_row = Gauss_transform_matrix(matrix, mxsize, i);

                if (matrix[leading_el_row][i] != 0.0) {
                    MX_TYPE temp = f[i];
                    f[i] = f[leading_el_row];
                    f[leading_el_row] = temp;

                    sgn = -sgn;
                } else {
                    break;
                }
            }

                //вычисление коэффициента
            double coef = matrix[j][i] / matrix[i][i];
            ++divG;

                //вычитание из j-той строчки i-той строки, умноженной на коэффициент
            for (int k = i; k < mxsize; ++k) {
                matrix[j][k] -= coef * matrix[i][k];

                ++addG;
                ++multG;
            }

            f[j] -= f[i] * coef;

            ++addG;
            ++multG;
        }
    }

    for (int i = 0; i < mxsize; ++i) {
        matrix_det *= matrix[i][i];
    }
    matrix_det *= sgn;

    return matrix_det;
}


// обратный ход метода Гаусса
void
Gauss_back_substitution(MX_TYPE* *matrix, int mxsize, MX_TYPE *f, double *res)
{
        //обратный ход метода Гаусса
    if (matrix[mxsize - 1][mxsize - 1] == 0) {
        res[mxsize - 1] = 0;
    } else {
        res[mxsize - 1] = f[mxsize - 1] / matrix[mxsize - 1][mxsize - 1];
    }

    ++divG;

    for (int i = mxsize - 2; i >= 0; --i) {
        MX_TYPE temp = f[i];

        for (int j = mxsize - 1; j > i; --j) {
            temp -= matrix[i][j] * res[j];

            ++addG;
            ++multG;
        }

        if (temp == 0) {
            res[i] = 0.0;
        } else {
            res[i] = (double) temp / matrix[i][i];
        }

        ++divG;
    }
}

//Метод Гаусса
double
Gauss(MX_TYPE* *matrix, int mxsize, MX_TYPE *f, double *res)
{
    MX_TYPE* *matrix_copy = init(mxsize, mxsize);
    copy_matrix_values(matrix, mxsize, mxsize, matrix_copy);

    MX_TYPE* *f_copy = init(1, mxsize);
    copy_matrix_values(&f, 1, mxsize, f_copy);

    double matrix_det = Gauss_forward_elimination(matrix_copy, mxsize, *f_copy);
    if (matrix_det != 0) {
        Gauss_back_substitution(matrix_copy, mxsize, *f_copy, res);
    }

    free_matrix_space(matrix_copy, mxsize);
    free_matrix_space(f_copy, 1);

    return matrix_det;
}

// поиск максимального по модулю элемента в строке с номером cur_row и возврат номера столбца
int
get_idx_row_max(MX_TYPE* *matrix, int col_cnt, int cur_row)
{
    MX_TYPE max = fabs(matrix[cur_row][cur_row]);

    int col_idx = cur_row;

    for (int k = cur_row + 1; k < col_cnt; ++k) {
        if (fabs(matrix[cur_row][k]) > max) {
            col_idx = k;
            max = fabs(matrix[cur_row][k]);
        }
    }

    return col_idx;
}

// переставить местами столбцы с номерами col1 и col2
void
swap_matrix_columns(MX_TYPE* *matrix, int row_cnt, int col1, int col2)
{
    for (int k = 0; k < row_cnt; ++k) {
        MX_TYPE temp = matrix[k][col1];
        matrix[k][col1] = matrix[k][col2];
        matrix[k][col2] = temp;
    }
}

//Метод Гаусса с выбором главного элемента
void
Gauss_select(MX_TYPE* *matrix, int mxsize, MX_TYPE *f, double *res)
{
    MX_TYPE* *matrix_copy = init(mxsize, mxsize);
    copy_matrix_values(matrix, mxsize, mxsize, matrix_copy);

    MX_TYPE* *f_copy = init(1, mxsize);
    copy_matrix_values(&f, 1, mxsize, f_copy);

    double temp;

    int *res_num = malloc(mxsize * sizeof(int));
    for (int i = 0; i < mxsize; ++i) {
        res_num[i] = i;
    }

        //прямой ход метода Гаусса
    for (int i = 0; i < mxsize - 1; ++i) {
            //поиск максимального по модулю элемента в строке
        int index = get_idx_row_max(matrix_copy, mxsize, i);

            //меняем столбцы с номерами i и index местами
        if (i != index) {
            swap_matrix_columns(matrix_copy, mxsize, i, index);
        }

        int buf = res_num[i];
        res_num[i] = res_num[index];
        res_num[index] = buf;

        for (int j = i + 1; j < mxsize; j++) {
                //вычисление коэффициента
            temp = matrix_copy[j][i] / matrix_copy[i][i];

                //вычитание из j-той строки i-той строки, умноженной на коэффициент
            for (int k = i; k < mxsize; k++) {
                matrix_copy[j][k] -= temp * matrix_copy[i][k];
            }

            (*f_copy)[j] -= temp * (*f_copy)[i];
        }
    }

    //обратный ход метода Гаусса
    Gauss_back_substitution(matrix_copy, mxsize, *f_copy, res);

    for (int i = 0; i < mxsize; ++i) {
        if (res_num[i] != i) {
            temp = res[i];
            res[i] = res[res_num[i]];
            res[res_num[i]] = temp;

            int buf = res_num[i];
            res_num[i] = res_num[buf];
            res_num[buf] = buf;
        }
    }

    free_matrix_space(matrix_copy, mxsize);
    free_matrix_space(f_copy, 1);
}

//создание матрицы по примеру из приложения 2
void
create_matrix(MX_TYPE* *matrix, MX_TYPE *f)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                matrix[i][j] = (i + j) / (m + n);
            } else {
                matrix[i][j] = n + pow(m, 2) + j/m + i/n;
            }
        }

        f[i] = i * i - n;
    }
}

// Вычисление нормы разности 2х векторов
long double
vector_diff(double *prev_x, double *curr_x, int mxsize)
{
    long double diff = 0.0;

    for (int i = 0; i < mxsize; ++i) {
        diff += pow((curr_x[i] - prev_x[i]), 2);
    }

    return sqrt(diff);
}

// вычисление очередной итерации в методе верхней релаксации
long double
SOR_next_iter(MX_TYPE* *matrix, int mxsize, MX_TYPE *f, double *prev_x, double *curr_x, double omega)
{
    double sum;
    int i, j;

    for (i = 0; i < mxsize; ++i) {
        sum = 0.0;

        for (j = 0; j < i; ++j) {
            sum += -(matrix[i][j] * curr_x[j]);

            ++addR;
            ++multR;

        }

        for (j = i; j < mxsize; ++j) {
            sum += -(matrix[i][j] * prev_x[j]);
            ++addR;
            ++multR;
        }

        sum += f[i];
        sum *= omega;
        sum /= matrix[i][i];
        curr_x[i] = prev_x[i] + sum;

        addR += 2;
        multR += 2;
        ++divR;
    }

    return vector_diff(prev_x, curr_x, mxsize);
}

// x имеет конечную норму?
int
is_finite(double *x, int x_size)
{
    for (int i = 0; i < x_size; ++i) {
        if (!isfinite(x[i])) {
                // x[i] - бесконечность или не число
            return 0;
        }
    }

    return 1;
}

// метод верхней релаксации
int
SOR(MX_TYPE* *matrix, int mxsize, MX_TYPE *f, double *fast_x, double eps)
{
    int SOR_converges = 0;
    int converge; // флаг для сходимости

    double omega;
    double conv_max_omega = 2.0; // максимальное значение омега для сходимости метода
    double start_omega = 0.1;
    double omega_step = 0.1;

    double *prev_x = malloc(mxsize * sizeof(double));
    double *curr_x = malloc(mxsize * sizeof(double));
    double *tmp = NULL;

    for (omega = start_omega; omega <= conv_max_omega; omega +=omega_step){ //смотрим для 0 < omega <= 2 сходимость
        for (int i = 0; i < mxsize; ++i) {
            prev_x[i] = 0.0; // принимаем за начальное приближение нулевой вектор
        }

        addR = 0;
        multR = 0;
        divR = 0;

        converge = 1;
        iter_cnt = 0;

        while (SOR_next_iter(matrix, mxsize, f, prev_x, curr_x, omega) > eps) {
            ++iter_cnt;
            if (iter_cnt > max_iter_cnt) {
                converge = 0;
                break;
            }

            tmp = prev_x;
            prev_x = curr_x;
            curr_x = tmp;
        }

        if (converge) {
            converge = is_finite(curr_x, mxsize);
        }

        if (converge) {
            printf("omega = %f %lld iterations\n", omega, iter_cnt);
            if ((omega < start_omega * 1.1) || (iter_cnt < min_iter_cnt)) {
                min_iter_cnt = iter_cnt;
                SOR_converges = 1;

                for (int i = 0; i < mxsize; ++i) {
                    fast_x[i] = curr_x[i];
                }

                addRfast = addR;
                multRfast = multR;
                divRfast = divR;

                fast_omega = omega;
            }
        } else {
            printf("omega = %f method diverges\n", omega);
        }
    }

    return SOR_converges;
}

int
main(int argc, char *argv[])
{
    int i, j;

    int var;
    do {
        printf("How do you want to input data?\n");
        printf("Using input file - input %d\n", INPUT_FILE);
        printf("Using formulas - input %d\n", INPUT_FORMULAS);
        scanf("%d", &var);
    } while (var != INPUT_FILE && var != INPUT_FORMULAS);

    int N; // размер матрицы
    if (var == INPUT_FILE) {
        do {
            printf("Input matix size:\n");
            scanf("%d", &N);
        } while (N < 1);
    } else {
        N = n;
    }

    MX_TYPE* *matrix = init(N, N); // исходная матрица
    MX_TYPE* *invmatrix = init(N, N); // обратная матрица

    MX_TYPE *f = malloc(N * sizeof(MX_TYPE)); // вектор значений
    MX_TYPE *res = malloc(N * sizeof(MX_TYPE)); // вектор решений

    if (var == INPUT_FILE) {
        printf("INPUT FORMAT:\n");
        printf("a11 a12 ... a1n f1\na21 a22 ... a2n f2\n... ... ... ... ...\nan1 an2 ... ann fn\n\n");
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                scanf("%lf", &matrix[i][j]);
            }

            scanf("%lf", &f[i]);
        }
    } else {
        if (var == INPUT_FORMULAS) {
            create_matrix(matrix, f);
        }
    }

    double matrix_det = Gauss(matrix, N, f, res);
    printf("Determinant of matrix: %.10g\n", matrix_det);

           if (matrix_det != 0.0) {
        inverse_matrix(matrix, N, invmatrix); // вычисляем обратную матрицу
        printf("Inverse matrix:\n");
        print_matrix(invmatrix, N, N);
        free_matrix_space(invmatrix, N);

        printf("Roots of SLAE — Gauss' method:\n"); // решаем систему методом Гаусса
        for (int i = 0; i < N; ++i) {
            printf("x%d = %f\n", i + 1, res[i]);
        }

        long long int add_cnt = addG, mult_cnt = multG, div_cnt = divG;

        Gauss_select(matrix, N, f, res);
        printf("Roots of SLAE — Gauss' method with selection:\n"); // решаем систему методом Гаусса с выбором главного элемента
        for (int i = 0; i < N; ++i) {
            printf("x%d = %f\n", i + 1, res[i]);
        }

        printf("Operations needed for Gauss' method and Gauss' method with selection:\n");
        printf(" addition %lld\n multiplication %lld\n division %lld\n", add_cnt, mult_cnt, div_cnt);

            // Метод верхней релаксации
        double eps;
        printf("Input epsilon:\n");
        do {
            scanf("%lf", &eps);
        } while (eps <= 0.0);

        int converge = SOR(matrix, N, f, res, eps);
        if (converge) {
            printf("method converges the most quickly with omega = %f\n", fast_omega);
            printf("Roots of SLAE — successive over-relaxation method:\n");
            for (i = 0; i < N; ++i) {
                printf("x%d = %f\n", i + 1, res[i]);
            }
            printf("Operations needed for successive over-relaxation method: \n");

            printf("%lld iterations\n\%lld addition\n%lld multiplication\n%lld division\n", min_iter_cnt, addRfast, multRfast, divRfast);
        } else {
            printf("%s\n", "successive over-relaxation method diverges");
        }
           } else {
               printf("The matrix is degenerate! There is no the one solution!\n");
           }

    free_matrix_space(matrix, N);
    free_matrix_space(&f, 1);
    free_matrix_space(&res, 1);

    return 0;
}
