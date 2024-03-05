//
// Created by Artemii Kazakov, ITMO.
//

#include "errors.h"
#include "return_codes.h"

#include <malloc.h>
#include <math.h>
#include <stdio.h>

// ------ CHANGE TYPE OF COMPUTING: 1 - float 2 - double 3 - long double --------

#define COMPUTE_T_NUM 1
#define ACCURACY 65

#if COMPUTE_T_NUM == 1
#define COMPUTE_T float
#define EPS 1e-6
#define EPS_OUT 1e-2
#define FABS(NUMBER) fabsf(NUMBER)
#define SQRT(NUMBER) sqrtf(NUMBER)
#elif COMPUTE_T_NUM == 2
#define COMPUTE_T double
#define EPS 1e-6
#define EPS_OUT 1e-4
#define FABS(NUMBER) fabs(NUMBER)
#define SQRT(NUMBER) sqrt(NUMBER)
#elif COMPUTE_T_NUM == 3
#define COMPUTE_T long double
#define EPS 1e-7
#define EPS_OUT 1e-4
#define FABS(NUMBER) fabsl(NUMBER)
#define SQRT(NUMBER) sqrtl(NUMBER)
#endif

// ------ PROTOTYPES OF FUNCTIONS ----------

// - MAJOR -

int check_qr_end(int n, COMPUTE_T **mat, int *result, int *retry, COMPUTE_T *accur, int *cycle);

void qr_iteration(int n, COMPUTE_T ***mat, COMPUTE_T *x, COMPUTE_T *u, COMPUTE_T **p, COMPUTE_T **tmp_matrix, COMPUTE_T **q);

void compute_u(COMPUTE_T *x, int dm, COMPUTE_T **res);

void compute_p(int n, COMPUTE_T *u, int dm, COMPUTE_T ***res);

int read_matrix_from_file(FILE *, int *, COMPUTE_T ***);

void write_matrix_output_file(FILE *output, COMPUTE_T **matrix, int n);

// - MINOR -

COMPUTE_T scalar_product(int, const COMPUTE_T *, const COMPUTE_T *);

void multiply_matrix(int, COMPUTE_T **, COMPUTE_T **, COMPUTE_T ***);

void get_sub_matrix(int n, COMPUTE_T **mat, int k, COMPUTE_T **res);

// - UTILS -

void copy_matrix(int n, COMPUTE_T ***, COMPUTE_T **);

int allocate_vector(int, COMPUTE_T **);

int allocate_matrix(int, COMPUTE_T ***);

void free_matrix(int, COMPUTE_T **);

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Number of arguments is %d, but must be 2 (... input_file output_file).", argc - 1);
		return ERROR_PARAMETER_INVALID;
	}

	FILE *input;
	if ((input = fopen(argv[1], "r")) == NULL)
	{
		ERROR_MESSAGE_CANNOT_OPEN_FILE(argv[1], "r", ERROR_CANNOT_OPEN_FILE)
	}

	int n;
	COMPUTE_T **matrix;
	CHECK_ERROR_WITH_FREE(SUCCESS, read_matrix_from_file(input, &n, &matrix), input_read, fclose(input);)

	fclose(input);

	if (n == 1)
	{
		goto one;
	}

	int context = 0;
	int error_qr = SUCCESS;
	COMPUTE_T **q;
	ERROR_GOTO_WITH_CONTEXT(SUCCESS, allocate_matrix(n, &q), error_qr, ERROR_OUT_OF_MEMORY, error_qr_label, context)
	COMPUTE_T *x;
	ERROR_GOTO_WITH_CONTEXT(SUCCESS, allocate_vector(n, &x), error_qr, ERROR_OUT_OF_MEMORY, error_qr_label, context)
	COMPUTE_T **tmp_matrix;
	ERROR_GOTO_WITH_CONTEXT(SUCCESS, allocate_matrix(n, &tmp_matrix), error_qr, ERROR_OUT_OF_MEMORY, error_qr_label, context)
	COMPUTE_T *u;
	ERROR_GOTO_WITH_CONTEXT(SUCCESS, allocate_vector(n, &u), error_qr, ERROR_OUT_OF_MEMORY, error_qr_label, context)
	COMPUTE_T **p;
	ERROR_GOTO_WITH_CONTEXT(SUCCESS, allocate_matrix(n, &p), error_qr, ERROR_OUT_OF_MEMORY, error_qr_label, context)

	int result = 0;
	COMPUTE_T accur = 0;
	int retry = 0;
	int cycle = 0;

	while (!result)
	{
		qr_iteration(n, &matrix, x, u, p, tmp_matrix, q);
		int check_check_qr_end = check_qr_end(n, matrix, &result, &retry, &accur, &cycle);
		if (check_check_qr_end != SUCCESS)
		{
			error_qr = check_check_qr_end;
			goto error_qr_label;
		}
	}

error_qr_label:
	for (int i = 0; i < context; i++)
	{
		switch (i)
		{
		case 0:
			free_matrix(n, q);
			break;
		case 1:
			free(x);
			break;
		case 2:
			free_matrix(n, tmp_matrix);
			break;
		case 3:
			free(u);
			break;
		case 4:
			free_matrix(n, p);
			break;
		}
	}
	if (error_qr != SUCCESS)
	{
		free_matrix(n, matrix);
		return error_qr;
	}

one:;

	FILE *output;
	if ((output = fopen(argv[2], "w")) == NULL)
	{
		free_matrix(n, matrix);
		ERROR_MESSAGE_CANNOT_OPEN_FILE(argv[2], "w", ERROR_CANNOT_OPEN_FILE)
	}

	write_matrix_output_file(output, matrix, n);

	fclose(output);

	free_matrix(n, matrix);
	return SUCCESS;
}

// - MAJOR -

int check_qr_end(int n, COMPUTE_T **mat, int *result, int *retry, COMPUTE_T *accur, int *cycle)
{
	int is_zero = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (mat[i][j] != mat[i][j])
			{
				fprintf(stderr, "One ore more numbers in computing equals nan.\n");
				return ERROR_UNKNOWN;
			}
			if (FABS(mat[i][j]) < EPS)
			{
				is_zero++;
			}
		}
	}
	COMPUTE_T proc = (COMPUTE_T)is_zero / (COMPUTE_T)(n * n) * 100;
	if (proc - ACCURACY >= EPS)
	{
		*result = 1;
	}
	else
	{
		*result = 0;
	}
	if (proc - (*accur) <= EPS)
	{
		if (*retry > n * n * n * n)
		{
			*result = 1;
		}
		else
		{
			(*retry)++;
		}
	}
	else
	{
		if (*cycle > n * n)
		{
			*result = 1;
		}
		*retry = 0;
		(*cycle)++;
	}
	*accur = proc;

	return SUCCESS;
}

void qr_iteration(int n, COMPUTE_T ***mat, COMPUTE_T *x, COMPUTE_T *u, COMPUTE_T **p, COMPUTE_T **tmp_matrix, COMPUTE_T **q)
{
	for (int i = 0; i < n - 1; i++)
	{
		get_sub_matrix(n, *mat, i, &x);

		int dm = n - i;

		compute_u(x, dm, &u);

		compute_p(n, u, dm, &p);

		multiply_matrix(n, p, *mat, &tmp_matrix);
		copy_matrix(n, mat, tmp_matrix);
		if (i == 0)
		{
			copy_matrix(n, &q, p);
		}
		else
		{
			multiply_matrix(n, q, p, &tmp_matrix);
			copy_matrix(n, &q, tmp_matrix);
		}
	}

	multiply_matrix(n, *mat, q, &tmp_matrix);
	copy_matrix(n, mat, tmp_matrix);
}

void compute_u(COMPUTE_T *x, int dm, COMPUTE_T **res)
{
	COMPUTE_T scalar = SQRT(scalar_product(dm, x, x));
	for (int i = 0; i < dm; i++)
	{
		(*res)[i] = x[i];
	}
	(*res)[0] = x[0] - scalar;
}

void compute_p(int n, COMPUTE_T *u, int dm, COMPUTE_T ***res)
{
	COMPUTE_T scalar = scalar_product(dm, u, u);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			(*res)[i][j] = (i == j ? 1.0f : 0.0f);
		}
	}

	for (int i = 0; i < dm; i++)
	{
		for (int j = 0; j < dm; j++)
		{
			(*res)[i + (n - dm)][j + (n - dm)] -= 2 * (u[i] * u[j]) / scalar;
		}
	}
}

int read_matrix_from_file(FILE *input, int *n, COMPUTE_T ***matrix)
{
	if (!fscanf(input, "%d\n", n))
	{
		fprintf(stderr, "Invalid input file format.");
		return ERROR_CANNOT_OPEN_FILE;
	}
	CHECK_ERROR(SUCCESS, allocate_matrix(*n, matrix), res_al_matrix)
	for (int i = 0; i < *n; i++)
	{
		for (int j = 0; j < *n; j++)
		{
#if COMPUTE_T_NUM == 1
			if (!fscanf(input, "%f", &((*matrix)[i][j])))
#elif COMPUTE_T_NUM == 2
			if (!fscanf(input, "%lf", &((*matrix)[i][j])))
#elif COMPUTE_T_NUM == 3
			if (!fscanf(input, "%Lf", &((*matrix)[i][j])))
#endif
			{
				free_matrix(*n, *matrix);
				ERROR_MESSAGE_CANNOT_OPEN_FILE("input", "r", ERROR_CANNOT_OPEN_FILE)
			}
		}
	}
	return SUCCESS;
}

void write_matrix_output_file(FILE *output, COMPUTE_T **matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		if (i + 1 < n && FABS(matrix[i][i + 1]) > EPS_OUT && FABS(matrix[i + 1][i]) > EPS_OUT &&
			FABS(matrix[i][i + 1]) / FABS(matrix[i + 1][i]) <= 1.5)
		{
			COMPUTE_T tmp = matrix[i][i] + matrix[i + 1][i + 1];
			COMPUTE_T res = FABS((matrix[i][i] * matrix[i + 1][i + 1] - matrix[i][i + 1] * matrix[i + 1][i]) / (tmp ? tmp : 1));

			fprintf(
				output,
#if COMPUTE_T_NUM == 1
				"%g +%gi\n%g -%gi\n",
#elif COMPUTE_T_NUM == 2
				"%lg +%lgi\n%lg -%lgi\n",
#elif COMPUTE_T_NUM == 3
				"%Lg +%Lgi\n%Lg -%Lgi\n",
#endif
				tmp / 2,
				res,
				tmp / 2,
				res);

			i++;
			continue;
		}
		fprintf(output,
#if COMPUTE_T_NUM == 1
				"%g\n",
#elif COMPUTE_T_NUM == 2
				"%lg\n",
#elif COMPUTE_T_NUM == 3
				"%Lg\n",
#endif
				matrix[i][i]);
	}
}

// - MINOR -

COMPUTE_T scalar_product(int n, const COMPUTE_T *vector1, const COMPUTE_T *vector2)
{
	COMPUTE_T scalar_result = 0;
	for (int i = 0; i < n; i++)
	{
		scalar_result += (vector1[i] * vector2[i]);
	}
	return scalar_result;
}

void multiply_matrix(int n, COMPUTE_T **mat1, COMPUTE_T **mat2, COMPUTE_T ***result)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			(*result)[i][j] = 0;
			for (int k = 0; k < n; k++)
			{
				(*result)[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
	}
}

void get_sub_matrix(int n, COMPUTE_T **mat, int k, COMPUTE_T **res)
{
	for (int i = k; i < n; i++)
	{
		(*res)[i - k] = mat[i][k];
	}
}

// - UTILS -

void copy_matrix(int n, COMPUTE_T ***matrix, COMPUTE_T **copy_matrix)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			(*matrix)[i][j] = copy_matrix[i][j];
		}
	}
}

int allocate_vector(int n, COMPUTE_T **vector)
{
	*vector = malloc(sizeof(COMPUTE_T) * n);
	if (*vector == NULL)
	{
		ERROR_MESSAGE_OUT_OF_MEMORY("vector", ERROR_OUT_OF_MEMORY)
	}
	return SUCCESS;
}

int allocate_matrix(int n, COMPUTE_T ***matrix)
{
	*matrix = malloc(sizeof(COMPUTE_T *) * n);
	if (*matrix == NULL)
	{
		ERROR_MESSAGE_OUT_OF_MEMORY("matrix", ERROR_OUT_OF_MEMORY)
	}
	for (int i = 0; i < n; i++)
	{
		if (allocate_vector(n, &((*matrix)[i])) != SUCCESS)
		{
			free_matrix(i, *matrix);
			ERROR_MESSAGE_OUT_OF_MEMORY("vector in matrix", ERROR_OUT_OF_MEMORY)
		}
	}
	return SUCCESS;
}

void free_matrix(int n, COMPUTE_T **matrix)
{
	for (int i = 0; i < n; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}
