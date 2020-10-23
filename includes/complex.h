#ifndef COMPLEX_H
#define COMPLEX_H

/* const double PI */
#define PI 3.141592653589793

/*
Input: complex number c, complex number c2
Result: complex number c = c * c2
*/
#define mul_complex_self(c, c2)       \
    do {                              \
        double ca = c.a;              \
        c.a = ca * c2.a - c.b * c2.b; \
        c.b = c.b * c2.a + ca * c2.b; \
    } while (0)

/* Complex number: z = a + bi */
typedef struct complex {
    double a; /* real part */
    double b; /* imaginary part */
} complex;

/*
Input: real part a, imaginary part b
Output: complex number z = a + bi
*/
complex create_complex(double a, double b);

/*
Input: complex number c1, complex number c2
Output: complex number z = c1 + c2
*/
complex add_complex(complex c1, complex c2);

/*
Input: complex number c1, complex number c2
Output: complex number z = c1 - c2
*/
complex sub_complex(complex c1, complex c2);

/*
Input: complex number c1, complex number c2
Output: complex number z = c1 * c2
*/
complex mul_complex(complex c1, complex c2);

/*
Input: real x
Output: complex number z = e^ix = cos(x) + i sin(x)
*/
complex euler_formula(double x);

#endif