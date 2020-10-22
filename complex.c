#include "./includes/complex.h"

#include <math.h>

/*
Input: real part a, imaginary part b
Output: complex number z = a + bi
*/
complex create_complex(double a, double b) {
    complex z = {a, b};
    return z;
}

/*
Input: complex number c1, complex number c2
Output: complex number z = c1 + c2
*/
complex add_complex(complex c1, complex c2) {
    complex z = {.a = c1.a + c2.a, .b = c1.b + c2.b};
    return z;
}

/*
Input: complex number c1, complex number c2
Output: complex number z = c1 - c2
*/
complex sub_complex(complex c1, complex c2) {
    complex z = {.a = c1.a - c2.a, .b = c1.b - c2.b};
    return z;
}

/*
Input: complex number c1, complex number c2
Output: complex number z = c1 * c2
*/
complex mul_complex(complex c1, complex c2) {
    complex z;
    z.a = c1.a * c2.a - c1.b * c2.b;
    z.b = c1.b * c2.a + c1.a * c2.b;
    return z;
}

/*
Input: real x
Output: complex number z = e^ix = cos(x) + i sin(x)
*/
complex euler_formula(double x) {
    complex z = {.a = cos(x), .b = sin(x)};
    return z;
}