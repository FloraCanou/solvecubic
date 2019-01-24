/* Copyright 2019 Flora Canou | V. C0-1.0.0 | Cubic Equation Solver
 * This Source Code Form is licensed under the Mozilla Public License, v. 2.0. 
 * If a copy of the MPL was not distributed with this file, you can obtain one at https://mozilla.org/MPL/2.0/. 
 * The program solves cubic equation in trigonometric/hyperbolic method. 
 */
 
 #include <math.h>
 #include <float.h>
 #define pi	3.14159265359
 
int cubic_equation_solver (double a, double b, double c, double d, double *x)
{
	int n = 0; //number of real solutions
	if (fabs(a) < DBL_EPSILON)
	{
		if (fabs(b) < DBL_EPSILON)
		{
			if (fabs(c) < DBL_EPSILON)
			{
				if (fabs(d) < DBL_EPSILON) // zero
					n = -1; //infinitely many solutions
				else // non-zero constant
					n = 0;
			}
			else // linear
			{
				n = 1;
				x[0] = -d/c;
			}
		}
		else //quadratic
		{
			double discriminant = c*c - 4*b*d;
			if (discriminant >= 0)
			{
				n = 2;
				x[0] = (-c + sqrt (discriminant)) / (2*b);
				x[1] = (-c - sqrt (discriminant)) / (2*b);
			}
			else
				n = 0;
		}
	}
	else //cubic
	{
		double p, q, discriminant, u;
		p = (3*a*c - b*b) / (3*a*a);
		q = (2*b*b*b - 9*a*b*c + 27*a*a*d) / (27*a*a*a);
		discriminant = 4*p*p*p + 27*q*q;
		if (discriminant <= 0) // three real solutions
		{
			n = 3;
			if (fabs(p) < DBL_EPSILON)
			{
				x[0] = -b / (3*a);
				x[1] = x[0];
				x[2] = x[0];
			}
			else
			{
				u = sqrt (-4*p/3); 
				x[0] = u * cos (acos (-4*q/(u*u*u)) / 3) - b / (3*a);
				x[1] = u * cos ((acos (-4*q/(u*u*u)) - 2*pi) / 3) - b / (3*a);
				x[2] = u * cos ((acos (-4*q/(u*u*u)) + 2*pi) / 3) - b / (3*a);
			}
		}
		else // one real solution
		{
			n = 1;
			if (fabs(p) < DBL_EPSILON)
				x[0] = -cbrt(q) - b / (3*a);
			else if (p < 0)
			{
				u = copysign (sqrt (-4*p/3), -q);
				x[0] = u * cosh (acosh (-4*q/(u*u*u)) / 3) - b / (3*a);
			}
			else
			{
				u = -sqrt (4*p/3);
				x[0] = u * sinh (asinh (-4*q/(u*u*u)) / 3) - b / (3*a);
			}
		}
	}
	return n;
}

