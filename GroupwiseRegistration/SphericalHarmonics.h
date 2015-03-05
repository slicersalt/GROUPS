#pragma once
#include "Geom.h"

class SphericalHarmonics
{
public:
	static void basis(int degree, float *p, float *Y)
	{
		// real spherical harmonics basis functions
		// polar coordinate
		double phi, theta;
		double dp[3] = {p[0], p[1], p[2]};
		Coordinate::cart2sph(dp, &phi, &theta);
		theta = PI / 2 - theta;  // convert to interval [0, PI]
		double *Pm = new double[degree + 1];

		// square root of 2
		double sqr2 = sqrt(2.0);

		for (int l = 0; l <= degree; l++)
		{
			// legendre part
			Series::legendre(l, cos(theta), Pm);
			float lconstant = sqrt((2 * l + 1) / (4 * PI));

			int center = (l + 1) * (l + 1) - l - 1;

			Y[center] = (float)(lconstant * Pm[0]);

			for (int m = 1; m <= l; m++)
			{
				double precoeff = lconstant * sqrt(1 / Series::factorial(l + m, l - m + 1));

				if (m % 2 == 1) precoeff = -precoeff;
				Y[center + m] = (float)(sqr2 * precoeff * Pm[m] * cos(m * phi));
				Y[center - m] = (float)(sqr2 * precoeff * Pm[m] * sin(m * phi));
			}
		}

		delete [] Pm;
	}
};

