/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2016  <copyright holder> <email>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "test_common.h"

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

const double values[y_size][x_size] = 
{
	{ 1, 1, 1, 1, 1 },
	{ 2, 2, 2, 2, 2 },
	{ 3, 3, 3, 3, 3 },
	{ 4, 4, 4, 4, 4 },
	{ 5, 6, 7, 8, 9 }
};

const double errs_u[y_size][x_size] = 
{
	{ 0.1, 0.1, 0.1, 0.1, 0.1 },
	{ 0.2, 0.2, 0.2, 0.2, 0.2 },
	{ 0.3, 0.3, 0.3, 0.3, 0.3 },
	{ 0.4, 0.4, 0.4, 0.4, 0.4 },
	{ 0.5, 0.6, 0.7, 0.8, 0.9 }
};

const double errs_l[y_size][x_size] = 
{
	{ 1.1, 1.1, 1.1, 1.1, 1.1 },
	{ 1.2, 1.2, 1.2, 1.2, 1.2 },
	{ 1.3, 1.3, 1.3, 1.3, 1.3 },
	{ 1.4, 1.4, 1.4, 1.4, 1.4 },
	{ 1.5, 1.6, 1.7, 1.8, 1.9 }
};

void print_default(int selection)
{
	for (uint y = 0; y < y_size; ++y)
	{
		for (uint x = 0; x < x_size; ++x)
		{
			if (selection == 0)
				printf("  %6g", values[y][x]);
			else if (selection > 0)
				printf("  %6g", errs_u[y][x]);
			else
				printf("  %6g", errs_l[y][x]);
		}
		putchar('\n');
	}
}

void fill_default(TH2DA * h)
{
	for (uint x = 0; x < x_size; ++x)
	{
		for (uint y = 0; y < y_size; ++y)
		{
			h->SetBinContent(x+1, y+1, values[y][x]);
			h->SetBinErrorH(x+1, y+1, errs_u[y][x]);
			h->SetBinErrorL(x+1, y+1, errs_l[y][x]);
		}
	}
}

void defaults_integral(double & value, double & err_u, double & err_l,
	size_t x_sta, size_t x_sto, size_t y_sta, size_t y_sto)
{
	value = err_l = err_u = 0;

	for (uint x = x_sta; x < x_sto; ++x)
	{
		for (uint y = y_sta; y < y_sto; ++y)
		{
			value += values[y][x];
			err_u += errs_u[y][x] * errs_u[y][x];
			err_l += errs_l[y][x] * errs_l[y][x];
		}
	}

	err_l = sqrt(err_l);
	err_u = sqrt(err_u);
}

void fill_random(TH2 * h)
{
	int xnum = h->GetXaxis()->GetNbins();
	int ynum = h->GetYaxis()->GetNbins();

	for (int i = 1; i <= xnum; ++i)
	{
		for (int j = 1; j <= ynum; ++j)
		{
			h->Fill(h->GetXaxis()->GetBinCenter(i), h->GetYaxis()->GetBinCenter(j), rand() % 1000);
		}
	}
}

void fill_flat(TH2 * h)
{
	int xnum = h->GetXaxis()->GetNbins();
	int ynum = h->GetYaxis()->GetNbins();

	for (int i = 1; i <= xnum; ++i)
	{
		for (int j = 1; j <= ynum; ++j)
		{
			h->Fill(h->GetXaxis()->GetBinCenter(i), h->GetYaxis()->GetBinCenter(j), ynum*xnum);
		}
	}
}

void fill_rising(TH2 * h)
{
	int xnum = h->GetXaxis()->GetNbins();
	int ynum = h->GetYaxis()->GetNbins();

	for (int j = 1; j <= ynum; ++j)
	{
		for (int i = 1; i <= xnum; ++i)
		{
			h->Fill(h->GetXaxis()->GetBinCenter(i), h->GetYaxis()->GetBinCenter(j), (j-1)*xnum + i-1);
		}
	}
}

void PrintErrors(TH2D * h)
{
	int xnum = h->GetXaxis()->GetNbins();
	int ynum = h->GetYaxis()->GetNbins();

	printf("\t");
	for (int i = 1; i <= xnum; ++i)
	{
		printf("%d\t", i-1);
	}
	printf("\n\n");

	for (int j = 1; j <= ynum; ++j)
	{
		printf("%d:\t", j-1);

		for (int i = 1; i <= xnum; ++i)
		{
			printf("%g\t", h->GetBinError(i, j));
		}
		printf("\n\t");
		for (int i = 1; i <= xnum; ++i)
		{
			printf("%g\t", h->GetBinContent(i, j));
		}
		printf("\n");
	}
	printf("\n\n");
}