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

#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <TApplication.h>
#include <TROOT.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>

#include <TH2DA.h>

static const size_t x_size = 5;
static const size_t y_size = 5;

extern const double values[y_size][x_size];
extern const double errs_u[y_size][x_size];
extern const double errs_l[y_size][x_size];

void print_default(int selection = 0);
void fill_default(TH2DA * h);
void defaults_integral(double & value, double & err_u, double & err_l,
					   size_t x_sta = 0, size_t x_sto = x_size, size_t y_sta = 0, size_t y_sto = y_size);

void fill_random(TH2 * h);
void fill_flat(TH2 * h);
void fill_rising(TH2 * h);
void PrintErrors(TH2D * h);


#endif // TEST_COMMON_H
