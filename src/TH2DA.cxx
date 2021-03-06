/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2013  Rafał Lalik <rafal.lalik@ph.tum.de>
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

#include "TH2DA.h"

#include <TEfficiency.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

#include <iostream>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

// Internal exceptions for the CheckConsistency method
class DifferentDimension: public std::exception {};
class DifferentNumberOfBins: public std::exception {};
class DifferentAxisLimits: public std::exception {};
class DifferentBinLimits: public std::exception {};
class DifferentLabels: public std::exception {};

// ClassImp(TH2DA)

TH2DA::TH2DA() : TH2D()
{
	init();
}

TH2DA::TH2DA(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup,Int_t nbinsy,Double_t ylow,Double_t yup) : TH2D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup)
{
	init();
}

TH2DA::TH2DA(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins,Int_t nbinsy,Double_t ylow,Double_t yup) : TH2D(name, title, nbinsx, xbins, nbinsy, ylow, yup)
{
	init();
}

TH2DA::TH2DA(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup ,Int_t nbinsy,const Double_t *ybins) : TH2D(name, title, nbinsx, xlow, xup, nbinsy, ybins)
{
	init();
}

TH2DA::TH2DA(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins,Int_t nbinsy,const Double_t *ybins) : TH2D(name, title, nbinsx, xbins, nbinsy, ybins)
{
	init();
}

TH2DA::TH2DA(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins,Int_t nbinsy,const Float_t *ybins) : TH2D(name, title, nbinsx, xbins, nbinsy, ybins)
{
	init();
}

TH2DA::TH2DA(const TMatrixDBase &m) : TH2D(m)
{
	init();
}

TH2DA::TH2DA(const TH2D &h2d) : TH2D(h2d)
{
	init();
}

TH2DA::~TH2DA()
{
}

void TH2DA::Copy(TObject& obj) const
{
	TH2D::Copy(obj);
	((TH2DA&)obj).fSumw2h = fSumw2h;
	((TH2DA&)obj).fSumw2l = fSumw2l;
}

void TH2DA::init()
{
	Sumw2();
	fSumw2h.Set(fNcells);
	fSumw2l.Set(fNcells);
}


Double_t TH2DA::GetBinErrorH(Int_t bin) const
{
	if (bin < 0) bin = 0;
	if (bin >= fNcells) bin = fNcells-1;
	if (fSumw2h.fN) return TMath::Sqrt(fSumw2h.fArray[bin]);

	return GetBinError(bin);
}

Double_t TH2DA::GetBinErrorL(Int_t bin) const
{
	if (bin < 0) bin = 0;
	if (bin >= fNcells) bin = fNcells-1;
	if (fSumw2l.fN) return TMath::Sqrt(fSumw2l.fArray[bin]);

	return GetBinError(bin);
}

void TH2DA::SetBinError(Int_t binx, Int_t biny, Double_t error)
{
	int bin = GetBin(binx,biny);
	_SetBinError(fSumw2h, bin, error);
	_SetBinError(fSumw2l, bin, error);

	TH2D::SetBinError(binx, biny, error);
}

void TH2DA::SetBinError(Int_t binx, Int_t biny, Double_t error_u, Double_t error_l)
{
	int bin = GetBin(binx,biny);
	_SetBinError(fSumw2h, bin, error_u);
	_SetBinError(fSumw2l, bin, error_l);
}

void TH2DA::SetBinError(Int_t binx, Int_t biny, Int_t binz, Double_t error)
{
	int bin = GetBin(binx,biny,binz);
	_SetBinError(fSumw2h, bin, error);
	_SetBinError(fSumw2l, bin, error);
}

void TH2DA::SetBinErrorH(Int_t bin, Double_t error)
{
	_SetBinError(fSumw2h, bin, error);
}

void TH2DA::SetBinErrorL(Int_t bin, Double_t error)
{
	_SetBinError(fSumw2l, bin, error);
}

void TH2DA::SetBinErrorH(Int_t binx, Int_t biny, Double_t error_u)
{
	_SetBinError(fSumw2h, GetBin(binx, biny), error_u);
}

void TH2DA::SetBinErrorL(Int_t binx, Int_t biny, Double_t error_l)
{
	_SetBinError(fSumw2l, GetBin(binx, biny), error_l);
}

void TH2DA::PrintErrors()
{
	int n_num = this->GetXaxis()->GetNbins();
	int p_num = this->GetYaxis()->GetNbins();

	printf("\t");
	for (int i = 1; i <= n_num; ++i)
	{
		printf("%10d\t", i-1);
	}
	printf("\n\n");

	for (int j = 1; j <= p_num; ++j)
	{
		printf("%d:\t", j-1);

		for (int i = 1; i <= n_num; ++i)
		{
			printf("%10g\t", GetBinErrorH(i, j));
		}
		printf("\n\t");
		for (int i = 1; i <= n_num; ++i)
		{
			printf("%10g\t", GetBinContent(i, j));
		}
		printf("\n\t");
		for (int i = 1; i <= n_num; ++i)
		{
			printf("%10g\t", GetBinErrorL(i, j));
		}
		printf("\n\t");
		for (int i = 1; i <= n_num; ++i)
		{
			printf("%10g\t", GetBinError(i, j));
		}
		printf("\n");
	}
	printf("\n\n");
}

Int_t TH2DA::Fill(Double_t x, Double_t y)
{
	Int_t binx, biny, bin;
	binx = fXaxis.FindBin(x);
	biny = fYaxis.FindBin(y);
	bin = GetBin(binx,biny);

	TH2D::Fill(x, y);
	++fSumw2h[bin];
	++fSumw2l[bin];

	return bin;
}

Int_t TH2DA::Fill(Double_t x, Double_t y, Double_t w)
{
	Int_t binx, biny, bin;
	binx = fXaxis.FindBin(x);
	biny = fYaxis.FindBin(y);
	bin = GetBin(binx,biny);

	TH2D::Fill(x, y, w);
	fSumw2h[bin] += w*w;
	fSumw2l[bin] += w*w;

	return bin;
}

void TH2DA::Scale(Double_t c1, Option_t * option)
{
	TH2D::Scale(c1, option);

	int n_num = this->GetXaxis()->GetNbins();
	int p_num = this->GetYaxis()->GetNbins();

	for (int i = 1; i <= n_num; ++i)
	{
		for (int j = 1; j <= p_num; ++j)
		{
			int bin = GetBin(i,j);
			fSumw2h[bin] *= c1*c1;
			fSumw2l[bin] *= c1*c1;
		}
	}
}

// ** protected, general functions ** //

Double_t TH2DA::GetBinError(const TArrayD & arr, Int_t bin) const
{PR(1);
	if (bin < 0) bin = 0;PR(1);
	if (bin >= fNcells) bin = fNcells-1;PR(1);
	
	if (arr.fN) return TMath::Sqrt(arr.fArray[bin]);PR(1);

	return GetBinError(bin);
}

void TH2DA::_SetBinError(TArrayD & arr, Int_t bin, Double_t error)
{
	if (bin < 0 || bin>= arr.fN) return;
	arr.fArray[bin] = error * error;
}

Bool_t TH2DA::Add(const TH1 *h2, Double_t c1)
{
	TH1::Add(h2, c1);
	//   - Loop on bins (including underflows/overflows)
// 	Int_t bin, binx, biny, binz;
// 	Double_t c0,c1,w;
// 	for (binz=0;binz<=nbinsz+1;binz++) {
// 		for (biny=0;biny<=nbinsy+1;biny++) {
// 			for (binx=0;binx<=nbinsx+1;binx++) {
// 				bin = GetBin(binx,biny,binz);
// 				c0  = GetBinContent(bin);
// 				c1  = h1->GetBinContent(bin);
// 				w   = c0*c1;
// 				SetBinContent(bin,w);
// 				if (fSumw2.fN)
// 				{
// 					Double_t e0 = GetBinError(bin);
// 					Double_t e1 = h1->GetBinError(bin);
// 
// 					Double_t e0u = GetBinErrorH(bin);
// 					Double_t e0l = GetBinErrorL(bin);
// 
// 					Double_t e1u = 0.;
// 					Double_t e1l = 0.;
// 
// 					if (is_assymetric)
// 					{
// 						e1u = ((TH2DA*)h1)->GetBinErrorH(bin);
// 						e1l = ((TH2DA*)h1)->GetBinErrorL(bin);
// 					}
// 					else
// 					{
// 						e1u = h1->GetBinError(bin);
// 						e1l = e1u;
// 					}
// 
// 					fSumw2.fArray[bin] = (e0*e0*c1*c1 + e1*e1*c0*c0);
// 					fSumw2h.fArray[bin] = (e0u*e0u*c1*c1 + e1u*e1u*c0*c0);
// 					fSumw2l.fArray[bin] = (e0l*e0l*c1*c1 + e1l*e1l*c0*c0);
// 				}
// 			}
// 		}
// 	}
// 	ResetStats();

	if (!h2->InheritsFrom("TH2DA"))
		return kTRUE;

	TH2DA * _h2 = (TH2DA*)h2;
	for (int i = 0; i < fNcells; ++i)
	{
// 		double err_b = fSumw2h.fArray[i];
// 		double err_p = _h2->GetBinErrorH(i);
		if (fSumw2h.fN)
			fSumw2h.fArray[i] += c1*c1 * _h2->GetBinErrorH(i)*_h2->GetBinErrorH(i);
// 		double err_a = fSumw2h.fArray[i];
// 		printf("%g ^2 + %g ^2 = %g ^2\n", err_b, err_p, err_a);
		if (fSumw2l.fN)
			fSumw2l.fArray[i] += c1*c1 * _h2->GetBinErrorL(i)*_h2->GetBinErrorL(i);
	}
	return kTRUE;
}

Bool_t TH2DA::Divide(TF1 *f1, Double_t c1)
{
	if (!f1)
	{
		Error("Add","Attempt to divide by a non-existing function");
		return kFALSE;
	}

	// delete buffer if it is there since it will become invalid
	if (fBuffer) BufferEmpty(1);

	Int_t nbinsx = GetNbinsX();
	Int_t nbinsy = GetNbinsY();
	Int_t nbinsz = GetNbinsZ();
	if (fDimension < 2) nbinsy = -1;
	if (fDimension < 3) nbinsz = -1;

	SetMinimum();
	SetMaximum();

	//    Reset the kCanRebin option. Otherwise SetBinContent on the overflow bin
	//    would resize the axis limits!
	ResetBit(kCanRebin);

	//   - Loop on bins (including underflows/overflows)
	Int_t bin, binx, biny, binz;
	Double_t cu,w;
	Double_t xx[3];
	Double_t *params = 0;
	f1->InitArgs(xx,params);
	for (binz=0;binz<=nbinsz+1;binz++)
	{
		xx[2] = fZaxis.GetBinCenter(binz);
		for (biny=0;biny<=nbinsy+1;biny++)
		{
			xx[1] = fYaxis.GetBinCenter(biny);
			for (binx=0;binx<=nbinsx+1;binx++)
			{
				xx[0] = fXaxis.GetBinCenter(binx);
				if (!f1->IsInside(xx)) continue;
				TF1::RejectPoint(kFALSE);
				bin = binx +(nbinsx+2)*(biny + (nbinsy+2)*binz);
				Double_t error1 = GetBinError(bin);
				Double_t error1u = GetBinErrorH(bin);
				Double_t error1l = GetBinErrorL(bin);
				cu  = c1*f1->EvalPar(xx);
				if (TF1::RejectedPoint()) continue;
				if (cu) w = GetBinContent(bin)/cu;
				else    w = 0;
				SetBinContent(bin,w);
				if (fSumw2.fN)
				{
					if (cu != 0)
					{
						fSumw2.fArray[bin] = error1*error1/(cu*cu);
						fSumw2h.fArray[bin] = error1u*error1u/(cu*cu);
						fSumw2l.fArray[bin] = error1l*error1l/(cu*cu);
					}
					else
					{
						fSumw2.fArray[bin] = 0;
						fSumw2h.fArray[bin] = 0;
						fSumw2l.fArray[bin] = 0;
					}
				}
			}
		}
	}
	ResetStats();
	return kTRUE;
}

Bool_t TH2DA::Divide(const TH1 *h1)
{
	if (!h1)
	{
		Error("Divide","Attempt to divide by a non-existing histogram");
		return kFALSE;
	}

	// delete buffer if it is there since it will become invalid
	if (fBuffer) BufferEmpty(1);

	Int_t nbinsx = GetNbinsX();
	Int_t nbinsy = GetNbinsY();
	Int_t nbinsz = GetNbinsZ();

	try {
		CheckConsistency(this,h1);
	} catch(DifferentNumberOfBins&) {
		Error("Divide","Attempt to divide histograms with different number of bins");
		return kFALSE;
	} catch(DifferentAxisLimits&) {
		Warning("Divide","Attempt to divide histograms with different axis limits");
	} catch(DifferentBinLimits&) {
		Warning("Divide","Attempt to divide histograms with different bin limits");
	} catch(DifferentLabels&) {
		Warning("Divide","Attempt to divide histograms with different labels");
	}


	if (fDimension < 2) nbinsy = -1;
	if (fDimension < 3) nbinsz = -1;

	//    Create Sumw2 if h1 has Sumw2 set
	if (fSumw2.fN == 0 && h1->GetSumw2N() != 0) Sumw2();


	//    Reset the kCanRebin option. Otherwise SetBinContent on the overflow bin
	//    would resize the axis limits!
	ResetBit(kCanRebin);

	Bool_t is_assymetric = h1->InheritsFrom("TH2DA");
	//   - Loop on bins (including underflows/overflows)
	Int_t bin, binx, biny, binz;
	Double_t c0,c1,w;
	for (binz=0;binz<=nbinsz+1;binz++)
	{
		for (biny=0;biny<=nbinsy+1;biny++)
		{
			for (binx=0;binx<=nbinsx+1;binx++)
			{
				bin = GetBin(binx,biny,binz);
				c0  = GetBinContent(bin);
				c1  = h1->GetBinContent(bin);
				if (c1) w = c0/c1;
				else    w = 0;
				SetBinContent(bin,w);
				if (fSumw2.fN)
				{

					Double_t e0 = GetBinError(bin);
					Double_t e1 = h1->GetBinError(bin);

					Double_t e0u = GetBinErrorH(bin);
					Double_t e0l = GetBinErrorL(bin);

					Double_t e1u = 0.;
					Double_t e1l = 0.;
					if (is_assymetric)
					{
						e1u = ((TH2DA*)h1)->GetBinErrorH(bin);
						e1l = ((TH2DA*)h1)->GetBinErrorL(bin);
					}
					else
					{
						e1u = h1->GetBinError(bin);
						e1l = e1u;
					}

					Double_t c12= c1*c1;
					if (!c1)
					{
						fSumw2.fArray[bin] = 0;
						fSumw2h.fArray[bin] = 0;
						fSumw2l.fArray[bin] = 0;
						continue;
					}
					fSumw2.fArray[bin] = (e0*e0*c1*c1 + e1*e1*c0*c0)/(c12*c12);
					fSumw2h.fArray[bin] = (e0u*e0u*c1*c1 + e1u*e1u*c0*c0)/(c12*c12);
					fSumw2l.fArray[bin] = (e0l*e0l*c1*c1 + e1l*e1l*c0*c0)/(c12*c12);
				}
			}
		}
	}
	ResetStats();
	return kTRUE;
}

Bool_t TH2DA::Divide(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2, Option_t *option)
{
	TString opt = option;
	opt.ToLower();
	Bool_t binomial = kFALSE;
	if (opt.Contains("b")) binomial = kTRUE;
	if (!h1 || !h2) {
		Error("Divide","Attempt to divide by a non-existing histogram");
		return kFALSE;
	}

	// delete buffer if it is there since it will become invalid
	if (fBuffer) BufferEmpty(1);

	Int_t nbinsx = GetNbinsX();
	Int_t nbinsy = GetNbinsY();
	Int_t nbinsz = GetNbinsZ();

	try {
		CheckConsistency(h1,h2);
		CheckConsistency(this,h1);
	} catch(DifferentNumberOfBins&) {
		Error("Divide","Attempt to divide histograms with different number of bins");
		return kFALSE;
	} catch(DifferentAxisLimits&) {
		Warning("Divide","Attempt to divide histograms with different axis limits");
	} catch(DifferentBinLimits&) {
		Warning("Divide","Attempt to divide histograms with different bin limits");
	}  catch(DifferentLabels&) {
		Warning("Divide","Attempt to divide histograms with different labels");
	}


	if (!c2) {
		Error("Divide","Coefficient of dividing histogram cannot be zero");
		return kFALSE;
	}

	if (fDimension < 2) nbinsy = -1;
	if (fDimension < 3) nbinsz = -1;

	//    Create Sumw2 if h1 or h2 have Sumw2 set
	if (fSumw2.fN == 0 && (h1->GetSumw2N() != 0 || h2->GetSumw2N() != 0)) Sumw2();

	SetMinimum();
	SetMaximum();

	//    Reset the kCanRebin option. Otherwise SetBinContent on the overflow bin
	//    would resize the axis limits!
	ResetBit(kCanRebin);

	Bool_t is_assymetric1 = h1->InheritsFrom("TH2DA");
	Bool_t is_assymetric2 = h2->InheritsFrom("TH2DA");

	//   - Loop on bins (including underflows/overflows)
	Int_t bin, binx, biny, binz;
	Double_t b1,b2,w,d1,d2;
	d1 = c1*c1;
	d2 = c2*c2;
	for (binz=0;binz<=nbinsz+1;binz++)
	{
		for (biny=0;biny<=nbinsy+1;biny++)
		{
			for (binx=0;binx<=nbinsx+1;binx++)
			{
				bin = binx +(nbinsx+2)*(biny + (nbinsy+2)*binz);
				b1  = h1->GetBinContent(bin);
				b2  = h2->GetBinContent(bin);
				if (b2) w = c1*b1/(c2*b2);
				else    w = 0;
				SetBinContent(bin,w);
				if (fSumw2.fN)
				{
					Double_t e1 = h1->GetBinError(bin);
					Double_t e2 = h2->GetBinError(bin);

					Double_t e1u = 0.0;
					Double_t e1l = 0.0;

					Double_t e2u = 0.0;
					Double_t e2l = 0.0;

					if (is_assymetric1)
					{
						e1u = ((TH2DA*)h1)->GetBinErrorH(bin);
						e1l = ((TH2DA*)h1)->GetBinErrorL(bin);
					}
					else
					{
						e1u = h1->GetBinError(bin);
						e1l = e1u;
					}

					if (is_assymetric2)
					{
						e2u = ((TH2DA*)h2)->GetBinErrorH(bin);
						e2l = ((TH2DA*)h2)->GetBinErrorL(bin);
					}
					else
					{
						e2u = h2->GetBinError(bin);
						e2l = e2u;
					}

					Double_t b22= b2*b2*d2;
					if (!b2) { fSumw2.fArray[bin] = 0; continue;}
					if (binomial) {
						if (b1 != b2) {
							// in the case of binomial statistics c1 and c2 must be 1 otherwise it does not make sense
							w = b1/b2;    // c1 and c2 are ignored
							//fSumw2.fArray[bin] = TMath::Abs(w*(1-w)/(c2*b2));//this is the formula in Hbook/Hoper1
							//fSumw2.fArray[bin] = TMath::Abs(w*(1-w)/b2);     // old formula from G. Flucke
							// formula which works also for weighted histogram (see http://root.cern.ch/phpBB2/viewtopic.php?t=3753 )
							fSumw2.fArray[bin] = TMath::Abs( ( (1.-2.*w)*e1*e1 + w*w*e2*e2 )/(b2*b2) );
							fSumw2h.fArray[bin] = TMath::Abs( ( (1.-2.*w)*e1u*e1u + w*w*e2u*e2u )/(b2*b2) );
							fSumw2l.fArray[bin] = TMath::Abs( ( (1.-2.*w)*e1l*e1l + w*w*e2l*e2l )/(b2*b2) );
						} else {
							//in case b1=b2 error is zero
							//use  TGraphAsymmErrors::BayesDivide for getting the asymmetric error not equal to zero
							fSumw2.fArray[bin] = 0;
							fSumw2h.fArray[bin] = 0;
							fSumw2l.fArray[bin] = 0;
						}
					} else {
						fSumw2.fArray[bin] = d1*d2*(e1*e1*b2*b2 + e2*e2*b1*b1)/(b22*b22);
						fSumw2h.fArray[bin] = d1*d2*(e1u*e1u*b2*b2 + e2u*e2u*b1*b1)/(b22*b22);
						fSumw2l.fArray[bin] = d1*d2*(e1l*e1l*b2*b2 + e2l*e2l*b1*b1)/(b22*b22);
					}
				}
			}
		}
	}
	ResetStats();
	if (binomial)
		// in case of binomial division use denominator for number of entries
		SetEntries ( h2->GetEntries() );

	return kTRUE;
}

Bool_t TH2DA::Multiply(TF1 *f1, Double_t c1)
{
	if (!f1)
	{
		Error("Add","Attempt to multiply by a non-existing function");
		return kFALSE;
	}

	// delete buffer if it is there since it will become invalid
	if (fBuffer) BufferEmpty(1);

	Int_t nbinsx = GetNbinsX();
	Int_t nbinsy = GetNbinsY();
	Int_t nbinsz = GetNbinsZ();
	if (fDimension < 2) nbinsy = -1;
	if (fDimension < 3) nbinsz = -1;

	// reset min-maximum
	SetMinimum();
	SetMaximum();

	//    Reset the kCanRebin option. Otherwise SetBinContent on the overflow bin
	//    would resize the axis limits!
	ResetBit(kCanRebin);

	//   - Loop on bins (including underflows/overflows)
	Int_t bin, binx, biny, binz;
	Double_t cu,w;
	Double_t xx[3];
	Double_t *params = 0;
	f1->InitArgs(xx,params);
	for (binz=0;binz<=nbinsz+1;binz++)
	{
		xx[2] = fZaxis.GetBinCenter(binz);
		for (biny=0;biny<=nbinsy+1;biny++)
		{
			xx[1] = fYaxis.GetBinCenter(biny);
			for (binx=0;binx<=nbinsx+1;binx++)
			{
				xx[0] = fXaxis.GetBinCenter(binx);
				if (!f1->IsInside(xx)) continue;
				TF1::RejectPoint(kFALSE);
				bin = binx +(nbinsx+2)*(biny + (nbinsy+2)*binz);
				Double_t error1 = GetBinError(bin);
				Double_t error1u = GetBinErrorH(bin);
				Double_t error1l = GetBinErrorL(bin);
				cu  = c1*f1->EvalPar(xx);
				if (TF1::RejectedPoint()) continue;
				w = GetBinContent(bin)*cu;
				SetBinContent(bin,w);
				if (fSumw2.fN)
				{
					fSumw2.fArray[bin] = cu*cu*error1*error1;
					fSumw2h.fArray[bin] = cu*cu*error1u*error1u;
					fSumw2l.fArray[bin] = cu*cu*error1l*error1l;
				}
			}
		}
	}
	ResetStats();
	return kTRUE;
}

Bool_t TH2DA::Multiply(const TH1 *h1)
{
	if (!h1)
	{
		Error("Multiply","Attempt to multiply by a non-existing histogram");
		return kFALSE;
	}

	Int_t nbinsx = GetNbinsX();
	Int_t nbinsy = GetNbinsY();
	Int_t nbinsz = GetNbinsZ();

	// delete buffer if it is there since it will become invalid
	if (fBuffer) BufferEmpty(1);

	try {
		CheckConsistency(this,h1);
	} catch(DifferentNumberOfBins&) {
		Error("Multiply","Attempt to multiply histograms with different number of bins");
		return kFALSE;
	} catch(DifferentAxisLimits&) {
		Warning("Multiply","Attempt to multiply histograms with different axis limits");
	} catch(DifferentBinLimits&) {
		Warning("Multiply","Attempt to multiply histograms with different bin limits");
	} catch(DifferentLabels&) {
		Warning("Multiply","Attempt to multiply histograms with different labels");
	}


	if (fDimension < 2) nbinsy = -1;
	if (fDimension < 3) nbinsz = -1;

	//    Create Sumw2 if h1 has Sumw2 set
	if (fSumw2.fN == 0 && h1->GetSumw2N() != 0) Sumw2();

	//   - Reset min-  maximum
	SetMinimum();
	SetMaximum();

	//    Reset the kCanRebin option. Otherwise SetBinContent on the overflow bin
	//    would resize the axis limits!
	ResetBit(kCanRebin);

	Bool_t is_assymetric = h1->InheritsFrom("TH2DA");

	//   - Loop on bins (including underflows/overflows)
	Int_t bin, binx, biny, binz;
	Double_t c0,c1,w;
	for (binz=0;binz<=nbinsz+1;binz++) {
		for (biny=0;biny<=nbinsy+1;biny++) {
			for (binx=0;binx<=nbinsx+1;binx++) {
				bin = GetBin(binx,biny,binz);
				c0  = GetBinContent(bin);
				c1  = h1->GetBinContent(bin);
				w   = c0*c1;
				SetBinContent(bin,w);
				if (fSumw2.fN)
				{
					Double_t e0 = GetBinError(bin);
					Double_t e1 = h1->GetBinError(bin);

					Double_t e0u = GetBinErrorH(bin);
					Double_t e0l = GetBinErrorL(bin);

					Double_t e1u = 0.;
					Double_t e1l = 0.;

					if (is_assymetric)
					{
						e1u = ((TH2DA*)h1)->GetBinErrorH(bin);
						e1l = ((TH2DA*)h1)->GetBinErrorL(bin);
					}
					else
					{
						e1u = h1->GetBinError(bin);
						e1l = e1u;
					}

					fSumw2.fArray[bin] = (e0*e0*c1*c1 + e1*e1*c0*c0);
					fSumw2h.fArray[bin] = (e0u*e0u*c1*c1 + e1u*e1u*c0*c0);
					fSumw2l.fArray[bin] = (e0l*e0l*c1*c1 + e1l*e1l*c0*c0);
				}
			}
		}
	}
	ResetStats();
	return kTRUE;
}

Bool_t TH2DA::Multiply(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2, Option_t *option)
{
	TString opt = option;
	opt.ToLower();
	//   Bool_t binomial = kFALSE;
	//   if (opt.Contains("b")) binomial = kTRUE;
	if (!h1 || !h2) {
		Error("Multiply","Attempt to multiply by a non-existing histogram");
		return kFALSE;
	}

	// delete buffer if it is there since it will become invalid
	if (fBuffer) BufferEmpty(1);

	Int_t nbinsx = GetNbinsX();
	Int_t nbinsy = GetNbinsY();
	Int_t nbinsz = GetNbinsZ();

	try {
		CheckConsistency(h1,h2);
		CheckConsistency(this,h1);
	} catch(DifferentNumberOfBins&) {
		Error("Multiply","Attempt to multiply histograms with different number of bins");
		return kFALSE;
	} catch(DifferentAxisLimits&) {
		Warning("Multiply","Attempt to multiply histograms with different axis limits");
	} catch(DifferentBinLimits&) {
		Warning("Multiply","Attempt to multiply histograms with different bin limits");
	} catch(DifferentLabels&) {
		Warning("Multiply","Attempt to multiply histograms with different labels");
	}

	if (fDimension < 2) nbinsy = -1;
	if (fDimension < 3) nbinsz = -1;

	//    Create Sumw2 if h1 or h2 have Sumw2 set
	if (fSumw2.fN == 0 && (h1->GetSumw2N() != 0 || h2->GetSumw2N() != 0)) Sumw2();

	//   - Reset min - maximum
	SetMinimum();
	SetMaximum();

	//    Reset the kCanRebin option. Otherwise SetBinContent on the overflow bin
	//    would resize the axis limits!
	ResetBit(kCanRebin);

	Bool_t is_assymetric1 = h1->InheritsFrom("TH2DA");
	Bool_t is_assymetric2 = h2->InheritsFrom("TH2DA");

	//   - Loop on bins (including underflows/overflows)
	Int_t bin, binx, biny, binz;
	Double_t b1,b2,w,d1,d2;
	d1 = c1*c1;
	d2 = c2*c2;
	for (binz=0;binz<=nbinsz+1;binz++)
	{
		for (biny=0;biny<=nbinsy+1;biny++)
		{
			for (binx=0;binx<=nbinsx+1;binx++)
			{
				bin = binx +(nbinsx+2)*(biny + (nbinsy+2)*binz);
				b1  = h1->GetBinContent(bin);
				b2  = h2->GetBinContent(bin);
				w   = (c1*b1)*(c2*b2);
				SetBinContent(bin,w);
				if (fSumw2.fN)
				{
					Double_t e1 = h1->GetBinError(bin);
					Double_t e2 = h2->GetBinError(bin);

					Double_t e1u = 0.0;
					Double_t e1l = 0.0;

					Double_t e2u = 0.0;
					Double_t e2l = 0.0;

					if (is_assymetric1)
					{
						e1u = ((TH2DA*)h1)->GetBinErrorH(bin);
						e1l = ((TH2DA*)h1)->GetBinErrorL(bin);
					}
					else
					{
						e1u = h1->GetBinError(bin);
						e1l = e1u;
					}

					if (is_assymetric2)
					{
						e2u = ((TH2DA*)h2)->GetBinErrorH(bin);
						e2l = ((TH2DA*)h2)->GetBinErrorL(bin);
					}
					else
					{
						e2u = h2->GetBinError(bin);
						e2l = e2u;
					}

					fSumw2.fArray[bin] = d1*d2*(e1*e1*b2*b2 + e2*e2*b1*b1);
					fSumw2h.fArray[bin] = d1*d2*(e1u*e1u*b2*b2 + e2u*e2u*b1*b1);
					fSumw2l.fArray[bin] = d1*d2*(e1l*e1l*b2*b2 + e2l*e2l*b1*b1);
				}
			}
		}
	}
	ResetStats();
	return kTRUE;
}

void TH2DA::Reset(Option_t* option)
{
	TH2D::Reset(option);
	fSumw2h.Reset();
	fSumw2l.Reset();
}

TH2DA * TH2DA::Efficiency(const char * name, const TH2DA* total, const TH2DA* pass)
{
	static const Double_t conf = 0.682689492137;

	try {
		CheckConsistency(total,pass);
	} catch(DifferentNumberOfBins&) {
		total->Error("Divide","Attempt to divide histograms with different number of bins");
		return 0;
	} catch(DifferentAxisLimits&) {
		total->Warning("Divide","Attempt to divide histograms with different axis limits");
	} catch(DifferentBinLimits&) {
		total->Warning("Divide","Attempt to divide histograms with different bin limits");
	} catch(DifferentLabels&) {
		total->Warning("Divide","Attempt to divide histograms with different labels");
	}

	TH2DA * heff = (TH2DA *)pass->Clone(name);
	heff->Reset();

	int n_num = total->GetXaxis()->GetNbins();
	int p_num = total->GetYaxis()->GetNbins();

	for (int i = 1; i <= n_num; ++i)
	{
		for (int j = 1; j <= p_num; ++j)
		{
			int bin = total->GetBin(i,j);

			Double_t bc_total = total->GetBinContent(bin);
			Double_t bc_pass = pass->GetBinContent(bin);
			Double_t bUpper = TEfficiency::ClopperPearson(bc_total, bc_pass, conf, true);
			Double_t bLower = TEfficiency::ClopperPearson(bc_total, bc_pass, conf, false);

			Double_t p = 0.0;

			if (bc_total)
				p = bc_pass/bc_total;

			heff->SetBinContent(bin, p);

			heff->SetBinErrorH(bin, bUpper - p);
			heff->SetBinErrorL(bin, p - bLower);

			Double_t q = 1.0 - p;
			heff->SetBinError(bin, TMath::Sqrt(p * q * bc_total)/bc_total);
		}
	}

	return heff;
}

Double_t TH2DA::GetTotalErrorH() const
{
	return GetTotalError(fSumw2h);
}

Double_t TH2DA::GetTotalErrorL() const
{
	return GetTotalError(fSumw2l);
}

Double_t TH2DA::GetTotalError(const TArrayD& arr) const
{
	Double_t total_err = 0.0;
	int n_num = this->GetXaxis()->GetNbins();
	int p_num = this->GetYaxis()->GetNbins();

	for (int i = 1; i <= n_num; ++i)
	{
		for (int j = 1; j <= p_num; ++j)
		{
			int bin = GetBin(i,j);

			total_err += arr.fArray[bin];
		}
	}

	return TMath::Sqrt(total_err);
}

TGraphAsymmErrors * TH2DA::BuildAsymmErrorsGraph(Size_t points, Double_t* x, Double_t* xe, Double_t* y, Double_t* yl, Double_t* yu) const
{
	TGraphAsymmErrors * graph = new TGraphAsymmErrors(points);

	for (int i = 0; i < points; ++i)
	{
		graph->SetPoint(i, x[i], y[i]);
		graph->SetPointError(i, xe[i], xe[i], yl[i], yu[i]);
	}

	return graph;
}

TGraphAsymmErrors * TH2DA::BuildErrorsProjection(bool use_xaxis, int sta, int sto) const
{
	// proj ection axis, is the one on which sta-sto range is evaluated
	// normal axis the other one
	TAxis * axis_proj = nullptr;
	TAxis * axis_norm = nullptr;

	// use this counters in the loop
	int n = 0;		// iterate over n axis
	int p = 0;		// iterate over p axis

	// assign axes according to projection (use_xaxis) selection
	int * real_proj_axis = nullptr;
	int * real_norm_axis = nullptr;

	if (use_xaxis)
	{
		axis_proj = this->GetXaxis();
		axis_norm = this->GetYaxis();

		real_proj_axis = &p;
			real_norm_axis = &n;
	}
	else
	{
		axis_proj = this->GetYaxis();
		axis_norm = this->GetXaxis();

		real_proj_axis = &n;
		real_norm_axis = &p;
	}

	// from now on, work in n-p space, where
	// n is treated as x-axis and p as y-axis
	int n_num = axis_norm->GetNbins();
	int p_num = axis_proj->GetNbins();

	if (sto < sta)
	{
		int temp = sta;
		sta = sto;
		sto = temp;
	}

	if (sta < 0)
		sta = 0;

	if (sto > p_num+1)
		sto = p_num+1;

	Double_t * x = new Double_t[n_num];
	Double_t * xe = new Double_t[n_num];
	Double_t * y = new Double_t[n_num];
	Double_t * yu = new Double_t[n_num];
	Double_t * yl = new Double_t[n_num];

	size_t real_n = 0;
	for (n = 0; n < n_num; ++n)
	{
		double tmp_v = 0.0;
		double tmp_u = 0.0;
		double tmp_l = 0.0;

		for (p = sta; p <= sto; ++p)
		{
			Int_t bin = GetBin(*real_proj_axis+1, *real_norm_axis+1);
			tmp_v += GetBinContent(bin);
			tmp_u += GetBinErrorH(bin)*GetBinErrorH(bin);
			tmp_l += GetBinErrorL(bin)*GetBinErrorL(bin);
		}

		if (tmp_v == 0.0 and tmp_u == 0.0 and tmp_l == 0.0)
			continue;

		x[real_n] = axis_norm->GetBinCenter(n+1);
		xe[real_n] = axis_norm->GetBinWidth(n+1)/2.0;

		y[real_n] = tmp_v;
		yu[real_n] = sqrt(tmp_u);
		yl[real_n] = sqrt(tmp_l);

		++real_n;
	}

	TGraphAsymmErrors * graph = BuildAsymmErrorsGraph(real_n, x, xe, y, yl, yu);

	delete [] x;
	delete [] xe;
	delete [] y;
	delete [] yu;
	delete [] yl;

	return graph;
}

TGraphAsymmErrors * TH2DA::ErrorsProjectionX(int sta, int sto) const
{
	return BuildErrorsProjection(false, sta, sto);
}

TGraphAsymmErrors * TH2DA::ErrorsProjectionY(int sta, int sto) const
{
	return BuildErrorsProjection(true, sta, sto);
}