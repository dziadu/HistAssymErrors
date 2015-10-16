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

#ifndef TH2DA_H
#define TH2DA_H

#include <TArrayD.h>
#include <TH2.h>

class TF1;
class TGraphAsymmErrors;

class TH2DA : public TH2D
{
public:
	TH2DA();
	TH2DA(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup,Int_t nbinsy,Double_t ylow,Double_t yup);
	TH2DA(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins,Int_t nbinsy,Double_t ylow,Double_t yup);
	TH2DA(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup ,Int_t nbinsy,const Double_t *ybins);
	TH2DA(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins,Int_t nbinsy,const Double_t *ybins);
	TH2DA(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins,Int_t nbinsy,const Float_t *ybins);
	TH2DA(const TMatrixDBase &m);
	TH2DA(const TH2D &h2d);
	virtual ~TH2DA();

	using TH2D::Fill;
	using TH2D::GetBinError;
	using TH2D::SetBinContent;
	using TH2D::SetBinError;

	virtual void Copy(TObject & obj) const;

	virtual Bool_t   Divide(TF1 *f1, Double_t c1=1);
	virtual Bool_t   Divide(const TH1 *h1);
	virtual Bool_t   Divide(const TH1 *h1, const TH1 *h2, Double_t c1=1, Double_t c2=1, Option_t *option="");

	virtual Int_t    Fill(Double_t x, Double_t y);
	virtual Int_t    Fill(Double_t x, Double_t y, Double_t w);

	virtual Double_t GetBinErrorU(Int_t bin) const;
	virtual Double_t GetBinErrorL(Int_t bin) const;

	virtual Double_t GetBinErrorU(Int_t binx, Int_t biny) const { return GetBinErrorU(GetBin(binx, biny)); }
	virtual Double_t GetBinErrorL(Int_t binx, Int_t biny) const { return GetBinErrorL(GetBin(binx, biny)); }

	virtual Bool_t   Multiply(TF1 *f1, Double_t c1=1);
	virtual Bool_t   Multiply(const TH1 *h1);
	virtual Bool_t   Multiply(const TH1 *h1, const TH1 *h2, Double_t c1=1, Double_t c2=1, Option_t *option="");

	virtual void     Reset(Option_t *option="");

	virtual void     Scale(Double_t c1, Option_t * option);

	virtual void     SetBinError(Int_t binx, Int_t biny, Double_t error);
	virtual void     SetBinError(Int_t binx, Int_t biny, Double_t error_u, Double_t error_l);
	virtual void     SetBinError(Int_t binx, Int_t biny, Int_t binz, Double_t error);
	
	virtual void     SetBinErrorU(Int_t bin, Double_t error);
	virtual void     SetBinErrorL(Int_t bin, Double_t error);
	virtual void     SetBinErrorU(Int_t binx, Int_t biny, Double_t error_u);
	virtual void     SetBinErrorL(Int_t binx, Int_t biny, Double_t error_l);

	virtual void     PrintErrors();

	virtual Double_t GetTotalErrorU() const;
	virtual Double_t GetTotalErrorL() const;

	virtual TGraphAsymmErrors * GetAsymErrorsGraphX(int xcol) const;
	virtual TGraphAsymmErrors * GetAsymErrorsGraphY(int yrow) const;

	static  TH2DA   *Efficiency(const char * name, const TH2DA * total,  const TH2DA * pass);

	ClassDef(TH2DA,1);

protected:
	void init();

	virtual Double_t GetBinError(const TArrayD & arr, Int_t bin) const;
	virtual void     SetBinError(TArrayD & arr, Int_t bin, Double_t error);

	virtual Double_t GetTotalError(const TArrayD & arr) const;

	virtual TGraphAsymmErrors * BuildAsymmErrorsGraph(Size_t points, Double_t *x, Double_t *xe, Double_t *y, Double_t *yl, Double_t *yu) const;

protected:
	TArrayD fSumw2u;
	TArrayD fSumw2l;

};

#endif // TH2DA_H
