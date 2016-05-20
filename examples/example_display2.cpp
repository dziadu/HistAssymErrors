#include <cstdlib>

#undef NDEBUG
#include <assert.h>

#include <iostream>
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH2DA.h"
#include "TFile.h"
#include <TMath.h>

#include <TApplication.h>
#include <TROOT.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

int main(int argc, char ** argv)
{
	TROOT AnalysisDST_Cal1("TreeAnalysis","compiled analysisDST macros");
	TApplication app("treeanal", NULL, NULL, NULL, 0);

	size_t bins_x = 10;
	size_t bins_y = 2;
	TH2DA * h = new TH2DA("hist1", "Histogram 1", bins_x, -5, 5, bins_y, -5, 5);

	for (uint i = 0; i < bins_x; ++i)
	{
		// linear
		h->SetBinContent(1+i, 1, 3-0.05*i);
		h->SetBinErrorH(1+i, 1, 1.+0.1*i);
		h->SetBinErrorL(1+i, 1, 2.-0.1*i);

		// gaussian
		h->SetBinContent(1+i, 2, 2);
		h->SetBinErrorH(1+i, 2, 1.+ TMath::Gaus(-4.5 + i, 0, 2.)*2.);
		h->SetBinErrorL(1+i, 2, 2.- TMath::Gaus(-4.5 + i, 0, 2.));
	}

	uint xbins = h->GetXaxis()->GetNbins();
	uint ybins = h->GetYaxis()->GetNbins();

	TCanvas * can_x = new TCanvas("can_x", "can_x", 1000, 1000);
	can_x->DivideSquare(ybins);

	char buff[100];
	for (uint i = 0; i < ybins; ++i)
	{
		sprintf(buff, "h_p_y_%d", i);

		// Get projection of the histogram
		TH1D * h_p = h->ProjectionX(buff, i+1, i+1);
		can_x->cd(i+1);
		h_p->Draw("l");

		// Get projection of the errors
		// Pay attention that bin==1 in the hist is bin==0 in the errors
		TGraphAsymmErrors * err = h->ErrorsProjectionX(i, i);
		err->SetFillColor(7);
		err->SetFillStyle(3001);
		err->Draw("same,3");
		h_p->Draw("same,l");

		double range_max = err->GetX()[0];
		double range_min = err->GetY()[0];

		for (uint i = 0; i < xbins; ++i)
		{
			range_max = std::max<double>(range_max, err->GetY()[i] + err->GetErrorYhigh(i));
			range_min = std::min<double>(range_min, err->GetY()[i] - err->GetErrorYlow(i));
		}
		h_p->GetYaxis()->SetRangeUser(
			range_min - 0.1 * fabs(range_max-range_min),
			range_max + 0.1 * fabs(range_max-range_min));
	}

	can_x->Draw();
	// file must be closed by user
	can_x->Print();

	h->PrintErrors();

// 	while(true) {}
	app.Run();
}
