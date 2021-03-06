#include <cstdlib>

#undef NDEBUG
#include <assert.h>

#include <iostream>
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH2DA.h"
#include "TFile.h"

#include <TApplication.h>
#include <TROOT.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

int main(int argc, char ** argv)
{
	TROOT AnalysisDST_Cal1("TreeAnalysis","compiled analysisDST macros");
	TApplication app("treeanal", NULL, NULL, NULL, 0);

	if (argc != 3)
		exit(EXIT_FAILURE);

	TFile * f = TFile::Open(argv[1], "READ");

	TH2DA * h = nullptr;

	f->GetObject(argv[2], h);

	if (!h)
	{
		std::cout << "No object f" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!h->InheritsFrom("TH2DA"))
	{
		std::cout << "Object not of type TH2DA" << std::endl;
		exit(EXIT_FAILURE);
	}

	uint xbins = h->GetXaxis()->GetNbins();
	uint ybins = h->GetYaxis()->GetNbins();

	TCanvas * can_x = new TCanvas("can_x", "can_x", 1000, 1000);
	can_x->DivideSquare(ybins);

	char buff[100];
	for (uint i = 0; i < ybins; ++i)
	{
		sprintf(buff, "h_p_y_%d", i);
		TH1D * h_p = h->ProjectionX(buff, i+1, i+1);
// 		PR(h->GetBinErrorH(i+1,3));
		can_x->cd(i+1);
		h_p->Draw("l");
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

// 	while(true) {}
	app.Run();

	f->Close();
}
