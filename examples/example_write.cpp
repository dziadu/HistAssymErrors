#include "SmartFactory.h"
#include "TH2DA.h"

#include <TMath.h>

int main()
{
	// create factory
	
	SmartFactory * fac = new SmartFactory("factory1");

	size_t bins_x = 10;
	size_t bins_y = 2;

	TH2DA * hist1 = fac->RegTH2<TH2DA>("hist1", "Histogram 1", bins_x, -5, 5, bins_y, -5, 5);

	for (uint i = 0; i < bins_x; ++i)
	{
		// linear
		hist1->SetBinContent(1+i, 1, 3-0.05*i);
		hist1->SetBinErrorH(1+i, 1, 1.+0.1*i);
		hist1->SetBinErrorL(1+i, 1, 2.-0.1*i);

		// gaussian
		hist1->SetBinContent(1+i, 2, 2);
		hist1->SetBinErrorH(1+i, 2, 1.+ TMath::Gaus(-4.5 + i, 0, 2.)*2.);
		hist1->SetBinErrorL(1+i, 2, 2.- TMath::Gaus(-4.5 + i, 0, 2.));
	}

	fac->rename("renamed_factory");
	fac->chdir("renamed_directory");

	// list objects
	fac->listRegisteredObjects();

	// export factory to file
	fac->exportStructure("example.root", true);
}
