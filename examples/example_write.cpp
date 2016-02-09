#include "SmartFactory.h"
#include "TH2DA.h"

int main()
{
	// create factory
	
	SmartFactory * fac = new SmartFactory("factory1");

	TH2DA * hist1 = fac->RegTH2<TH2DA>("dir1/hist1", "Histogram 1", 10, -5, 5, 10, -5, 5);

	hist1->SetBinContent(3,3, 11);
	hist1->SetBinErrorH(3,3, 1.1);
	hist1->SetBinErrorL(3,3, 1.2);

	fac->rename("renamed_factory");
	fac->chdir("renamed_directory");

	// list objects
	fac->listRegisteredObjects();

	// export factory to file
	fac->exportStructure("example.root", true);
}