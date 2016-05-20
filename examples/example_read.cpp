#include <cstdlib>

#undef NDEBUG
#include <assert.h>

#include "SmartFactory.h"
#include "TH2DA.h"

#include "TFile.h"

int main()
{
	// create factory
	SmartFactory * fac = new SmartFactory("factory1");

	// import from file and register in the factory
	// data will be stored in memory, file remains open
	TFile * f = fac->importStructure("example.root");

	// list of registered objects 
	fac->listRegisteredObjects();

	// you can fetch specific object by its name
	TH2DA * h1 = (TH2DA*)fac->getObject("hist1");

	// if failed, then objects are not read from file
	assert(h1 != nullptr);

	// file must be closed by user
	f->Close();
}
