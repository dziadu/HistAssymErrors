#include <cppunit/extensions/HelperMacros.h>

#include <TApplication.h>
#include <TROOT.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>

#include <TH2DA.h>

#include "test_common.h"

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

class ProjectionsCase_tests : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( ProjectionsCase_tests );
	CPPUNIT_TEST( MyTest );
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();

protected:
	void MyTest();

	TH2DA * h;
	TH2D * hd;
};

CPPUNIT_TEST_SUITE_REGISTRATION( ProjectionsCase_tests );

void ProjectionsCase_tests::setUp()
{
	h = new TH2DA("h", "h", 10, -5, 5, 10, -5, 5);
	fill_random(h);

	hd = new TH2D("hd", "hd", 10, -5, 5, 10, -5, 5);
	hd->Sumw2();

	fill_random(hd);
}

void ProjectionsCase_tests::MyTest()
{
	TH2DA * h_1 = new TH2DA("h_1", "h_1", x_size, 0, 10, y_size, 0, 10);
	fill_default(h_1);

	double def_int;
	double def_err_u;
	double def_err_l;

	// test proper integral calculation
	defaults_integral(def_int, def_err_u, def_err_l);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( def_int, h_1->Integral(), 0.0001 );

	TH1D * h_1_p01 = h_1->ProjectionY("h_1_p01",0, 2);
	defaults_integral(def_int, def_err_u, def_err_l, 0, 2, 0, y_size);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( def_int, h_1_p01->Integral() , 0.0001 );

	print_default(1);
	TGraphAsymmErrors * assym_errors = nullptr;
	for (uint i = 0; i < y_size; ++i)
	{
		assym_errors = h_1->ErrorsProjectionY(0, i+1);
		for (uint j = 0; j < x_size; ++j)
		{
			defaults_integral(def_int, def_err_u, def_err_l, j, j+1, 0, i+1);
// 			printf("testing dy=%d, x=%d\n", i+1, j);

			CPPUNIT_ASSERT_DOUBLES_EQUAL( def_err_u, assym_errors->GetErrorYhigh(j), 0.0000000001 );
			CPPUNIT_ASSERT_DOUBLES_EQUAL( def_err_l, assym_errors->GetErrorYlow(j), 0.0000000001 );
		}
	}

	for (uint i = 0; i < x_size; ++i)
	{
		assym_errors = h_1->ErrorsProjectionX(0, i+1);
		for (uint j = 0; j < y_size; ++j)
		{
			defaults_integral(def_int, def_err_u, def_err_l, 0, i+1, j, j+1);
// 			printf("testing dx=%d, y=%d\n", i+1, j);

			CPPUNIT_ASSERT_DOUBLES_EQUAL( def_err_u, assym_errors->GetErrorYhigh(j), 0.0000000001 );
			CPPUNIT_ASSERT_DOUBLES_EQUAL( def_err_l, assym_errors->GetErrorYlow(j), 0.0000000001 );
		}
	}
}
