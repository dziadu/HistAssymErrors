#include <cppunit/extensions/HelperMacros.h>

#include <TApplication.h>
#include <TROOT.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>

#include <TH2DA.h>

#include "test_common.h"

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

class BasicCase_tests : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( BasicCase_tests );
	CPPUNIT_TEST( MyTest );
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();

protected:
	void MyTest();

	TH2DA * h;
	TH2D * hd;
};

CPPUNIT_TEST_SUITE_REGISTRATION( BasicCase_tests );

void BasicCase_tests::setUp()
{
	h = new TH2DA("h", "h", 10, -5, 5, 10, -5, 5);
	fill_random(h);

	hd = new TH2D("hd", "hd", 10, -5, 5, 10, -5, 5);
	hd->Sumw2();

	fill_random(hd);
}

void BasicCase_tests::MyTest()
{
	float fnum = 2.00001f;
// 	CPPUNIT_FAIL("zxczc");
	CPPUNIT_ASSERT_DOUBLES_EQUAL( fnum, 2.0f, 0.0005 );

	std::string pattern_string("%%d pattern");
	std::string test_string("test pattern");
	std::string replace_string("test");

	TROOT troot("TreeAnalysis","compiled analysisDST macros");
	TApplication app();

// 	h->Draw();
// 	h->Print();
// 	hd->Print();

// 	h->PrintErrors();
// 	PrintErrors(hd);

	TH2DA * h2 = (TH2DA *)h->Clone("new h");
// 	h2->PrintErrors();

	h2->Divide(h);
// 	h2->PrintErrors();
	CPPUNIT_ASSERT_DOUBLES_EQUAL( sqrt(2), h2->GetBinErrorH(1, 1), 0.0005 );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, h2->GetBinContent(1, 1), 0.0005 );

	TH2D * hd2 = (TH2D*)hd->Clone("new hd");
	hd2->Divide(hd);
// 	PrintErrors(hd2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( sqrt(2), hd2->GetBinError(1, 1), 0.0005 );

	TH2DA * h3 = (TH2DA *)h->Clone("h3: new h");
	TH2DA * h4 = (TH2DA *)h->Clone("h4: new h");

	h3->Reset();
	h4->Reset();

	fill_flat(h3);
	fill_rising(h4);

// 	h3->PrintErrors();
// 	h4->PrintErrors();

	TH2DA * h5 = TH2DA::Efficiency("h5_eff", h3, h4);
// 	h5->PrintErrors();

// 	PR(h5->GetTotalErrorU());
// 	PR(h5->Integral());
// 	PR(h5->GetTotalErrorL());

	CPPUNIT_ASSERT_DOUBLES_EQUAL( (h5->GetTotalErrorH() - h5->GetTotalErrorL()) == 0.0, false, 0.0005 );

// 	TFile * f = new TFile("/tmp/test.root", "RECREATE");
// 	if (f->IsOpen())
// 	{
// 		f->cd();
// 		h5->Write();
// 		h5->GetAsymErrorsGraphX(4)->Write("grx");
// 		h5->GetAsymErrorsGraphY(4)->Write("gry");
// 		f->Close();
// 	}
}
