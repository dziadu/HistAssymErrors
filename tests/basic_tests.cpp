#include <cppunit/extensions/HelperMacros.h>

#include <TApplication.h>
#include <TROOT.h>
#include <TF2.h>
#include <TFile.h>

#include <TH2DA.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

class BasicCase : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( BasicCase );
	CPPUNIT_TEST( MyTest );
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();

protected:
	void MyTest();

	TH2DA * h;
	TH2D * hd;
};

CPPUNIT_TEST_SUITE_REGISTRATION( BasicCase );

void fill_random(TH2 * h)
{
	int xnum = h->GetXaxis()->GetNbins();
	int ynum = h->GetYaxis()->GetNbins();

	for (int i = 1; i <= xnum; ++i)
	{
		for (int j = 1; j <= ynum; ++j)
		{
			h->Fill(h->GetXaxis()->GetBinCenter(i), h->GetYaxis()->GetBinCenter(j), rand() % 1000);
		}
	}
}

void fill_flat(TH2 * h)
{
	int xnum = h->GetXaxis()->GetNbins();
	int ynum = h->GetYaxis()->GetNbins();

	for (int i = 1; i <= xnum; ++i)
	{
		for (int j = 1; j <= ynum; ++j)
		{
			h->Fill(h->GetXaxis()->GetBinCenter(i), h->GetYaxis()->GetBinCenter(j), ynum*xnum);
		}
	}
}

void fill_rising(TH2 * h)
{
	int xnum = h->GetXaxis()->GetNbins();
	int ynum = h->GetYaxis()->GetNbins();

	for (int j = 1; j <= ynum; ++j)
	{
		for (int i = 1; i <= xnum; ++i)
		{
			h->Fill(h->GetXaxis()->GetBinCenter(i), h->GetYaxis()->GetBinCenter(j), (j-1)*xnum + i-1);
		}
	}
}

void PrintErrors(TH2D * h)
{
	int xnum = h->GetXaxis()->GetNbins();
	int ynum = h->GetYaxis()->GetNbins();

	printf("\t");
	for (int i = 1; i <= xnum; ++i)
	{
		printf("%d\t", i-1);
	}
	printf("\n\n");

	for (int j = 1; j <= ynum; ++j)
	{
		printf("%d:\t", j-1);

		for (int i = 1; i <= xnum; ++i)
		{
			printf("%g\t", h->GetBinError(i, j));
		}
		printf("\n\t");
		for (int i = 1; i <= xnum; ++i)
		{
			printf("%g\t", h->GetBinContent(i, j));
		}
		printf("\n");
	}
	printf("\n\n");
}

void BasicCase::setUp()
{
	h = new TH2DA("h", "h", 10, -5, 5, 10, -5, 5);
	fill_random(h);

	hd = new TH2D("hd", "hd", 10, -5, 5, 10, -5, 5);
	hd->Sumw2();

	fill_random(hd);
}

void BasicCase::MyTest()
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
	CPPUNIT_ASSERT_DOUBLES_EQUAL( sqrt(2), h2->GetBinErrorU(1, 1), 0.0005 );
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
	h5->PrintErrors();

	PR(h5->GetTotalErrorU());
	PR(h5->Integral());
	PR(h5->GetTotalErrorL());

// 	TFile * f = new TFile("/tmp/test.root", "RECREATE");
// 	if (f->IsOpen())
// 	{
// 		f->cd();
// 		h5->Write();
// 		f->Close();
// 	}
}
