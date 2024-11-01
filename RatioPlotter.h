#ifndef RATIOPLOTTER
#define RATIOPLOTTER
#include "TH1.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"

//TODO: Maybe take ownership of the histograms instead of cloning?
//TODO: Check for empty histograms
//TODO: Upperpad not getting setlogy from canvas...

class RatioPlotter{
public:
	RatioPlotter(TH1 *h1, TH1 *h2, Option_t *option="");
	~RatioPlotter();
	
	//Accessors
	TPad* GetUpperPad();
	TPad* GetLowerPad();
	TAxis* GetXaxis();
	TAxis* GetUpperYaxis();
	TAxis* GetLowerYaxis();

	const TH1* GetH3();
	
	//Modifiers
	void SetTitle(const char *title);
	void SetNormalize(Bool_t normalize=true);
	void SetAxisTitles(std::string xtitle, std::string ytitle);
	void SetAxisTitlesWithUnits(std::string xtitle, std::string units);
	TLegend* BuildLegend(Double_t x1=0.3, Double_t y1=0.21, Double_t x2=0.3, Double_t y2=0.21, const char *title="", Option_t *option="");
	
	void SetH1LineColor(Color_t color);
	void SetH2LineColor(Color_t color);
	void SetH3LineColor(Color_t color);
	void SetH1DrawOpt(Option_t *options);
	void SetH2DrawOpt(Option_t *options);
	void SetH3DrawOpt(Option_t *options);
	
	void SetUpperYaxisTitle(const char *title);
	void SetUpperYaxisTitleOffset(Float_t offset);
	void SetLowerXaxisTitle(const char *title);
	void SetLowerXaxisTitleOffset(Float_t offset);
	void SetLowerYaxisTitle(const char *title);
	void SetLowerYaxisTitleOffset(Float_t offset);
	
	//Methods
	void Draw(TCanvas *c);
	
private:
	void SetupPads();
	void SetupHists();
	Bool_t BuildLowerHist();
	void FinalizePads();
	
	enum CalculationMode{
		kDivideHist = 1, 	//TH1::Divide
		kRelDiff = 2,		//Relative Difference
	};
	
	//RatioPlotter Options
	Bool_t _normalize = false;
	Bool_t _drawn = false;
	CalculationMode _mode = kDivideHist;
	
	//Histogram Options
	Color_t _h1Color = gStyle->GetHistLineColor(), _h2Color = kRed;
	Color_t _h3Color = _h1Color;
	Bool_t _h3Overridden = false;
	TString _h1DrawOpt = "hist", _h2DrawOpt = "hist", _h3DrawOpt = "hist";
	
	//Pad Options
	Float_t _InsetWidth = 0.0025;
	Float_t _UpTopMargin = 0.1, _UpBottomMargin = 0.05;
	Float_t _LowTopMargin = 0.05, _LowBottomMargin = 0.3;
	Float_t _LeftMargin = 0.1, _RightMargin = 0.1;
	Float_t _SplitFraction = 0.3;
	
	//Histograms
	TH1 *_h1=nullptr, *_h2=nullptr;
	TH1 *_h3=nullptr;
	TH1 *_hzero=nullptr;
	const char *_h1XaxisTitle;
	
	//Pads
	TPad *_UpperPad=nullptr, *_LowerPad=nullptr, *_TopPad=nullptr;
};

#endif//RATIOPLOTTER
