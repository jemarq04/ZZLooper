#include "RatioPlotter.h"
#include "TClass.h"
#include "TCanvas.h"
#include <iostream>

RatioPlotter::RatioPlotter(TH1 *h1, TH1 *h2, Option_t *option){
	if (!h1 || !h2){
		std::cerr << "RatioPlotter needs two histograms." << std::endl;
		return;
	}
	
	if (!h1->IsA()->InheritsFrom(TH1::Class()) 
			|| !h2->IsA()->InheritsFrom(TH1::Class())){
		std::cerr << "RatioPlotter needs two histograms derived from TH2 or TH3" << std::endl;
		return;
	}
	
	_h1 = (TH1*)h1->Clone("__h1");
	_h2 = (TH1*)h2->Clone("__h2");
	_hzero = (TH1*)_h1->Clone("__hzero");
	_hzero->Add(_hzero, -1);
	_hzero->SetLineColor(kBlack);
	_hzero->SetLineStyle(kDashed);
	
	_h1XaxisTitle = _h1->GetXaxis()->GetTitle();
	_h1->GetXaxis()->SetTitle("");
	
	TString optionString = TString(option);
	if (optionString.Contains("div")){
		optionString.ReplaceAll("div", "");
		_mode = kDivideHist;
	}
	else if (optionString.Contains("reldiff")){
		optionString.ReplaceAll("reldiff", "");
		_mode = kRelDiff;
	}
	else{
		_mode = kDivideHist;
	}
	
	SetupPads();
	if (!BuildLowerHist()) 
		std::cerr << "Fatal issue building lower histogram" << std::endl;
}

RatioPlotter::~RatioPlotter(){
	if (!_drawn){
		if (_UpperPad) delete _UpperPad;
		if (_LowerPad) delete _LowerPad;
		if (_TopPad) delete _TopPad;
	}
}

TPad* RatioPlotter::GetUpperPad(){return _UpperPad;}
TPad* RatioPlotter::GetLowerPad(){return _LowerPad;}
TAxis* RatioPlotter::GetXaxis(){return _h3->GetXaxis();}
TAxis* RatioPlotter::GetUpperYaxis(){return _h1->GetYaxis();}
TAxis* RatioPlotter::GetLowerYaxis(){return _h3->GetYaxis();}
const TH1* RatioPlotter::GetH3(){return _h3;}

void RatioPlotter::SetTitle(const char *title){
	if (!_drawn) _h1->SetTitle(title);
}
void RatioPlotter::SetNormalize(Bool_t normalize){
	_normalize = normalize;
}
void RatioPlotter::SetAxisTitles(std::string xtitle, std::string ytitle){
	_h3->GetXaxis()->SetTitle(xtitle.c_str());
	_h1->GetYaxis()->SetTitle(ytitle.c_str());
}
void RatioPlotter::SetAxisTitlesWithUnits(std::string xtitle, std::string units){
	_h1->BufferEmpty(-1);
	char y[100];
	Double_t bin_width = _h1->GetBinWidth(_h1->GetMaximumBin());
	if (units == "")
		snprintf(y, 100, "Number of Entries/%.2f", bin_width);
	else{
		xtitle += " (" + units + ")";
		snprintf(y, 100, "Number of Entries/%.2f (%s^{-1})", bin_width, units.c_str());
	}
	
	_h3->GetXaxis()->SetTitle(xtitle.c_str());
	_h1->GetYaxis()->SetTitle(y);
}
TLegend* RatioPlotter::BuildLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, const char *title, Option_t *option){
	return _UpperPad->BuildLegend(x1, y1, x2, y2, title, option);
}

void RatioPlotter::SetH1LineColor(Color_t color){
	_h1Color = color;
}
void RatioPlotter::SetH2LineColor(Color_t color){
	_h2Color = color;
}
void RatioPlotter::SetH3LineColor(Color_t color){
	_h3Color = color;
}
void RatioPlotter::SetH1DrawOpt(Option_t *options){
	_h1DrawOpt = options;
}
void RatioPlotter::SetH2DrawOpt(Option_t *options){
	_h2DrawOpt = options;
}
void RatioPlotter::SetH3DrawOpt(Option_t *options){
	_h3DrawOpt = options;
}
void RatioPlotter::SetUpperYaxisTitle(const char *title){
	_h1->GetYaxis()->SetTitle(title);
}
void RatioPlotter::SetUpperYaxisTitleOffset(Float_t offset){
	_h1->GetYaxis()->SetTitleOffset(offset);
}
void RatioPlotter::SetLowerXaxisTitle(const char *title){
	_h3->GetXaxis()->SetTitle(title);
}
void RatioPlotter::SetLowerXaxisTitleOffset(Float_t offset){
	_h3->GetXaxis()->SetTitleOffset(offset);
}
void RatioPlotter::SetLowerYaxisTitle(const char *title){
	_h3->GetYaxis()->SetTitle(title);
}
void RatioPlotter::SetLowerYaxisTitleOffset(Float_t offset){
	_h3->GetYaxis()->SetTitleOffset(offset);
}

void RatioPlotter::SetupPads(){
	if (_UpperPad){
		delete _UpperPad;
		_UpperPad = nullptr;
	}
	if (_LowerPad){
		delete _LowerPad;
		_LowerPad = nullptr;
	}
	if (!gPad){
		std::cerr << "Canvas needs to be created before first" << std::endl;
		return;
	}
	double pm = _InsetWidth;
	double width = gPad->GetWNDC();
	double height = gPad->GetHNDC();
	double f = height/width;
	_UpperPad = new TPad("upper_pad", "", pm*f, _SplitFraction, 1.-pm*f, 1.-pm);
	_LowerPad = new TPad("lower_pad", "", pm*f, pm, 1.-pm*f, _SplitFraction);
	
	_UpperPad->SetTopMargin(_UpTopMargin);
	_UpperPad->SetBottomMargin(_UpBottomMargin);
	_UpperPad->SetLeftMargin(_LeftMargin);
	_UpperPad->SetRightMargin(_RightMargin);
	_LowerPad->SetTopMargin(_LowTopMargin);
	_LowerPad->SetBottomMargin(_LowBottomMargin);
	_LowerPad->SetLeftMargin(_LeftMargin);
	_LowerPad->SetRightMargin(_RightMargin);
	
	if (_TopPad){
		delete _TopPad;
		_TopPad = nullptr;
	}
	_TopPad = new TPad("top_pad", "", pm*f, pm, 1-pm*f, 1-pm);
	_TopPad->SetBit(kCannotPick);
}

void RatioPlotter::SetupHists(){
	Double_t max1 = _h1->GetMaximum()*1.1, max2 = _h2->GetMaximum()*1.1;
	Double_t max = (max1 > max2)? max1 : max2;
	_h1->GetYaxis()->SetRangeUser(0, max);
	_h2->GetYaxis()->SetRangeUser(0, max);

	if (_normalize){
		_h1->Scale(1./_h1->Integral());
		_h2->Scale(1./_h2->Integral());
		BuildLowerHist();
	}

	_h1->SetLineColor(_h1Color);
	_h2->SetLineColor(_h2Color);
}

Bool_t RatioPlotter::BuildLowerHist(){
	if (_h3){
		delete _h3;
		_h3 = nullptr;
	}
	if (_mode == kDivideHist){
		_h3 = (TH1*)_h1->Clone("__h3");
		_h3->Reset();
		_h3->Divide(_h1, _h2);
	}
	else if (_mode == kRelDiff){
		_h3 = (TH1*)_h1->Clone("__h3");
		_h3->Add(_h2, -1);
		_h3->Divide(_h2);
		_h3->GetYaxis()->SetRangeUser(-1, 1);
		_h3->SetTitle("");
		_h3->GetXaxis()->SetTitle(_h1XaxisTitle);
		_h3->GetYaxis()->SetTitle("");
		_h3->GetYaxis()->SetNdivisions(8);
	}
	else return false;
	return true;
}

void RatioPlotter::FinalizePads(){
	_h1->GetXaxis()->SetLabelSize(0);
	_h3->SetStats(0);
	_h3->GetXaxis()->SetLabelSize(0.1);
	_h3->GetXaxis()->SetTitleOffset(1);
	_h3->GetXaxis()->SetTitleSize(0.1);
	_h3->GetYaxis()->SetTitleSize(0.1);
	_h3->GetYaxis()->SetTitleOffset(0.4);
	_h3->GetYaxis()->SetLabelSize(0.08);
	double ratio = ( (_UpperPad->GetBottomMargin()-(1-_UpperPad->GetTopMargin()))*(1-0.3) ) 
		/ ( (_LowerPad->GetBottomMargin()-(1-_LowerPad->GetTopMargin()))*0.3 );
	_h3->GetXaxis()->SetTickLength(_h1->GetYaxis()->GetTickLength() * ratio);
}

void RatioPlotter::Draw(TCanvas *c){
	_drawn = true;
	SetupHists();
	
	_UpperPad->SetLogy(c->GetLogy());
	
	c->cd();
	_TopPad->Draw();
	_UpperPad->Draw();
	_TopPad->cd();
	_LowerPad->Draw();
	
	_UpperPad->cd();
	_h1->Draw(_h1DrawOpt);
	_h2->Draw(_h2DrawOpt + ",same");

	_LowerPad->cd();
	_h3->Draw(_h3DrawOpt);
	_hzero->Draw("hist,same");
	
	FinalizePads();
	
	c->Modified();
	c->Update();
}
