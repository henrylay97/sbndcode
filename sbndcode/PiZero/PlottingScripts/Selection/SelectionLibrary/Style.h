#include "TStyle.h"

namespace Style 
{
  void SetStyle()
  {
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(4);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetStatColor(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.04);

    gStyle->SetPaperSize(20,26);
    gStyle->SetCanvasDefH(1000);
    gStyle->SetCanvasDefW(1400);
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadRightMargin(0.06);
    gStyle->SetPadBottomMargin(0.2);
    gStyle->SetPadLeftMargin(0.2);

    gStyle->SetTextFont(62);
    gStyle->SetTextSize(0.09);
    gStyle->SetLabelFont(62,"x");
    gStyle->SetLabelFont(62,"y");
    gStyle->SetLabelFont(62,"z");
    gStyle->SetLabelSize(0.07,"x");
    gStyle->SetTitleSize(0.07,"x");
    gStyle->SetLabelSize(0.07,"y");
    gStyle->SetTitleSize(0.07,"y");
    gStyle->SetLabelSize(0.07,"z");
    gStyle->SetTitleSize(0.07,"z");
    gStyle->SetLabelFont(62,"t");
    gStyle->SetTitleFont(62,"x");
    gStyle->SetTitleFont(62,"y");
    gStyle->SetTitleFont(62,"z");
    gStyle->SetTitleFont(62,"t"); 
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(1,"z");
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetTitleFont(62,"pad");
    gStyle->SetTitleBorderSize(0);

    gStyle->SetMarkerStyle(20);
    gStyle->SetHistLineWidth(5);
    gStyle->SetLineStyleString(2,"[12 12]");

    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetNdivisions(508);
    gStyle->SetPalette(kBlueRedYellow);

    gStyle->SetPalette(1,0);
    const int NRGBs = 5;
    const int NCont = 255;
  
    double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    double red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
				     NCont);
    gStyle->SetNumberContours(NCont); 
  }
}
