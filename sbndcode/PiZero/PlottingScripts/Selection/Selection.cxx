#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.h"
#include "SelectionLibrary/Selections.h"
#include "SelectionLibrary/CommonUtils.h"
#include "SelectionLibrary/Style.h"
#include "LatexHeaders.h"
#include "Plots.h"

#include "TSystem.h"
#include "TROOT.h"

#include <fstream>

void Selection(const TString productionVersion, const SelectionParams &selection, std::vector<Plot> &plots = selection_plots);

void ProduceCutTable(const TString &saveDir, std::vector<Sample> &samples, const SelectionParams &selection);

int main()
{
  Selection("NCPiZeroAv11", Selections::ncpizero_incl, no_plots);
}

void Selection(const TString productionVersion, const SelectionParams &selection, std::vector<Plot> &plots)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/selection/" + selection.name;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_intime.root";

  Style::SetStyle();
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *intimeEvents = new TChain("ncpizeroana/events");
  intimeEvents->Add(intimeFile);

  TChain *rockboxsubruns = new TChain("ncpizeroana/subruns");
  rockboxsubruns->Add(rockboxFile);
  TChain *intimesubruns = new TChain("ncpizeroana/subruns");
  intimesubruns->Add(intimeFile);

  double rockboxScaling, intimeScaling;
  std::cout << rockboxsubruns->GetEntries() << std::endl;
  CommonUtils::GetScaling(rockboxsubruns, intimesubruns, rockboxScaling, intimeScaling);

  std::vector<Sample> samples = { { "rockbox", rockboxEvents, rockboxScaling },
                                          { "intime", intimeEvents, intimeScaling }
  };

  ProduceCutTable(saveDir, samples, selection);

  TCut currentCut = "";

  for(auto cut : selection.cuts)
    {
      currentCut += cut.cut;
      cut.cut = currentCut;
      
      if(plots.size() != 0)
        gSystem->Exec("mkdir -p " + saveDir + "/" + cut.name);
      for(auto plot : plots)
        {
          TCanvas *canvas = new TCanvas("c_" + plot.name + "_" + cut.name,
                                        "c_" + plot.name + "_" + cut.name);
          canvas->cd();

          plot.axes_labels += CommonUtils::POTString();

	  Plotting::MakeStackedPlot(canvas, samples, plot, cut, selection.categories, {.25, .8, .8, .87}, 4);

          canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
          canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");

          delete canvas;
        }
    }

  gSystem->Exec("pdflatex -output-directory " + saveDir + " " + saveDir + "/cut_table.tex");
}

void ProduceCutTable(const TString &saveDir, std::vector<Sample> &samples, const SelectionParams &selection)
{
  std::ofstream texFile;
  texFile.open(saveDir + "/cut_table.tex");

  double totalSignal = 0, totalSignalSlices = 0, totalBackSlices = 0;

  for(auto const& sample : samples)
    {
      totalSignal       += sample.scaling * sample.tree->Draw("", selection.true_category);
      totalSignalSlices += sample.scaling * sample.tree->Draw("", selection.categories[0].cut);
      totalBackSlices   += sample.scaling * sample.tree->Draw("", !selection.categories[0].cut);
    }

  texFile << docStart;

  texFile << '\n'
          << "Total Signal: " << totalSignal << "\\\\ \n"
          << "Total Signal Slices: " << totalSignalSlices << "\\\\ \n"
          << "Total Background Slices: " << totalBackSlices << "\\\\ \n" << std::endl;

  texFile << tableStart 
          << "\\hline\n"
          << "Cut Name & $\\epsilon$ (\\%) & $\\rho$ (\\%) & $\\epsilon\\rho$ & Selection $\\epsilon$ (\\%) & Selection $\\epsilon\\rho$ & BR (\\%) \\\\ \\hline" 
          << std::endl;

  TCut currentCut = "";

  for(unsigned i = 0; i < selection.cuts.size(); ++i)
    {
      Cut cut = selection.cuts[i];
      currentCut += cut.cut;
      cut.cut = currentCut;

      double sigSlices = 0., backSlices = 0.;
      for(unsigned j = 0; j < selection.categories.size(); ++j)
        {
          for(auto const& sample : samples)
            {
              if(j == 0)
                sigSlices += sample.scaling * sample.tree->Draw("", cut.cut + selection.categories[j].cut);
              else
                backSlices += sample.scaling * sample.tree->Draw("", cut.cut + selection.categories[j].cut);
            }
        }

      const double eff     = sigSlices * 100. / totalSignal;
      const double selEff  = sigSlices * 100. / totalSignalSlices;
      const double pur     = sigSlices * 100./ (sigSlices + backSlices);
      const double backRej = 100. - 100. * (backSlices / totalBackSlices);
      
      texFile << cut.printed_name << " & " << Form("%.2f", eff) << " & " << Form("%.2f", pur)
              << " & " << Form("%.2f", (eff * pur) / 100.)
              << " & " << Form("%.2f", selEff) << " & " << Form("%.2f", (selEff * pur) / 100.)
              << " & " << Form("%.2f", backRej) << "\\\\ \\hline" << std::endl;
    }

  texFile << tableEnd << docEnd;
}
