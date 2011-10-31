#include "tdrStyle.C"

void rsLimit7TeV(){

  setTDRStyle();

//=========Macro generated from canvas: cLimit/Limit
//=========  (Mon Feb 22 22:44:48 2010) by ROOT version5.18/00a
//   TCanvas *cLimit = new TCanvas("cLimit", "Limit",450,40,800,550);
//  TCanvas *cLimit = new TCanvas("cLimit", "Limit",100,122,600,600);
  cLimit = new TCanvas("cLimit","cLimit",800,600);

   gStyle->SetOptStat(0);
   //   cLimit->Range(595.5973,-0.03483694,1345.283,0.2539571);
   cLimit->SetFillColor(0);
   cLimit->SetBorderMode(0);
   cLimit->SetBorderSize(2);
   cLimit->SetLeftMargin(0.139262);
   cLimit->SetRightMargin(0.0604027);
   cLimit->SetTopMargin(0.0804196);
   cLimit->SetBottomMargin(0.120629);
   cLimit->SetFrameBorderMode(0);
   cLimit->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(11);
   graph->SetName("Graph");
   graph->SetTitle("");
   graph->SetFillColor(1);


   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff0000");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.0);
   graph->SetPoint(0,  863,0.01);
   graph->SetPoint(1, 1132,0.02);
   graph->SetPoint(2, 1274,0.03);
   graph->SetPoint(3, 1395,0.04);
   graph->SetPoint(4, 1503,0.05);
   graph->SetPoint(5, 1597,0.06);
   graph->SetPoint(6, 1676,0.07);
   graph->SetPoint(7, 1742,0.08);
   graph->SetPoint(8, 1801,0.09);
   graph->SetPoint(9, 1844,0.1);
   graph->SetPoint(10,1881,0.11);
   
   cout << " Don't forget to change this..." << endl;
   TH1 *Graph6 = new TH1F("Graph6","",100,863,1881);
   Graph6->SetMinimum(0);
   Graph6->SetMaximum(0.12);
   Graph6->SetDirectory(0);
   Graph6->SetStats(0);
   Graph6->GetXaxis()->SetTitle("M_{1} [GeV]");
   Graph6->GetYaxis()->SetTitle("Coupling k/#bar{M}_{Pl}");
   Graph6->GetXaxis()->SetLabelFont(42);
   Graph6->GetYaxis()->SetLabelFont(42);
   Graph6->GetYaxis()->SetTitleOffset(1.8);
   graph->SetHistogram(Graph6);
   
   graph->GetYaxis()->SetTitleOffset(1.19);
   graph->Draw("al");


   TGraph* graph2 = new TGraph(11);
   graph2->SetMarkerColor(ci);
   graph2->SetMarkerStyle(20);
   graph2->SetMarkerSize(0.0);
   graph2->SetPoint(0,  845,0.01);
   graph2->SetPoint(1, 1133,0.02);
   graph2->SetPoint(2, 1275,0.03);
   graph2->SetPoint(3, 1396,0.04);
   graph2->SetPoint(4, 1504,0.05);
   graph2->SetPoint(5, 1598,0.06);
   graph2->SetPoint(6, 1677,0.07);
   graph2->SetPoint(7, 1743,0.08);
   graph2->SetPoint(8, 1801,0.09);
   graph2->SetPoint(9, 1844,0.1);
   graph2->SetPoint(10,1881,0.11);
   graph2->SetLineStyle(kDashed);
   graph2->SetLineColor(ci);
   graph2->SetLineWidth(3);
   graph2->GetXaxis()->SetLabelFont(42);
   graph2->GetYaxis()->SetLabelFont(42);
   graph2->Draw("plsame");


//    graph = new TGraph(3);
//    graph->SetName("Graph");
//    graph->SetTitle("");
//    graph->SetFillColor(1);

//    ci = TColor::GetColor("#0000ff");
//    graph->SetLineColor(ci);
//    graph->SetLineWidth(3);

//    ci = TColor::GetColor("#0000ff");
//    graph->SetMarkerColor(ci);
//    graph->SetMarkerStyle(22);
//    graph->SetMarkerSize(1.4);
//    graph->SetPoint(0,750,0.03355478);
//    graph->SetPoint(1,1000,0.07617687);
//    graph->SetPoint(2,1250,0.1542326);
   
//    TH1 *Graph7 = new TH1F("Graph7","",100,700,1300);
//    Graph7->SetMinimum(0.021487);
//    Graph7->SetMaximum(0.1663004);
//    Graph7->SetDirectory(0);
//    Graph7->SetStats(0);
//    Graph7->GetXaxis()->SetTitle("Graviton Mass (GeV)");
//    Graph7->GetYaxis()->SetTitle("Coupling k/#bar{M}_{Pl}");
//    graph->SetHistogram(Graph7);
   
//    graph->Draw("pc");
   
//    graph = new TGraph(3);
//    graph->SetName("Graph");
//    graph->SetTitle("");
//    graph->SetFillColor(1);

//    ci = TColor::GetColor("#00ff00");
//    graph->SetLineColor(ci);
//    graph->SetLineWidth(3);

//    ci = TColor::GetColor("#00ff00");
//    graph->SetMarkerColor(ci);
//    graph->SetMarkerStyle(21);
//    graph->SetMarkerSize(1.3);
//    graph->SetPoint(0,750,0.02431275);
//    graph->SetPoint(1,1000,0.05519538);
//    graph->SetPoint(2,1250,0.1117521);
   
//    TH1 *Graph8 = new TH1F("Graph8","",100,700,1300);
//    Graph8->SetMinimum(0.01556881);
//    Graph8->SetMaximum(0.1204961);
//    Graph8->SetDirectory(0);
//    Graph8->SetStats(0);
//    Graph8->GetXaxis()->SetTitle("Graviton Mass (GeV/c^{2})");
//    Graph8->GetYaxis()->SetTitle("Coupling k/#bar{M}_{Pl}");
//    graph->SetHistogram(Graph8);
   
//    graph->Draw("pc");
   
//    TF1 *LambdaPi = new TF1("LambdaPi","pol1",500,2500);
//    LambdaPi->SetFillColor(15);
//    LambdaPi->SetFillStyle(3004);
//    LambdaPi->SetLineColor(15);
//    LambdaPi->SetLineWidth(1);
//    LambdaPi->SetParameter(0,0);
//    LambdaPi->SetParError(0,0);
//    LambdaPi->SetParLimits(0,0,0);
//    LambdaPi->SetParameter(1,2.61097e-05);
//    LambdaPi->SetParError(1,0);
//    LambdaPi->SetParLimits(1,0,0);
//    LambdaPi->Draw("same");
// 
    TF1 *LambdaPi = new TF1("LambdaPi","pol1",250,2500);
    LambdaPi->SetFillColor(15);
    LambdaPi->SetFillStyle(3004);
    LambdaPi->SetLineColor(15);
    LambdaPi->SetLineWidth(1);
    LambdaPi->SetParameter(0,0);
    LambdaPi->SetParError(0,0);
    LambdaPi->SetParLimits(0,0,0);
    LambdaPi->SetParameter(1,2.61097e-05);
    LambdaPi->SetParError(1,0);
    LambdaPi->SetParLimits(1,0,0);
    LambdaPi->GetXaxis()->SetLabelFont(42);
    LambdaPi->GetYaxis()->SetLabelFont(42);
    LambdaPi->Draw("same");
   
    graph = new TGraph(27);
    graph->SetName("Graph");
    graph->SetTitle("Graph");
    graph->SetFillColor(1);
    graph->SetFillStyle(3004);    
    graph->SetLineStyle(5);
    graph->SetLineWidth(3);
    graph->SetPoint(0,180,0.1071);
    graph->SetPoint(1,183,0.1062);
    graph->SetPoint(2,190,0.1043);
    graph->SetPoint(3,200,0.1016);
    graph->SetPoint(4,210,0.0989);
    graph->SetPoint(5,220,0.0963);
    graph->SetPoint(6,230,0.0938);
    graph->SetPoint(7,240,0.0913);
    graph->SetPoint(8,250,0.0889);
    graph->SetPoint(9,260,0.0866);
    graph->SetPoint(10,270,0.0843);
    graph->SetPoint(11,280,0.0821);
    graph->SetPoint(12,290,0.0799);
    graph->SetPoint(13,300,0.0778);
    graph->SetPoint(14,400,0.0603);
    graph->SetPoint(15,500,0.0481);
    graph->SetPoint(16,600,0.0397);
    graph->SetPoint(17,700,0.0337);
    graph->SetPoint(18,800,0.0292);
    graph->SetPoint(19,900,0.0258);
    graph->SetPoint(20,1000,0.0231);
    graph->SetPoint(21,1100,0.0209);
    graph->SetPoint(22,1200,0.0191);
    graph->SetPoint(23,1300,0.0176);
    graph->SetPoint(24,1400,0.0163);
    graph->SetPoint(25,1500,0.0152);
    graph->SetPoint(26,2200,0.01);
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetLabelFont(42);
    graph->SetFillColor(kBlack);
    graph->SetFillStyle(3004);
    graph->SetLineWidth(1);
    graph->SetFillStyle(3002);
    graph->SetFillColor(kBlack);
    graph->SetLineWidth(-10000);

    graph->Draw("c");
    //    graph->Draw("same");

//    TPaveText *pt = new TPaveText(0.194631,0.748252,0.436242,0.816434,"blNDC");
//    pt->SetName("95% CL Limit");
//    pt->SetFillColor(0);
//    pt->SetBorderSize(1);
//    pt->SetLineColor(0);
//    pt->SetTextSize(0.04);
//    TText *text = pt->AddText("95% CL Limit");
//    pt->Draw();
  
   TLegend *leg = new TLegend(0.2072864,0.4667832,0.4007538,0.7517483,NULL,"brNDC");
   //    TLegend *leg = new TLegend(0.2072864,0.5332168,0.4007538,0.7517483,NULL,"brNDC"); 
    //   TLegend *leg = new TLegend(0.2110553,0.4248252,0.4120603,0.6433566,NULL,"brNDC");
   //   TLegend *leg = new TLegend(0.2181208,0.5297203,0.4194631,0.7482517,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetTextSize(0.05);
   leg->SetLineColor(0);
   leg->SetLineStyle(0);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   //   TLegendEntry *entry=leg->AddEntry("graph","50/pb","pl");

//    ci = TColor::GetColor("#ff0000");
//    entry->SetLineColor(ci);
//    entry->SetLineStyle(1);
//    entry->SetLineWidth(3);

    entry=leg->AddEntry("graph","Electroweak Limits","lf");
    entry->SetLineColor(1);
    entry->SetLineStyle(5);
    entry->SetLineWidth(3);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1);

    ci = TColor::GetColor("#ff0000");
    entry->SetMarkerColor(ci);
    entry->SetMarkerStyle(20);
    entry->SetMarkerSize(1.3);
    entry=leg->AddEntry("Graph","95% CL Limit","l");


    leg->AddEntry(graph2,"Expected Limit","l");
    ci = TColor::GetColor("#00ff00");
    entry->SetMarkerColor(ci);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1.3);
    //    entry=leg->AddEntry("LambdaPi","#Lambda_{#pi}> 10TeV","lf");
    entry=leg->AddEntry("LambdaPi","M_{D} > 10TeV","lf");

//     entry->SetFillColor(15);
//     entry->SetFillStyle(3004);
//     entry->SetLineColor(15);
//     entry->SetLineStyle(1);
//     entry->SetLineWidth(1);
//     entry->SetMarkerColor(1);
//     entry->SetMarkerStyle(1);
//     entry->SetMarkerSize(1);

    leg->Draw();
   
   TPaveText *pt = new TPaveText(0.2236181,0.7884615,0.4736181,0.8583916,"blNDC");
   //   TPaveText *pt = new TPaveText(0.1959799,0.7954545,0.4459799,0.8653846,"blNDC");
   //   TPaveText *pt = new TPaveText(0.2130872,0.8339161,0.4630872,0.9038462,"blNDC");
   pt->SetName("CMS Preliminary");
   pt->SetBorderSize(1);
   pt->SetLineColor(0);
   pt->SetFillColor(0);
   //   pt->SetTextFont(72);
   pt->SetTextSize(0.06);
   text = pt->AddText("CMS Preliminary");
   pt->Draw();
   
   pt = new TPaveText(0.629397,0.798951,0.8002513,0.8653846,"blNDC");
   //   pt = new TPaveText(0.6005025,0.8059441,0.7713568,0.8723776,"blNDC");
   //   pt = new TPaveText(0.2738693,0.7027972,0.4447236,0.7692308,"blNDC");
   //   pt = new TPaveText(0.2416107,0.7727273,0.4127517,0.8391608,"blNDC");
   pt->SetFillColor(0);
   pt->SetBorderSize(1);
   pt->SetLineColor(0);
   pt->SetTextSize(0.06);
   text = pt->AddText("2.2 fb^{-1} at 7 TeV");
   pt->Draw();

   //   pt = new TPaveText(0.6879195,0.8129371,0.7885906,0.9108392,"blNDC");
   //   pt->SetFillColor(0);
   //   pt->SetBorderSize(1);
   //   pt->SetLineColor(0);
   //   pt->SetTextSize(0.0454545);
   //   text = pt->AddText("L = 1091 pb^{-1}");
   //  pt->Draw();


   cLimit->Modified();
   cLimit->cd();
   cLimit->SetSelected(cLimit);
}
