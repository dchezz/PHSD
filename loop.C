
//Macro used to analyse PHSD trees
//Framework: ROOT
//Author: Dario Chaires-Arciniega



#include <TString.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>



Int_t loop(){

 gROOT->Reset();
  TChain mychain("phsd");
  mychain.Add("/home/dario/Desktop/PHSD/100000_events/phsd_100000.root");


// Declaration of leaf types
   Int_t           isub; // external loop (run over different values for the impact parameter)
   Int_t           irun; //internal loop (parallel events)
   Float_t         b;
   Float_t         ibw;
   Float_t         psi[4];
   Float_t         epsilon[4];
   Int_t           np;//participants
   Int_t           n;//total particles
   Int_t           id[9000];   //[n]
   Short_t         q[90000];   //[n]
   Float_t         e[9000];   //[n]
   Float_t         px[9000];   //[n]
   Float_t         py[9000];   //[n]
   Float_t         pz[9000];   //[n]
   Int_t           code1[9000];   //[n]
   Int_t           code2[9000];   //[n]

   // List of branches
   TBranch        *b_isub;   //!
   TBranch        *b_irun;   //!
   TBranch        *b_b;   //!
   TBranch        *b_ibw;   //!
   TBranch        *b_psi;   //!
   TBranch        *b_epsilon;   //!
   TBranch        *b_np;   //!
   TBranch        *b_n;   //!
   TBranch        *b_id;   //!
   TBranch        *b_q;   //!
   TBranch        *b_e;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_code1;   //!
   TBranch        *b_code2;   //!

////////////All particles, protons, neutrons, and all charged particles////////////////////////////

//Pseudorapidity distributions

TH1F *heta1 = new TH1F("h1b","PHSD, Au + Au, 11GeV, 100 min. bias events;#eta;# entries",317,-8,8);
TH1F *heta2 = new TH1F("h2b","eta",317,-8,8);
TH1F *heta3 = new TH1F("h3b","eta",317,-8,8);
TH1F *heta4 = new TH1F("h4b","eta",317,-8,8);

//Theta distributions
TH1F *htheta1 = new TH1F("h1a","PHSD, Au + Au, 11GeV, 100 min. bias events;#theta;# entries",317,-0.1,3.15);
TH1F *htheta2 = new TH1F("h2a","theta",317,-0.1,3.15);
TH1F *htheta3 = new TH1F("h3a","theta",317,-0.1,3.15);
TH1F *htheta4 = new TH1F("h4a","theta",317,-0.1,3.15);


//Pions
TH1F *pheta1 = new TH1F("h1c","PHSD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#eta;# entries",317,-8,8);
TH1F *pheta2 = new TH1F("h2c","eta",317,-8,8);

TH1F *phtheta1 = new TH1F("h1d","PHSD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#theta(rad);# entries",317,-0.1,3.15);
TH1F *phtheta2= new TH1F("h2d","theta",317,-0.1,3.15);

TH1F *hImpact = new TH1F("hb_100000", "PHSD Au+Au, 11 GeV, 100000 events, bmin=0 fm, bmax=14 fm;b; entries",100,-0.1,14.2); //Impact parameter

TF1 *f = new TF1("f","[0]*x",-0.1,14.2); //Sample fit

TF1 *g = new TF1("g","2*x/(196)",-0.1,14.2); //Min.Bias distribution
 
//Setting Branch Adress
   mychain.SetBranchAddress("isub", &isub, &b_isub);
   mychain.SetBranchAddress("irun", &irun, &b_irun);
   mychain.SetBranchAddress("b", &b, &b_b);
   mychain.SetBranchAddress("ibw", &ibw, &b_ibw);
   mychain.SetBranchAddress("psi", psi, &b_psi);
   mychain.SetBranchAddress("epsilon", epsilon, &b_epsilon);
   mychain.SetBranchAddress("np", &np, &b_np);
   mychain.SetBranchAddress("n", &n, &b_n);
   mychain.SetBranchAddress("id", id, &b_id);
   mychain.SetBranchAddress("q", q, &b_q);
   mychain.SetBranchAddress("e", e, &b_e);
   mychain.SetBranchAddress("px", px, &b_px);
   mychain.SetBranchAddress("py", py, &b_py);
   mychain.SetBranchAddress("pz", pz, &b_pz);
   mychain.SetBranchAddress("code1", code1, &b_code1);
   mychain.SetBranchAddress("code2", code2, &b_code2);
  

Int_t nevent = mychain.GetEntries();

for (Int_t jentry=0; jentry < nevent; jentry++) //Loop over events
    {
mychain.GetEvent(jentry);
hImpact->Fill(b); //Fill b dist.


for(Int_t j=0; j < n; j++) {  //Loop over particles

Float_t pt[9000];
Float_t theta[9000];
Float_t eta[9000];
pt[j] = sqrt(px[j]*px[j]+py[j]*py[j]);
theta[j] = atan2(pt[j],pz[j]);
eta[j] = -log(tan(theta[j]/2)); 

///Filling histograms/////////////////////////////////////////////////////////////////////////////////
    heta1->Fill(eta[j]);
    htheta1->Fill(theta[j]);             
if(id[j]==2212) {
    heta2->Fill(eta[j]);
    htheta2->Fill(theta[j]); 
        }
if(id[j]==2112) {
     heta3->Fill(eta[j]);
     htheta3->Fill(theta[j]);    
        }
if(q[j]!=0) {
     heta4->Fill(eta[j]);
     htheta4->Fill(theta[j]); 
        }
if(id[j]==-211)  {
  pheta1->Fill(eta[j]);
  phtheta1->Fill(theta[j]);
     }
if(id[j]==211)  {
  pheta2->Fill(eta[j]);
  phtheta2->Fill(theta[j]);
     }

      }//for j
      
   }//for i

/////////////////////////////////////// A E S T H E T I C S ///////////////////////////////////////////

TCanvas *c1=new TCanvas("c1","c1",1080,720);
c1->cd();

htheta1->SetMarkerColor(kBlue);
htheta1->SetMarkerStyle(kFullCross);
htheta1->Scale(1./(htheta1->GetEntries()));
htheta1->Draw();
htheta2->SetMarkerColor(kGreen);
htheta2->SetMarkerStyle(kFullStar);
htheta2->Scale(1./(htheta1->GetEntries()));
htheta2->Draw("Same");
htheta3->SetMarkerColor(kMagenta);
htheta3->SetMarkerStyle(kCircle);
htheta3->Scale(1./(htheta1->GetEntries()));
htheta3->Draw("Same");
htheta4->SetMarkerColor(kBlack);
htheta4->SetMarkerStyle(kFullDotLarge);
htheta4->Scale(1./(htheta1->GetEntries()));
htheta4->Draw("Same");

TLegend *legend=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend->SetHeader("");
legend->AddEntry(htheta1,"all particles","p");
legend->AddEntry(htheta2,"protons","p");
legend->AddEntry(htheta3,"neutrons","p");
legend->AddEntry(htheta4,"charged particles","p");
legend->Draw();


TCanvas *c2=new TCanvas("c2","c2",1080,720);
c2->cd();

heta1->SetMarkerColor(kBlue);
heta1->SetMarkerStyle(kFullCross);
heta1->Scale(1./(htheta1->GetEntries()));
heta1->Draw();
heta2->SetMarkerColor(kGreen);
heta2->SetMarkerStyle(kFullStar);
heta2->Scale(1./(htheta1->GetEntries()));
heta2->Draw("Same");
heta3->SetMarkerColor(kMagenta);
heta3->SetMarkerStyle(kCircle);
heta3->Scale(1./(htheta1->GetEntries()));
heta3->Draw("Same");
heta4->SetMarkerColor(kBlack);
heta4->SetMarkerStyle(kFullDotLarge);
heta4->Scale(1./(htheta1->GetEntries()));
heta4->Draw("Same");

TLegend *legend2=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend2->SetHeader("");
legend2->AddEntry(heta1,"all particles","p");
legend2->AddEntry(heta2,"protons","p");
legend2->AddEntry(heta3,"neutrons","p");
legend2->AddEntry(heta4,"charged particles","p");
legend2->Draw();


TCanvas *c3=new TCanvas("c3","c3",1080,720);
c3->cd();

pheta1->SetMarkerColor(kBlue);
pheta1->SetMarkerStyle(kFullDotLarge);
pheta1->Scale(1./(htheta1->GetEntries()));
pheta1->Draw();
pheta2->SetMarkerColor(kBlack);
pheta2->SetMarkerStyle(kFullStar);
pheta2->Scale(1./(htheta1->GetEntries()));
pheta2->Draw("Same");

TLegend *legend3=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend3->SetHeader("");
legend3->AddEntry(pheta1,"#pi^{-}","p");
legend3->AddEntry(pheta2,"#pi^{+}","p");
legend3->Draw();


TCanvas *c4=new TCanvas("c4","c4",1080,720);
c4->cd();

phtheta1->SetMarkerColor(kBlue);
phtheta1->SetMarkerStyle(kFullDotLarge);
phtheta1->Scale(1./(htheta1->GetEntries()));
phtheta1->Draw();
phtheta2->SetMarkerColor(kBlack);
phtheta2->SetMarkerStyle(kFullStar);
phtheta2->Scale(1./(htheta1->GetEntries()));
phtheta2->Draw("Same");

TLegend *legend4=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend4->SetHeader("");
legend4->AddEntry(phtheta1,"#pi^{-}","p");
legend4->AddEntry(phtheta2,"#pi^{+}","p");
legend4->Draw();

TCanvas *c5=new TCanvas("c5","c5",1080,720);
c5->cd();

hImpact->SetLineColor(kBlack);
hImpact->SetLineWidth(1);
Double_t scale = hImpact->GetXaxis()->GetBinWidth(5)/(hImpact->Integral());
hImpact->Scale(scale);

hImpact->Fit(f);
hImpact->Draw();
g->SetLineColor(kBlue);
g->Draw("Same");


return 0;

}//loop
