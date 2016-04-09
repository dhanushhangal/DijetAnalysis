#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <vector>
#include "TROOT.h"
#include "TChain.h"
#include <cmath>
#include "TStyle.h"
#include "TCanvas.h"

using namespace std;

void MyClass::Loop()
{
   if (fChain == 0) return;


const double ptmaxcut = 300. ;
const double ptmincut = 120. ;  
double leadingjetcut = 120. ;
double subleadingjetcut = 50. ;
double subsubleadingjetcut = 29.;
int third_highest_idx_gen = -1;
int second_highest_idx_gen = -1 ;
int highest_idx_gen = -1 ;
int inclusive_jets = 0;
const double trketamaxcut = 2.4;
const double etacut = 1.6; 
int low_Aj_events = 0;
int high_Aj_events = 0;
int dijet = 0;
int thirdjetevents = 0;
int lowAjthirdjet = 0;
int highAjthirdjet = 0;
const double Ajcut = 0.22 ;
int thirdjetwithetacut = 0;
int lowAjthirdjetwithetacut = 0;
int highAjthirdjetwithetacut = 0;
 int jetsml = 0;
 int jetlrg = 0;
 int jetsmlsub = 0;
 int jetlrgsub = 0;


TCanvas *cEta = new TCanvas("cEta","phi dist",600,400);
//TCanvas *cpT = new TCanvas("cpT","pT distribution",600,400);

TH1F *h1 = new TH1F("h1", "phi distribution", 180, 2, 4);
TH1F *h10 = new TH1F("h10", "jet size distribution", 6, 0, 5);
TH1F *h11 = new TH1F("h11", "jet size distribution in acc", 6, 0, 5);
TH2F *h5 = new TH2F("h5", "2D phi distribution", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
TH2F *h7 = new TH2F("h7", "2D d eta distribution", 100, -2, 2, 100, -2, 2);
TH2F *h8 = new TH2F("h8", "13 eta-phi correlation", 100, -5, 5, 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
TH2F *h9 = new TH2F("h9", "23 eta-phi correlation", 100, -5, 5, 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
TH1F *h2 = new TH1F("h2", "Aj distribution for all dijets", 6, 0, 0.66);
TH1F *h3 = new TH1F("h3", "pT dist for leading jet", 75, 25, 250);
TH1F *h4 = new TH1F("h4", "pT dist for subleading jet", 75, 25, 250);
TH1F *h6 = new TH1F("h6", "pT dist for subsubleading jet", 75, 25, 250);
TH1F *h12 = new TH1F("h12", "Aj dist for dijets with third jet", 6, 0, 0.66);
TH1F *h13 = new TH1F("h13", "Aj dist for dijets with third jet in eta cut", 6, 0, 0.66);
TH1F *hmixedeta12 = new TH1F("hmixedeta12","d eta12 dist for mixed events", 50, -3, 3);
TH1F *hmixedeta23 = new TH1F("hmixedeta23","d eta23 dist for mixed events", 50, -3, 3);
TH1F *hmixedeta13 = new TH1F("hmixedeta13","d eta13 dist for mixed events", 50, -3, 3);
TH2F *hmixedeta1213 = new TH2F("hmixedeta1213","2D d eta dist for mixed events", 100, -2, 2, 100, -2, 2);
TH1F *hmixedphi = new TH1F("hmixedphi","d phi23 dist for mixed events", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
TH1F *heta1 = new TH1F("heta1","eta distfor 1", 100, -5, 5);
TH1F *heta2 = new TH1F("heta2","eta distfor 2", 100, -5, 5);
TH1F *heta3 = new TH1F("heta3","eta distfor 3", 100, -5, 5);
TH1F *hphi1 = new TH1F("hphi1","phi distfor 1", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
TH1F *hphi2 = new TH1F("hphi2","phi distfor 2", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
TH1F *hpih3 = new TH1F("hphi3","phi distfor 3", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
TH1F *hdeta12 = new TH1F("hdeta12","d eta12 dist", 50, -3, 3);
TH1F *hdeta13 = new TH1F("hdeta13","d eta13 dist", 50, -3, 3);
TH1F *hdeta23 = new TH1F("hdeta23","d eta23 dist", 50, -3, 3);
TH1F *htrkphisml = new TH1F("htrkphisml","|eta_jet| < 0.5", 50,-TMath::Pi()/2, TMath::Pi()/2);
TH1F *htrkphilrg = new TH1F("htrkphilrg","0.5 < |eta_jet| < 1", 50,-TMath::Pi()/2, TMath::Pi()/2);
TH1F *htrkdphi = new TH1F("htrkphi","|d eta| < 0.5", 50,-TMath::Pi()/2, TMath::Pi()/2);
TH1F *htrkphisub = new TH1F("htrkphisub","|d eta| < 0.5", 50,-TMath::Pi()/2, TMath::Pi()/2);
TH2F *hjettrk = new TH2F("hjettrk","2D histo leading_jet-track", 50, -TMath::Pi()/2, TMath::Pi()/2, 50, -3, 3);
TH1F *htrkphismlsub = new TH1F("htrkphismlsub","|eta_jet| < 0.5", 50,-TMath::Pi()/2, TMath::Pi()/2);
TH1F *htrkphilrgsub = new TH1F("htrkphilrgsub","0.5 < |eta_jet| < 1", 50,-TMath::Pi()/2, TMath::Pi()/2);
TH1F *htrkphi = new TH1F("htrkphi","|d eta| < 0.5", 50,-TMath::Pi()/2, TMath::Pi()/2);
TH1F *htrkmixdphi = new TH1F("htrkmixdphi","|d eta| < 0.5", 50,-TMath::Pi()/2, TMath::Pi()/2);
TH1F *htrketa = new TH1F("htrketa","|d eta| < 0.5", 50, -3, 3);
TH1F *htrkmixdeta = new TH1F("htrkmixdeta","|d eta| < 0.5", 50, -3, 3);

  
Long64_t nentries = fChain->GetEntriesFast();

Long64_t nbytes = 0, nb = 0;
 


  for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);

      double lead_gen_pt = 0. ;
      double sublead_gen_pt = 0. ;
      double subsublead_gen_pt = 0. ;

      // if (jentry%1000 == 0){
      //  cout<<"event number: "<<jentry<<endl;
      //	}

      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // std::cout<<"  test "<<nbytes<<std::endl;      
  
      int jets = genpt->size();
      h10 -> Fill(jets); 

      //cout<<jets<<endl;
      //if (jets>100) break;      

      //search for leading jet
      for(int j = 0; j < jets ; j++) {
        
        double jet_pt= genpt->at(j);
        double eta = geneta->at(j);

        if(jet_pt <= leadingjetcut || abs(eta) >= etacut) continue ;
        if(jet_pt > lead_gen_pt){
          lead_gen_pt = jet_pt;
          highest_idx_gen = j;
          //cout<<"leadgen pT: "<<lead_gen_pt<<endl;
          //cout<<"lead jet: "<<highest_idx_gen<<endl;
          //cout<<"lead genphi: "<<genphi->at(j)<<endl;        
        }
	if(jet_pt >= leadingjetcut && abs(eta) <= etacut){
          inclusive_jets ++;
          if(fabs(eta) <= 0.5) jetsml++;
          if(fabs(eta) > 0.5 && abs(eta) <= 1) jetlrg++ ;
	}
      }
          
      //search for subleading jet
      for(int ijet = 0 ; ijet < jets ; ijet++){
        if(highest_idx_gen < 0) continue;
        if(ijet==highest_idx_gen) continue ;
        if(genpt->at(ijet) <= subleadingjetcut || abs(geneta->at(ijet)) >= etacut) continue ;
        if(genpt->at(ijet) > sublead_gen_pt){
          sublead_gen_pt = genpt->at(ijet);
          second_highest_idx_gen = ijet;
          //cout<<"sublead jet: "<<second_highest_idx_gen<<endl;
          //cout<<"subleadgen pT: "<<sublead_gen_pt<<endl;
          //cout<<"sub genphi: "<<genphi->at(ijet)<<endl;
          //heta2 -> Fill(geneta->at(second_highest_idx_gen));
        }
     
        if(sublead_gen_pt >= subleadingjetcut && abs(eta) <= etacut){
          if(fabs(eta) <= 0.5) jetsmlsub++;
          if(fabs(eta) > 0.5 && abs(eta) <= 1) jetlrgsub++ ;
	}

      }
      
      // jet track corr

      //leading jet and tracks
      if(lead_gen_pt >120){
        int tracks = trkPt->size();
        //cout<<"tracks in this event: "<<tracks<<endl;
        for(int trk = 0 ; trk < tracks ; trk++){
	  //cout<<trkPt->at(trk)<<endl;
	  if (trkPt->at(trk) <= 1 || trkPt->at(trk) >= 2) continue;
          if (fabs(trkPhi->at(trk) - genphi->at(highest_idx_gen)) >= TMath::Pi()/2) continue;
          //cout<<"track eta: "<< trkEta->at(trk) <<endl;

          if(fabs(trkEta->at(trk) - geneta->at(highest_idx_gen)) <= 0.5){
	    htrketa -> Fill(trkEta->at(trk));
	    
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            htrkmixdeta -> Fill(geneta->at(highest_idx_gen) - htrketa -> GetRandom());
            
            htrkphi -> Fill(trkPhi->at(trk));
            /*
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            htrkmixdphi -> Fill(genphi->at(highest_idx_gen) - htrkphi -> GetRandom());
            */
            hjettrk -> Fill(trkPhi->at(trk) - genphi->at(highest_idx_gen),trkEta->at(trk) - geneta->at(highest_idx_gen));
            htrkdphi -> Fill(trkPhi->at(trk) - genphi->at(highest_idx_gen));
            if (fabs(geneta->at(highest_idx_gen)) <= 0.5){
              htrkphisml -> Fill(trkPhi->at(trk) - genphi->at(highest_idx_gen));
            }
            if (fabs(geneta->at(highest_idx_gen)) > 0.5 && fabs(geneta->at(highest_idx_gen)) <= 1){
              htrkphilrg -> Fill(trkPhi->at(trk) - genphi->at(highest_idx_gen));
            }
          }
	}      
      }

      //subleading jet and tracks      
      if(lead_gen_pt >120 && sublead_gen_pt >50){
        int tracks = trkPt->size();
        //cout<<"tracks in this event: "<<tracks<<endl;
        for(int trk = 0 ; trk < tracks ; trk++){
	  //cout<<trkPt->at(trk)<<endl;
	  if (trkPt->at(trk) <= 1 || trkPt->at(trk) >= 2) continue;
          if (fabs(trkPhi->at(trk) - genphi->at(second_highest_idx_gen)) >= TMath::Pi()/2) continue;
          //cout<<"track eta: "<< trkEta->at(trk) <<endl;
          if(fabs(trkEta->at(trk) - geneta->at(second_highest_idx_gen)) <= 0.5){
	    //hjettrk -> Fill(trkPhi->at(trk) - genphi->at(highest_idx_gen),trkEta->at(trk) - geneta->at(highest_idx_gen));
            htrkphisub -> Fill(trkPhi->at(trk) - genphi->at(second_highest_idx_gen));
            if (fabs(geneta->at(second_highest_idx_gen)) <= 0.5){
              htrkphismlsub -> Fill(trkPhi->at(trk) - genphi->at(second_highest_idx_gen));
            }
            if (fabs(geneta->at(second_highest_idx_gen)) > 0.5 && fabs(geneta->at(second_highest_idx_gen)) <= 1){
              htrkphilrgsub -> Fill(trkPhi->at(trk) - genphi->at(second_highest_idx_gen));
            }
          }
	}      
      }

     
      //condition for two jets
      if(lead_gen_pt < 120 || sublead_gen_pt < 50){

        highest_idx_gen = -1;
        second_highest_idx_gen = -1;

      }
      
      
      if(highest_idx_gen> -1 && second_highest_idx_gen> -1 ){

        //defining phi and eta
        double dphi12 = genphi->at(highest_idx_gen) - genphi->at(second_highest_idx_gen);
        if(dphi12 >= 3*TMath::Pi()/2) dphi12 = dphi12 - 2*TMath::Pi();
        if(dphi12 <= -TMath::Pi()/2) dphi12 = dphi12 + 2*TMath::Pi();
        double deta12 = geneta->at(highest_idx_gen) - geneta->at(second_highest_idx_gen); 

        //applying phi cut      
        if (abs(dphi12) >= (5*(TMath::Pi())/6) && abs(dphi12 <= (7*(TMath::Pi())/6))){    
	  //cout<<"event has dijet"<<endl;
          dijet ++;
        
          //defining Aj
          double Aj = (lead_gen_pt-sublead_gen_pt)/(lead_gen_pt+sublead_gen_pt);
          h2 -> Fill(Aj);
          //cAj->Update();

          if(Aj > Ajcut) high_Aj_events ++;
          else low_Aj_events ++;

          //checking for third jet
          if (jets > 2){
            thirdjetevents ++;
            h12 -> Fill(Aj);
            if(Aj > Ajcut) highAjthirdjet ++;
            else lowAjthirdjet ++;

            for(int l = 0; l < jets ; l++) {
              if(l==highest_idx_gen) continue ;
              if(l==second_highest_idx_gen) continue ;
              //cout<<"geneta :"<<geneta->at(l)<<endl;
              if(fabs(geneta->at(l)) >= etacut) continue ;
            
              if(genpt->at(l) > subsublead_gen_pt){
                subsublead_gen_pt = genpt->at(l);
                third_highest_idx_gen = l;
                //cout<<"geneta :"<<geneta->at(l)<<endl;
                //cout<<"subsubleadgen pT: "<<subsublead_gen_pt<<endl;
                //cout<<"subsub lead genphi: "<<genphi->at(l)<<endl;
	      }
           
	    }

	    if(subsublead_gen_pt > 0){ 
	       thirdjetwithetacut ++;
	       if(Aj < Ajcut) lowAjthirdjetwithetacut++;
	       else highAjthirdjetwithetacut++;
               double dphi13 = genphi->at(highest_idx_gen) - genphi->at(third_highest_idx_gen);
               if(dphi13 >= 3*TMath::Pi()/2) dphi13 = dphi13 - 2*TMath::Pi();
               if(dphi13 <= -TMath::Pi()/2) dphi13 = dphi13 + 2*TMath::Pi();
               double deta13 = geneta->at(highest_idx_gen) - geneta->at(third_highest_idx_gen);
               double dphi23 = genphi->at(second_highest_idx_gen) - genphi->at(third_highest_idx_gen);
               if(dphi23 >= 3*TMath::Pi()/2) dphi23 = dphi23 - 2*TMath::Pi();
               if(dphi23 <= -TMath::Pi()/2) dphi23 = dphi23 + 2*TMath::Pi();
               double deta23 = geneta->at(second_highest_idx_gen) - geneta->at(third_highest_idx_gen);

               hdeta12 -> Fill(geneta->at(highest_idx_gen) - geneta->at(second_highest_idx_gen));
               hdeta13 -> Fill(geneta->at(highest_idx_gen) - geneta->at(third_highest_idx_gen));
               hdeta23 -> Fill(geneta->at(second_highest_idx_gen) - geneta->at(third_highest_idx_gen));
              
               heta1 -> Fill(geneta->at(highest_idx_gen));
               heta2 -> Fill(geneta->at(second_highest_idx_gen));
               heta3 -> Fill(geneta->at(third_highest_idx_gen));
               
               hphi1 -> Fill(genphi->at(highest_idx_gen));
               hphi2 -> Fill(genphi->at(second_highest_idx_gen));
               hphi3 -> Fill(genphi->at(third_highest_idx_gen));  
               
               h13 -> Fill(Aj);
               h5 -> Fill(dphi12,dphi13);
               h7 -> Fill(deta12,deta13);
               h8 -> Fill(deta13,dphi13);
               h9 -> Fill(deta23,dphi23);
	      }
      
	  }
      
        if((genpt->at(highest_idx_gen)<= leadingjetcut )||
           //(genpt->at(third_highest_idx_gen)<= subsubleadingjetcut )||
	   //(genpt->at(highest_idx_gen)>= ptmaxcut ) ||
	   //(genpt->at(highest_idx_gen)<= ptmincut ) ||
	   (genpt->at(second_highest_idx_gen)<= subleadingjetcut ))
           //(TMath::Abs(geneta->at(highest_idx_gen)) >= etacut ) ||
	   //(TMath::Abs(geneta->at(second_highest_idx_gen)) >= etacut ))
        {
	  highest_idx_gen = -1;  
	  second_highest_idx_gen = -1; 
	}                	
	
        h1 -> Fill(dphi12);
        //cpT -> Update();
        h3 -> Fill(lead_gen_pt);
        h4 -> Fill(sublead_gen_pt);
        h6 -> Fill(subsublead_gen_pt);
	}
    
      }
  } 

  //for (Int_t i=0 ; i<1000000 ; i++) {
  
  //hmixedeta1213 -> Fill(heta1->GetRandom() - heta2->GetRandom(), heta1->GetRandom() - heta3->GetRandom());  
    //hmixedeta23 -> Fill(heta2->GetRandom() - heta3->GetRandom());
    
  //}
  
  h3 -> Sumw2();
  h4 -> Sumw2();
  hdeta23 -> Sumw2();
  hmixedeta23 -> Sumw2();
  htrkphisub -> Sumw2();
  htrkphismlsub -> Sumw2();
  htrkphilrgsub -> Sumw2();
  //h1 -> Draw();
  //gStyle -> SetOptStat(0);
  //h2 -> SetFillColorAlpha (33, 0.35);
  //h2 -> Draw("");
  h3 -> SetLineColor(kRed);
    //h3 -> Draw();
  h4 -> SetLineColor(kBlue);
    //h4 -> Draw("same");     
  //h5 -> GetXaxis() -> SetTitle("dphi12")
  //h5 -> Draw("COLZ");
  h6 -> SetLineColor(kGreen);
    //h6 -> Draw("same");
  //h7 -> Draw("COLZ");
  //h8 -> Draw("COLZ");
  //h9 -> Draw("COLZ");
  //h10 -> Draw();
  h12 -> SetFillColorAlpha(46, 0.35);  
  //h12 -> SetLineColor(kRed);
  //h12 -> Scale(scale);
  //h12 -> Draw("same");
  h13 -> SetFillColorAlpha(41, 0.35);
  //h13 -> Draw("same");
  //TH1D * projh2X = h8->ProjectionX();
  //projh2X -> Draw();
  //hmixedeta13 -> Scale(1./38925);
  //hmixedeta23 -> Scale(1./38719);
  //hmixedeta23 -> Draw();
  //heta2 ->Draw();
  //hmixedeta1213 -> Scale(1./1150);
  //h7 -> Divide(hmixedeta1213);
  //hdeta23 -> Draw();
  //h7 -> Draw("COLZ");
  //TH1D * hprjxcorreta = h7->ProjectionX();
  //TH1D * hprjycorreta = h7->ProjectionY();
  //hprjxcorreta -> Sumw2();
  //hprjxcorreta -> Draw();
  htrkphismlsub -> Scale(1./(jetsmlsub));
  htrkphilrgsub -> Scale(1./(jetlrgsub));
  //cEta -> Update();
  //htrkphismlsub -> Draw();
  //htrkphilrgsub -> Draw("same");
  //htrkphisub->Draw();
  /*  TAxis *xaxis = hjettrk->GetXaxis();
  TAxis *yaxis = hjettrk->GetYaxis();
  Int_t biny1 = yaxis->FindBin(-0.5);
  Int_t biny2 = yaxis->FindBin(0.5);
  Int_t binx1 = xaxis->FindBin(-TMath::Pi()/2);
  Int_t binx2 = xaxis->FindBin(TMath::Pi()/2);*/
  //cout<<biny1<<" and " <<biny2<<endl;
  //TH1D * hjettrkx = hjettrk->ProjectionX(" ",biny1,biny2,"o");
  //TH1D * hjettrky = hjettrk->ProjectionY(" ",binx1,binx2,"o");
  //hjettrky -> Sumw2();
  //hjettrky -> Draw();
  //hjettrk -> Draw("COLZ");
  //htrkphi -> Draw();
  //htrketa -> Draw();
  htrkmixdeta -> Sumw2();
  htrkmixdeta -> Scale(1./20340); 
  //hjettrky -> Divide(htrkmixdeta);  
  //hjettrky -> Sumw2();
  //hjettrky -> Draw();
  //htrkmixdphi -> Sumw2();
  //htrkmixdphi -> Scale(1./20340);
  //htrkmixdphi -> Draw();


//cout<< nentries <<endl;      
cout<<"inclusive jets: "<<inclusive_jets<<endl;
cout<<"dijets: "<<dijet<<endl;
cout<<"low_Aj_events: "<<low_Aj_events<<endl;
cout<<"high_Aj_events: "<<high_Aj_events<<endl;
cout<<"third jet events: "<<thirdjetevents<<endl;
cout<<"low Aj third jet events: "<<lowAjthirdjet<<endl;
cout<<"high Aj third jet events: "<<highAjthirdjet<<endl;
cout<<"third jet with eta cut:"<<thirdjetwithetacut<<endl;
cout<<"lowAjthirdjetwithetacut:"<<lowAjthirdjetwithetacut<<endl;
cout<<"highAjthirdjetwithetacut:"<<highAjthirdjetwithetacut<<endl;


}
