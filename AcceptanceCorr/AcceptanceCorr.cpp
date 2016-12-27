#include <iostream>
#include <string>
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TAxis.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TCut.h"
#include "TEventList.h"

TString dataLoc;
TString fDataExt;
TString fSimuExt;
TString elecExt;
TString pionExt;
TString metals[6] = {"C", "Fe", "Pb", "DC", "DFe", "DPb"};
TString Metal;

Int_t nSimuFiles;
Int_t met = 6;

Double_t Q2_MIN;
Double_t Q2_MAX;
Double_t XB_MIN;
Double_t XB_MAX;
Double_t NU_MIN;
Double_t NU_MAX;
Double_t ZH_MIN;
Double_t ZH_MAX;
Double_t PT_MIN;
Double_t PT_MAX;
Double_t PHI_MIN;
Double_t PHI_MAX;

Int_t N_Q2;
Int_t N_XB;
Int_t N_NU;
Int_t N_ZH;
Int_t N_PT;
Int_t N_PHI;

Double_t delta_Q2;
Double_t delta_XB;
Double_t delta_NU;
Double_t delta_ZH;
Double_t delta_PT;
Double_t delta_PHI;

Double_t *v_Q2;
Double_t *v_XB;
Double_t *v_NU;
Double_t *v_ZH;
Double_t *v_PT;
Double_t *v_PHI;

void getBinning(char **argv);
void genAcceptance(TFile *f, TString Metal);
TCut getCutData(TString Metal, Int_t q2i, Int_t xbi, Int_t pti, Int_t zhi, Int_t phii);
TCut getCutSimul(Int_t q2i, Int_t xbi, Int_t pti, Int_t zhi, Int_t phii);

// Create a TNtuple with data, accepted, thrown and acceptance in the binning specified 1 dimensional
// Also possible to make 1 dimensional with Nu instead of Xb
int main(int argc, char **argv){
	if(argc < 25){
		std::cout << "The number of arguments is incorrect" << std::endl;
		std::cout << argc << std::endl;
		return 0;
	}
		
	getBinning(argv);
	
	dataLoc = (TString) argv[19];
	fDataExt = (TString) argv[20];
	fSimuExt = (TString) argv[21];
	nSimuFiles = (Int_t) std::stoi(argv[22]);
	elecExt = (TString) argv[23];
	pionExt = (TString) argv[24];
	
	TFile *f = new TFile("out.root", "RECREATE");
    
    for(Int_t i = 0; i < met; i++){
    	Metal = metals[i];
    	genAcceptance(f, Metal);
    }
    
    f->Close();
	delete f;
	
	delete v_Q2;
	delete v_XB;
	delete v_ZH;
	delete v_PT;
	delete v_PHI;
	delete v_NU;

	return 0;
}

void getBinning(char **argv){
	Q2_MIN = (Double_t) std::stod(argv[1]);
	Q2_MAX = (Double_t) std::stod(argv[2]);
	N_Q2 = (Int_t) std::stoi(argv[3]);
	XB_MIN = (Double_t) std::stod(argv[4]);
	XB_MAX = (Double_t) std::stod(argv[5]);
	N_XB = (Int_t) std::stoi(argv[6]);
	NU_MIN = (Double_t) std::stod(argv[7]);
	NU_MAX = (Double_t) std::stod(argv[8]);
	N_NU = (Int_t) std::stoi(argv[9]);
	ZH_MIN = (Double_t) std::stod(argv[10]);
	ZH_MAX = (Double_t) std::stod(argv[11]);
	N_ZH = (Int_t) std::stoi(argv[12]); 
	PT_MIN = (Double_t) std::stod(argv[13]);
	PT_MAX = (Double_t) std::stod(argv[14]);
	N_PT = (Int_t) std::stoi(argv[15]);
	PHI_MIN = (Double_t) std::stod(argv[16]);
	PHI_MAX = (Double_t) std::stod(argv[17]);
	N_PHI = (Int_t) std::stoi(argv[18]);
	
	delta_Q2 = (Q2_MAX-Q2_MIN)/N_Q2;
	delta_XB = (XB_MAX-XB_MIN)/N_XB;
	delta_ZH = (ZH_MAX-ZH_MIN)/N_ZH;
	delta_PT = (PT_MAX-PT_MIN)/N_PT;
	delta_PHI = (PHI_MAX-PHI_MIN)/N_PHI;
	delta_NU = (NU_MAX-NU_MIN)/N_NU;
	v_Q2 = new Double_t[N_Q2+1];
	v_NU = new Double_t[N_NU+1];
	v_ZH = new Double_t[N_ZH+1];
	v_XB = new Double_t[N_XB+1];
	v_PT = new Double_t[N_PT+1];
	v_PHI = new Double_t[N_PHI+1];
	for(Int_t i = 0; i < N_Q2+1; i++){
        if(i == 0) v_Q2[i] = Q2_MIN;
        else v_Q2[i] = v_Q2[i-1] + delta_Q2;
    }
    for(Int_t i = 0; i < N_XB+1; i++){    
        if(i == 0) v_XB[i] = XB_MIN;
        else v_XB[i] = v_XB[i-1] + delta_XB;    
    }
    for(Int_t i = 0; i < N_ZH+1; i++){
        if(i == 0) v_ZH[i] = ZH_MIN;
        else v_ZH[i] = v_ZH[i-1] + delta_ZH;
    }
    for(Int_t i = 0; i < N_PT+1; i++){
        if(i == 0) v_PT[i] = PT_MIN;
        else v_PT[i] = v_PT[i-1] + delta_PT;
    }    
    for(Int_t i = 0; i < N_PHI+1; i++){
        if(i == 0) v_PHI[i] = PHI_MIN;
        else v_PHI[i] = v_PHI[i-1] + delta_PHI;
    }
}

void genAcceptance(TFile *f, TString Metal){
	TString simuMetal;
    Int_t BinTot = N_Q2*N_XB*N_PT*N_ZH*N_PHI;
    TCut cut;
    
    TFile *fPion = new TFile(dataLoc + Metal + fDataExt + pionExt);
    TNtuple *ntuplePion = (TNtuple*) fPion->Get("data_pion");
    
    if(Metal.Contains("D")) simuMetal = "D";
    else simuMetal = Metal;
    TChain *accept = new TChain("accept_pion");
	for(Int_t q = 0; q < nSimuFiles; q++)
		accept->Add(dataLoc + simuMetal + std::to_string(q+1) + fSimuExt + pionExt);
	accept->SetEstimate(accept->GetEntries());
		
	TChain *thrown = new TChain("thrown_pion");
	for(Int_t q = 0; q < nSimuFiles; q++)
		thrown->Add(dataLoc + simuMetal + std::to_string(q+1) + fSimuExt + pionExt);
	thrown->SetEstimate(thrown->GetEntries());
    
    TH1F *hData = new TH1F("hData", "", BinTot, 0, BinTot);
    TH1F *hAcc = new TH1F("hAcc", "", BinTot, 0, BinTot);
    TH1F *hThr = new TH1F("hThr", "", BinTot, 0, BinTot);
    TH1F *hAcceptance = new TH1F("hAcceptance", "", BinTot, 0, BinTot);
    TEventList *el;
    Int_t idx;
    
    // The sequence of the binning is the following
    // 1) Q2 2) Xb 3) Pt 4) Zh 5) Phi
    
    for(Int_t q2i = 0; q2i < N_Q2; q2i++){
    	for(Int_t xbi = 0; xbi < N_XB; xbi++){
    		for(Int_t pti = 0; pti < N_PT; pti++){
    			for(Int_t zhi = 0; zhi < N_ZH; zhi++){
    				for(Int_t phii = 0; phii < N_PHI; phii++){
    					idx = q2i*N_XB*N_ZH*N_PT*N_PHI+xbi*N_ZH*N_PT*N_PHI+pti*N_PT*N_PHI+zhi*N_PHI+phii;
    					
    					cut = getCutData(Metal, q2i, xbi, pti, zhi, phii);
    					ntuplePion->Draw(">>listdata", cut, "goff");
						el = (TEventList*) gDirectory->Get("listdata");
    					hData->Fill(idx+0.5, el->GetN());
    					
    					cut = getCutSimul(q2i, xbi, pti, zhi, phii);
    					accept->Draw(">>listacc", cut, "goff");
    					el = (TEventList*) gDirectory->Get("listacc");
    					hAcc->Fill(idx+0.5, el->GetN());
    					thrown->Draw(">>listthr", cut, "goff");
    					el = (TEventList*) gDirectory->Get("listthr");
    					hThr->Fill(idx+0.5, el->GetN());
    				}
    			}
    		}
    	}
    }
    
    hData->Sumw2();
    hAcc->Sumw2();
    hThr->Sumw2();
    hAcceptance->Divide(hAcc, hThr, 1, 1);
    hAcceptance->Divide(hData, hAcceptance, 1, 1);
    
    f->cd();
    hData->Write(Form("hData%s", (const char*) Metal));
    hAcc->Write(Form("hAcc%s", (const char*) Metal));
    hThr->Write(Form("hThr%s", (const char*) Metal));
    
    return;
}

TCut getCutData(TString Metal, Int_t q2i, Int_t xbi, Int_t pti, Int_t zhi, Int_t phii){
	TCut cut;
	TCut liquid = "TargType==1";
    TCut solid = "TargType==2";
    TCut q2cut, xbcut, ptcut, zhcut, phicut;
    if(Metal == "C" || Metal == "Fe" || Metal == "Pb") cut = solid;
    else cut = liquid;
    q2cut = Form("Q2>%f && Q2<%f", v_Q2[q2i], v_Q2[q2i+1]);
    xbcut = Form("Xb>%f && Xb<%f", v_XB[xbi], v_XB[xbi+1]);
    ptcut = Form("Pt>%f && Pt<%f", v_PT[pti], v_PT[pti+1]);
    zhcut = Form("Zh>%f && Zh<%f", v_ZH[zhi], v_ZH[zhi+1]);
    phicut = Form("Phi>%f && Phi<%f", v_PHI[phii], v_PHI[phii+1]);
    cut = cut && q2cut && xbcut && ptcut && zhcut && phicut;
    
    return cut;   
}

TCut getCutSimul(Int_t q2i, Int_t xbi, Int_t pti, Int_t zhi, Int_t phii){
	TCut cut;
	TCut q2cut, xbcut, ptcut, zhcut, phicut;
	q2cut = Form("Q2>%f && Q2<%f", v_Q2[q2i], v_Q2[q2i+1]);
    xbcut = Form("Xb>%f && Xb<%f", v_XB[xbi], v_XB[xbi+1]);
    ptcut = Form("Pt>%f && Pt<%f", v_PT[pti], v_PT[pti+1]);
    zhcut = Form("Zh>%f && Zh<%f", v_ZH[zhi], v_ZH[zhi+1]);
    phicut = Form("Phi>%f && Phi<%f", v_PHI[phii], v_PHI[phii+1]);
    cut = cut && q2cut && xbcut && ptcut && zhcut && phicut;
    
    return cut;   
}

