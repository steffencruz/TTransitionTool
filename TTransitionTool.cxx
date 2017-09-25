#include "TTransitionTool.h"

#include <TFile.h>
#include <TLegend.h>
#include<TLegendEntry.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <TArrow.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TVirtualFitter.h>

#include<fstream>
#include <cmath>
#include <algorithm>    // std::count
#include <vector>     

#include <TTigress.h>  


ClassImp(TTransitionTool)

std::string TTransitionTool::nndcfile = "";
	
////////////////////////////////////////////////////////////////////////////////

Bool_t TTransitionTool::verbose = true;
TList *TTransitionTool::list = 0;

Double_t TTransitionTool::ExcSig = 1.0;
Double_t TTransitionTool::GamSig = 1.0; //   new

Double_t TTransitionTool::NSig = 1.5;

std::vector<double> TTransitionTool::energies;
std::vector<int> TTransitionTool::states;
std::map<int,std::vector<int> > TTransitionTool::cascades;
TH1D *TTransitionTool::hint = NULL;
TH2F *TTransitionTool::htrans = NULL;
TH2F *TTransitionTool::hgams = NULL;
TH2F *TTransitionTool::hseq = NULL;
Int_t TTransitionTool::nlines = 0;
Int_t TTransitionTool::nsequences = 0;
Int_t TTransitionTool::ncascades = 0;
Int_t TTransitionTool::nstates = 0;
TF1 *TTransitionTool::fTigSigma = NULL;
TF1 *TTransitionTool::fTigEff = NULL;
TGraphErrors *TTransitionTool::gTigRes = NULL;
TGraphErrors *TTransitionTool::gTigEff = NULL;


TTransitionTool::TTransitionTool()	{	
}

TTransitionTool::~TTransitionTool()	{	}

void TTransitionTool::Print(Option_t *opt) {

	printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
	printf("\n\n\t____TTransitionTool____");
  printf("\n\t NNDC File    :- %s",nndcfile.c_str());
  printf("\n\t Verbose      :- %s",verbose?"TRUE":"FALSE");

  int from_state = atoi(opt);
  if(from_state && from_state<nstates){
   // PrintDecays(from_state,true);
    PrintCascades(from_state,0.0,true);
  }

	printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n\n");
		
}


std::vector<int> TTransitionTool::PrintStates(Double_t emin, Double_t emax){

	std::vector<int> vals;
	if(!energies.size()){
		printf("\n\n No States Have Been Loaded !!\n\n");
		return vals;
	}
	if(!hint) SetIntensity(0,1);
	
	printf("\n\n List of States:-");

	for(int i=0; i<(Int_t)energies.size(); i++){
		if(energies.at(i)<emin)
			continue;
					
		if(energies.at(i)>emax)
			break;
						
		vals.push_back(i+1);	
		printf("\n\tState %2i =   %6.1f keV",i+1,energies.at(i));
		if(verbose) printf("\t\t[ Pop. Strength = %.3E ]",hint->GetBinContent(i));
	}
	printf("\n\n");
	return vals;
}

//std::vector<int>  TTransitionTool::PrintDecays(Int_t from_state, Double_t egam, Bool_t printeng)

std::vector<int>  TTransitionTool::PrintCascades(Int_t from_state, Double_t egam, Bool_t printeng){

	std::string msgam = "";
	if(egam) msgam = Form("with %.1f keV gamma gate ",egam);
	printf("\n\n List of Decays from state %2i %s%s:-",from_state,msgam.c_str(),printeng?"{With Gamma Energy}":"");
	std::vector<int> indx = GetCascadeIndex(from_state,egam);
	Int_t ncas = indx.size();
	int j=0;
	
	for(int k=0; k<ncas; k++){

		j = indx.at(k);
		printf("\n\tCascade %2i = ",k+1);
		
		if(printeng) 
			printf(" %6.1f ",energies.at(cascades[j].at(0)-1));
		else 
			printf("  %2i",cascades[j].at(0));
			
		for(int i=1; i<cascades[j].size(); i++){
			if(printeng) {
				printf(" -{ %6.1f }->  %6.1f ",	energies.at(cascades[j].at(i-1)-1)-energies.at(cascades[j].at(i)-1)
																			,	energies.at(cascades[j].at(i)-1));
			} else
				printf("  -> %2i",cascades[j].at(i));
		}
		if(printeng)printf(" keV * * *");
	}
	
	printf("\n\n");
	return indx;
}


// State strength function gives the strength of the state after a gamma gate. It is normalized.
// CalcStates()
// Use resolution and make realistic

void TTransitionTool::ClearVars(Option_t *opt){
	
	nndcfile.clear();
	list = new TList();
	hseq = 0;
	hint = 0;
	nsequences = 0;

	states.clear();	
	energies.clear();
	htrans = 0;	
	hgams  = 0;
	nlines = 0;
	nstates = 0;	

}


Bool_t TTransitionTool::Init(std::string fname, Int_t nmax, Bool_t verb){
//	gStyle->SetOptStat(0);
//	gStyle->SetTitleXOffset(1.3);
//	gStyle->SetTitleYOffset(1.5);	
	
	SetVerbose(verb);
	Bool_t success = LoadLevelsGammas(fname.c_str(),nmax);
	
	if(!fTigSigma) SetResolutionCurve();
  if(!fTigEff)   SetEfficiencyCurve();

	return success;
}

Bool_t TTransitionTool::LoadLevelsGammas(std::string fname, int nmax){

	ClearVars();
	
	std::ifstream infile(fname.c_str());	
	if(!infile.is_open()){ // try looking in this directory for it
	  fname = Form("%s/%s",TIGDIR,fname.c_str());
	  infile.open(fname.c_str());
  }
	if(!infile.is_open()){
		printf("\nFAILED TO OPEN FILE :  %s\n\n",fname.c_str());
		return false;
	}
	
	nndcfile.assign(fname);
	
	nlines = std::count(std::istreambuf_iterator<char>(infile),std::istreambuf_iterator<char>(), '\n')+1;
	infile.close();
	infile.open(fname.c_str());	
	
	GetRid("TmpMat");
	GetRid("TmpMat2");
	TH2F *htmp = new TH2F("TmpMat","TmpMat",nlines,0,nlines,nlines,0,nlines);
	TH2F *htmp2 = new TH2F("TmpMat2","TmpMat2",nlines,0,nlines,10000,0,10000);

	Double_t initial, final, intensity, egamma, maxgam=0;
	int statenum;
	nstates = 0;
	ncascades = 0;
	
	printf("\n Reading Input File ' %s ' :\n",fname.c_str());
	if(verbose) printf("\n\t   State Energy\tGamma Energy\tIntensity[%%]\tFinal State [num]");
	while(infile.good()){
		
		infile >> initial >> egamma >> intensity >> final;

		GetStateIndex(initial,true); // add new state, but don't duplicate		

		statenum = GetStateIndex(final)+1;
		
		if(verbose) printf("\t\t%5.1f\t%8.1f\t%9.1f\t%7.1f   [%i]\n",initial,egamma,intensity,final,statenum);
		
		if(statenum==0){ // final state doesn't match any state
			printf("\n\t! Error :  Final state %.1f not valid !\n\n",final);
			continue;
		}
				
		htmp->SetBinContent(nstates,statenum,intensity);	// fill intensity-transition matrix	
		if(maxgam<egamma)
			maxgam = egamma;
			
		htmp2->SetBinContent(nstates,egamma,intensity);  // fill gammas histogram
		
		if(nstates==nmax){
			printf("\n..............................................................\n\t\tSTOPPED READING [reached %i states] ",nstates);			
			break;
		}		
	}
	printf("\n\t --- Complete! [eof = %i  fail = %i  bad = %i]\n\n",infile.eof(),infile.fail(),infile.bad());
	infile.close();
	
	GetRid("TransitionMat");
	htrans = new TH2F("TransitionMat","",nstates,0,nstates,nstates,0,nstates);
	htrans->SetTitle("Transition Intensity Matrix; Initial State; Final State");
	
	for(int i=1; i<=nstates; i++) // go to the maximum row
		for(int j=1; j<=i; j++) // go up to the diagonal
			htrans->SetBinContent(i,j,htmp->GetBinContent(i,j));
	
	htrans->Scale(1./htrans->GetMaximum()); // normalize to max strength 1

	GetRid("GammasMat");
	hgams = new TH2F("GammasMat","",nstates,0,nstates,6000.0,0,6000.0);
	hgams->SetTitle("Gammas Emitted Matrix; State; E_{#gamma} [keV]");

	for(int i=1; i<=nstates; i++){ // go to the maximum row
		TH1D *hprojy = htmp2->ProjectionY("",i,i);
		for(int j=hprojy->FindFirstBinAbove(); j<=hprojy->FindLastBinAbove(); j++){ // go up to the diagonal
			if(htmp2->GetBinContent(i,j))
				hgams->SetBinContent(i,j,htmp2->GetBinContent(i,j));
		}
	}	
	hgams->Scale(1./hgams->GetMaximum()); // normalize all transitions to max strength 1
		
	for(int i=1; i<=nstates; i++)
		BuildDecayScheme(i);
				
	return true;
}

TF1 *TTransitionTool::SetEfficiencyCurve(double eff, double eng){

  fTigEff = new TF1("func_eff","[0]*pow(10.,[1]*log10(x)+[2]*pow(log10(x),2.)+[3]*pow(1/x,2.))",0,4000);
  fTigEff->SetNpx(2000);
  fTigEff->SetLineWidth(2);
  fTigEff->SetParameters(313.754,-0.681343,0.0418179,-3154.89);	
  fTigEff->SetParameter(0,fTigEff->GetParameter(0)*eff/fTigEff->Eval(eng));
  printf("\n Efficiency curve has been set so that Eff @ %.1f keV = %.2f %%\n\n",eng,eff*100.0);
  
  fTigEff->SetNameTitle("EffCurve","Gamma Efficiency Curve; E_{#gamma} [keV]; #varepsilon_{#gamma}");
       
  TF1 *func = (TF1*)fTigEff->Clone("Efficiency");   
  return func;
}

Double_t TTransitionTool::Efficiency(Double_t eng){

	if(!fTigEff)
    SetEfficiencyCurve();
	
	return fTigEff->Eval(eng)*0.01; 
}

TF1 *TTransitionTool::SetResolutionCurve(double sig, double eng){

	fTigSigma = new TF1("func_sig","[0]*pow(x,0.5)",0,4000); // power law?
	fTigSigma->SetNpx(2000);
  fTigSigma->SetLineWidth(2);
	fTigSigma->SetParameter(0,1.0); // from fitting singles data	
  fTigSigma->SetParameter(0,fTigSigma->GetParameter(0)*sig/fTigSigma->Eval(eng)); 
	fTigSigma->SetNameTitle("ResCurve","Gamma Resolution Curve; E_{#gamma} [keV]; #sigma_{#gamma} [keV]");
  printf("\n Resolution curve has been set so that Sigma @ %.1f keV = %.2f keV\n",eng,sig);

  TF1 *func = (TF1*)fTigSigma->Clone("Resolution");     
  return func;
}

Double_t TTransitionTool::Resolution(Double_t eng, Bool_t FWHM){

	if(!fTigSigma)
	  SetResolutionCurve();
	  
	Double_t sig = fTigSigma->Eval(eng);
	if(FWHM)
	  sig*=2.355;
  return sig;
}


////////////////////////////////////////////////////////////////////////////////
// NOTES

// Gate on exc energy range produce gamma spectrum. Compare to data and eliminate what we dont see experimentally
// normalize theory to specific peak (direct to ground if possible)

// ACCOUNTS FOR THE FACT THAT THE REAL DATA HAS GAMMAS FROM STATES OUTSIDE OF THE SELECTED RANGE
// BY INCLUDING ALL STATES THAT ARE UP TO ±1.5 SIGMA BEYOND THE SPECIFIED EXC RANGE AND
// WEIGHTING THEIR GAMMA SPECTRA BY THE COUNTS EXPECTED WITHIN THE EXC RANGE [MIN. ~7%]

////////////////////////////////////////////////////////////////////////////////

Double_t TTransitionTool::BranchingRatio(Double_t state_eng, Double_t egam){
// This allows gates on transitions such as B->C in the cascade A->B->C

	Int_t from_state = GetStateIndex(state_eng,false);
	if(!BuildDecayScheme(from_state+1))
		return 0;
	
	Double_t eng1 = energies.at(from_state);	
  TH1D *hg = CalcGammas(from_state+1); // uses gamma spectrum
  //hg->Draw();
  Double_t val = hg->Integral((int)egam-1,(int)egam+1);
	if(!val){
	  printf("\n\t Error :  The transition %.1f keV -> %.1f keV was not found in file.\n\n",state_eng,state_eng-egam);
	  return 0;
  }  
  
	if(verbose)printf("\n\t %.1f keV gamma is produced from %.1f keV state cascades with total intensity = %.2e\n",egam,eng1,val);
  
  return val;  
}		

Bool_t TTransitionTool::SetIntensity(Int_t state, Double_t strength){

  if(!energies.size()){
    printf("\n\t! Error : No states have been loaded yet. Intensities cannot be set now.\n\n"); 
    return false;
  }

	if(!hint){	
		hint = new TH1D("StateIntensites","",nstates,0,nstates);
		hint->SetTitle("State Intensities Set Using Data; State Number; Populated Strength");
		for(int i=1; i<=nstates; i++)
			hint->SetBinContent(i,strength);
	} else if(!state){ // fill all bins with given strength
		for(int i=1; i<=nstates; i++)
			hint->SetBinContent(i,strength);
	}
	
	hint->SetBinContent(state,strength);
	return true;
}	

Bool_t TTransitionTool::ReadIntensities(const char *fname){

	std::ifstream infile(fname);
	if(!infile.is_open()){
		printf("\n\t! Error : Could not locate file ' %s '\n\n",fname);
		return false;	
	}
	
	Bool_t success = SetIntensity(0,0);// initiate the histogram
	if(!success)
	  return false;
	  
	Int_t bin;
	Double_t eng, val;
	while(infile.good()){
		infile >> bin >> eng >> val;
		hint->SetBinContent(bin,val);
	}
	printf("\n\n\tRead in file ' %s ' and set TH1D*'hint' intensities.\n\n",fname);
	return true;
}

void TTransitionTool::WriteIntensities(const char *fname){

	if(!hint){
		printf("\n\t Error :  At least one intensity must be set!");
		return;
	}

	std::ofstream ofile(fname);
	
	for(int i=1; i<hint->GetNbinsX(); i++)
		ofile << i << "\t" << energies.at(i-1) << "\t" << hint->GetBinContent(i) << "\n";

	printf("\n\n\tWrote Contents of TH1D*'hint' to file!\n\n");
	ofile.close();
}

TH1D *TTransitionTool::CalcGammas(Int_t from_state, Double_t egam){
	
	TH1D *hgcalc;
	if(!nndcfile.size()){
	  printf("\n Error :  Load a file first.\n\n");
	  return hgcalc;
	}
	
	if(from_state>energies.size()){
		printf("\n State ' %i ' out of range [max = %lu].\n",from_state,energies.size());	
	  return 0;
  }
	const char *name = Form("GammaSpectrum_State%i%s",from_state,egam>0?Form("_GamGate%.1f",egam):"");
	hgcalc = (TH1D*)list->FindObject(name);
	if(hgcalc){
		if(verbose)printf("\n %s already exists. Using this histogram.\n",hgcalc->GetName());
		return hgcalc;		
	}
	
	if(verbose)printf("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * *");	
	if(verbose)printf("\nDrawing gammas from state %i [%.1f keV] :-\n",from_state,energies.at(from_state-1));	
				
	if(!BuildDecayScheme(from_state))
		return 0;
	
	GetRid(name);
	hgcalc = new TH1D(name,"",4000.0,0,4000.0);
	const char *gnamefull = Form("With %.1f keV Gamma Gate",egam);	
	hgcalc->SetTitle(Form("Gammas Emitted From %.1f keV State [%i] %s; E_{#gamma}^{thry}[keV]; Intensity [# of #gamma's emitted per state];",energies.at(from_state-1),from_state,egam>0?gnamefull:""));
	hgcalc->SetLineColor(kRed);
	hgcalc->SetLineWidth(2);
	list->Add(hgcalc);
	
	std::vector<int> indx = GetCascadeIndex(from_state,egam);
	Int_t ncas = indx.size(), j=0, initial, final;
	Double_t gam, strength, intensity, scale, nbranches;
			
	GetRid("tmp");
	TH1D *htmpi, *htmpf, *htmpg = new TH1D("tmp","",4000.0,0,4000.0);
	
	for(int k=0; k<ncas; k++){

		j	= indx.at(k);
		intensity = 1.0;
		intensity /= htrans->GetBinContent(cascades[j].at(0),cascades[j].at(1));
		if(verbose) printf("\n Drawing Cascade %2i :-",k+1);
		for(int i=0; i<(Int_t)cascades[j].size()-1; i++){

			initial = cascades[j].at(i);
			final = cascades[j].at(i+1);
			if(verbose) printf("\n\tDecay %2i -> %2i : ",initial,final);
			
			// remember that vector index = bin number -1
			gam = energies.at(initial-1)-energies.at(final-1);			
			// strength of transition [expressed as a ratio to 100%]
			strength = htrans->GetBinContent(initial,final);
			// project out all transition strengths from this state so that we can sum over all strengths
			htmpi = htrans->ProjectionY("transtmpi",initial,initial);				
			// sum all transition strengths to normalize and convert this transition to a probability
			strength /= htmpi->Integral(1,initial);
			
			// the transition strength is shared equally between each cascade
			htmpf = htrans->ProjectionY("transtmpf",final,final);		
			nbranches = 0.;
			for(int i=1; i<final; i++)
			  if(htmpf->GetBinContent(i)>0)
			    nbranches+=1.0;
			// we want to make sure that the sum of all branches from a state are 1.    
			if(nbranches>1.0)  
        strength /= nbranches;
	//		if(htmpf->GetEntries())
	//		  strength /= htmpf->GetEntries(); // reduce the strength of a specific cascade
			
			intensity*=strength;
			if(verbose) printf(" Gam. Eng. = %4.1f\tInt. = %5.4f \t Tot. Int. = %5.4f",gam,strength,intensity);
	
			htmpg->SetBinContent(gam,intensity);
		}
		scale = 1.0;
		if(egam>0){
			Int_t bin = floor(egam);
			scale = htmpg->GetBinContent(bin); // rescale by gamma intensity
			if(verbose){
				printf("\n\t---- Gamma Spectrum Contents ----");
				for(int m=htmpg->FindFirstBinAbove(); m<=htmpg->FindLastBinAbove(); m++){
					if(htmpg->GetBinContent(m))
						printf("\n\t  Bin = %i, Content = %.4f",m,htmpg->GetBinContent(m));
				}		
				printf("\n\tGAMMA GATE:  Bin = %i, Content = %.2f\n",bin,scale);		
			}
			htmpg->SetBinContent(bin,0.0); // we don't see gated transition
		}			
		hgcalc->Add(htmpg,scale);
		htmpg->Reset();
	}
	
	if(verbose){
		printf("\n\t---- TOTAL Spectrum Contents ----");
		Double_t tot = 0.0;
		for(int m=hgcalc->FindFirstBinAbove(); m<=hgcalc->FindLastBinAbove(); m++){
				if(hgcalc->GetBinContent(m)){
					printf("\n\t  Bin = %4i, Content = %.4f",m,hgcalc->GetBinContent(m));
			    tot+=hgcalc->GetBinContent(m);
        }
			}	
    printf("\n\t Total Content = %.4f",tot);			
		printf("\n\t --- Complete!\n\n");
	}		
	
	return hgcalc;
}

TH1D *TTransitionTool::CalcGammasRng(Double_t emin, Double_t emax, Double_t egam){

	TH1D *h;
	if(!nndcfile.size()){
	  printf("\n Error :  Load a file first.\n\n");
	  return h;
	}
	const char *name = Form("GammaSpectra%s%s",emin>=0&&emax>emin?Form("_Exc%.1fTo%.1f",emin,emax):"",egam>0?Form("_Gam%.1f",egam):"");
	
	h = (TH1D*)list->FindObject(name);
	
	if(h){
		printf("\n %s already exists. Using this histogram.\n",h->GetName());
		return h;		
	}
	fflush(stdout);
// excitation energy resolution means that states outside of emin-emax range will still be 
// in this window so the theory should do the same in order to reproduce what we see.	
// Include states which are outside of emin-emax range, but only counts their contribution in this range	
// a peak at ±1.5 sigma from the emin-emax range limits would contribute a maximum of ~6.7%
	Double_t emin_ext=emin, emax_ext=emax;	
	Bool_t extend_excrng = true;	
	
	if(emin<0 && emax<0){
		emin_ext = 0.0;
		emax_ext = energies.back()+500.0;
		extend_excrng = false;
	}
	if(extend_excrng){ 
		emin_ext -= NSig*ExcSig;
		emax_ext += NSig*ExcSig;
	}
	
	printf("\n\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");	
	printf("\n\nDrawing gammas from states in%s range %.1f - %.1f keV :-",extend_excrng?Form(" *EXTENDED [+/- %.1f keV]*",NSig*ExcSig):"",emin_ext,emax_ext);	
		
	states = PrintStates(emin_ext,emax_ext);
	if(!states.size()){
		printf("\n\t! Error : No states in this range\n\n");
		return 0;
	}
	
	TH1D *htmp;
	THStack *hstack = new THStack(Form("GammaStacked_%.1fTo%.1f%s",emin,emax,extend_excrng?"_Extended":""),"");
	hstack->SetTitle(Form("Gamma Spectrum From Selected States %s; E_{#gamma}^{thry} [keV]; Intensity [# of #gamma's per state]",egam>0?Form("With %.1f keV Gam Gate",egam):""));
	
	h = new TH1D(name,"",4000.0,0,4000.0);
	h->SetTitle(Form("Gammas Emitted From %5.1f-%5.1f keV State [%i-%i], %s; E_{#gamma}^{thry} [keV]; Intensity [# of #gamma's per state];",
			energies.at(states.at(0)-1),energies.at(states.back()-1),states.at(0),states.back(),egam>0?Form("With %.1f keV Gam Gate",egam):""));
	h->SetLineColor(kRed);
	h->SetLineWidth(2);
	
	
	TF1 *func = new TF1("gaus","gaus",0,8000.0);
	func->SetParameters(1.0/(sqrt(2*3.1415926)*ExcSig),0.0,ExcSig);		
	
	Double_t strength;
	TLegend *leg = new TLegend(0.6,0.4,0.98,0.92);
	leg->SetFillStyle(1001);
	
	TH1D *hgcalc;
	Int_t n=1;
	for(int i=0; i<(Int_t)states.size(); i++){
		hgcalc = CalcGammas(states.at(i),egam);
		if(!hgcalc)
			continue;

		// calculate the contribution from each state is its integral intensity over excitation gate range 
		// take each state to be a normalized gaussian 
		func->SetParameter(1,energies.at(states.at(i)-1));					
		// rescale this spectrum to account for expected contribution in emin-emax range
		strength = func->Integral(emin,emax);
		
		GetRid(Form("%s_Strength",hgcalc->GetName()));
		htmp = (TH1D*) hgcalc->Clone(Form("%s_Strength",hgcalc->GetName()));
		htmp->Scale(strength);
		list->Add(htmp);

		// keep track of the strength 
		htmp->SetTitle(Form("%s @ Strength %4.1f%%",htmp->GetTitle(),strength*100));
		
		// all low [out of range] state gammas are black and all high [o.o.r.] gammas are blue
		if(energies.at(states.at(i)-1)<emin || energies.at(states.at(i)-1)>emax)
			htmp->SetLineColorAlpha(kBlack,strength);		
		else{
		 	if(n==9)
		 		n=19;
			htmp->SetLineColor(++n);
		}	
		
		h->Add(htmp);
		hstack->Add(htmp);
	//	leg->AddEntry(htmp,Form("State %2i, %6.1f keV [%2.f%%]",states.at(i),energies.at(states.at(i)-1),strength*100),"lp");		
		leg->AddEntry(htmp,Form("State %2i [%5.1f keV @ %2.f%%]",states.at(i),energies.at(states.at(i)-1),strength*100.0),"l");
	}

	GetRid(h->GetName(),false);
	list->Add(h);	

	GetRid(hstack->GetName(),false);
	list->Add(hstack);	

	GetRid("CanvasEnergyRange");
	TCanvas *c = new TCanvas("CanvasEnergyRange","Canvas",1100,500);
	c->Divide(2,1);
	c->cd(1);
	GetRid(Form("%s_Clone",hstack->GetName()));	
	THStack *hstack2 = (THStack*)hstack->Clone(Form("%s_Clone",hstack->GetName()));	
	hstack2->Draw();
	leg->Draw();
	
	c->cd(2);
	GetRid(Form("%s_Clone",h->GetName()));		
	TH1D *hh = (TH1D*)h->Clone(Form("%s_Clone",h->GetName()));	
	hh->DrawCopy();
	
	return h;
}
/*
TH1D *TTransitionTool::CalcExc(double egam, bool use_int=true){

  states = PrintStates(egam);
  TF1 *func;
  TH1D *h, *hsum;
	for(int i=0; i<(Int_t)states.size(); i++){
    func = new TF1(Form("state_%i_exc",states.at(i)),"gaus",0,emax);
    func->SetNPX((int)emax);
    h = func->GetHistogram();
    h->Scale(strength/h->Integral()); // normalize hist to contain population strength

    hsum->Add(h,StateStrength(egam)); // recale with gamma gating strength, if applicable
  }
}
*/
TH2F *TTransitionTool::CalcExcGam(Double_t egam, Bool_t use_int){


	TH2F *h;
	if(!nndcfile.size()){
	  printf("\n Error :  Load a file first.\n\n");
	  return h;
	}

	const char *name = Form("ExcGam%s%s",egam>0?Form("_Gam%.1f",egam):"",use_int?"_UsingSetIntensities":"");
  h = (TH2F*)list->FindObject(name);
	if(h){
		printf("\n %s already exists. Using this histogram.\n",h->GetName());
		return h;		
	}

	printf("\n\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");	
	printf("\n\nDrawing excitation energy verus gamma energy :-");	
	
	Double_t emax = energies.back()+500.0;	
	states = PrintStates(0.0,emax);
	if(!states.size()){
		printf("\n\t! Error : No states in this range\n\n");
		return 0;
	}

	h = new TH2F(name,"",2000,0,4000.0,emax/40.0,0.0,emax);
	h->SetTitle(Form("Excitation Versus Gamma %s; E_{#gamma}^{thry} [keV]; E_{exc} [keV];",
	                egam>0?Form("With %.1f keV Gam Gate",egam):""));		

	TAxis *yax = h->GetYaxis();
	Double_t binexwid = yax->GetBinWidth(0), exeng, energy, scale=1.0, val=1.0;	
	Int_t binexlo, binexhi;
	TF1 *func = new TF1("gaus","gaus",0,8000.0);
	func->SetParameters(1.0/(sqrt(2*3.1415926)*ExcSig),0.0,ExcSig);		
		
	TH1D *hgcalc;
	printf("\n\t Progress...\r");
	for(int i=0; i<(Int_t)states.size(); i++){
    printf("\t Progress... %4.2f %%\r",(Double_t)i/(Double_t)states.size()*100.0);     
	  fflush(stdout);
		hgcalc = CalcGammas(states.at(i),egam);
		if(!hgcalc)
			continue;
    energy = energies.at(states.at(i)-1);
    hgcalc = RealisticCalcGammas(hgcalc); // include realistic gamma peak width

    if(use_int && hint)
      scale = hint->GetBinContent(states.at(i));

		// take each state to be a normalized gaussian 
		func->SetParameter(1,energy);	     
    binexlo = yax->FindBin(energy-2.5*ExcSig);  
    binexhi = yax->FindBin(energy+2.5*ExcSig);
    for(int k=binexlo; k<=binexhi; k++){
      exeng = yax->GetBinCenter(k);
      for(int j=hgcalc->FindFirstBinAbove(); j<=hgcalc->FindLastBinAbove(); j++){
        val = hgcalc->GetBinContent(j)*func->Eval(exeng)*binexwid;
        if(val)
          h->Fill(hgcalc->GetBinCenter(j),exeng,val*scale);
      }
    }
	}
  printf("\t Progress... 100.00%%     COMPLETE\n\n");
	list->Add(h);	

	return h;
}

TCanvas *TTransitionTool::DrawDecayMats(Int_t from_state, Bool_t engaxis){
	
	TCanvas *c = new TCanvas("Canvas","Canvas",1100,500);
	
	if(!htrans || !cascades.size()){
		printf("\n\tError :  Read Input File First!\n");
		return c;
	}
	
	c->Divide(2,1);
	c->cd(1);

	TH2F *hint2 = (TH2F*)htrans->Clone(Form("%s_Copy",htrans->GetName()));
	hint2->GetXaxis()->SetNdivisions(from_state);
	hint2->GetXaxis()->SetRangeUser(0,from_state);
	hint2->GetYaxis()->SetNdivisions(from_state);
	hint2->GetYaxis()->SetRangeUser(0,from_state);	
	gPad->SetGrid();
	hint2->Draw("colz");
	
	TF1 *ff = new TF1("pol1","pol1",0,from_state);
	ff->SetParameters(0,1);
	ff->SetLineColor(1);
	ff->SetLineWidth(1);	
	ff->Draw("same");
	
	c->cd(2);
	//CalcGammas(from_state)->Draw();
	THStack *hdecay = GetDecayScheme(from_state,engaxis);
	hdecay->Draw("nostack"); // first draw is necessary
	hdecay->SetTitle(Form("Cascades From State %i : %.1f keV ; Cascade Number; State %s",from_state,energies.at(from_state-1),engaxis?"E_{#gamma}^{thry} [keV]":"Number"));
	hdecay->GetYaxis()->SetTitleOffset(1.5);
	hdecay->GetXaxis()->SetNdivisions(nsequences);
	gPad->SetGridx();	
	if(!engaxis)hdecay->GetYaxis()->SetNdivisions(from_state+1);	
	hdecay->GetYaxis()->SetRangeUser(1,from_state);
	gPad->Modified();
	gPad->Update();
	
	DrawTransitions(from_state,engaxis);
	
	return c;
}

Int_t TTransitionTool::GetStateIndex(Double_t val, bool add_element){ // compares a value to contents of vector to see if it exists as an element
	for(int i=0; i<(Int_t)energies.size(); i++){
		if(fabs(val-energies.at(i))<1.0) // within 1 keV means same state
			return i;
	}
	
	// if element hasn't been found we can add it in
 	if(add_element){
 		energies.push_back(val);
		nstates++; 		
 		if(verbose) printf("\n State %2i :  \r",nstates);
 		return energies.size()-1;
 	}
	return -1;
}

std::vector<int> TTransitionTool::GetCascadeIndex(Int_t from_state, Double_t egam){

	std::vector<int> indx;
	int s1 = 0, s2 = 0;
	if(egam>0){
		TH1D *htmp;
		Double_t binval;
		Int_t bin = round(egam);
		for(int i=1; i<=from_state; i++){
			htmp = hgams->ProjectionY("",i,i);
			/*
			for(int k=htmp->FindFirstBinAbove(); k<=htmp->FindLastBinAbove(); k++){
				binval = htmp->GetBinContent(k);
				if(binval) printf("\n\t state %i @ %.1f keV [bin %i] = %.1f",i,htmp->GetBinCenter(k),k,binval);
			}
			*/
			binval = htmp->GetBinContent((int)egam);
		
			if(binval){
				s2 = i;
				s1 = GetStateIndex(energies.at(s2-1)-egam)+1; // indx + 1 = state number 
				//printf("\n\t State %i @ %.1f keV [bin %i] = %.1f",i,egam,htmp->GetXaxis()->FindBin(egam),binval);				
				if(verbose)printf("\n\t >--{ Found transition : State %i -> %i }--> ",s2,s1);
			}
		}
		if(!s2 || s1<0){
			if(verbose)printf("\n\t No Gamma @ %.1f keV found in cascades from state %i!\n",egam,from_state);
			return indx;
		} 
	}
	
	for(int i=0; i<cascades.size(); i++){
		if(cascades[i].at(0) == from_state){
			if(s1 && s2){
				// search through cascade for s2->s1 transition
				for(int j=0; j<(Int_t)cascades[i].size()-1; j++){
					if(cascades[i].at(j) == s2){ // if s2 is found then s1 MUST be next
					 if(cascades[i].at(j+1) ==s1)
					 		indx.push_back(i);
					 else 
					 	break;
					}
				}					 	
			} else{
				indx.push_back(i);
			}
		}
	}
	
	return indx;		
}

Int_t TTransitionTool::BuildDecayScheme(int from_state){

	hseq = 0;
	nsequences = 0;	
	
	if(!htrans){
		printf("\n\t! Error : Read Input File First !\n\n");
		return 0;
	} else if(from_state>nstates){
		printf("\n\t! Error: Maximum allowed state = % i!\n\n",nstates);
		return 0;
	}
	
	// make a sequence hist for this state and see if there are already cascades built
	if(MakeSequenceHist(from_state)) // if cascade is already built then just return 
		return nsequences;
	
	std::vector<int> cas;  // index to show the current column on each row
	cas.push_back(from_state);
	int m = from_state, i=0;
	
	if(verbose) printf("\n Building Cascades from state %i..",from_state);
	for(i=1; i<=m; i++){
			
		if(i==from_state)
			break;
		
		if(i==m){ // end of loop for row m, so go back to previous row
			i = cas.at(cas.size()-1); // s[end] is previous column [the -1 is immediately incremented to zero..]
			m	= cas.at(cas.size()-2); // s[end-1] is previous row
			cas.pop_back(); // remove last element
			continue;
		}
					
		if(htrans->GetBinContent(m,i)>0){
			cas.push_back(i);
			if(i==1) // ground state has been reached
				AddCascade(cas);
			m = i; // set current row
			i = 0; // restart loop on next row	(0 will be incremented immediately to 1)		
			continue;	
		}	
	}
	
	MakeSequenceHist(from_state);
	if(verbose) printf("\n\t --- Complete!\n");
	
	return nsequences;
}

void TTransitionTool::AddCascade(std::vector<int> states){

	cascades[ncascades] = states;
	ncascades++;
	
	if(verbose){
		printf("\n\t* Added Cascade %3i :  [ ",ncascades);
		for(int i=0; i<(Int_t)states.size(); i++)
			printf( "%i ",states.at(i));
		printf("]");
	}
	return;
}

Int_t TTransitionTool::MakeSequenceHist(Int_t from_state){
	
	std::vector<int> cas = GetCascadeIndex(from_state); 
	nsequences = cas.size();
	if(!nsequences)
		return 0;
	if(verbose) printf("\n\t State %i has %i cascades.\n",from_state,nsequences);	
	
	Int_t maxnstates=0;
	for(int i=0; i<nsequences; i++){
		///printf("\n\t Cascade %i :  ",cas.at(i));
		for(int j=0; j<cascades[cas.at(i)].size(); j++)
//printf(" -> %i",cascades[cas.at(i)].at(j));
		if(cascades[cas.at(i)].size()>maxnstates)
			maxnstates = cascades[cas.at(i)].size();
	}//printf("\n maxnstates = %i\n\n",maxnstates);
	
	const char *name = Form("CascadeMat_State%i",from_state);	
	GetRid(name);
	hseq = new TH2F(name,"",maxnstates,0,maxnstates,nsequences,0,nsequences);
	hseq->SetTitle("Cascade To Ground State; Series of States; Cascade Number");	
	hseq->SetLineWidth(2);
	
	Int_t j=0;
	for(int k=0; k<nsequences; k++){
		j = cas.at(k);
		for(int i=0; i<(Int_t)cascades[j].size(); i++){
		//	printf("\n\t Bin [ %i , %i ] = %i",i+1,k+1,cascades[j].at(i));
			hseq->SetBinContent(i+1,k+1,cascades[j].at(i));
			hseq->SetBinError(i+1,k+1,0.0001);
		}
	}
	
	return nsequences;
}

THStack *TTransitionTool::GetDecayScheme(Int_t from_state, Bool_t engaxis){
	
	MakeSequenceHist(from_state);	
	
	TH1D *hproj;
	THStack *hdecay = new THStack(Form("DecayScheme_State%i",from_state),Form("Cascades From State %i",from_state));
	
	Int_t binmax = hseq->ProjectionX()->FindLastBinAbove();
	for(int i=1; i<=binmax; i++){
		hproj = hseq->ProjectionY(Form("Step%i",i),i,i);
		
		if(engaxis) ConvertToEnergyAxis(hproj);
		
		if(i==1){
			hproj->SetLineWidth(4);		
			hproj->SetLineColor(kBlue);
		}
		else{
			hproj->SetLineWidth(2);
			hproj->SetLineColor(kBlack);
		}

		hdecay->Add(hproj);
	}

	return hdecay;
}

void TTransitionTool::ConvertToEnergyAxis(TH1D *h){

	for(int j=1; j<=h->FindLastBinAbove();j++){
		if(h->GetBinContent(j)==0)
			continue;
		h->SetBinContent(j,energies.at((Int_t)h->GetBinContent(j)-1));
	}
	return;
}

void TTransitionTool::DrawTransitions(int from_state, Bool_t engaxis){
		
	TArrow *arrow;
	TPaveText *pt;
	TH1D *h;
	Double_t dh=0.1, ypos;
	if(engaxis)
		dh=50;	
	
	for(int i=1; i<=nsequences; i++){
		h = hseq->ProjectionX("",i,i);
		if(engaxis) ConvertToEnergyAxis(h);		
			
		for(int j=1; j<=h->FindLastBinAbove(); j++){
			if(!engaxis && h->GetBinContent(j+1)==0)
				continue;
			arrow = new TArrow(i-0.5,h->GetBinContent(j)-dh,i-0.5,h->GetBinContent(j+1)+dh,0.01,"|>");
			arrow->SetAngle(40);
			arrow->SetFillColor(kRed);
			arrow->SetFillStyle(3001);
			arrow->SetLineColorAlpha(kRed,0.25);
			arrow->Draw();
			if(engaxis){
				ypos = 0.5*(h->GetBinContent(j)+h->GetBinContent(j+1));
				pt = new TPaveText(i-0.95,ypos-dh,i-0.05,ypos+dh,"nb");
				TText *text = pt->AddText(Form("%6.1f keV",h->GetBinContent(j)-h->GetBinContent(j+1)));
				pt->SetFillColor(0);
				pt->SetTextColorAlpha(kBlack,0.5);
				text->SetTextAngle(25.);
		//		pt->SetTextSize(2);
				pt->Draw();
			}
		}
	}
	
	return;
}

TH1D *TTransitionTool::RealisticCalcGammas(int from_state, Double_t egam, Int_t color){

  if(!energies.size()){
    verbose = false;
    LoadLevelsGammas();
  }  
  
   if(!fTigEff){
    printf("\n\n Warning : Efficiency has been set to default values.\n");
    SetEfficiencyCurve();
  } 

  if(!fTigSigma){
    printf("\n\n Warning : Energy resolution has been set to default values.\n");
    SetResolutionCurve();
  }  
  
  TH1D *hgcalc = CalcGammas(from_state,egam);
  if(!hgcalc)
    return hgcalc;
  
	// converts delta function spectrum into sigma(E) and efficiency corrected spectrum
	GetRid(Form("%s_Realistic",hgcalc->GetName()));
	TH1D *htmp = new TH1D(Form("%s_Realistic",hgcalc->GetName()),hgcalc->GetTitle(),hgcalc->GetNbinsX(),0,4000);		
	htmp->GetXaxis()->SetTitle(hgcalc->GetXaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(hgcalc->GetYaxis()->GetTitle());
	htmp->SetLineColor(color);	
	
	Double_t x, val, sig, eff;
	TF1 *func = new TF1("gauspeak","gaus",0,4000);
	func->SetNpx(4000);
		
	// integrate the counts under the data peaks 
	for(int i=1; i<hgcalc->GetNbinsX(); i++){
		
		val = hgcalc->GetBinContent(i);
		if(!val)
			continue;
		
		// define a gaussian centered at bin center with sigma = sig and integral = hcalc->GetBinContent(i)		
		x = hgcalc->GetBinCenter(i);
		sig = fTigSigma->Eval(x);
		eff = Efficiency(x);		
		
		func->SetParameters(eff*val/(sqrt(2*3.14159)*sig),x,sig);
		// now get this histogram and use it as a filter over the data
		htmp->Add((TH1D*)func->GetHistogram());
	}
	
	return htmp;
}

TH1D *TTransitionTool::RealisticCalcGammas(TH1D *hgcalc, Bool_t abseff){

  if(!fTigSigma){
    printf("\n\n Error : Energy resolution must be set first.\n\n");
    return 0;  
  }
  
	// converts delta function spectrum into sigma(E) and efficiency corrected spectrum
	GetRid(Form("%s_Realistic",hgcalc->GetName()));
	TH1D *htmp = new TH1D(Form("%s_Realistic",hgcalc->GetName()),hgcalc->GetTitle(),hgcalc->GetNbinsX(),0,4000);		
	htmp->GetXaxis()->SetTitle(hgcalc->GetXaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(hgcalc->GetYaxis()->GetTitle());	
	
	Double_t x, val, sig, eff;
	TF1 *func = new TF1("gauspeak","gaus",0,4000);
	func->SetNpx(4000);
		
	// integrate the counts under the data peaks 
	for(int i=1; i<hgcalc->GetNbinsX(); i++){
		
		val = hgcalc->GetBinContent(i);
		if(!val)
			continue;
		
		// define a gaussian centered at bin center with sigma = sig and integral = hcalc->GetBinContent(i)		
		x = hgcalc->GetBinCenter(i);
		sig = fTigSigma->Eval(x);
		eff = Efficiency(x);		
		
		func->SetParameters(eff*val/(sqrt(2*3.14159)*sig),x,sig);
		// now get this histogram and use it as a filter over the data
		htmp->Add((TH1D*)func->GetHistogram());
	}
	
	return htmp;
}

void TTransitionTool::GetRid(const char *name, Bool_t delete_all){
	TObject *obj;

	obj = list->FindObject(name);
	if(obj) list->Remove(obj);

	if(delete_all){
		obj = gROOT->FindObjectAny(name);
		if(!obj) return;
		/*if(strcmp(obj->ClassName(),"TCanvas")==0){
			TCanvas *c = (TCanvas*)obj;
			c->Close();
		}
		else */if(obj) obj->Delete();
	}
	
}

