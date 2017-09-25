#ifndef TTRANSITIONTOOL_H
#define TTRANSITIONTOOL_H

#include<TObject.h>
#include<Rtypes.h>
#include<TList.h>

#include<TH3S.h>
#include<TH2F.h>
#include<TH1D.h>
#include<TF1.h>
#include<TGraphErrors.h>
#include<THStack.h>
#include<TCanvas.h>

#include <string>
#include <map>

#ifndef TIGDIR
#define TIGDIR "$PROGDIR/TTransitionTool"
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// _____________ TTransitionTool : a gamma-particle analysis class ____________
// 
//  TTransitionTool is used to analyze data by applying [in almost any combination] :-
//   1. Gamma-gates
//   2. Particle-gates
//   3. Angular-gates
// 
//
// _____ Optional additions _____
//  Same as above, but for specific regions of SHARC [UQQQ, UBOX and DBOX].
//  Names and conventions will be taken care of by MakeExcGamThetaMats.C
//
//
// _____ Known Bugs _____
// * Background subtraction must include a region both above and below the peak 
// Intensities in calculated spectra may be incorrect.  
// * If a cascade begins with only 1 transition which later branches out the first transition
// will be ncascades too strong because it will use strength==1 every time.
//
//    Made by Steffen Cruz [2016]
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


class TTransitionTool 	{
	
	public:

		TTransitionTool(void);
		~TTransitionTool();
	
		static void Print(Option_t * = "");		
	  static std::vector<int> PrintStates(Double_t emin=0.0, Double_t emax=10000.0);
	  static std::vector<int> PrintCascades(Int_t from_state, Double_t egam=0.0, Bool_t printeng=false);
			
		static void ClearVars(Option_t * = "");		

	  static void SetVerbose(bool verb) { verbose = verb;}
		
		static Bool_t Init(std::string fname="Sr96_LevelsGammas.txt", Int_t nmax=100, Bool_t verb=true);
    static Bool_t LoadLevelsGammas(std::string fname="Sr96_LevelsGammas.txt", int nmax=100);

	  static TCanvas *DrawDecayMats(Int_t from_state, Bool_t engaxis = true);

//	  static void FixIntensities(Double_t exmin=0.0, Double_t exmax=10000.0, Double_t egam=0.0, Int_t level=1, Double_t strength=50.0);
	  static Bool_t SetIntensity(Int_t state, Double_t strength);
	  static Bool_t ReadIntensities(const char *fname = "SavedIntensities.txt");
	  static void WriteIntensities(const char *fname = "SavedIntensities.txt");

	  static TH1D *CalcGammas(Int_t from_state, Double_t egam=0.0);
		static TH1D *CalcGammasRng(Double_t exmin=-1.0, Double_t exmax=-1.0, Double_t egam=0.0);
    static TH2F *CalcExcGam(Double_t egam=-1.0, Bool_t use_int=false);	

	  static TH1D *RealisticCalcGammas(int from_state, Double_t egam=0.0, Int_t color=1);
	  static TH1D *RealisticCalcGammas(TH1D *hgam, Bool_t abseff = false);
				
		static TH1D *DrawCalcIntensitites(){ return hint; }				
		static Double_t BranchingRatio(Double_t state_eng, Double_t egam);			
														
    static TF1 *SetResolutionCurve(double sig=2.0, double eng=1000.0);    
    static Double_t Resolution(Double_t eng, Bool_t FWHM=false);

    static TF1 *SetEfficiencyCurve(double eff=0.1, double eng=1000.0);		
		static Double_t Efficiency(Double_t eng);		
		
	private: 
     
		static Int_t GetStateIndex(Double_t val, bool add_element=false);
		static std::vector<int> GetCascadeIndex(Int_t from_state, Double_t egam=0.0);
	
	  static Int_t BuildDecayScheme(Int_t from_state);	  	
	  static void AddCascade(std::vector<int> states);
		static Int_t MakeSequenceHist(Int_t from_state);
	  static THStack *GetDecayScheme(Int_t from_state, Bool_t engaxis = true);

	  static void ConvertToEnergyAxis(TH1D *h);
	  static void DrawTransitions(Int_t from_state, Bool_t engaxis);	  
	
	  static void GetRid(const char *name, Bool_t delete_all = true); // removes object

		static std::string nndcfile;
		
		static TF1 *fTigSigma, *fTigEff, *fTigEffLow, *fTigEffUp;
		static TGraphErrors *gTigEff, *gTigRes;			
		
		static Bool_t verbose;
		static TList *list;

		static Int_t nlines; 
		static Int_t nsequences;
		static Int_t ncascades;
		static Int_t nstates;		
		static Double_t ExcSig;
		static Double_t GamSig;
		static Double_t NSig;			
		static std::vector<double> energies;
		static std::vector<int> states;
		static std::map<int,std::vector<int> > cascades;
		
		static TH2F *htrans, *hgams, *hseq;
		static TH1D *hint;
		        
	ClassDef(TTransitionTool,0)
};


#endif
