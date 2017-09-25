#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include <ctype.h>


struct levelscheme{

  int ntransitions;
  std::vector<double> state_initial; 
  std::vector<double> gamma_energy;
  std::vector<double> gamma_strength;
  std::vector<double> state_final;
  
  int AddTransition(double state_i, double gam_e, double gam_i, double state_f){
    if(fabs(gam_e+state_f-state_i)>5.0){
      printf("\n Error {levelscheme} : values Ei=%.1f Eg=%.1f Ig=%.1f Ef=%.1f violate energy conservation. Ignored. \n",state_i,gam_e,gam_i,state_f);
      return -1;
    }
    state_initial.push_back(state_i);
    gamma_energy.push_back(gam_e);
    gamma_strength.push_back(gam_i);
    state_final.push_back(state_f);
    ntransitions++;
    
    return (int)state_initial.size()-1; // returns index of new transition
  }
  
  int Print(int indx_beg=0, int indx_end=-1 ){
       
    if(indx_beg>ntransitions-1 || indx_beg<0)
      return 0;
    if(indx_end>ntransitions-1 || indx_end<indx_beg)
      indx_end=ntransitions-1;

    printf("\n\n\t\t _____LEVEL_SCHEME_[%i-%i]______\n",indx_beg,indx_end);  
  printf("\nentry# : State [keV] : Gamma energy [keV] : Gamma strength [%%] : Final State [keV]");
    for(int i=indx_beg; i<=indx_end; i++){
      printf("\n%3i %13.1f %20.1f %20.1f %18.1f",i,state_initial.at(i),gamma_energy.at(i),gamma_strength.at(i),state_final.at(i));
      if(i<ntransitions-1 && state_initial.at(i)!=state_initial.at(i+1)) printf("\n");
    }
    
    return indx_end-indx_beg;
  }
  
  void Write(std::string fname){
    std::ofstream ofile(fname.c_str());
    int n=0;
    
    for(int i=0; i<(int)state_initial.size(); i++){
       ofile << state_initial.at(i) << "\t\t" << gamma_energy.at(i) << "\t\t" << gamma_strength.at(i) << "\t\t" << state_final.at(i) << "\n";
      if(i<(int)state_initial.size()-1 && state_initial.at(i)!=state_initial.at(i+1)){
        n++;
        ofile << "\n";
      }
    }
    ofile.close();
    printf("\n\n Wrote %i states and %lu transitions to '%s'.\n\n",n,state_final.size(),fname.c_str());
  }
};

levelscheme LVLZ;

static const unsigned long npos = std::string::npos;	

std::vector<std::string> GetFullLines(std::string fname);
std::vector<std::string> GetFields(std::string);
std::vector<double> GetNumericalValues(std::string);
int CheckNumericalValues(double , std::vector<double> &, std::vector<double> &, std::vector<double> &);


void ParseCSV(std::string fname = "test.csv", int state_i_col=0, int gam_eng_col=1, int gam_int_col=2, int state_f_col=3){

  std::vector<std::string> fulllines = GetFullLines(fname);
  std::vector<std::string> strvec;
  std::vector<double> numvals;
  std::string line, numstr;
  
  // state info
  double state_eng;
  int ntrans, indx;
  std::vector<double> state_initial, gam_eng, gam_int, state_final;  
  
  for(int i=0; i<(int)fulllines.size(); i++){
    line.assign(fulllines.at(i));  
    printf("\n\n------------------------------------------------------------------------------------------------------\n");    
    printf("\nPARSED LINE: %s\n",line.c_str());
    strvec = GetFields(line);      
    
    state_initial.clear();
    gam_eng.clear();
    gam_int.clear();
    state_final.clear();
    
    // look over fields and pull out numerical values within each field
    for(int j=0; j<(int)strvec.size(); j++){
  //    printf("\n%i.\tstrvec.at(%i) = '%s'",j,j,strvec.at(j).c_str());
    
      numvals = GetNumericalValues(strvec.at(j));
      if(!numvals.size())
        continue;

      // now assign the numerical values to their respective vectors
      if(j==state_i_col){ // first number in state energy column
        state_eng = numvals.at(0);      
      }if(j==gam_eng_col){ // store numerical values as gamma energies 
        gam_eng = numvals;
      }if(j==gam_int_col){ // store numerical values as gamma strengths 
        gam_int = numvals;
      }if(j==state_f_col){ // store numerical values as final state energies
        state_final = numvals;
      }
    }
    
    ntrans = CheckNumericalValues(state_eng,gam_eng,gam_int,state_final);
    state_initial.assign(ntrans,state_eng); 
    
    std::vector<int> indx;
    for(int i=0; i<ntrans; i++)
      indx.push_back(LVLZ.AddTransition(state_initial.at(i),gam_eng.at(i),gam_int.at(i),state_final.at(i)));
    
    if(indx.size())
      LVLZ.Print(indx.front(),indx.back());
  }    

  printf("\n\n\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");    
  printf("\n* * * * * * * * * * * * * * * * * COMPLETE * * * * * * * * * * * * * * * * * * * \n");    

  LVLZ.Print();
  LVLZ.Write("LVLZ.txt");
}

int CheckNumericalValues(double state_eng, std::vector<double> &gam_eng, std::vector<double> &gam_int, std::vector<double> &state_final){

  int ntrans = state_final.size();
  // there are 
  int ngameng, ngamint;
  ngameng = gam_eng.size();
  ngamint = gam_int.size();


  if(ntrans>1 && ngameng!=ngamint)
    printf("\n-> Warning: There are a different number of gamma energy entries [%i] than gamma intensity entries [%i]. The 0th, 2nd, 4th... elements will be used but this may not be correct. Check the csv file and add additional uncertainties to make these fields equal.",ngameng,ngamint);
  if(ngameng!=ntrans || ngamint!=ntrans)
    printf("\n-> Warning: The number of gammas does not match the number of final states [initial=%i, gam_eng=%lu, gam_int=%lu, final=%lu] ...",ntrans,gam_eng.size(),gam_int.size(),state_final.size());

  bool match;
  std::vector<int> indx;
  double etol = 5.0; // allow up to etol energy discrepancy between initial-gam_eng-final
  
  for(int i=0; i<(int)gam_eng.size(); i++){ // loop over all gamma energies and look for energy matches
    match = false;
    for(int j=0; j<ntrans; j++){
    
      double etot = gam_eng.at(i)+state_final.at(j);
      if(fabs(state_eng - etot)<5.0){
     //   printf("\nENERGY MATCH:  gamma %i [%.1f] + final state %i [%.1f] = %.1f keV [inital state = %.1f keV]",
     //     i,gam_eng.at(i),j,final_state.at(j),etot,state_eng);
        match = true;
        indx.push_back(i);
        break;
      }
    }
    if(!match){ // Couldn't find a final state for gamma i so remove it
  //    printf("\n Couldn't find a final state for gamma %i [%.1f keV], so deleting...",i,gam_eng.at(i));
      gam_eng.erase(gam_eng.begin()+i);
      i--;
    }
  }
  
  if((int)gam_int.size()>ntrans){
    for(int i=(int)gam_int.size()-1; i>0; i--)
      if(i%2) gam_int.erase(gam_int.begin()+i);
  }

  if(ngameng != (int)gam_eng.size() || ngamint != (int)gam_int.size())
    printf(" ntrans=%i : ngameng now %lu, and ngamint now %lu\n",ntrans,gam_eng.size(),gam_int.size()); 

  return ntrans;
} 

std::vector<double> GetNumericalValues(std::string field){

  std::vector<std::string> strvec;
  std::vector<double> numvals;
  std::string line, numstr;
  
  unsigned long pos_e=field.find("e+"), pos_E=field.find("E+"), pos_exp=npos;
  int expval;
  if(pos_e!=npos){
    pos_exp = pos_e;
    expval = (int)field.at(pos_e+2)-48;
  }if(pos_E!=npos){
    pos_exp = pos_E;
    expval = (int)field.at(pos_E+2)-48;
  }
  
  if(pos_exp!=npos && expval>=0){
    printf("\n-> Warning: Scientific notation has been detected [ %s == 10^%i ] in this field '%s'",field.substr(pos_exp,3).c_str(),expval,field.c_str());
    field.replace(pos_exp,3,expval,'0');
    printf(" ... Automatically changed to '%s'",field.c_str());
  }
  // step through field string and extract numerical values
  for(int i=0; i<(int)field.length(); i++){
   //   printf("\n\"%c\"",strval.at(i));

    if(isdigit(field.at(i)) || field.at(i)=='.'){
      numstr = numstr + field.at(i);
    //  printf(" is a number or decimal point - \"%s\"",numstr.c_str());              
    } else {
  //    printf(" is not useful");
      if(numstr.length()){
        numvals.push_back(atof(numstr.data()));
  //      printf("\n I found a number = %s [ as a number %f ]\n",numstr.c_str(),vals.back());         
      }
      numstr.clear();
      continue;   
    } 
  }
  if(numstr.length()){
    numvals.push_back(atof(numstr.data()));    
//    printf("\n I found a number = %s [ as a number %f ]\n",numstr.c_str(),atof(numstr.data()));
  }    
  //  printf("\n");

  return numvals;
}

std::vector<std::string> GetFields(std::string line){

  std::vector<std::string> strvec;

  int n=0;
  std::string strval;
  unsigned long pos_curr=0, pos_next=0, posq_beg=0, posq_end=0;  
  
  pos_next = line.find(',',pos_curr);   // line begins with a comma, skip
  while(pos_next!=npos){
  //  printf("\n\n\tQUOTES AND COMMAS [pos_curr=%i posq_beg=%i pos1_end=%i pos_next=%i]",(int)pos_curr,(int)posq_beg,(int)posq_end,(int)pos_next);
  
    //  extract individual fields using comma separation
    pos_curr = pos_next;                // next time search from position of previous comma 
    posq_beg = line.find('"',pos_curr+1);   // search for a quotation mark from prev comma
    posq_end = line.find('"',posq_beg+1); // search for a closing quotation mark after first q.m.
    pos_next = line.find(',',pos_curr+1);   // search for next comma from position of previous comma + 1
    
    // if a quotation mark is found before a comma, make a string with entire quotation contents
    if(posq_beg<pos_next){
      strval.assign(line.substr(posq_beg+1,posq_end-posq_beg-1).data());
      pos_next = posq_end+1;
    } else 
      strval.assign(line.substr(pos_curr,pos_next-pos_curr).data());    

    strvec.push_back(strval);
   // printf("\nField %i -> \"%s\"",n,strval.c_str());    
  }

  return strvec;

}

std::vector<std::string> GetFullLines(std::string fname){

  
  std::string nextline, thisline, prevline;
  std::ifstream infile(fname.c_str());
  std::vector<std::string> fulllines;
  
  std::getline(infile,thisline);
  while(std::getline(infile,nextline)){
    
    if(nextline.find(",")!=0){ // this line is a continuation
      thisline = thisline + nextline;
    //  printf("\n %s \n",thisline.c_str());
      continue;
    } 
    
 //   printf("\n-> %s\n",thisline.c_str());
    fulllines.push_back(thisline);
    thisline = nextline;
    
  }

  return fulllines;
}


