#include "Domain.h"

//#include <CompactNSearch> //NEW WAY

using namespace std;
namespace SPH {
  
#ifdef NONLOCK_TEST
void Domain::CheckParticlePairs(const int &i){
  cout << "Particle i: "<<i<<endl;
  cout << "Nb Count "<<Particles[i]->Nb<<endl;
  cout << "ipairs: "<<ipair_SM[i]<<endl;
  cout << "jpairs: "<<jpair_SM[i]<<endl;
  
  cout << "First pair pos"<<ipl_SM[i]<<endl;
  cout << "Nb list"<<endl;
  for (int n=0;n< MAX_NB_PER_PART;n++){
    cout << Anei[i][n] <<", ";
  } 
  cout  << endl<<"Ref list "<<endl;
  for (int n=0;n< MAX_NB_PER_PART;n++)
    cout << Aref[i][n]<<", ";
  cout << endl;
  cout << "Pair detail "<<endl;
  cout << "i<j pairs"<<endl;
  for (int n=0;n< ipair_SM[i];n++)
    cout << pair_test[Aref[i][n]].first<<", "<<pair_test[Aref[i][n]].second<<endl;
  cout <<"j>i pairs"<<endl;
  for (int n = 0;n< jpair_SM[i];n++)
    cout <<pair_test[Aref[i][MAX_NB_PER_PART-1-n]].first<<", "<<pair_test[Aref[i][MAX_NB_PER_PART-1-n]].second<<endl;
  
  cout << "Done. "<<endl;
}
#endif

inline void Domain::CalcPairPosList(){                             //Calculate position list for every particle
  
  #ifdef NONLOCK_TEST
    pair_test.clear();
  #endif
  first_pair_perproc[0] = 0;
  pair_count = 0;
  for (int p=0;p<Nproc;p++)first_pair_perproc[p] = 0;
  int count = 0;
  for (int p=0;p<Nproc;p++) {
    if (p > 0 ){
      for (int j=0;j<p;j++)
        first_pair_perproc[p]+=SMPairs[j].Size();
    }
    pair_count += SMPairs[p].Size();
  }
  // cout << "Pairs per proc"<<endl;
  // for (int p=0;p<Nproc;p++)
    // cout << SMPairs[p].Size()<<endl;
   // cout << "First index per proc "<<endl;
  // for (int p=0;p<Nproc;p++) cout << first_pair_perproc[p]<<endl;

  pair_force.resize(pair_count);
  temp_force.resize(pair_count);
  pair_StrainRate.resize(pair_count);
  pair_RotRate.resize(pair_count);
  pair_densinc.resize(pair_count);
  //cout << "Pair Count: " << pair_count << endl;

  #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i = 0;i<Particles.Size();i++){
      ipair_SM[i]=0;jpair_SM[i]=0;
      //ipl_SM[i]=0;  
      for (int n=0;n<MAX_NB_PER_PART;n++){
        Anei[i][n]=0;
        Aref[i][n]=0;
      }
    }
 
  //cout << "calc tables "<<endl;
  //TODO: parallelize?
  //#pragma omp parallel for schedule (static) num_threads(Nproc)
  for (int k=0;k<Nproc;k++){
    for (int pp=0;pp<SMPairs[k].Size();pp++){
      //int p = first_pair_perproc[k] + pp;

      int i = std::min(SMPairs[k][pp].first,SMPairs[k][pp].second);
      int j = std::max(SMPairs[k][pp].first,SMPairs[k][pp].second);
      #ifdef NONLOCK_TEST
        pair_test.push_back(std::make_pair(i,j)); //ONLY FOR TEST
      #endif
      // if (i==2000||j==2000){ 
        // cout << "Pair i j "<< i << ", "<<j<<", "<<first_pair_perproc[k]+pp<<endl;
        // cout << "k pp "<<k <<", "<<pp<<endl;
      // }
      Anei[i][ipair_SM[i]] = j; //Only stores j>i
      Anei[j][MAX_NB_PER_PART - 1 - jpair_SM[j]] = i; //Only stores j>i
      Aref[i][ipair_SM[i]] = first_pair_perproc[k]+pp;
      Aref[j][MAX_NB_PER_PART - 1 - jpair_SM[j]] = first_pair_perproc[k]+pp;
      ipair_SM[i]++;            //ngji in 
      jpair_SM[j]++;            //njli, pairs in which j has particles with index smaller than it
    }//Pairs
  }//Proc
  
  // ipl_SM[0] = 0;
  // #pragma omp parallel for schedule (static) num_threads(Nproc)
  // for (int i=0; i<Particles.Size();i++){
    // for (int j=0; j<i;j++)
      // ipl_SM[i] += ipair_SM[j];//Nisimura 2015 Eqn 6
  // }
  
  int max_nb=0;
  for (int i=0; i<Particles.Size();i++){
    if (ipair_SM[i]+jpair_SM[i]>max_nb)
      max_nb = ipair_SM[i]+jpair_SM[i];
  }
  //cout << "Max nb"<<max_nb<<endl;
}

inline void Domain::CalcRefTable(){

  #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i = 0;i<Particles.Size();i++){
      for (int n=0;n < ipair_SM[i];n++){ //neighbour i < j count 
        int pair = ipl_SM[i] + n; //According with ipl_SM[i]
        // if (pair >= pair_count) {
          // cout << "ERROR pair " <<pair << "> pair_count "<< pair_count << endl;
          // cout << "i n "<<i<<", "<<n<<endl;
          // throw new Fatal("ERROR");
        // }
        Aref [i][n] = pair;
        int j = Anei[i][n];
        for (int k=0; k < jpair_SM[j];k++){
          if ( Anei[j][MAX_NB_PER_PART-1-k] == i)
            //Aref [j][ipair_SM[j]+k]= pair;
            Aref [j][MAX_NB_PER_PART-1-k]= pair;
        }
      }//nb
    }//particle
    
  //cout << "aref "<<Aref[0][0]<<endl;
}

};