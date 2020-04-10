/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "WaterpKwBase.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/Matrix.h"
#include "ActionRegister.h"
#include <string>
#include <cmath>
#include <iostream>
#include "time.h"
#include <iomanip>


using namespace std;

namespace PLMD{
namespace colvar{

void WaterpKwBase::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","GROUPA","First Acid/Base group");
  keys.add("atoms","GROUPB","List of Hydrogen Atoms");
  keys.add("compulsory","LAMBDA","1","The lambda parameter of the sum_exp function; 0 implies 1");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.addOutputComponent("sd","default","Acid-base distance order parameter");
  keys.addOutputComponent("tc","default","Total charge order parameter");
}

WaterpKwBase::WaterpKwBase(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  list_a = ga_lista;
  list_b = gb_lista;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parse("D_0",d0);
  parse("LAMBDA",lambda);

// pair stuff
  bool dopair=false;
  parseFlag("PAIR",dopair);

// neighbor list stuff
  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parse("NL_CUTOFF",nl_cut);
   if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
   parse("NL_STRIDE",nl_st);
   if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }
  
  //TODO: add neighbor list for gc_lista
  addComponentWithDerivatives("sd"); componentIsNotPeriodic("sd");
  addComponentWithDerivatives("tc"); componentIsNotPeriodic("tc");
  if(gb_lista.size()>0){
    if(doneigh)  nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc(),nl_cut,nl_st);
    else         nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc());
  } else {
    if(doneigh)  nl= new NeighborList(ga_lista,pbc,getPbc(),nl_cut,nl_st);
    else         nl= new NeighborList(ga_lista,pbc,getPbc());
  }

//  requestAtoms(nl->getFullAtomList());
 
  std::vector<AtomNumber> atoms;
  atoms.insert(atoms.end(),list_a.begin(),list_a.end());
  atoms.insert(atoms.end(),list_b.begin(),list_b.end());
  requestAtoms(atoms);

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh){
   log.printf("  using neighbor lists with\n");
   log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }

}

WaterpKwBase::~WaterpKwBase(){
  //cout << "##########Crazy1###########";
  delete nl;
  //cout << "##########Crazy2###########";
}

void WaterpKwBase::prepare(){
  if(nl->getStride()>0){
    if(firsttime || (getStep()%nl->getStride()==0)){
//      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    }else{
//      requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void WaterpKwBase::calculate()
{

 //MCA: clock stuff
 //clock_t t0,tf;
 //double total_time;

//t0 = clock();

 //cout << "Tag -2 ###########################";
 //The 2 scalar CVs
 double IonDistance=0.0;
 double TotalCharge=0.0;
 //Length of acids groups, with and without the H atoms
 unsigned len_acids = list_a.size();
 unsigned len_acids_hyd = len_acids + list_b.size();

 //Parameter controlling the smotthness of the |x| function
 double alpha=0.0001;
 double rcut=3.0; //MCA: neighbor list

 vector<double> sum_exp(len_acids_hyd);
 fill(sum_exp.begin(),sum_exp.end(),0.);

 Tensor virial;
 Tensor virial_dist; //MCA: IonDistance CV
 std::vector<Vector> deriv_dist(len_acids_hyd);
 std::vector<Vector> deriv_tc(len_acids_hyd);
 Vector zeros;
 zeros.zero();
 fill(deriv_dist.begin(), deriv_dist.end(), zeros);
 fill(deriv_tc.begin(), deriv_tc.end(), zeros);

 if(nl->getStride()>0 && invalidateList){
   nl->update(getPositions());
 }

unsigned nt=OpenMP::getNumThreads();

const unsigned nn=nl->size();

//if(nt*stride*10>nn) nt=nn/stride/10;
if(nt==0)nt=1;

//cout << "Tag -1 " << endl;

 std::vector<Vector> omp_deriv_dist(len_acids_hyd);
 std::vector<Vector> omp_deriv_tc(len_acids_hyd);
 fill(omp_deriv_dist.begin(), omp_deriv_dist.end(), zeros);
 fill(omp_deriv_tc.begin(), omp_deriv_tc.end(), zeros);
 Tensor omp_virial;

 Matrix<double> c(len_acids_hyd,len_acids_hyd);
 //Matrix<double> dist(len_acids_hyd,len_acids_hyd);
 //std::vector<std::vector<Vector>> dist;
 vector<vector<Vector>> dist(len_acids, vector<Vector>(len_acids_hyd));
 Matrix<double> distmod(len_acids,len_acids_hyd);
 vector<double> coord(len_acids);
 vector<double> charge(len_acids);
 fill(coord.begin(),coord.end(),0.);

 std::vector<double> d;
 d.insert(d.end(),list_a.size(),d0);

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "00 Initial allocation: "  << fixed << setprecision(5) << total_time << endl;
//
//t0 = clock();

#pragma omp parallel for
for(unsigned int i=0;i<len_acids;i++) {   
  for(unsigned int j=i+1;j<len_acids;j++) {   
     if(pbc){
        dist[i][j]=pbcDistance(getPosition(i),getPosition(j));
     } else {
        dist[i][j]=delta(getPosition(i),getPosition(j));
     }
     dist[j][i] = -dist[i][j];
     distmod[i][j] = dist[i][j].modulo();
     distmod[j][i] = distmod[i][j];
  }
  for(unsigned int j=len_acids;j<len_acids_hyd;j++) {   
     if(pbc){
        dist[i][j]=pbcDistance(getPosition(i),getPosition(j));
     } else {
        dist[i][j]=delta(getPosition(i),getPosition(j));
     }
     distmod[i][j] = dist[i][j].modulo();
  }
  distmod[i][i]=0.;
}

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "01 Dist compute: "  << fixed << setprecision(5) << total_time << endl;
//
//t0 = clock();

for(unsigned int j=len_acids;j<len_acids_hyd;j++) {   
   double sum_tmp = 0.0;
   #pragma omp parallel for reduction(+:sum_tmp)
   for(unsigned int i=0;i<len_acids;i++) {   
      sum_tmp += exp(lambda * distmod[i][j]);
   }
   sum_exp[j] = sum_tmp;
}

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "02 sum_exp compute: "  << fixed << setprecision(5) << total_time << endl;
//
//t0 = clock();

//MCA: Ion distance CV
 for(unsigned int i=0;i<len_acids;i++) {   
   double sum_tmp = 0.0;
   #pragma omp parallel for reduction(+:sum_tmp)
   for(unsigned int j=len_acids;j<len_acids_hyd;j++) {   

      c[i][j] = exp( lambda * distmod[i][j] ) / sum_exp[j];
      //coord[i] += c[i][j];
      sum_tmp += c[i][j];

   }
   coord[i] = sum_tmp;
   charge[i] = coord[i] - d[i];
}

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "03 c_ij compute: "  << fixed << setprecision(5) << total_time << endl;
//
//t0 = clock();

 std::vector<vector<Vector>> ompdfunc_delta(len_acids, vector<Vector>(len_acids_hyd));
 std::vector<vector<Vector>> dfunc_delta(len_acids, vector<Vector>(len_acids_hyd));

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "04 dfunc_delta allocation: "  << fixed << setprecision(5) << total_time << endl;
//
////cout << "Tag 1 " << dist[0][0].modulo() << endl;
//
//t0 = clock();

#pragma omp parallel for 
for(unsigned i=0;i<len_acids;i++) {
  for(unsigned j=0;j<len_acids_hyd;j++) {
    dfunc_delta[i][j]=zeros;
    ompdfunc_delta[i][j]=zeros;
  }
}

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "05 dfun_delta compute: "  << fixed << setprecision(5) << total_time << endl;

//cout << "Tag 2 " << endl;
double dfunc_coord;

//t0 = clock();

#pragma omp parallel for private(dfunc_coord)
 for(unsigned int i=0;i<len_acids;i++) {
   for(unsigned int j=len_acids;j<len_acids_hyd;j++){
     if(distmod[i][j]<rcut) {
        for(unsigned int k=0;k<len_acids;k++) {
          dfunc_coord = -lambda *  c[k][j] * c[i][j];
          ompdfunc_delta[i][k] += dfunc_coord * dist[k][j]/distmod[k][j];
          ompdfunc_delta[i][j] -= dfunc_coord * dist[k][j]/distmod[k][j];
        }       
        dfunc_coord = lambda *  c[i][j];
        ompdfunc_delta[i][i] += dfunc_coord * dist[i][j]/distmod[i][j];
        ompdfunc_delta[i][j] -= dfunc_coord * dist[i][j]/distmod[i][j];
     }
   }
 }

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "06 ompdfun_delta compute: "  << fixed << setprecision(5) << total_time << endl;

//cout << "Tag 3 " << endl;
//

//delete[] dfunc_coord;

//t0 = clock();

#pragma omp critical
for(unsigned i=0;i<len_acids;i++) {
  for(unsigned j=0;j<len_acids_hyd;j++) {
    dfunc_delta[i][j]=ompdfunc_delta[i][j];
  }
}

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "07 dfunc_delta compute: "  << fixed << setprecision(5) << total_time << endl;

//delete[] ompdfunc_delta;

 std::vector<double> dfunc_theta(len_acids);

//cout << "Tag 4 " << endl;

//t0 = clock();

#pragma omp parallel for reduction(+:TotalCharge)
 for(unsigned int i=0;i<len_acids;i++) {
     TotalCharge    += sqrt(pow(charge[i],2)+alpha);
     dfunc_theta[i]  = charge[i]/sqrt(charge[i]*charge[i]+alpha);
 }

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "08 Total charge compute: "  << fixed << setprecision(5) << total_time << endl;

TotalCharge -= len_acids * sqrt(alpha);

//cout << "Tag 5 " << endl;
//t0 = clock();

//MCA: Adding the Distace CV here
#pragma omp parallel for reduction(+:IonDistance)
 for(unsigned int i=0;i<len_acids;i++) {
   for(unsigned int k=i+1;k<len_acids;k++) {
     IonDistance -= distmod[i][k] * charge[i] * charge[k];
   }
 }

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "09 Ion distance compute: "  << fixed << setprecision(5) << total_time << endl;
//
////cout << "Tag 6 " << endl;
//t0 = clock();

//MCA: derivatives for the WaterpKw CV
#pragma omp parallel for
 for(unsigned int m=0;m<len_acids_hyd;m++) {
   for(unsigned int i=0;i<len_acids;i++) {   
     omp_deriv_tc[m] -= dfunc_theta[i] * dfunc_delta[i][m];
  }
}

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "10 omp_deriv_tc compute: "  << fixed << setprecision(5) << total_time << endl;
//
//t0 = clock();
//cout << "Tag 7 " << endl;

//MCA: deriv_distatives for the IonDistance CV
#pragma omp parallel for
for(unsigned int m=0;m<len_acids_hyd;m++) {
   for( unsigned int n=0;n<len_acids;n++) {
      if((m<len_acids)and(m!=n)) {
        omp_deriv_dist[m] += charge[m] * charge[n] * dist[m][n]/distmod[m][n];
      }
      if(distmod[n][m]<rcut) {
         for( unsigned int k=0;k<len_acids;k++) {
           omp_deriv_dist[m] += distmod[k][n] 
                    * charge[k] * dfunc_delta[n][m]; 
         }
      }
   }
}

//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "11 omp_deriv_dist compute: "  << fixed << setprecision(5) << total_time << endl;
//
////cout << "Tag 8 " << endl;
//t0 = clock();

#pragma omp critical
for(unsigned i=0;i<len_acids_hyd;i++) deriv_dist[i]+=omp_deriv_dist[i];
#pragma omp critical
for(unsigned i=0;i<len_acids_hyd;i++) deriv_tc[i]+=omp_deriv_tc[i];

 Value* vsd=getPntrToComponent("sd");
 Value* vtc=getPntrToComponent("tc");

 for(unsigned i=0;i<deriv_dist.size();++i) setAtomsDerivatives(vsd,i,deriv_dist[i]);
 //setValue           (vsd,IonDistance);
 vsd->set(IonDistance);
 //setBoxDerivatives  (vsd,virial_dist);
 
 for(unsigned i=0;i<deriv_tc.size();++i) setAtomsDerivatives(vtc,i,deriv_tc[i]);
 //setValue           (vsd,IonDistance);
 vtc->set(TotalCharge);
 //setBoxDerivatives  (vsd,virial_dist);
 //
//tf = clock();
//total_time = double(tf-t0)/double(CLOCKS_PER_SEC);
//cout << "12 critical compute: "  << fixed << setprecision(5) << total_time << endl;

 
 }
}
}
