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
#include "SolvationShellBase.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/Matrix.h"
#include "ActionRegister.h"
#include <string>
#include <cmath>
#include <iostream>


using namespace std;

namespace PLMD{
namespace colvar{

void SolvationShellBase::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","GROUPA","First Acid/Base group");
  keys.add("atoms","GROUPB","Second Acid/Base group)");
  keys.add("atoms","GROUPC","Third Acid/Base group");
  keys.add("atoms","GROUPD","List of Hydrogen Atoms");
  keys.add("compulsory","LAMBDA","1","The lambda parameter of the sum_exp function; 0 implies 1");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","D_1","0.0","The d_1 parameter of the switching function");
  keys.add("compulsory","D_2","0.0","The d_2 parameter of the switching function");
  keys.addOutputComponent("sp","default","Protonation state order parameter");
  keys.addOutputComponent("sd","default","Acid-base distance order parameter");
}

SolvationShellBase::SolvationShellBase(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista,gb_lista,gc_lista,gd_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);
  parseAtomList("GROUPC",gc_lista);
  parseAtomList("GROUPD",gd_lista);

  list_a = ga_lista;
  list_b = gb_lista;
  list_c = gc_lista;
  list_d = gd_lista;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parse("D_0",d0);
  parse("D_1",d1);
  parse("D_2",d2);
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
  addComponentWithDerivatives("sp"); componentIsNotPeriodic("sp");
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
  atoms.insert(atoms.end(),list_c.begin(),list_c.end());
  atoms.insert(atoms.end(),list_d.begin(),list_d.end());
  requestAtoms(atoms);


  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0;i<gb_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n  third group:\n");
  for(unsigned int i=0;i<gc_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", gc_lista[i].serial());
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

SolvationShellBase::~SolvationShellBase(){
  delete nl;
}

void SolvationShellBase::prepare(){
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
void SolvationShellBase::calculate()
{

 //The 3 scalar CVs
 double SolvationShell=0.0;
 double IonDistance=0.0;
 double TotalCharge=0.0;
 //Length of acids groups, with and without the H atoms
 unsigned len_acids = list_a.size()+list_b.size()+list_c.size();
 unsigned len_acids_hyd = len_acids + list_d.size();

 vector<double> sum_exp(len_acids_hyd);
 fill(sum_exp.begin(),sum_exp.end(),0.);

 Tensor virial;
 Tensor virial_dist; //MCA: IonDistance CV
 vector<Vector> deriv(len_acids_hyd);
 vector<Vector> deriv_dist(len_acids_hyd);
 vector<Vector> deriv_tc(len_acids_hyd);
 Vector zeros;
 zeros.zero();
 fill(deriv.begin(), deriv.end(), zeros);
 fill(deriv_dist.begin(), deriv_dist.end(), zeros);

 if(nl->getStride()>0 && invalidateList){
   nl->update(getPositions());
 }

unsigned nt=OpenMP::getNumThreads();

const unsigned nn=nl->size();

//if(nt*stride*10>nn) nt=nn/stride/10;
if(nt==0)nt=1;

 std::vector<Vector> omp_deriv(len_acids_hyd);
 std::vector<Vector> omp_deriv_dist(len_acids_hyd);
 fill(omp_deriv.begin(), omp_deriv.end(), zeros);
 fill(omp_deriv_dist.begin(), omp_deriv_dist.end(), zeros);
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
 d.insert(d.end(),list_a.size(),d0/list_a.size());
 d.insert(d.end(),list_b.size(),d1/list_b.size());
 d.insert(d.end(),list_c.size(),d2/list_c.size());

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
}

for(unsigned int j=len_acids;j<len_acids_hyd;j++) {   
   for(unsigned int i=0;i<len_acids;i++) {   
      sum_exp[j] += exp(lambda * distmod[i][j]);
   }
}

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

 //MCA: double check this vector assignment. It was (len_acids_hyd)**3 before 
 vector<vector<vector<double> > > dfunc_coord(len_acids, vector<vector<double> >(len_acids_hyd, vector<double>(len_acids)));

#pragma omp parallel for
 for(unsigned int i=0;i<len_acids;i++) {
   for(unsigned int j=len_acids;j<len_acids_hyd;j++){

     dfunc_coord[i][j][i] = lambda *  c[i][j] * (1 - c[i][j]);
     for(unsigned int n=i+1;n<len_acids;n++) {

       dfunc_coord[i][j][n] = -lambda *  c[n][j] * c[i][j];
       dfunc_coord[n][j][i] = dfunc_coord[i][j][n];
     }
   }
 }


vector<int> acid_index(len_acids);
 for(unsigned int i=0;i<len_acids;i++) {
     if(i<list_a.size()) { 
       acid_index[i]=0; 
     } else if(i<list_a.size()+list_b.size()) {
       acid_index[i]=1; 
     } else {
       acid_index[i]=2; 
     }
 }

 vector<double> square(3);

 square[0] = pow(2,0);
 square[1] = pow(2,1);
 square[2] = pow(2,2);

#pragma omp parallel for reduction(+:SolvationShell,TotalCharge)
 for(unsigned int i=0;i<len_acids;i++) {
     SolvationShell += square(acid_index[i])*charge[i];
     TotalCharge    += sqrt(pow(charge[i],2)+alpha)
     dfunc_theta[i]  = charge[i]/sqrt(pow(charge[i],2)+alpha);
 }

//MCA: Adding the Distace CV here
#pragma omp parallel for reduction(+:IonDistance)
 for(unsigned int i=0;i<len_acids;i++) {
   for(unsigned int k=i+1;k<len_acids;k++) {
     if(acid_index[i]!=acid_index[k]) {
       IonDistance -= distmod[i][k] * charge[i] * charge[k];
     }
   }
 }

//MCA: derivatives for the SolvationShell CV
#pragma omp parallel for
 for(unsigned int m=0;m<len_acids;m++) {
   for(unsigned int i=0;i<len_acids;i++) {   
     for(unsigned int j=len_acids;j<len_acids_hyd;j++) {
       omp_deriv[m] -= square[acid_index[i]] * dfunc_coord[i][j][m] 
                 * dist[m][j]/distmod[m][j];
     } 
  }
}
#pragma omp parallel for
 for(unsigned int m=len_acids;m<len_acids_hyd;m++) {
   for(unsigned int i=0;i<len_acids;i++) {   
     for(unsigned int n=0;n<len_acids;n++) {   
        omp_deriv[m] += square[acid_index[i]] * dfunc_coord[i][m][n] 
                  * dist[n][m]/distmod[n][m];
      } 
   }
}

//MCA: deriv_distatives for the IonDistance CV
#pragma omp parallel for
for(unsigned int m=0;m<len_acids;m++) {
   for( unsigned int n=0;n<len_acids;n++) {
      if(acid_index[m]!=acid_index[n]) {
        omp_deriv_dist[m] += charge[m] * charge[n] * dist[m][n]/distmod[m][n];
      }
      for( unsigned int k=n+1;k<len_acids;k++) {
         if(acid_index[n]!=acid_index[k]) {
            for(unsigned int h=len_acids;h<len_acids_hyd;h++) {

               //MCA: double check the i, k indexes
               omp_deriv_dist[m] += distmod[k][n] 
                        * ( charge[k] * dfunc_coord[n][h][m] 
                        +   charge[n] * dfunc_coord[k][h][m] ) 
                        * dist[m][h]/distmod[m][h];
            }
         }    
      }
   }
}
#pragma omp parallel for
for(unsigned int m=len_acids;m<len_acids_hyd;m++) {
   for(unsigned int i=0;i<len_acids;i++) {
      for(unsigned int k=i+1;k<len_acids;k++) { 
         //MCA: double check the i, k indexes
         if(acid_index[i]!=acid_index[k]){
            for( unsigned int n=0;n<len_acids;n++) {
               omp_deriv_dist[m] -= distmod[i][k] 
                              * ( charge[k] * dfunc_coord[i][m][n] 
                              +   charge[i] * dfunc_coord[k][m][n] ) 
                              * dist[n][m]/distmod[n][m];
            }          
         }    
      }
   }
}

#pragma omp critical
for(unsigned i=0;i<len_acids_hyd;i++) deriv[i]+=omp_deriv[i];
#pragma omp critical
for(unsigned i=0;i<len_acids_hyd;i++) deriv_dist[i]+=omp_deriv_dist[i];
//virial+=omp_virial;

 Value* vsp=getPntrToComponent("sp");
 Value* vsd=getPntrToComponent("sd");
 Value* vtc=getPntrToComponent("tc");

 for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(vsp,i,deriv[i]);
 //setValue           (vsp,SolvationShell);
 vsp->set(SolvationShell);
 //setBoxDerivatives  (vsp,virial);

 for(unsigned i=0;i<deriv_dist.size();++i) setAtomsDerivatives(vsd,i,deriv_dist[i]);
 //setValue           (vsd,IonDistance);
 vsd->set(IonDistance);
 //setBoxDerivatives  (vsd,virial_dist);
 
 for(unsigned i=0;i<deriv_tc.size();++i) setAtomsDerivatives(vtc,i,deriv_tc[i]);
 //setValue           (vsd,IonDistance);
 vsd->set(TotalCharge);
 //setBoxDerivatives  (vsd,virial_dist);
 
 }
}
}
