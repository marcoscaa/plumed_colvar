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
#include "IonDistanceBase.h"
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

void IonDistanceBase::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
  keys.add("atoms","GROUPC","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
//  keys.add("compulsory","NN_n","6","The n parameter of the switching function ");
//  keys.add("compulsory","MM_n","0","The m parameter of the switching function; 0 implies 2*NN");
//  keys.add("compulsory","NN_m","6","The n parameter of the switching function ");
//  keys.add("compulsory","MM_m","0","The m parameter of the switching function; 0 implies 2*NN");
//  keys.add("compulsory","NN_q","6","The n parameter of the switching function ");
//  keys.add("compulsory","MM_q","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","LAMBDA","1","The lambda parameter of the sum_exp function; 0 implies 1");
//  keys.add("compulsory","LESS_THAN","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","D_1","0.0","The d_0 parameter of the switching function");
//  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
//  keys.add("compulsory","N_0","The n_0 parameter of the switching function");
//  keys.add("compulsory","EQ_SHELL","Number of atoms in the shell at equilibrium");
//  keys.add("compulsory","SWITCH_SIGN","Sign of the switching function for coordination");
}

IonDistanceBase::IonDistanceBase(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista,gb_lista,gc_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);
  parseAtomList("GROUPC",gc_lista);


  list_a = ga_lista;
  list_b = gb_lista;
  list_c = gc_lista;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

//  parse("NN_n",pn);
//  parse("MM_n",qn);
//  parse("NN_m",pm);
//  parse("MM_m",qm);
//  parse("NN_q",pq);
//  parse("MM_q",qq);
  parse("D_0",d0);
  parse("D_1",d1);
//  parse("R_0",r0);
  parse("LAMBDA",lambda);
//  parse("LESS_THAN",less_than);

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
  
  addValueWithDerivatives(); setNotPeriodic();
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
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh){
   log.printf("  using neighbor lists with\n");
   log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }
}

IonDistanceBase::~IonDistanceBase(){
  delete nl;
}

void IonDistanceBase::prepare(){
  if(nl->getStride()>0){
    if(firsttime || (getStep()%nl->getStride()==0)){
      //requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    }else{
      //requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void IonDistanceBase::calculate()
{

 //double qsolv=0.;
 double IonDistance=0.0;
 //int lista_size = list_a.size();
 vector<double> sum_exp(getNumberOfAtoms());
 fill(sum_exp.begin(),sum_exp.end(),0.);

 Tensor virial;
 vector<Vector> deriv(getNumberOfAtoms());
 Vector zeros;
 zeros.zero();
 fill(deriv.begin(), deriv.end(), zeros);
 //deriv.resize(getPositions().size());

 if(nl->getStride()>0 && invalidateList){
   nl->update(getPositions());
 }

// unsigned stride=comm.Get_size();
// unsigned rank=comm.Get_rank();
// if(serial){
//   stride=1;
//   rank=0;
// }else{
//   stride=comm.Get_size();
//   rank=comm.Get_rank();
// }

unsigned nt=OpenMP::getNumThreads();

const unsigned nn=nl->size();

//if(nt*stride*10>nn) nt=nn/stride/10;
if(nt==0)nt=1;

#pragma omp parallel num_threads(nt)
{
 std::vector<Vector> omp_deriv(getPositions().size());
 Tensor omp_virial;

 Matrix<double> c(getNumberOfAtoms(),getNumberOfAtoms());
 vector<double> coord(list_a.size()+list_b.size());
 vector<double> charge(list_a.size()+list_b.size());
 fill(coord.begin(),coord.end(),0.);

 std::vector<double> d;
 d.insert(d.end(),list_a.size(),d0);
 d.insert(d.end(),list_b.size(),d1);

//#pragma omp for reduction(+:voronoi) nowait
// for(unsigned int i=rank;i<nn;i+=stride) {   
 for(unsigned int i=0;i<list_a.size()+list_b.size();i++) {
    for(unsigned int j=list_a.size()+list_b.size();j<list_a.size()+list_b.size()+list_c.size();j++) {

       Vector distance;

       if(pbc){
          distance=pbcDistance(getPosition(i),getPosition(j));
       } else {
          distance=delta(getPosition(i),getPosition(j));
       }

       sum_exp[j] += exp(lambda * distance.modulo());
    }
 }

 for(unsigned int i=0;i<list_a.size()+list_b.size();i++) {
    for(unsigned int j=list_a.size()+list_b.size();j<list_a.size()+list_b.size()+list_c.size();j++) {

       Vector distance;

       if(pbc){
        distance=pbcDistance(getPosition(i),getPosition(j));
       } else {
        distance=delta(getPosition(i),getPosition(j));
       }

       c[i][j] = exp( lambda * distance.modulo()) / sum_exp[j];
       coord[i] += c[i][j];
    }
    charge[i] = coord[i] - d[i];
 }

// if(!serial){
//    comm.Sum(voronoi);
//    comm.Sum(dfunc_vor);
// }

 vector<vector<vector<double> > > dfunc_coord(getNumberOfAtoms(), vector<vector<double> >(getNumberOfAtoms(), vector<double>(getNumberOfAtoms())));


// for(unsigned int i=rank;i<list_a.size();i+=stride) {
 for(unsigned int i=0;i<list_a.size()+list_b.size();i++) {
   for(unsigned int j=list_a.size()+list_b.size();j<getNumberOfAtoms();j++){
     for(unsigned int n=0;n<list_a.size()+list_b.size();n++) {
       int d_in;

       if (i == n) d_in = 1;
       if (i != n) d_in = 0;

       dfunc_coord[i][j][n] = lambda * c[n][j] * (d_in - c[i][j]);
     }
   }
 }

// if(!serial){
//    comm.Sum(abs_coord_tot);
//    comm.Sum(dfunc_abs_coord);
// }

 //Vector distance;

 for(unsigned int i=0;i<list_a.size();i++) {
   for(unsigned int k=list_a.size();k<list_a.size()+list_b.size();k++) {
     Vector distance_ik;

       if(pbc){
         distance_ik=pbcDistance(getPosition(i),getPosition(k));
       } else {
         distance_ik=delta(getPosition(i),getPosition(k));
       }

     IonDistance -= distance_ik.modulo() * charge[i] * charge[k];
//cout << " i = " << i << " k = " << k << " charge[i] = " << charge[i] << " charge[k] = " << charge[k] << " d_ik = " << distance_ik.modulo() << " tot = " << IonDistance << endl;
   }
 }

// if(!serial){
//    comm.Sum(IonDistance);
// }

 //Tensor virial;


 //for(unsigned int m=rank;m<getNumberOfAtoms();m+=stride) {
 for(unsigned int m=0;m<getNumberOfAtoms();m++) {

  Vector distance_nj;
  Vector distance_ik;

//#pragma omp for reduction(+:ncoord) nowait
      for(unsigned int i=0;i<list_a.size();i++) {   
        if (m == i){
          for(unsigned int k=list_a.size();k<list_a.size()+list_b.size();k++) {   
            if(pbc){
              distance_ik=pbcDistance(getPosition(i),getPosition(k));
            } else {
              distance_ik=delta(getPosition(i),getPosition(k));
            }

            deriv[m] += charge[i] * charge[k] * distance_ik/distance_ik.modulo();
          }
        continue;
        }

        for(unsigned int k=list_a.size();k<list_a.size()+list_b.size();k++) {
          if (m == k){
  
              if(pbc){
                distance_ik=pbcDistance(getPosition(i),getPosition(k));
              } else {
                distance_ik=delta(getPosition(i),getPosition(k));
              }
  
              deriv[m] -= charge[i] * charge[k] * distance_ik/distance_ik.modulo();

              continue;
          }
        }
      }

      for(unsigned int j=list_a.size()+list_b.size();j<getNumberOfAtoms();j++) {   
         if (m == j){
            for(unsigned int n=0;n<list_a.size()+list_b.size();n++) {   

               if(pbc){
                  distance_nj=pbcDistance(getPosition(n),getPosition(j));
               } else {
                  distance_nj=delta(getPosition(n),getPosition(j));
               }
   
               for(unsigned int i=0;i<list_a.size();i++) {
                  for(unsigned int k=list_a.size();k<list_a.size()+list_b.size();k++) { 

                     if(pbc){
                        distance_ik=pbcDistance(getPosition(i),getPosition(k));
                     } else {
                        distance_ik=delta(getPosition(i),getPosition(k));
                     }

                     deriv[m] -= distance_ik.modulo() * (charge[k] * dfunc_coord[i][j][n] + charge[i] * dfunc_coord[k][j][n]) * distance_nj/distance_nj.modulo();
               
                  }
               }
            }
            continue;
            }

         for(unsigned int n=0;n<list_a.size()+list_b.size();n++){   
            if (m == n){

               if(pbc){
                 distance_nj=pbcDistance(getPosition(n),getPosition(j));
               } else {
                 distance_nj=delta(getPosition(n),getPosition(j));
               }
   
               for(unsigned int i=0;i<list_a.size();i++) {
                  for(unsigned int k=list_a.size();k<list_a.size()+list_b.size();k++) { 
   
                     if(pbc){
                       distance_ik=pbcDistance(getPosition(i),getPosition(k));
                     } else {
                       distance_ik=delta(getPosition(i),getPosition(k));
                     }
     
     
                     deriv[m] += distance_ik.modulo() * (charge[k] * dfunc_coord[i][j][n] + charge[i] * dfunc_coord[k][j][n]) * distance_nj/distance_nj.modulo();
               
                  }
               }
               continue;
            }
         }
      }
  }
#pragma omp critical
 if(nt>1){
  for(int i=0;i<getPositions().size();i++) deriv[i]+=omp_deriv[i];
  virial+=omp_virial;
 }
}

// if(!serial){
//   comm.Sum(IonDistance);
//   if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
//   comm.Sum(virial);
// }

 for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
 setValue           (IonDistance);
 setBoxDerivatives  (virial);

 }
}
}
