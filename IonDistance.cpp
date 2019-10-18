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
#include "tools/SwitchingFunction.h"
#include "ActionRegister.h"
#include <cmath> 

#include <string>
#include <iostream>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR COORDINATION
/*
Calculate coordination numbers.

This keyword can be used to calculate the number of contacts between two groups of atoms
and is defined as
\f[
\sum_{i\in A} \sum_{i\in B} s_{ij}
\f]
where \f$s_{ij}\f$ is 1 if the contact between atoms \f$i\f$ and \f$j\f$ is formed,
zero otherwise.
In practise, \f$s_{ij}\f$ is replaced with a switching function to make it differentiable.
The default switching function is:
\f[
s_{ij} = \frac{ 1 - \left(\frac{{\bf r}_{ij}-d_0}{r_0}\right)^n } { 1 - \left(\frac{{\bf r}_{ij}-d_0}{r_0}\right)^m }
\f]
but it can be changed using the optional SWITCH option.

To make your calculation faster you can use a neighbor list, which makes it that only a
relevant subset of the pairwise distance are calculated at every step.

If GROUPB is empty, it will sum the \f$\frac{N(N-1)}{2}\f$ pairs in GROUPA. This avoids computing 
twice permuted indexes (e.g. pair (i,j) and (j,i)) thus running at twice the speed.

Notice that if there are common atoms between GROUPA and GROUPB the switching function should be
equal to one. These "self contacts" are discarded by plumed (since version 2.1),
so that they actually count as "zero".


\par Examples

The following example instructs plumed to calculate the total coordination number of the atoms in group 1-10 with the atoms in group 20-100.  For atoms 1-10 coordination numbers are calculated that count the number of atoms from the second group that are within 0.3 nm of the central atom.  A neighbour list is used to make this calculation faster, this neighbour list is updated every 100 steps.
\verbatim
COORDINATION GROUPA=1-10 GROUPB=20-100 R_0=0.3 NLIST NL_CUTOFF=0.5 NL_STRIDE=100 
\endverbatim

The following is a dummy example which should compute the value 0 because the self interaction
of atom 1 is skipped. Notice that in plumed 2.0 "self interactions" were not skipped, and the
same calculation should return 1.
\verbatim
c: COORDINATION GROUPA=1 GROUPB=1 R_0=0.3
PRINT ARG=c STRIDE=10
\endverbatim

\verbatim
c1: COORDINATION GROUPA=1-10 GROUPB=1-10 R_0=0.3
x: COORDINATION GROUPA=1-10 R_0=0.3
c2: COMBINE ARG=x COEFFICIENTS=2
# the two variables c1 and c2 should be identical, but the calculation of c2 is twice faster
# since it runs on half of the pairs. Notice that to get the same result you
# should double it
PRINT ARG=c1,c2 STRIDE=10
\endverbatim
See also \ref PRINT and \ref COMBINE



*/
//+ENDPLUMEDOC
   
class IonDistance : public IonDistanceBase{
  SwitchingFunction switchingFunction;

public:
  explicit IonDistance(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
//  virtual double pairing(double distance,double&dfunc,unsigned i,unsigned j)const;
  //virtual double pairing(double arg,double&dfunc, double sw_cutoff, double sw_shift, int p, int q);
  virtual double sw_func(double arg, double sw_cutoff, double sw_shift, int p, int q);
  virtual double dsw_func(double arg, double sw_cutoff, double sw_shift, int p, int q);
};

PLUMED_REGISTER_ACTION(IonDistance,"IONDISTANCE")

void IonDistance::registerKeywords( Keywords& keys ){
  IonDistanceBase::registerKeywords(keys);
//  keys.add("compulsory","NN","6","The n parameter of the switching function ");
//  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
//  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
//  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
//  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
//                               "The following provides information on the \\ref switchingfunction that are available. " 
//                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords."); 
}

IonDistance::IonDistance(const ActionOptions&ao):
Action(ao),
IonDistanceBase(ao) {}
//{
//
//  string sw,errors;
//  parse("SWITCH",sw);
//  if(sw.length()>0){
//    switchingFunction.set(sw,errors);
//    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
//  } else {
//    int nn=6;
//    int mm=0;
//    double d0=0.0;
//    double r0=0.0;
//    parse("R_0",r0);
//    if(r0<=0.0) error("R_0 should be explicitly specified and positive");
//    parse("D_0",d0);
//    parse("NN",nn);
//    parse("MM",mm);
//    switchingFunction.set(nn,mm,r0,d0);
//  }
//  
//  checkRead();
//
//  log<<"  contacts are counted with cutoff "<<switchingFunction.description()<<"\n";
//}

double IonDistance::sw_func(double arg, double sw_cutoff, double sw_shift, int p, int q){
 double sw;
 if(abs(arg-sw_shift)/sw_cutoff>(1.-100.0*epsilon) && abs(arg-sw_shift)/sw_cutoff<(1+100.0*epsilon)){
   sw=(double)p/(double)q;
 }else{
   sw = (1-pow (((arg-sw_shift)/sw_cutoff),p))/(1-pow (((arg-sw_shift)/sw_cutoff),q));
 }
 return(sw);
}

double IonDistance::dsw_func(double arg, double sw_cutoff, double sw_shift, int p, int q){
 double d_sw;
 if(abs(arg-sw_shift)/sw_cutoff>(1.-100.0*epsilon) && abs(arg-sw_shift)/sw_cutoff<(1+100.0*epsilon)){
   d_sw=(arg-sw_shift)/abs(arg-sw_shift)*(double)p*((double)p-(double)q)/(double)q;
 }else{
   d_sw = -((p*pow((arg-sw_shift)/sw_cutoff,-1+p))/(sw_cutoff*(1-pow((arg-sw_shift)/sw_cutoff,q)))) + (q*(1-pow((arg-sw_shift)/sw_cutoff,p))*pow((arg-sw_shift)/sw_cutoff,-1+q))/(sw_cutoff*pow(1-pow((arg-sw_shift)/sw_cutoff,q),2));
 }
 return(d_sw);
}




//double IonDistance::pairing(double arg,double&dfunc, double sw_cutoff, double sw_shift, int p, int q){
// sw = (1-pow(((arg-sw_shift)/sw_cutoff),p))/(1-pow(((arg-sw_shift)/sw_cutoff),q));
// d_sw = -(p*pow(((arg-sw_shift)/sw_cutoff),p-1))/((1-pow(((arg-sw_shift)/sw_cutoff),q))*sw_cutoff)+(q*(1-pow(((arg-sw_shift)/sw_cutoff),p))*pow(((arg-sw_shift)/sw_cutoff),q-1))/((pow(1-pow(((arg-sw_shift)/sw_cutoff),q)),2)*sw_cutoff);
//
// switchingFunction.set(p,q,sw_cutoff,sw_shift);
// return switchingFunction.calculate(arg,dfunc);
//}


}

}
