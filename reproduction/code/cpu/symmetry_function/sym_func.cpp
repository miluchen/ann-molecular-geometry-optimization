/*
 * Generate the atomic environment vector (AEV) given xyz coordinates
 */

 #include <cmath>
 #include <vector>
 #include <string>
 #include <cstring>
 #include <fstream>
 #include <cassert>
 #include <iostream>
 #include <algorithm>
 #include <unordered_map>
 
 #include "sym_func.h"
 
 #define ATOMTYPE_NUM 4
 
 using namespace std;
 
 /************************** Program parameters *********************************/
//  const bool _isDebug  = true;	// debug mode
 // the number below should be a long double, a double = 3.141592653589793, so as to compare with the CUDA results
 const double PI = 3.14159265358979323846264338328; 
 
 /************************** Helper functions ***********************************/
 
 /*************************** Member function for Sym_Param **********************/
 Sym_Param::Sym_Param(): atomtype_num(0), rrc(0), arc(0), rad_eta_num(0), rad_rs_num(0), ang_zeta_num(0), ang_theta_num(0),
                         ang_eta_num(0), ang_rs_num(0), rad_etas(0), rad_rss(0), ang_zetas(0), ang_thetas(0),
                         ang_etas(0), ang_rss(0) {}
 
 Sym_Param::~Sym_Param() {
     delete [] rad_etas;
     delete [] rad_rss;
     delete [] ang_zetas;
     delete [] ang_thetas;
     delete [] ang_etas;
     delete [] ang_rss;
 }
 
 // create the sym_param structure from a file, it doesn't check the validity of the input parameter file TODO
 void Sym_Param::create(string fname) {
     ifstream paramf (fname);
     // open the parameter file and initialize sym_param
     if (paramf.is_open()) {
         string line;
         paramf >> atomtype_num;		getline(paramf, line);
         paramf >> rrc;				getline(paramf, line);
         paramf >> arc;				getline(paramf, line);
         paramf >> rad_eta_num;		getline(paramf, line);	
         paramf >> rad_rs_num;		getline(paramf, line);
         paramf >> ang_zeta_num;		getline(paramf, line);
         paramf >> ang_theta_num;	getline(paramf, line);
         paramf >> ang_eta_num;		getline(paramf, line);
         paramf >> ang_rs_num;		getline(paramf, line);
         // allocate space
         rad_etas   = new double[rad_eta_num];
         rad_rss    = new double[rad_rs_num];
         ang_zetas  = new double[ang_zeta_num];
         ang_thetas = new double[ang_theta_num];
         ang_etas   = new double[ang_eta_num];
         ang_rss    = new double[ang_rs_num];
         // read in more parameters
         for (int i=0; i<rad_eta_num; i++)
             paramf >> rad_etas[i];
         getline(paramf, line);
         for (int i=0; i<rad_rs_num; i++)
             paramf >> rad_rss[i];
         getline(paramf, line);
         for (int i=0; i<ang_zeta_num; i++)
             paramf >> ang_zetas[i];
         getline(paramf, line);
         for (int i=0; i<ang_theta_num; i++)
             paramf >> ang_thetas[i];
         getline(paramf, line);
         for (int i=0; i<ang_eta_num; i++)
             paramf >> ang_etas[i];
         getline(paramf, line);
         for (int i=0; i<ang_rs_num; i++)
             paramf >> ang_rss[i];
     }
     else {
         cout << fname << ": can not open" << endl;
         assert(false);
     }
 }
 
 // print out the parameters
 void Sym_Param::print(void) {
     cout << atomtype_num	<< endl;
     cout << rrc				<< endl;
     cout << arc				<< endl;
     cout << rad_eta_num	<< endl;
     cout << rad_rs_num	<< endl;
     cout << ang_zeta_num	<< endl;
     cout << ang_theta_num	<< endl;
     cout << ang_eta_num		<< endl;
     cout << ang_rs_num		<< endl;
     for (int i=0; i<rad_eta_num; i++)
         cout << rad_etas[i] << " ";
     cout << endl;
     for (int i=0; i<rad_rs_num; i++)
         cout << rad_rss[i] << " ";
     cout << endl;
     for (int i=0; i<ang_zeta_num; i++)
         cout << ang_zetas[i] << " ";
     cout << endl;
     for (int i=0; i<ang_theta_num; i++)
         cout << ang_thetas[i] << " ";
     cout << endl;
     for (int i=0; i<ang_eta_num; i++)
         cout << ang_etas[i] << " ";
     cout << endl;
     for (int i=0; i<ang_rs_num; i++)
         cout << ang_rss[i] << " ";
     cout << endl;
 }
 
 /**************************** Member functions for AEV class ***************************/
 AEV::AEV(): rad_dim(0), ang_dim(0), tot_dim(0), params() {}
 
 AEV::~AEV() {}
 
 void AEV::create(string fname) {
     params.create(fname);
 
     int n = ATOMTYPE_NUM;
     int N = ATOMTYPE_NUM + ATOMTYPE_NUM*(ATOMTYPE_NUM+1)/2;
     // initialize the dimentions of the radial and angular vectors
     rad_dim = params.rad_eta_num * params.rad_rs_num;
     ang_dim = params.ang_eta_num * params.ang_theta_num * params.ang_eta_num * params.ang_rs_num;
     tot_dim = rad_dim * n + ang_dim * (N-n);
 
     // initialize the indexmap
     atypes = {"C", "H", "N", "O", "CC", "CH", "CN", "CO", "HH", "HN", "HO", "NN", "NO", "OO"}; // sorted array
     for (int i=0; i<n; i++)
         indexmap[atypes[i]] = i * rad_dim;
     for (int i=n; i<N; i++)
         indexmap[atypes[i]] = n * rad_dim + (i-n) * ang_dim;
 }
 
 double** AEV::buildaev(double *coords, char *atoms, int num) { // char *atoms may not be enough, need string TODO
     // need to build a distance matrix first
     double **dists = new double*[num * sizeof(double *)];
     for (int i=0; i<num; i++) {
         dists[i] = new double[num * sizeof(double)];
     }
     for (int i=0; i<num-1; i++) {
         for (int j=i+1; j<num; j++) {
             dists[i][j] = getDist(coords, i, j);
             dists[j][i] = dists[i][j];
         }
     }
     // build the atomic environment vector
     double **aevs = new double*[num * sizeof(double *)];
     for (int i=0; i<num; i++) {
         aevs[i] = new double[tot_dim * sizeof(double)];
         memset(aevs[i], 0, tot_dim * sizeof(double));
     }
 
     buildrad(aevs, dists, atoms, num); // this is not efficient TODO
     buildang(aevs, dists, atoms, num); // this is not efficient TODO
 
     for (int i=0; i<num; i++)
         delete [] dists[i];
     delete [] dists;
 
     return aevs;	// assume the caller free the space TODO
 }
 
 void AEV::buildrad(double **radialaev, double **dists, char *atoms, int num) {
     for (int i=0; i<num-1; i++) {	// maybe do it differently with GPU TODO
         for (int j=i+1; j<num; j++) {
             if (dists[i][j] >= params.rrc)
                 continue;
             string key1 = string() + atoms[j];
             string key2 = string() + atoms[i];
             int index1 = indexmap[key1];
             int index2 = indexmap[key2];
             for (int k=0; k<params.rad_eta_num; k++) {
                 for (int l=0; l<params.rad_rs_num; l++) {
                     double tmp = rcutoff(dists[i][j]) * exp(-params.rad_etas[k] * (dists[i][j] - params.rad_rss[l]) * (dists[i][j] - params.rad_rss[l]));
                     radialaev[i][index1++] += tmp;
                     radialaev[j][index2++] += tmp;
                 }
             }
         }
     }
 }
 
 void AEV::buildang(double **angularaev, double **dists, char *atoms, int num) { // not efficient! TODO no good
     for (int i=0; i<num; i++) {
         for (int j=0; j<num-1; j++) {
             // check the distance between atom and atom1
             if (i == j || dists[i][j] >= params.arc)
                 continue;
             for (int k=j+1; k<num; k++) {
                 // check the distance between atom and atom2
                 if (i == k || dists[i][k] >= params.arc)
                     continue;
                 // get the index
                 string key = atoms[j] > atoms[k] ? (string() + atoms[k] + atoms[j]) : (string() + atoms[j] + atoms[k]);
                 double theta = acos((dists[i][j] * dists[i][j] + dists[i][k] * dists[i][k] - dists[j][k] * dists[j][k]) / (2 * dists[i][j] * dists[i][k]));
                 int index = indexmap[key];
                 for (int ii=0; ii<params.ang_zeta_num; ii++) {
                     for (int jj=0; jj<params.ang_theta_num; jj++) {
                         for (int kk=0; kk<params.ang_eta_num; kk++) {
                             for (int ll=0; ll<params.ang_rs_num; ll++) {
                                 angularaev[i][index++] += acutoff(dists[i][j]) * acutoff(dists[i][k]) * pow(2, 1-params.ang_zetas[ii]) *
                                                             pow((1 + cos(theta - params.ang_thetas[jj])), params.ang_zetas[ii]) *
                                                             exp(-params.ang_etas[kk] * pow(((dists[i][j] + dists[i][k])/2 - params.ang_rss[ll]), 2));
                             }
                         }
                     }
                 }
             }
         }
     }
 }
 
 double AEV::rcutoff(double rij) {
     return (rij >= params.rrc) ? 0.0 : (0.5 * cos(PI*rij / params.rrc) + 0.5);
 }
 
 double AEV::acutoff(double rij) {
     return (rij >= params.arc) ? 0.0 : (0.5 * cos(PI*rij / params.arc) + 0.5);
 }
 
 double AEV::getDist(double *coords, int i, int j) {
     return sqrt((coords[3*i]-coords[3*j]) * (coords[3*i]-coords[3*j]) + 
                 (coords[3*i+1]-coords[3*j+1]) * (coords[3*i+1]-coords[3*j+1]) +
                 (coords[3*i+2]-coords[3*j+2]) * (coords[3*i+2]-coords[3*j+2]));
 }
 