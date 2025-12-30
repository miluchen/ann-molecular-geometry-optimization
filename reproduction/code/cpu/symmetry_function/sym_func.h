/*
 * Header file for sym_func.cpp
 */

 #ifndef SYM_H
 #define SYM_H
 
 #include <vector>
 #include <string>
 #include <unordered_map>
 
 /*
 struct ... {
     double *rad; int size1;
     double *ang; int size2;
     ...
 };TODO
 */
 
 // Sym_Param structure
 struct Sym_Param {
     int atomtype_num;		// the number of distinct elements
     double rrc;				// cutoff radius for radial environment
     double arc;				// cutoff radius for angular environment
     
     int rad_eta_num;		// number of parameter etas for radial
     int rad_rs_num;		// number of Rs for radial
     int ang_zeta_num;
     int ang_theta_num;
     int ang_eta_num;
     int ang_rs_num;
 
     double *rad_etas;
     double *rad_rss;
     double *ang_zetas;
     double *ang_thetas;
     double *ang_etas;
     double *ang_rss;
 
     Sym_Param();
     Sym_Param(const Sym_Param &other); // copy constructor
     ~Sym_Param();
 
     void create(std::string fname); // create the parameters from a file
     void print(void); // print out the parameters
 };
 
 // AEV structure
 class AEV {
     private:
         int rad_dim;	// dimension of radial vector
         int ang_dim;	// dimension of angular vector
 
         // 4 types of atoms, 4 + 10 = 14 types AEV: H, C, O, N, HH, HC, HO, HN, CC, CO, CN, OO, ON, NN
         // n types of atoms, n + n*(n+1)/2 = (n^2+3n)/2 types AEV
         std::array<std::string, 14> atypes;
         // a map structure to store the index mapping
         std::unordered_map<std::string, int> indexmap;
 
         double rcutoff(double r);
         double acutoff(double r);
         double getDist(double *coords, int i, int j);
     public:
         int tot_dim;	// dimension of the atom environment vector (aev)
         Sym_Param params;
 
         // store the results of corresponding elements TODO
         /*double *carbon;
         double *hydrogen;
         double *nitrogen;
         double *oxygen;*/
 
         // the resulting atome environment vector
         /*double **O_features;
         double **O_desiredoutput;
 
         double **aevs;*/
 
         AEV();
         ~AEV();
 
         // initialize parameters
         void create(std::string fname);
         // build the aevectors
         double **buildaev(double *coords, char *atoms, int num);
         // build the radial aev
         void buildrad(double **radialaev, double **dists, char *atoms, int num);
         // build the angular aev
         void buildang(double **angularaev, double **dists, char *atoms, int num);
         // get the index given the atom types
         // int getIndex();
         // read in an xyz file
         // void readxyz(std::string fname);
         // write/append the aev into a file
         // writeaev(std::string fname, bool append = true);
         // compute the AEV and build the AEV vector
 };
 
 #endif
 