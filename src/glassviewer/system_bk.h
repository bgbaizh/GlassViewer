#include <iostream>
#include <exception>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "atom.h"
#include <mutex> 
#include <wignerSymbols.h>
#include "fftw3.h"


namespace py = pybind11;
using namespace std;

struct cell{
  vector<int> members;
  vector<int> neighbor_cells;
};


class System{

    public:

        //-----------------------------------------------------
        // Constructor, Destructor and Access functions
        //-----------------------------------------------------
        System();
        ~System();
        int nop;
        int ghost_nop;
        int real_nop;

        //-----------------------------------------------------
        // Simulation box related methods
        //-----------------------------------------------------
        double rot[3][3];
        double rotinv[3][3];
        int triclinic;
        double boxx, boxy, boxz;//the length of the 3 egdes
        double boxdims[3][2];
        double box[3][3];
        void assign_triclinic_params(vector<vector<double>>, vector<vector<double>>);
        vector<vector<double>> get_triclinic_params();
        void sbox(vector<vector<double>>);
        vector<vector<double>> gbox();
        vector<double> remap_atom(vector<double>);

        //-----------------------------------------------------
        // Atom related methods
        //-----------------------------------------------------
        vector<Atom> atoms;
        void assign_particles( vector<Atom>);
        void read_particle_file(string);    // TBDep
        void set_atoms( vector<Atom>);
        vector<Atom> get_atoms();
        void add_atoms( vector<Atom>);
        vector<Atom> get_all_atoms();
        Atom gatom(int);
        void satom(Atom);


        //----------------------------------------------------
        // Neighbor methods
        //----------------------------------------------------
        int filter;
        int usecells;
        int nx, ny, nz;
        int total_cells;
        cell *cells;
        double neighbordistance;
        void get_all_neighbors_normal();
        void process_neighbor(int, int);
        int get_all_neighbors_sann(double);
        int get_all_neighbors_bynumber(double, int, int,vector<int> atomlist = vector<int>());
        int get_neighbors_from_temp(int);
        int get_all_neighbors_adaptive(double, int, double);
        void get_all_neighbors_voronoi();
        void reset_all_neighbors(vector<int> atomlist = vector<int>());
        void reset_main_neighbors();        
        double get_abs_distance(int,int,double&,double&,double&);
        double get_abs_distance(Atom , Atom );
        vector<double> get_distance_vector(Atom , Atom);
        //mutex pdfreslock;
        mutex pdfthreadflaglock;
        void set_neighbordistance(double);
        vector<int> get_pairdistances(double cut,bool partial,int centertype,int secondtype,int histnum,double histlow,int threadnum);
        class pdfpara{
            public:
            vector<int> res;
            vector<vector<int>> resthread;
            double deltacut;
            //double d_square,d;
            //double diffx,diffy,diffz;
            //int M[3]={0};
            //int N[3]={0};
            double Height[3]={0};
            double iCrossj[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            double iCrossjnorm[3]={0};
            double kdotiCrossj[3]={0};
            int index[3][3]={{0,1,2},{1,2,0},{2,0,1}};// 计算叉乘的时候，以角标k i j 为顺序 按照index数组的顺序进行计算
            double histlow_square;
            double cut_square;
            bool halftimes=false;
            double pointHeight[3]={0};
            double pdotiCrossj[3]={0};
            double cut;
            bool partial;
            int centertype;
            int secondtype;
            int histnum;
            double histlow;
            int threadnum;
            bool *threadflag;
        };
        static void pairditancethread(int atomsstart, int atomsfinish, int threadid, System* sys, pdfpara* s);
        bool pdf_halftimes;
        vector<int> get_pairangle(double histlow,double histhigh,int histnum);
        double get_angle(int,int,int);
        //variables for a filter
        void susecells(int);
        int gusecells();
        int cell_index(int, int, int);
        void set_up_cells();
        vector<int> cell_periodic(int, int, int);
        void get_all_neighbors_cells();
        void get_temp_neighbors_cells(vector<int> atomlist = vector<int>());
        void get_temp_neighbors_brute(vector<int> atomlist = vector<int>());
        void store_neighbor_info();
        void set_atom_cutoff(double);

        //---------------------------------------------------
        // Methods for q calculation
        //---------------------------------------------------
        int *reqdqs;
        int lenqs;
        int *reqdaqs;
        int lenaqs;
        vector<double> gqvals(int qq);
        vector<double> gaqvals(int qq);
        vector<int> rq_backup;
        void set_reqd_qs(vector<int>);
        void set_reqd_aqs(vector<int>);
        void calculate_w(vector <int>, vector <int> atomlist,bool);

        void calculate_q(vector <int>,vector <int> atomlist);
        void calculate_aq(vector <int>,vector <int> atomlist);
        double dfactorial(int ,int );
        void convert_to_spherical_coordinates(double , double , double , double &, double &, double &);
        double PLM(int, int, double);
        void YLM(int , int , double , double , double &, double &);
        void QLM(int ,int ,double ,double ,double &, double & );
        void calculate_complexQLM_6();
        double get_number_from_bond(int,int);
        double get_number_from_bond(Atom ,Atom );
        void calculate_frenkel_numbers();
        //disorder vars
        void calculate_disorder();
        void find_average_disorder();        
        //---------------------------------------------------
        // Methods for Global_Statics calculation
        //---------------------------------------------------
        vector<vector<double>>bondpos;
        vector<vector<double>>bondvec;
        vector<vector<vector<vector<double>>>> bondqlm;
        vector<vector<vector<double>>> global_Qlm;
        vector<double>global_Ql;
        vector<double>global_Wl;
        vector<double>global_Wlnorm;
        void GlobalBOO_Bond( vector <int> atomlist);
        void  GlobalBOO_Sum(vector <int>);
        vector<vector<double>> GlobalBOO_CF(vector <int>qs, double cut, int histnum, double histlow, bool norm, int n1, int n2, int n3,bool ffton);
        //-----------------------------------------------------
        // Solids and Clustering methods
        //-----------------------------------------------------
        double minfrenkel;
        double threshold;
        double avgthreshold;
        int maxclusterid;
        int solidq;
        int criteria;
        int comparecriteria;
        void find_solid_atoms();
        void find_clusters(double);
        void harvest_cluster(const int, const int);
        void find_clusters_recursive(double);
        int largest_cluster();
        void set_nucsize_parameters(double,double,double);
        void get_largest_cluster_atoms();

        //-----------------------------------------------------
        // Voronoi based methods
        //-----------------------------------------------------
        int alpha;
        void find_average_volume();
        int voronoiused;
        double face_cutoff;
        void set_face_cutoff(double);

        //-------------------------------------------------------
        // CNA parameters
        //-------------------------------------------------------
        double lattice_constant;
        void get_diamond_neighbors();
        vector<int> identify_diamond_structure();
        void identify_cndia();
        void get_cna_neighbors(int);
        void get_acna_neighbors(int);
        void get_common_neighbors(int);
        void get_common_bonds(int);
        void identify_cn12();
        void identify_cn14();
        vector<int> calculate_cna(int);

        //-------------------------------------------------------
        // Other order parameters
        //-------------------------------------------------------
        double switching_fn(double, double, int, int);
        void average_entropy();
        void average_entropy_switch(double, int, int);
        void entropy(double, double, double, double, double, double);
        void calculate_centrosymmetry_atom(int, int);
        void calculate_centrosymmetry(int);
        vector<double> get_centrosymmetry();

};
