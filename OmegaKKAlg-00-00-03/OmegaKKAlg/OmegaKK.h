#ifndef Physics_Analysis_OmegaKK_H
#define Physics_Analysis_OmegaKK_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include <vector>
#include "CLHEP/Vector/LorentzVector.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "CLHEP/Geometry/Point3D.h"
typedef HepGeom::Point3D<double> HepPoint3D;    

#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "VertexFit/WTrackParameter.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

//#include "VertexFit/ReadBeamParFromDb.h"

//using CLHEP::HepLorentzVector;

class OmegaKK : public Algorithm {

	public:
		OmegaKK(const std::string& name, ISvcLocator* pSvcLocator);
		StatusCode initialize();
		StatusCode execute();
		StatusCode finalize();  
	
		bool topology();

		void before_execute();

		bool find_charged_track();
		virtual bool check_charged_particle();

		bool find_gamma_track_and_assign_4mom();
		virtual bool check_gamma();

		bool pid_and_assign_4mom();
		virtual bool check_pid();


		virtual bool vertex_fit();
		VertexParameter get_vxpar();

		bool fit4c();
		
		bool fit5c();

		void write_ana_data();

        void corgen(HepMatrix &, HepVector &, int );
        void corset(HepSymMatrix &, HepMatrix &, int );
        void calibration(RecMdcKalTrack * , HepVector &, int );
	private:
		CLHEP::HepLorentzVector ecms;

		vector<int> iGood, ipGood, imGood, iGam;
	    vector<CLHEP::HepLorentzVector> vpGood, vmGood, m_Gam;
		int nGood, nCharge, m_nGam;
		vector<RecMdcTrack*> vGood;
		HepPoint3D IP;
		int m_runNo, m_event;

		vector<int> ipip, ipim, iKp, iKm;
		vector<CLHEP::HepLorentzVector> ppip, ppim, pKp, pKm;

		vector<CLHEP::HepLorentzVector> particle_four_mom;

        vector<CLHEP::HepLorentzVector> all_particles_4mom;

		vector<int> particle_type;

	    vector<HepLorentzVector> pGam;

		WTrackParameter wpip, wpim, wKp, wKm;
		KalmanKinematicFit *kmfit;
	
		int m_g1, m_g2;

//		vector<double> dedx_pid, tof1_pid, tof2_pid, prob_pid;
//		vector<double> ptrk_pid, cost_pid;
		vector<double> s_dang, s_eraw;

//		CLHEP::HepLorentzVector ecms;
		double m_vr0cut;
		double m_vz0cut;
		double m_energyThreshold;
		double m_gammaPhiCut;
		double m_gammaThetaCut;
		double m_gammaAngleCut;
		int m_test4C;
		int m_test5C;
		int m_checkDedx;
		int m_checkTof;
		double m_bar_energy_cut;
		double m_end_energy_cut;
		double m_dang_cut;                   // define Ntuples here
		double m_min_emctime;
		double m_max_emctime; 
		double m_bar_costheta_cut;
		double m_min_end_costheta_cut;
		double m_max_end_costheta_cut;
        double m_collider_energy;
        int m_save_4mom_for_mctruth;
        int m_5C_Calibration;

		NTuple::Tuple* m_tuple1;
		NTuple::Item<double> m_vx0;
		NTuple::Item<double> m_vy0;
		NTuple::Item<double> m_vz0;
		NTuple::Item<double>  m_vr0;
		NTuple::Item<double>  m_rvxy0;
		NTuple::Item<double>  m_rvz0;
		NTuple::Item<double>  m_rvphi0;

//		NTuple::Tuple*  m_tuple2;      // fake photon
//		NTuple::Item<double>  m_dthe;
//		NTuple::Item<double>  m_dphi;
//		NTuple::Item<double>  m_dang;
//		NTuple::Item<double>  m_eraw;


		NTuple::Tuple*  m_tuple_ana;     // rhopi 5C

//		NTuple::Tuple*  m_tuple4;     // rhopi 4C
		NTuple::Item<double>  m_chi1;
		NTuple::Item<double>  m_mpi0_4c;
        NTuple::Item<double>  m_chisq_1g;
        NTuple::Item<double>  m_chisq_3g;

		NTuple::Item<double> m_chi2;
		NTuple::Item<double> m_momega;
		NTuple::Item<double> m_mkaonp;
		NTuple::Item<double> m_mkaonm;

		NTuple::Item<double> m_pip;
		NTuple::Item<double> m_pim;

		NTuple::Item<double> m_mkaonmpip;
		NTuple::Item<double> m_mkaonppim;
		NTuple::Item<double> m_mpi0_5c;
		NTuple::Item<double> m_mkaon;
		NTuple::Item<double> m_mtot;
		NTuple::Item<double> m_fcos;

		NTuple::Item<long> m_particle_index;
		NTuple::Matrix<double> m4xyz;
		NTuple::Matrix<double> m4rft;

//		NTuple::Item<long> m_charged_index;
//		NTuple::Array<double> m_ptrk_pid;
//		NTuple::Array<double> m_cost_pid;
//		NTuple::Array<double> m_dedx_pid;
//		NTuple::Array<double> m_tof1_pid;
//		NTuple::Array<double> m_tof2_pid;
//		NTuple::Array<double> m_prob_pid;	

		NTuple::Item<long> m_gamma_index;
		NTuple::Array<double> m_dang;
		NTuple::Array<double> m_eraw;

		NTuple::Item<double>  m_m2gg;
		NTuple::Item<double>  m_etot;

		NTuple::Item<long> m_runnr;
		NTuple::Item<long> m_recnr;
		NTuple::Item<long> m_indexmc;
		NTuple::Array<long> m_pdgid;
		NTuple::Array<long> m_motheridx;
		
		NTuple::Tuple* m_tuple_truth;
		NTuple::Item<long> m_truth_index;
		NTuple::Matrix<double> m_truth_mom;
};

#endif 
