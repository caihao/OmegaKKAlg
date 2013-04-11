#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "McTruth/McParticle.h"

#include "DstEvent/TofHitStatus.h"

#include "TMath.h"
#include "TRandom.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
//typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "OmegaKKAlg/OmegaKK.h"

//#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"

#include <vector>

const double ME  = 0.000511;
const double MU  = 0.105658;
const double MPI = 0.139570;
const double MK  = 0.493677;
const double MP  = 0.938272;
//const double xmass[5] = {ME, MU, MPI, MK, MP};

//const double velc = 29.9792458;  tof_path unit in cm.
//const double velc = 299.792458;   // tof path unit in mm


static long m_cout_all(0);

int Ncut0, Ncut1, Ncut2, Ncut3, Ncut4, Ncut5, Ncut6;

/////////////////////////////////////////////////////////////////////////////

OmegaKK::OmegaKK(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) {

		//Declare the properties  
		declareProperty("Vr0cut",				m_vr0cut				= 1.0);
		declareProperty("Vz0cut",				m_vz0cut				= 10.0);
		declareProperty("EnergyThreshold",		m_energyThreshold		= 0.04);
		declareProperty("GammaPhiCut",			m_gammaPhiCut			= 20.0);
		declareProperty("GammaThetaCut",		m_gammaThetaCut			= 20.0);
		declareProperty("GammaAngleCut",		m_gammaAngleCut			= 20.0);
		declareProperty("Test4C",				m_test4C				= 1);
		declareProperty("Test5C",				m_test5C				= 1);
		declareProperty("CheckDedx",			m_checkDedx				= 0);
		declareProperty("CheckTof",				m_checkTof				= 0);

		declareProperty("BarEneCut",			m_bar_energy_cut		= 0.025);
		declareProperty("EndEneCut",			m_end_energy_cut		= 0.050);
		declareProperty("DangCut",				m_dang_cut				= 20);
		declareProperty("MinEstCut",			m_min_emctime			= -0.01);
		declareProperty("MaxEstCut",			m_max_emctime			= 14.01);
		declareProperty("BarCosthetaCut",		m_bar_costheta_cut		= 0.8);
		declareProperty("MinEndCosthetaCut",	m_min_end_costheta_cut	= 0.86);
		declareProperty("MaxEndCosthetaCut",	m_max_end_costheta_cut	= 0.92);
    
        declareProperty("EnergyOfCollider",     m_collider_energy       = 3.686);
        declareProperty("Save4MomForMCTruth",   m_save_4mom_for_mctruth = 0);
        declareProperty("Do5CCalibration",        m_5C_Calibration        = 0);

	}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode OmegaKK::initialize()
{
	StatusCode status;
	MsgStream log(msgSvc(), name());

	m_cout_all++;

	NTuplePtr nt1(ntupleSvc(), "FILE1/vxyz");
	if ( nt1 ) 
		m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("FILE1/vxyz", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if ( m_tuple1 ) {
			status = m_tuple1->addItem("vx0",		m_vx0);
			status = m_tuple1->addItem("vy0",		m_vy0);
			status = m_tuple1->addItem("vz0",		m_vz0);
			status = m_tuple1->addItem("vr0",		m_vr0);
			status = m_tuple1->addItem("rvxy0",		m_rvxy0);
			status = m_tuple1->addItem("rvz0",		m_rvz0);
			status = m_tuple1->addItem("rvphi0",	m_rvphi0);
		} else { 
			return StatusCode::FAILURE;
		}
	}
	NTuplePtr nt_truth(ntupleSvc(), "FILE1/truth");
	if ( nt_truth )
		m_tuple_truth = nt_truth;
	else {
		m_tuple_truth = ntupleSvc()->book("FILE1/truth", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if ( m_tuple_truth ) {
			status = m_tuple_truth->addItem("truth_index", m_truth_index, 0, 6);
			cout << "Begin truth!" << endl;
			if (m_save_4mom_for_mctruth != 0) 
                status = m_tuple_truth->addIndexedItem("truth_mom", m_truth_index, 4, m_truth_mom);
		} else {
			return StatusCode::FAILURE;
		}
	}

	NTuplePtr nt_ana(ntupleSvc(), "FILE1/ana");
	if ( nt_ana )
		m_tuple_ana = nt_ana;
	else {
		m_tuple_ana = ntupleSvc()->book("FILE1/ana", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if ( m_tuple_ana) {
			if (m_test4C == 1) {
				status = m_tuple_ana->addItem("chi1",   m_chi1);
				status = m_tuple_ana->addItem("mpi0_4c",   m_mpi0_4c);
                status = m_tuple_ana->addItem("chisq_1g", m_chisq_1g);
                status = m_tuple_ana->addItem("chisq_3g", m_chisq_3g);

			}
			if (m_test5C == 1) {
				status = m_tuple_ana->addItem("chi2",		m_chi2);
				status = m_tuple_ana->addItem("mpi0_5c",	m_mpi0_5c);
				status = m_tuple_ana->addItem("fcos",		m_fcos);
			}
			if (m_test5C ==1 || m_test4C ==1) {
				status = m_tuple_ana->addItem("momega",	m_momega);
				status = m_tuple_ana->addItem("mkaonp",	m_mkaonp);
				status = m_tuple_ana->addItem("mkaonm",	m_mkaonm);

				status = m_tuple_ana->addItem("mpip", m_pip);
				status = m_tuple_ana->addItem("mpim", m_pim);

				status = m_tuple_ana->addItem("mkaonppim",	m_mkaonppim);
				status = m_tuple_ana->addItem("mkaonmpip",	m_mkaonmpip);
				status = m_tuple_ana->addItem("mkaon",		m_mkaon);
				status = m_tuple_ana->addItem("mtot",		m_mtot);
			}
			status = m_tuple_ana->addItem("particle_index", m_particle_index, 0, 6);
			status = m_tuple_ana->addIndexedItem("m4xyz", m_particle_index, 4, m4xyz);
//			status = m_tuple_ana->addIndexedItem("m4rft", m_particle_index, 4, m4rft);

//			status = m_tuple_ana->addItem("charged_index", m_charged_index, 0, 4);
//			status = m_tuple_ana->addIndexedItem("ptrk_pid", m_charged_index,	m_ptrk_pid);
//			status = m_tuple_ana->addIndexedItem("cost_pid", m_charged_index,	m_cost_pid);
//			status = m_tuple_ana->addIndexedItem("dedx_pid", m_charged_index,	m_dedx_pid);
//			status = m_tuple_ana->addIndexedItem("tof1_pid", m_charged_index,	m_tof1_pid);
//			status = m_tuple_ana->addIndexedItem("tof2_pid", m_charged_index,	m_tof2_pid);
//			status = m_tuple_ana->addIndexedItem("prob_pid", m_charged_index,	m_prob_pid);

			status = m_tuple_ana->addItem("gamma_index", m_gamma_index, 0, 2);
			status = m_tuple_ana->addIndexedItem("dang", m_gamma_index, m_dang);
			status = m_tuple_ana->addIndexedItem("eraw", m_gamma_index, m_eraw);

			status = m_tuple_ana->addItem("m2gg",   m_m2gg); // invariant mass of two photons 
			status = m_tuple_ana->addItem("etot",   m_etot); // cut for e and mu

			status = m_tuple_ana->addItem("runnr",	m_runnr);
			status = m_tuple_ana->addItem("recnr",	m_recnr);
			status = m_tuple_ana->addItem("indexmc",	m_indexmc, 0, 100);
			status = m_tuple_ana->addIndexedItem("pdgid", m_indexmc, m_pdgid);
			status = m_tuple_ana->addIndexedItem("motheridx", m_indexmc, m_motheridx);
		} else {
			return StatusCode::FAILURE;
		}
	}
	ecms = HepLorentzVector(0.04, 0, 0, m_collider_energy);

	return StatusCode::SUCCESS;
}

bool OmegaKK::topology() {
	SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
	int m_numParticle = 0;
//	int m_truth_index = 0;

	if (!mcParticleCol) {
		return false;
	} else {
		bool psipDecay = false;
		int rootIndex = -1;
		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		for (; iter_mc != mcParticleCol->end(); iter_mc++) {
			if ((*iter_mc)->primaryParticle()) continue;
			if (!(*iter_mc)->decayFromGenerator()) continue;
			if ((*iter_mc)->particleProperty() == 100443) {
				psipDecay = true;
				rootIndex = (*iter_mc)->trackIndex();
			}
			if (!psipDecay) continue;
			int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
			int pdgid = (*iter_mc)->particleProperty();
			m_pdgid[m_numParticle] = pdgid;
			m_motheridx[m_numParticle] = mcidx;
			m_numParticle += 1;

            if (m_save_4mom_for_mctruth == 0) continue;
			if ((*iter_mc)->particleProperty() == 321) {
				HepLorentzVector vt = (*iter_mc)->initialFourMomentum();	
				m_truth_mom[0][0] = vt.e();
				m_truth_mom[0][1] = vt.vect().x();
				m_truth_mom[0][2] = vt.vect().y();
				m_truth_mom[0][3] = vt.vect().z();
				++m_truth_index;
//				cout << m_truth_index << endl;
			}
			if ((*iter_mc)->particleProperty() == -321) {
				HepLorentzVector vt = (*iter_mc)->initialFourMomentum();	
				m_truth_mom[1][0] = vt.e();
				m_truth_mom[1][1] = vt.vect().x();
				m_truth_mom[1][2] = vt.vect().y();
				m_truth_mom[1][3] = vt.vect().z();
				++m_truth_index;
//				cout << "2 = " << m_truth_index << endl;
			}
			if ((*iter_mc)->particleProperty() == 223) {
				HepLorentzVector vt = (*iter_mc)->initialFourMomentum();	
				m_truth_mom[2][0] = vt.e();
				m_truth_mom[2][1] = vt.vect().x();
				m_truth_mom[2][2] = vt.vect().y();
				m_truth_mom[2][3] = vt.vect().z();
//				cout << "3 = " << m_truth_index << endl;
				++m_truth_index;
			}
			if ((*iter_mc)->particleProperty() == 111) {
				HepLorentzVector vt = (*iter_mc)->initialFourMomentum();	
				m_truth_mom[3][0] = vt.e();
				m_truth_mom[3][1] = vt.vect().x();
				m_truth_mom[3][2] = vt.vect().y();
				m_truth_mom[3][3] = vt.vect().z();
				++m_truth_index;
			}
			if ((*iter_mc)->particleProperty() == 211) {
				HepLorentzVector vt = (*iter_mc)->initialFourMomentum();	
				m_truth_mom[4][0] = vt.e();
				m_truth_mom[4][1] = vt.vect().x();
				m_truth_mom[4][2] = vt.vect().y();
				m_truth_mom[4][3] = vt.vect().z();
				++m_truth_index;
			}
			if ((*iter_mc)->particleProperty() == -211) {
				HepLorentzVector vt = (*iter_mc)->initialFourMomentum();	
				m_truth_mom[5][0] = vt.e();
				m_truth_mom[5][1] = vt.vect().x();
				m_truth_mom[5][2] = vt.vect().y();
				m_truth_mom[5][3] = vt.vect().z();
				++m_truth_index;
			}

		}	
//		cout << "Truth ----> Write!" << endl;
	//	m_truth_index = m_truth_index;
		if (m_save_4mom_for_mctruth != 0) m_tuple_truth->write();
	}
	m_indexmc = m_numParticle;
	//		m_tuple12->write();
	return true;
}
void OmegaKK::before_execute() {
	iGood.clear(); ipGood.clear(); imGood.clear();	iGam.clear();
	vpGood.clear(); vmGood.clear(); m_Gam.clear();
	vGood.clear();

	ipip.clear(); ipim.clear(); iKp.clear(); iKm.clear();
	ppip.clear(); ppim.clear(); pKp.clear(); pKm.clear();

	particle_four_mom.clear();
    all_particles_4mom.clear();
	particle_type.clear();

	pGam.clear();

//	dedx_pid.clear(); tof1_pid.clear(); tof2_pid.clear(); prob_pid.clear();
//	ptrk_pid.clear(); cost_pid.clear();
	s_dang.clear(); s_eraw.clear();

	nCharge = 0;
}
bool OmegaKK::find_charged_track() {
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	Hep3Vector xorigin(0,0,0);
	//if (m_reader.isRunNumberValid(runNo)) {
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if (vtxsvc->isVertexValid()) {
		double* dbv = vtxsvc->PrimaryVertex(); 
		double* vv = vtxsvc->SigmaPrimaryVertex();  
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}
	IP = HepPoint3D(xorigin[0], xorigin[1], xorigin[2]); 
//	HepPoint3D point0(0. ,0. ,0.);   // the initial point for MDC recosntruction
	double xv = xorigin.x();
	double yv = xorigin.y();

	for(int i = 0; i < evtRecEvent->totalCharged(); i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
		if (!(*itTrk)->isMdcTrackValid()) continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		double pch = mdcTrk->p();
		double x0 = mdcTrk->x();
		double y0 = mdcTrk->y();
		double z0 = mdcTrk->z();
		double phi0 = mdcTrk->helix(1);
		double Rxy = (x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0);
		m_vx0 = x0;
		m_vy0 = y0;
		m_vz0 = z0;
		m_vr0 = Rxy;

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		//HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
		VFHelix helixip(HepPoint3D(0, 0, 0),a,Ea); 
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();

		m_rvxy0 = fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		m_rvz0 = vecipa[3];         //the nearest distance to IP in z direction
		m_rvphi0 = vecipa[1];

//		m_tuple1->write();

		if(fabs(m_rvz0) >= m_vz0cut) continue;
		if(fabs(m_rvxy0) >= m_vr0cut) continue;

		iGood.push_back(i);
		vGood.push_back(mdcTrk); // all tracks are stored here
		if (mdcTrk->charge() > 0) {
			ipGood.push_back(i);
		} else {
			imGood.push_back(i);
		}
		nCharge += mdcTrk->charge();
	}
	nGood = iGood.size();
	if (!check_charged_particle()) return false;
	return true;
}
bool OmegaKK::check_charged_particle() {
	return ((ipGood.size() == 2) && (imGood.size() == 2) && (nCharge == 0));
}

bool OmegaKK::find_gamma_track_and_assign_4mom() {
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();

		double m_abs_costhe(fabs(cos(emcTrk->theta())));
		if (m_abs_costhe > m_bar_costheta_cut && m_abs_costhe < m_min_end_costheta_cut) continue;
		if (m_abs_costhe > m_max_end_costheta_cut) continue;
		double eraw = emcTrk->energy();
//		log << MSG::DEBUG << "i = " << i << ";  ene = " << eraw << endreq;
		if (m_abs_costhe < m_bar_costheta_cut) {
			if (eraw < m_bar_energy_cut) continue;
		} else {
			if (eraw < m_end_energy_cut) continue;
		}
//		if (m_runNo > 0) {
			if(emcTrk->time() < m_min_emctime || emcTrk->time() > m_max_emctime) continue;	
//		}

		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
		// find the nearest charged track
		double dthe = 200.;
		double dphi = 200.;
		double dang = 200.; 
		for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition();
			//      double ctht = extpos.cosTheta(emcpos);
			double angd = extpos.angle(emcpos);
			double thed = extpos.theta() - emcpos.theta();
			double phid = extpos.deltaPhi(emcpos);
			thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			if (fabs(thed) < fabs(dthe)) dthe = thed;
			if (fabs(phid) < fabs(dphi)) dphi = phid;
			if(angd < dang){
				dang = angd;
			}
		}
		if(dang>=200) continue;
//		double eraw = emcTrk->energy();
		dthe = dthe * 180 / (CLHEP::pi);
		dphi = dphi * 180 / (CLHEP::pi);
		dang = dang * 180 / (CLHEP::pi);
//		m_dthe = dthe;
//		m_dphi = dphi;
//		m_dang = dang;
//		m_eraw = eraw;
//		m_tuple2->write();
//		if(eraw < m_energyThreshold) continue;
		//    if((fabs(dthe) < m_gammaThetaCut) && (fabs(dphi)<m_gammaPhiCut) ) continue;
		if(fabs(dang) < m_gammaAngleCut) continue;

		iGam.push_back(i);
		HepLorentzVector m_lv_temp(emcpos.unit() *= eraw, eraw);
		m_Gam.push_back(m_lv_temp);
		s_dang.push_back(dang);
		s_eraw.push_back(eraw);
	}
	m_nGam = iGam.size();
	if (!check_gamma()) return false;

	for(int i = 0; i < m_nGam; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i]; 
		RecEmcShower* emcTrk = (*itTrk)->emcShower();
		double eraw = emcTrk->energy();
		double phi = emcTrk->phi();
		double the = emcTrk->theta();
		HepLorentzVector ptrk;
		ptrk.setPx(eraw * sin(the) * cos(phi));
		ptrk.setPy(eraw * sin(the) * sin(phi));
		ptrk.setPz(eraw * cos(the));
		ptrk.setE(eraw);

		pGam.push_back(ptrk);
	}
	return true;
}
bool OmegaKK::check_gamma() {
	return m_nGam >= 2;
}

bool OmegaKK::pid_and_assign_4mom() {
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	ParticleID *pid = ParticleID::instance();
	for(int i = 0; i < nGood; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];

		//  if(pid) delete pid;
		pid->init();
		pid->setMethod(pid->methodProbability());
		//  pid->setMethod(pid->methodLikelihood());  //for Likelihood Method  

		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // use PID sub-system
		pid->identify(pid->onlyPion() | pid->onlyKaon());    // seperater Pion/Kaon
		//  pid->identify(pid->onlyPion());
		//  pid->identify(pid->onlyKaon());
		pid->calculate();
		if (!(pid->IsPidInfoValid())) continue;
		if (pid->probPion() < 0.001 && pid->probKaon() < 0.001) return false;

		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();

//		ptrk_pid.push_back(mdcTrk->p());
//		cost_pid.push_back(mdcTrk->theta());

		if (pid->probPion() > pid->probKaon()) {
//			dedx_pid.push_back(pid->chiDedx(2));
//			tof1_pid.push_back(pid->chiTof1(2));
//          tof2_pid.push_back(pid->chiTof2(2));
//			prob_pid.push_back(pid->probPion()); 

			particle_four_mom.push_back(vGood[i]->p4(MPI));

		ctivate.adobe.com	RecMdcKalTrack* mdcKalTrk_pi = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
			RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);//PID can set to electron, muon, pion, kaon and proton;

			HepLorentzVector ptrk;
			ptrk.setPx(mdcKalTrk_pi->px());
			ptrk.setPy(mdcKalTrk_pi->py());
			ptrk.setPz(mdcKalTrk_pi->pz());
			double p3 = ptrk.mag();

			ptrk.setE(sqrt(p3*p3+MPI*MPI));

			if(mdcKalTrk_pi->charge() >0) {
				ipip.push_back(iGood[i]);
				ppip.push_back(ptrk);

				particle_type.push_back(1); // 1 represents pi+
			} else {
				ipim.push_back(iGood[i]);
				ppim.push_back(ptrk);

				particle_type.push_back(2); // 2 represents pi-
			}
		} else {
//			dedx_pid.push_back(pid->chiDedx(3));
//			tof1_pid.push_back(pid->chiTof1(3));
//          tof2_pid.push_back(pid->chiTof2(3));
//			prob_pid.push_back(pid->probPion()); 

			particle_four_mom.push_back(vGood[i]->p4(MK));

			RecMdcKalTrack* mdcKalTrk_K = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
			RecMdcKalTrack::setPidType  (RecMdcKalTrack::kaon);//PID can set to electron, muon, pion, kaon and proton;The default setting is pion

			HepLorentzVector ptrk;
			ptrk.setPx(mdcKalTrk_K->px());
			ptrk.setPy(mdcKalTrk_K->py());
			ptrk.setPz(mdcKalTrk_K->pz());
			double p3 = ptrk.mag();

			ptrk.setE(sqrt(p3*p3+MK*MK));

			if(mdcKalTrk_K->charge() >0) {
				iKp.push_back(iGood[i]);
				pKp.push_back(ptrk);

				particle_type.push_back(3); // 3 represents K+
			} else {
				iKm.push_back(iGood[i]);
				pKm.push_back(ptrk);

				particle_type.push_back(4); // 4 represents K-
			}
		}
	}

//	cout << ipip.size() << ipim.size() << iKp.size() << iKm.size() << endl;
	if (!check_pid()) return false;
	return true;
}
bool OmegaKK::check_pid() {
	return (ipip.size() * ipim.size() * iKp.size() * iKm.size() == 1);
}


VertexParameter OmegaKK::get_vxpar() {
	VertexParameter tmp;
	HepPoint3D vx(0., 0., 0.);
	HepSymMatrix Evx(3, 0);
	double bx = 1E+6;
	double by = 1E+6;
	double bz = 1E+6;
	Evx[0][0] = bx * bx;
	Evx[1][1] = by * by;
	Evx[2][2] = bz * bz;

	tmp.setEvx(Evx);
	tmp.setVx(IP);  // fill the parameters obtained before!!!
	return tmp;
}


bool OmegaKK::vertex_fit() {
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack();
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack();
	RecMdcKalTrack *KpTrk = (*(evtRecTrkCol->begin() + iKp[0]))->mdcKalTrack();
	RecMdcKalTrack *KmTrk = (*(evtRecTrkCol->begin() + iKm[0]))->mdcKalTrack();

    WTrackParameter wvpipTrk, wvpimTrk;
    WTrackParameter wvKpTrk, wvKmTrk;

    if (m_5C_Calibration == 0) {
    
    	wvpipTrk = WTrackParameter(MPI, pipTrk->getZHelix(), pipTrk->getZError());
    	wvpimTrk = WTrackParameter(MPI, pimTrk->getZHelix(), pimTrk->getZError());
    
    	wvKpTrk = WTrackParameter(MK, KpTrk->getZHelix(), KpTrk->getZError());
    	wvKmTrk = WTrackParameter(MK, KmTrk->getZHelix(), KmTrk->getZError());
    } else {
        HepVector wpip_zHel(5,0);
        HepVector wpim_zHel(5,0);
        HepVector wkp_zHel(5,0);
        HepVector wkm_zHel(5,0);
    
        OmegaKK::calibration(pipTrk, wpip_zHel, 0);
        OmegaKK::calibration(pimTrk, wpim_zHel, 0);
        OmegaKK::calibration(KpTrk, wkp_zHel, 1);
        OmegaKK::calibration(KmTrk, wkm_zHel, 1);
    
    	wvpipTrk = WTrackParameter(MPI, wpip_zHel, pipTrk->getZError());
    	wvpimTrk = WTrackParameter(MPI, wpim_zHel, pimTrk->getZError());
    
    	wvKpTrk = WTrackParameter(MK, wkp_zHel, KpTrk->getZError());
    	wvKmTrk = WTrackParameter(MK, wkm_zHel, KmTrk->getZError());
        
    }
	VertexParameter vxpar = get_vxpar();

	VertexFit* vtxfit = VertexFit::instance();
	vtxfit->init();
	vtxfit->AddTrack(0,  wvpipTrk);
	vtxfit->AddTrack(1,  wvpimTrk);

	vtxfit->AddTrack(2,  wvKpTrk);
	vtxfit->AddTrack(3,  wvKmTrk);

	vtxfit->AddVertex(0, vxpar, 0, 1, 2, 3);
	if(!vtxfit->Fit(0)) return false;
	vtxfit->Swim(0);

	wpip = vtxfit->wtrk(0);
	wpim = vtxfit->wtrk(1);

	wKp = vtxfit->wtrk(2);
	wKm = vtxfit->wtrk(3);

	//KinematicFit * kmfit = KinematicFit::instance();
	kmfit = KalmanKinematicFit::instance();
	return true;
}
bool OmegaKK::fit4c() {
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
//	HepLorentzVector ecms(0.041, -0.001, 0.004, 3.686);
//	HepLorentzVector ecms(0.04, 0, 0, 3.65); // continue spectrum

	int ig1 = -1;
	int ig2 = -1;

    RecEmcShower *g1Trk = 0;
    RecEmcShower *g2Trk = 0;
    RecEmcShower *g3Trk = 0;

	double chisq = 9999.;
	for(int i = 0; i < m_nGam-1; i++) {
		g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
		for(int j = i+1; j < m_nGam; j++) {
			g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
			kmfit->init();
			kmfit->AddTrack(0, wpip);
			kmfit->AddTrack(1, wpim);
			kmfit->AddTrack(2, wKp);
			kmfit->AddTrack(3, wKm);
			kmfit->AddTrack(4, 0.0, g1Trk);
			kmfit->AddTrack(5, 0.0, g2Trk);
			kmfit->AddFourMomentum(0, ecms);
			bool oksq = kmfit->Fit();
			if(oksq) {
				double chi2 = kmfit->chisq();
				if(chi2 < chisq) {
					chisq = chi2;
					ig1 = iGam[i];
					ig2 = iGam[j];
					m_g1 = i;
					m_g2 = j;
				}
			}
		}
	}

	if(chisq > 200) return false;

    double chisq_3g = 9999.;
    if (m_nGam >= 3) {
        for (int i = 0; i < m_nGam - 2; ++i) {
            g1Trk = (*(evtRecTrkCol->begin() + iGam[i]))->emcShower();
            for (int j = i + 1; j < m_nGam - 1; ++j) {
                g2Trk = (*(evtRecTrkCol->begin() + iGam[j]))->emcShower();
                for (int k = j + 1; k < m_nGam; ++k) {
                    cout << "i = " << i << "; j = " << j << "; k = " << k << "; m_nGam = " << m_nGam << endl;
                    g3Trk = (*(evtRecTrkCol->begin() + iGam[k]))->emcShower();
                    kmfit->init();
                    kmfit->AddTrack(0, wpip);
                    kmfit->AddTrack(1, wpim);
                    kmfit->AddTrack(2, wKp);
                    kmfit->AddTrack(3, wKm);
                    kmfit->AddTrack(4, 0.0, g1Trk);
                    kmfit->AddTrack(5, 0.0, g2Trk);
                    kmfit->AddTrack(6, 0.0, g3Trk);
                    kmfit->AddFourMomentum(0, ecms);
                    bool oksq = kmfit->Fit();   
                    if (oksq) {
                        double chisq = kmfit->chisq();
                        if (chisq < chisq_3g) chisq_3g = chisq;                   
                    }           
                }
            }
        }
    }

    double chisq_1g = 9999.;
    for (int i = 0; i < m_nGam; ++i) {
        g1Trk = (*(evtRecTrkCol->begin() + iGam[i]))->emcShower();
        kmfit->init();
        kmfit->AddTrack(0, wpip);
        kmfit->AddTrack(1, wpim);
        kmfit->AddTrack(2, wKp);
        kmfit->AddTrack(3, wKm);
        kmfit->AddTrack(4, 0.0, g1Trk);
        kmfit->AddFourMomentum(0, ecms);
        bool oksq = kmfit->Fit();
        if (oksq) {
            double chisq = kmfit->chisq();
            if (chisq < chisq_1g) chisq_1g = chisq;
        }
    }


	g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
	g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();
	kmfit->init();
	kmfit->AddTrack(0, wpip);
	kmfit->AddTrack(1, wpim);
	kmfit->AddTrack(2, wKp);
	kmfit->AddTrack(3, wKm);
	kmfit->AddTrack(4, 0.0, g1Trk);
	kmfit->AddTrack(5, 0.0, g2Trk);
	kmfit->AddFourMomentum(0, ecms);
	if (!kmfit->Fit()) return false;
	HepLorentzVector ppi0 = kmfit->pfit(4) + kmfit->pfit(5);
	m_chi1 = kmfit->chisq();
	m_mpi0_4c = ppi0.m();
    m_chisq_1g = chisq_1g;
    m_chisq_3g = chisq_3g;
	return true;
}
bool OmegaKK::fit5c() {
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
//	HepLorentzVector ecms(0.04, -0.001, 0.004, 3.686);
//	HepLorentzVector ecms(0.04, 0, 0, 3.65); // continue spectrum
	double chisq = 9999.;
	int ig1 = -1;
	int ig2 = -1;
	for(int i = 0; i < m_nGam-1; i++) {
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
		for(int j = i+1; j < m_nGam; j++) {
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
			kmfit->init();
			kmfit->AddTrack(0, wpip);
			kmfit->AddTrack(1, wpim);
			kmfit->AddTrack(2, wKp);
			kmfit->AddTrack(3, wKm);
			kmfit->AddTrack(4, 0.0, g1Trk);
			kmfit->AddTrack(5, 0.0, g2Trk);
			kmfit->AddResonance(0, 0.135, 4, 5);
			kmfit->AddFourMomentum(1, ecms);
			if(!kmfit->Fit(0)) continue;
			if(!kmfit->Fit(1)) continue;
			if(kmfit->Fit()) {
				double chi2 = kmfit->chisq();
				if(chi2 < chisq) {
					chisq = chi2;
					ig1 = iGam[i];
					ig2 = iGam[j];
					m_g1 = i;
					m_g2 = j;
				}
			}
		}
	}

	if (chisq > 200) return false;
	RecEmcShower* g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
	RecEmcShower* g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();

	kmfit->init();
	kmfit->AddTrack(0, wpip);
	kmfit->AddTrack(1, wpim);
	kmfit->AddTrack(2, wKp);
	kmfit->AddTrack(3, wKm);
	kmfit->AddTrack(4, 0.0, g1Trk);
	kmfit->AddTrack(5, 0.0, g2Trk);
	kmfit->AddResonance(0, 0.135, 4, 5);
	kmfit->AddFourMomentum(1, ecms);
	//		kmfit->Fit(0);
	//		kmfit->Fit(1);
	if (!kmfit->Fit()) return false;
	HepLorentzVector ppi0 = kmfit->pfit(4) + kmfit->pfit(5);
	HepLorentzVector pomega = ppi0 + kmfit->pfit(0) + kmfit->pfit(1);
	HepLorentzVector pkaonp = kmfit->pfit(2);
	HepLorentzVector pkaonm = kmfit->pfit(3);
	HepLorentzVector pkaonmpip = pkaonm + kmfit->pfit(0);
	HepLorentzVector pkaonppim = pkaonp + kmfit->pfit(1);

	HepLorentzVector pkaon = pkaonp + pkaonm;
	HepLorentzVector ptot = pomega + pkaon;

	m_chi2  = kmfit->chisq();
	m_momega = pomega.m();
	m_mkaonp = pkaonp.m();
	m_mkaonm = pkaonm.m();
	m_mkaonmpip = pkaonmpip.m();
	m_mkaonppim = pkaonppim.m();
	m_mpi0_5c = ppi0.m();
	m_mkaon = pkaon.m();
	m_mtot = ptot.m();

	double eg1 = (kmfit->pfit(4)).e();
	double eg2 = (kmfit->pfit(5)).e();
	double fcos = abs(eg1-eg2)/ppi0.rho();
	m_fcos = fcos;

    all_particles_4mom.push_back(kmfit->pfit(0));
    all_particles_4mom.push_back(kmfit->pfit(1));
    all_particles_4mom.push_back(kmfit->pfit(2));
    all_particles_4mom.push_back(kmfit->pfit(3));
    all_particles_4mom.push_back(kmfit->pfit(4));
    all_particles_4mom.push_back(kmfit->pfit(5));

	return true;
}

void OmegaKK::write_ana_data() {
//	HepLorentzVector ecms(0.04, -0.001, 0.004, 3.686);
//	HepLorentzVector ecms(0.04, 0, 0, 3.65); // continue spectrum
	Hep3Vector m_boost_vec(-ecms.boostVector());

	HepLorentzVector pTot;
//	for(int i = 0; i < m_nGam - 1; i++){
//		for(int j = i+1; j < m_nGam; j++) {
			HepLorentzVector p2g = pGam[m_g1] + pGam[m_g2];
			pTot = ppip[0] + ppim[0] + pKp[0] + pKm[0];
			pTot += p2g;
			m_m2gg = p2g.m();
			m_etot = pTot.e();
//			m_tuple3 -> write();
//		}
//	}
	HepLorentzVector boosted_four_mom;
    m_particle_index = 0;
    
    for(int i = 0; i < 6; ++i) {
        boosted_four_mom = boostOf(all_particles_4mom[i], m_boost_vec);
        m4xyz[i][0] = boosted_four_mom.e();
        m4xyz[i][1] = boosted_four_mom.vect().x();
        m4xyz[i][2] = boosted_four_mom.vect().y();
        m4xyz[i][3] = boosted_four_mom.vect().z();    

//	    m4rft[i][0] = boosted_four_mom.e();
//	    m4rft[i][1] = boosted_four_mom.vect().mag();	
//	    m4rft[i][2] = boosted_four_mom.vect().phi();			
//	    m4rft[i][3] = boosted_four_mom.vect().theta();

        ++m_particle_index;
    }

//	m_particle_index = 0;
//	m_charged_index = 0;
//	for (int i = 0; i < nGood; ++i) {
//		for (int j = 0; j < nGood; ++j) {
//			if (particle_type[j] == i + 1) {
//				boosted_four_mom = boostOf(particle_four_mom[j], m_boost_vec);	
//				m4xyz[m_particle_index][0] = boosted_four_mom.e();
//				m4xyz[m_particle_index][1] = boosted_four_mom.vect().x();	
//				m4xyz[m_particle_index][2] = boosted_four_mom.vect().y();			
//				m4xyz[m_particle_index][3] = boosted_four_mom.vect().z();
//				m4rft[m_particle_index][0] = boosted_four_mom.e();
//				m4rft[m_particle_index][1] = boosted_four_mom.vect().mag();	
//				m4rft[m_particle_index][2] = boosted_four_mom.vect().phi();			
//				m4rft[m_particle_index][3] = boosted_four_mom.vect().theta();
//				++m_particle_index;
//				
////				m_ptrk_pid[m_charged_index] = ptrk_pid[j];
////				m_cost_pid[m_charged_index] = cost_pid[j];
////				m_dedx_pid[m_charged_index] = dedx_pid[j];
////				m_tof1_pid[m_charged_index] = tof1_pid[j];
////				m_tof2_pid[m_charged_index] = tof2_pid[j];
////				m_prob_pid[m_charged_index] = prob_pid[j];
////				++m_charged_index;
//			}
//		}
//	}
//
//	boosted_four_mom = boostOf(m_Gam[m_g1], m_boost_vec);	
//	m4xyz[nGood][0] = boosted_four_mom.e();
//	m4xyz[nGood][1] = boosted_four_mom.vect().x();
//	m4xyz[nGood][2] = boosted_four_mom.vect().y();
//	m4xyz[nGood][3] = boosted_four_mom.vect().z();
//	m4rft[nGood][0] = boosted_four_mom.e();
//	m4rft[nGood][1] = boosted_four_mom.vect().mag();
//	m4rft[nGood][2] = boosted_four_mom.vect().phi();
//	m4rft[nGood][3] = boosted_four_mom.vect().theta();
//	++m_particle_index;
//	boosted_four_mom = boostOf(m_Gam[m_g2], m_boost_vec);	
//	m4xyz[nGood + 1][0] = boosted_four_mom.e();
//	m4xyz[nGood + 1][1] = boosted_four_mom.vect().x();
//	m4xyz[nGood + 1][2] = boosted_four_mom.vect().y();
//	m4xyz[nGood + 1][3] = boosted_four_mom.vect().z();
//	m4rft[nGood + 1][0] = boosted_four_mom.e();
//	m4rft[nGood + 1][1] = boosted_four_mom.vect().mag();
//	m4rft[nGood + 1][2] = boosted_four_mom.vect().phi();
//	m4rft[nGood + 1][3] = boosted_four_mom.vect().theta();
//	++m_particle_index;
//
	m_gamma_index = 0;
	m_dang[0] = s_dang[m_g1];
	m_eraw[0] = s_eraw[m_g1];
	++m_gamma_index;
	m_dang[1] = s_dang[m_g2];
	m_eraw[1] = s_eraw[m_g2];
	++m_gamma_index;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode OmegaKK::execute() {
	MsgStream log(msgSvc(), name());
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
	m_runNo = eventHeader->runNumber();
	m_event = eventHeader->eventNumber();
    m_runnr = m_runNo;
    m_recnr = m_event;
	Ncut0++;

	before_execute();
	if (m_runNo < 0) {
		if (!topology()) return StatusCode::FAILURE;
	}
	if (!find_charged_track()) return StatusCode::SUCCESS;
	Ncut1++;
	if (!find_gamma_track_and_assign_4mom()) return StatusCode::SUCCESS;
	Ncut2++;
	if (!pid_and_assign_4mom()) return StatusCode::SUCCESS;
	Ncut3++;
	if (!vertex_fit()) return StatusCode::SUCCESS;

	if (m_test4C == 1) {
		if (!fit4c()) return StatusCode::SUCCESS;
	}
	Ncut4++;

	if (m_test5C == 1) {
		if (!fit5c()) return StatusCode::SUCCESS;	
	}
	Ncut5++;
	write_ana_data();
	m_tuple_ana->write();
	return StatusCode::SUCCESS;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode OmegaKK::finalize() {
	cout<<"total number:         "<<Ncut0<<endl;
	cout<<"nGood==4, nCharge==0: "<<Ncut1<<endl;
	cout<<"nGam>=2:              "<<Ncut2<<endl;
	cout<<"Pass Pid:             "<<Ncut3<<endl;
	cout<<"Pass 4C:              "<<Ncut4<<endl;
	cout<<"Pass 5C:              "<<Ncut5<<endl;
//	cout<<"J/psi->rho0 pi0:      "<<Ncut6<<endl;
	cout << "m_collider_energy = " << m_collider_energy << endl;
//	MsgStream log(msgSvc(), name());
//	log << MSG::INFO << "in finalize()" << endmsg;
		

	return StatusCode::SUCCESS;
}

//*************************************************************************
// *****************************************************************
// ** A macro to create correlated Gaussian-distributed variables **
// *****************************************************************

// corset(): sets up the generation by calculating C from V.
void OmegaKK::corset(HepSymMatrix &V, HepMatrix &C, int n)
{
//cout<<"v="<<V<<endl;
//cout<<"c="<<C<<endl;
  double sum;

  // Compute square root of matrix sigma
  for (int j=0; j<n; j++) {
    sum = 0;
    for (int k=0; k<j; k++) {
      sum = sum + C[j][k]*C[j][k];
      
    }
//cout<<"sum="<<sum<<endl;
//cout<<"v("<<j<<","<<j<<")="<<V[j][j]<<endl;
    C[j][j] = sqrt(abs(V[j][j] - sum));
//cout<<"c("<<j<<","<<j<<")="<<C[j][j]<<endl;
    // Off Diagonal terms
    for (int i=j+1; i<n; i++) {
      sum = 0;
      for (int k=0; k<j; k++) {
	sum = sum + C[i][k]*C[j][k];
      }
      C[i][j] = (V[i][j] - sum)/C[j][j];
//cout<<"C("<<i<<","<<j<<")="<<C[i][j]<<endl;
    }
  }
}

// corgen(): generates a set of n random numbers Gaussian-distributed with covariance
// matrix V (V = C*C') and mean values zero.
void OmegaKK::corgen(HepMatrix &C, HepVector &x, int n)
{
  int i,j;
  int nmax = 100;
  
  if (n > nmax ) {
    printf("Error in corgen: array overflown");
  }
  
  double tmp[3];
  for(int p = 0 ; p < n; p ++){
      tmp[p] = gRandom->Gaus(0,1);
//cout<<"tmp["<<p<<"]="<<tmp[p]<<endl;
  }
  for ( i=0; i<n; i++) {
    x[i] = 0.0;
    for (j=0; j<=i; j++) {
      x[i] = x[i]+C[i][j]*tmp[j];
    }
  }
}


void OmegaKK::calibration(RecMdcKalTrack *trk , HepVector &wtrk_zHel, int n )
{

   HepVector pip_calerr_d2(5,0);
   HepVector pim_calerr_d2(5,0);
   HepVector kp_calerr_d2(5,0);
   HepVector km_calerr_d2(5,0);

//   pip_calerr_d2[0] = 1.0;
//   pip_calerr_d2[1] = 1.189*1.050;
//   pip_calerr_d2[2] = 1.235*1.060;
//   pip_calerr_d2[3] = 1.0;
//   pip_calerr_d2[4] = 1.189*1.048;
//
//   pim_calerr_d2[0] = 1.0;
//   pim_calerr_d2[1] = 1.187*1.040;
//   pim_calerr_d2[2] = 1.194*1.061;
//   pim_calerr_d2[3] = 1.0;
//   pim_calerr_d2[4] = 1.168*1.036;
//
//   kp_calerr_d2[0] = 1.0;
//   kp_calerr_d2[1] = 1.164*1.019;
//   kp_calerr_d2[2] = 1.216*1.056;
//   kp_calerr_d2[3] = 1.0;
//   kp_calerr_d2[4] = 1.188*1.053;
//
//   km_calerr_d2[0] = 1.0;
//   km_calerr_d2[1] = 1.177*1.024;
//   km_calerr_d2[2] = 1.191*1.050;
//   km_calerr_d2[3] = 1.0;
//   km_calerr_d2[4] = 1.172*1.032;


// My first 5C calibration
//    pip_calerr_d2[0] = 1.0;
//    pip_calerr_d2[1] = 1.05*1.188;
//    pip_calerr_d2[2] = 1.042*1.195;
//    pip_calerr_d2[3] = 1.0;
//    pip_calerr_d2[4] = 1.0128*1.134;
//    
//    pim_calerr_d2[0] = 1.0;
//    pim_calerr_d2[1] = 1.038*1.173;
//    pim_calerr_d2[2] = 1.0364*1.182;
//    pim_calerr_d2[3] = 1.0;
//    pim_calerr_d2[4] = 1.004*1.129;
//    
//    kp_calerr_d2[0] = 1.0;
//    kp_calerr_d2[1] = 1.0094*1.165;
//    kp_calerr_d2[2] = 1.055*1.205;
//    kp_calerr_d2[3] = 1.0;
//    kp_calerr_d2[4] = 1.008*1.132;
//    
//    km_calerr_d2[0] = 1.0;
//    km_calerr_d2[1] = 1.02*1.167;
//    km_calerr_d2[2] = 1.044*1.188;
//    km_calerr_d2[3] = 1.0;
//    km_calerr_d2[4] = 0.9969*1.115;

//My squared 5C calibration
    pip_calerr_d2[0] = 1.0;
    pip_calerr_d2[1] = 1.05*1.188*1.05*1.188;
    pip_calerr_d2[2] = 1.042*1.195*1.042*1.195;
    pip_calerr_d2[3] = 1.0;
    pip_calerr_d2[4] = 1.0128*1.134*1.0128*1.134;
    
    pim_calerr_d2[0] = 1.0;
    pim_calerr_d2[1] = 1.038*1.173*1.038*1.173;
    pim_calerr_d2[2] = 1.0364*1.182*1.0364*1.182;
    pim_calerr_d2[3] = 1.0;
    pim_calerr_d2[4] = 1.004*1.129*1.004*1.129;
    
    kp_calerr_d2[0] = 1.0;
    kp_calerr_d2[1] = 1.0094*1.165*1.0094*1.165;
    kp_calerr_d2[2] = 1.055*1.205*1.055*1.205;
    kp_calerr_d2[3] = 1.0;
    kp_calerr_d2[4] = 1.008*1.132*1.008*1.132;
    
    km_calerr_d2[0] = 1.0;
    km_calerr_d2[1] = 1.02*1.167*1.02*1.167;
    km_calerr_d2[2] = 1.044*1.188*1.044*1.188;
    km_calerr_d2[3] = 1.0;
    km_calerr_d2[4] = 0.9969*1.115*0.9969*1.115;

   HepVector pip_calmean_d2(5,0);
   HepVector pim_calmean_d2(5,0);
   HepVector kp_calmean_d2(5,0);
   HepVector km_calmean_d2(5,0);
/*
   pip_calmean_d2[0] = 0;
   pip_calmean_d2[1] = -0.106-0.059;
   pip_calmean_d2[2] = 0.418-0.103;
   pip_calmean_d2[3] = 0;
   pip_calmean_d2[4] = 0.447/2.;
   pim_calmean_d2[0] = 0;
   pim_calmean_d2[1] = 0.160-0.016;
   pim_calmean_d2[2] = -0.391+0.104;
   pim_calmean_d2[3] = 0;
   pim_calmean_d2[4] = 0.490/2.;
   kp_calmean_d2[0] = 0;
   kp_calmean_d2[1] = 0.038-0.040;
   kp_calmean_d2[2] = 0.478-0.244;
   kp_calmean_d2[3] = 0;
   kp_calmean_d2[4] = 0.493/2.;
   km_calmean_d2[0] = 0;
   km_calmean_d2[1] = -0.134+0.057;
   km_calmean_d2[2] = -0.497+0.245;
   km_calmean_d2[3] = 0;
   km_calmean_d2[4] = 0.540/2.;
*/

   pip_calmean_d2[0] = 0;
   pip_calmean_d2[1] = 0;
   pip_calmean_d2[2] = 0;
   pip_calmean_d2[3] = 0;
   pip_calmean_d2[4] = 0;
   pim_calmean_d2[0] = 0;
   pim_calmean_d2[1] = 0;
   pim_calmean_d2[2] = 0;
   pim_calmean_d2[3] = 0;
   pim_calmean_d2[4] = 0;
   kp_calmean_d2[0] = 0;
   kp_calmean_d2[1] = 0;
   kp_calmean_d2[2] = 0;
   kp_calmean_d2[3] = 0;
   kp_calmean_d2[4] = 0;
   km_calmean_d2[0] = 0;
   km_calmean_d2[1] = 0;
   km_calmean_d2[2] = 0;
   km_calmean_d2[3] = 0;
   km_calmean_d2[4] = 0;

//cout<<"pip_calerr_d2="<<pip_calerr_d2<<endl;
//cout<<"pim_calerr_d2="<<pim_calerr_d2<<endl;
//cout<<"kp_calerr_d2="<<kp_calerr_d2<<endl;
//cout<<"km_calerr_d2="<<km_calerr_d2<<endl;

//cout<<"pip_calmean_d2="<<pip_calmean_d2<<endl;
//cout<<"pim_calmean_d2="<<pim_calmean_d2<<endl;
//cout<<"kp_calmean_d2="<<kp_calmean_d2<<endl;
//cout<<"km_calmean_d2="<<km_calmean_d2<<endl;

  if(trk->charge()>0 && n==0){
//pip
     HepSymMatrix wpip_zerr(5,0);
     wpip_zerr = trk->getZError();
//cout<<"wpip_zerr="<<wpip_zerr<<endl;
     HepSymMatrix wpip_zcal(3,0);
/*
     wpip_zcal[0][0] = (pip_calerr_d2[1]*pip_calerr_d2[1]-1)*wpip_zerr[1][1];
     wpip_zcal[0][1] = sqrt(pip_calerr_d2[1]*pip_calerr_d2[1]-1)*sqrt(pip_calerr_d2[2]*pip_calerr_d2[2]-1)*wpip_zerr[1][2];
     wpip_zcal[0][2] = sqrt(pip_calerr_d2[1]*pip_calerr_d2[1]-1)*sqrt(pip_calerr_d2[4]*pip_calerr_d2[4]-1)*wpip_zerr[1][4];
     wpip_zcal[1][1] = (pip_calerr_d2[2]*pip_calerr_d2[2]-1)*wpip_zerr[2][2];
     wpip_zcal[1][2] = sqrt(pip_calerr_d2[2]*pip_calerr_d2[2]-1)*sqrt(pip_calerr_d2[4]*pip_calerr_d2[4]-1)*wpip_zerr[2][4];
     wpip_zcal[2][2] = (pip_calerr_d2[4]*pip_calerr_d2[4]-1)*wpip_zerr[4][4];
*/

     wpip_zcal[0][0] = (pip_calerr_d2[1]*pip_calerr_d2[1]-1)*wpip_zerr[1][1];
     wpip_zcal[1][1] = (pip_calerr_d2[2]*pip_calerr_d2[2]-1)*wpip_zerr[2][2];
     wpip_zcal[2][2] = (pip_calerr_d2[4]*pip_calerr_d2[4]-1)*wpip_zerr[4][4];
/*
     wpip_zcal[0][0] = (pip_calerr_d2[1]*pip_calerr_d2[1]-1)*wpip_zerr[1][1];
     wpip_zcal[0][1] = sqrt(pip_calerr_d2[1]*pip_calerr_d2[1]-1)*sqrt(pip_calerr_d2[2]*pip_calerr_d2[2]-1)*sqrt(wpip_zerr[1][1]*wpip_zerr[2][2]);
     wpip_zcal[0][2] = sqrt(pip_calerr_d2[1]*pip_calerr_d2[1]-1)*sqrt(pip_calerr_d2[4]*pip_calerr_d2[4]-1)*sqrt(wpip_zerr[1][1]*wpip_zerr[4][4]);
     wpip_zcal[1][1] = (pip_calerr_d2[2]*pip_calerr_d2[2]-1)*wpip_zerr[2][2];
     wpip_zcal[1][2] = sqrt(pip_calerr_d2[2]*pip_calerr_d2[2]-1)*sqrt(pip_calerr_d2[4]*pip_calerr_d2[4]-1)*sqrt(wpip_zerr[2][2]*wpip_zerr[4][4]);
     wpip_zcal[2][2] = (pip_calerr_d2[4]*pip_calerr_d2[4]-1)*wpip_zerr[4][4];
*/
//cout<<"wpip_zcal="<<wpip_zcal<<endl;

     HepMatrix wpip_zerrc(3,3,0);
     OmegaKK::corset(wpip_zcal,wpip_zerrc,3);
     HepVector wpip_zgen(3,0);
     OmegaKK::corgen(wpip_zerrc,wpip_zgen,3);
     
     wtrk_zHel[0] = trk->getZHelix()[0];
     wtrk_zHel[1] = trk->getZHelix()[1]+pip_calmean_d2[1]*sqrt(wpip_zerr[1][1])+wpip_zgen[0];
     wtrk_zHel[2] = trk->getZHelix()[2]+pip_calmean_d2[2]*sqrt(wpip_zerr[2][2])+wpip_zgen[1];
     wtrk_zHel[3] = trk->getZHelix()[3];
     wtrk_zHel[4] = trk->getZHelix()[4]+pip_calmean_d2[4]*sqrt(wpip_zerr[4][4])+wpip_zgen[2];

//cout<<"wtrk_zHel="<<wtrk_zHel<<endl;
  }
  if(trk->charge()<0 && n==0)
  {
  //pim
     HepSymMatrix wpim_zerr(5,0);
     wpim_zerr = trk->getZError();

     HepSymMatrix wpim_zcal(3,0);
/*
     wpim_zcal[0][0] = (pim_calerr_d2[1]*pim_calerr_d2[1]-1)*wpim_zerr[1][1];
     wpim_zcal[0][1] = sqrt(pim_calerr_d2[1]*pim_calerr_d2[1]-1)*sqrt(pim_calerr_d2[2]*pim_calerr_d2[2]-1)*wpim_zerr[1][2];
     wpim_zcal[0][1] = sqrt(pim_calerr_d2[1]*pim_calerr_d2[1]-1)*sqrt(pim_calerr_d2[4]*pim_calerr_d2[4]-1)*wpim_zerr[1][4];
     wpim_zcal[1][1] = (pim_calerr_d2[2]*pim_calerr_d2[2]-1)*wpim_zerr[2][2];
     wpim_zcal[1][2] = sqrt(pim_calerr_d2[2]*pim_calerr_d2[2]-1)*sqrt(pim_calerr_d2[4]*pim_calerr_d2[4]-1)*wpim_zerr[2][4];
     wpim_zcal[2][2] = (pim_calerr_d2[4]*pim_calerr_d2[4]-1)*wpim_zerr[4][4];
*/

     wpim_zcal[0][0] = (pim_calerr_d2[1]*pim_calerr_d2[1]-1)*wpim_zerr[1][1];
     wpim_zcal[1][1] = (pim_calerr_d2[2]*pim_calerr_d2[2]-1)*wpim_zerr[2][2];
     wpim_zcal[2][2] = (pim_calerr_d2[4]*pim_calerr_d2[4]-1)*wpim_zerr[4][4];

/*
     wpim_zcal[0][0] = (pim_calerr_d2[1]*pim_calerr_d2[1]-1)*wpim_zerr[1][1];
     wpim_zcal[0][1] = sqrt(pim_calerr_d2[1]*pim_calerr_d2[1]-1)*sqrt(pim_calerr_d2[2]*pim_calerr_d2[2]-1)*sqrt(wpim_zerr[1][1]*wpim_zerr[2][2]);
     wpim_zcal[0][1] = sqrt(pim_calerr_d2[1]*pim_calerr_d2[1]-1)*sqrt(pim_calerr_d2[4]*pim_calerr_d2[4]-1)*sqrt(wpim_zerr[1][1]*wpim_zerr[4][4]);
     wpim_zcal[1][1] = (pim_calerr_d2[2]*pim_calerr_d2[2]-1)*wpim_zerr[2][2];
     wpim_zcal[1][2] = sqrt(pim_calerr_d2[2]*pim_calerr_d2[2]-1)*sqrt(pim_calerr_d2[4]*pim_calerr_d2[4]-1)*sqrt(wpim_zerr[2][2]*wpim_zerr[4][4]);
     wpim_zcal[2][2] = (pim_calerr_d2[4]*pim_calerr_d2[4]-1)*wpim_zerr[4][4];
*/
     HepMatrix wpim_zerrc(3,3,0);
     OmegaKK::corset(wpim_zcal,wpim_zerrc,3);
     HepVector wpim_zgen(3,0);
     OmegaKK::corgen(wpim_zerrc,wpim_zgen,3);


     wtrk_zHel[0] = trk->getZHelix()[0];
     wtrk_zHel[1] = trk->getZHelix()[1]+pim_calmean_d2[1]*sqrt(wpim_zerr[1][1])+wpim_zgen[0];
     wtrk_zHel[2] = trk->getZHelix()[2]+pim_calmean_d2[2]*sqrt(wpim_zerr[2][2])+wpim_zgen[1];
     wtrk_zHel[3] = trk->getZHelix()[3];
     wtrk_zHel[4] = trk->getZHelix()[4]+pim_calmean_d2[4]*sqrt(wpim_zerr[4][4])+wpim_zgen[2];

  }
  if(trk->charge()>0 && n==1)
  {
//kp
     HepSymMatrix wkp_zerr(5,0);
     wkp_zerr = trk->getZErrorK();

     HepSymMatrix wkp_zcal(3,0);
/*
     wkp_zcal[0][0] = (kp_calerr_d2[1]*kp_calerr_d2[1]-1)*wkp_zerr[1][1];
     wkp_zcal[0][1] = sqrt(kp_calerr_d2[1]*kp_calerr_d2[1]-1)*sqrt(kp_calerr_d2[2]*kp_calerr_d2[2]-1)*wkp_zerr[1][2];
     wkp_zcal[0][2] = sqrt(kp_calerr_d2[1]*kp_calerr_d2[1]-1)*sqrt(kp_calerr_d2[4]*kp_calerr_d2[4]-1)*wkp_zerr[1][4];
     wkp_zcal[1][1] = (kp_calerr_d2[2]*kp_calerr_d2[2]-1)*wkp_zerr[2][2];
     wkp_zcal[1][2] = sqrt(kp_calerr_d2[2]*kp_calerr_d2[2]-1)*sqrt(kp_calerr_d2[4]*kp_calerr_d2[4]-1)*wkp_zerr[2][4];
     wkp_zcal[2][2] = (kp_calerr_d2[4]*kp_calerr_d2[4]-1)*wkp_zerr[4][4];
*/

     wkp_zcal[0][0] = (kp_calerr_d2[1]*kp_calerr_d2[1]-1)*wkp_zerr[1][1];
     wkp_zcal[1][1] = (kp_calerr_d2[2]*kp_calerr_d2[2]-1)*wkp_zerr[2][2];
     wkp_zcal[2][2] = (kp_calerr_d2[4]*kp_calerr_d2[4]-1)*wkp_zerr[4][4];
/*
     wkp_zcal[0][0] = (kp_calerr_d2[1]*kp_calerr_d2[1]-1)*wkp_zerr[1][1];
     wkp_zcal[0][1] = sqrt(kp_calerr_d2[1]*kp_calerr_d2[1]-1)*sqrt(kp_calerr_d2[2]*kp_calerr_d2[2]-1)*sqrt(wkp_zerr[1][1]*wkp_zerr[2][2]);
     wkp_zcal[0][2] = sqrt(kp_calerr_d2[1]*kp_calerr_d2[1]-1)*sqrt(kp_calerr_d2[4]*kp_calerr_d2[4]-1)*sqrt(wkp_zerr[1][1]*wkp_zerr[4][4]);
     wkp_zcal[1][1] = (kp_calerr_d2[2]*kp_calerr_d2[2]-1)*wkp_zerr[2][2];
     wkp_zcal[1][2] = sqrt(kp_calerr_d2[2]*kp_calerr_d2[2]-1)*sqrt(kp_calerr_d2[4]*kp_calerr_d2[4]-1)*sqrt(wkp_zerr[2][2]*wkp_zerr[4][4]);
     wkp_zcal[2][2] = (kp_calerr_d2[4]*kp_calerr_d2[4]-1)*wkp_zerr[4][4];
*/
     HepMatrix wkp_zerrc(3,3,0);
     OmegaKK::corset(wkp_zcal,wkp_zerrc,3);
     HepVector wkp_zgen(3,0);
     OmegaKK::corgen(wkp_zerrc,wkp_zgen,3);

     wtrk_zHel[0] = trk->getZHelixK()[0];
     wtrk_zHel[1] = trk->getZHelixK()[1]+kp_calmean_d2[1]*sqrt(wkp_zerr[1][1])+wkp_zgen[0];
     wtrk_zHel[2] = trk->getZHelixK()[2]+kp_calmean_d2[2]*sqrt(wkp_zerr[2][2])+wkp_zgen[1];
     wtrk_zHel[3] = trk->getZHelixK()[3];
     wtrk_zHel[4] = trk->getZHelixK()[4]+kp_calmean_d2[4]*sqrt(wkp_zerr[4][4])+wkp_zgen[2];

  }

    if(trk->charge()<0 && n==1)
  {
//km
     HepSymMatrix wkm_zerr(5,0);
     wkm_zerr = trk->getZErrorK();

     HepSymMatrix wkm_zcal(3,0);
/*
     wkm_zcal[0][0] = (km_calerr_d2[1]*km_calerr_d2[1]-1)*wkm_zerr[1][1];
     wkm_zcal[0][1] = sqrt(km_calerr_d2[1]*km_calerr_d2[1]-1)*sqrt(km_calerr_d2[2]*km_calerr_d2[2]-1)*wkm_zerr[1][2];
     wkm_zcal[0][2] = sqrt(km_calerr_d2[1]*km_calerr_d2[1]-1)*sqrt(km_calerr_d2[4]*km_calerr_d2[4]-1)*wkm_zerr[1][4];
     wkm_zcal[1][1] = (km_calerr_d2[2]*km_calerr_d2[2]-1)*wkm_zerr[2][2];
     wkm_zcal[1][2] = sqrt(km_calerr_d2[2]*km_calerr_d2[2]-1)*sqrt(km_calerr_d2[4]*km_calerr_d2[4]-1)*wkm_zerr[2][4];
     wkm_zcal[2][2] = (km_calerr_d2[4]*km_calerr_d2[4]-1)*wkm_zerr[4][4];
*/

     wkm_zcal[0][0] = (km_calerr_d2[1]*km_calerr_d2[1]-1)*wkm_zerr[1][1];
     wkm_zcal[1][1] = (km_calerr_d2[2]*km_calerr_d2[2]-1)*wkm_zerr[2][2];
     wkm_zcal[2][2] = (km_calerr_d2[4]*km_calerr_d2[4]-1)*wkm_zerr[4][4];
/*
     wkm_zcal[0][0] = (km_calerr_d2[1]*km_calerr_d2[1]-1)*wkm_zerr[1][1];
     wkm_zcal[0][1] = sqrt(km_calerr_d2[1]*km_calerr_d2[1]-1)*sqrt(km_calerr_d2[2]*km_calerr_d2[2]-1)*sqrt(wkm_zerr[1][1]*wkm_zerr[2][2]);
     wkm_zcal[0][2] = sqrt(km_calerr_d2[1]*km_calerr_d2[1]-1)*sqrt(km_calerr_d2[4]*km_calerr_d2[4]-1)*sqrt(wkm_zerr[1][1]*wkm_zerr[4][4]);
     wkm_zcal[1][1] = (km_calerr_d2[2]*km_calerr_d2[2]-1)*wkm_zerr[2][2];
     wkm_zcal[1][2] = sqrt(km_calerr_d2[2]*km_calerr_d2[2]-1)*sqrt(km_calerr_d2[4]*km_calerr_d2[4]-1)*sqrt(wkm_zerr[2][2]*wkm_zerr[4][4]);
     wkm_zcal[2][2] = (km_calerr_d2[4]*km_calerr_d2[4]-1)*wkm_zerr[4][4];
*/
     HepMatrix wkm_zerrc(3,3,0);
     OmegaKK::corset(wkm_zcal,wkm_zerrc,3);
     HepVector wkm_zgen(3,0);
     OmegaKK::corgen(wkm_zerrc,wkm_zgen,3);

     wtrk_zHel[0] = trk->getZHelixK()[0];
     wtrk_zHel[1] = trk->getZHelixK()[1]+km_calmean_d2[1]*sqrt(wkm_zerr[1][1])+wkm_zgen[0];
     wtrk_zHel[2] = trk->getZHelixK()[2]+km_calmean_d2[2]*sqrt(wkm_zerr[2][2])+wkm_zgen[1];
     wtrk_zHel[3] = trk->getZHelixK()[3];
     wtrk_zHel[4] = trk->getZHelixK()[4]+km_calmean_d2[4]*sqrt(wkm_zerr[4][4])+wkm_zgen[2];

  }         


}

