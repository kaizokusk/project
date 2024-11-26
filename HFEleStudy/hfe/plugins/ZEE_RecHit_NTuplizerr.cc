// -*- C++ -*-
//NTuplelizer for working on miniAODs//////////////////////////////////////////
// Package:    Electron_GNN_Regression/ZEE_RecHit_NTuplizer
// Class:      ZEE_RecHit_NTuplizer
//
/*Skimmer to work on miniAODs
*/
//
// Original Author:  Rajdeep Mohan Chatterjee
//         Created:  Fri, 21 Feb 2020 11:38:58 GMT
//
// Modified by:  Avik Das
//



// system include files
#include <memory>
#include <iostream>
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TFile.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include <Math/Vector4D.h>

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
/*
#include "ZFinder/Event/interface/PDGID.h"  // PDGID enum (ELECTRON, POSITRON, etc.)
#include "ZFinder/Event/interface/TriggerList.h"  // ET_ET_TIGHT, ET_ET_DZ, ET_ET_LOOSE, ET_NT_ET_TIGHT, ET_HF_ET_TIGHT, ET_HF_ET_LOOSE, ET_HF_HF_TIGHT, ET_HF_HF_LOOSE, SINGLE_ELECTRON_TRIGGER, ALL_TRIGGERS
#include "ZFinder/Event/interface/PileupReweighting.h"  // RUN_2012_*_TRUE_PILEUP, SUMMER12_53X_MC_TRUE_PILEUP
*/

//
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace edm;
using namespace std;
using namespace reco;

using reco::TrackCollection;

double dR(const reco::Candidate* a1, const reco::GenParticle* a2){
	double eta1=a1->eta();
	double eta2=a2->eta();
	double phi1=a1->phi();
	double phi2=a2->phi();
	double dp=phi1-phi2;
	if(std::abs(dp)> M_PI ){
		dp=(2*M_PI )-std::abs(dp);// for crcting phi angle crction

	}

	
	double x = sqrt(pow((eta1-eta2),2)+pow((dp),2));
	return x;
}

double dR_gen_reco_jets(const reco::GenJet* gj1, const reco::PFJet* rj1){
        double eta1=gj1->eta();
        double eta2=rj1->eta();
        double phi1=gj1->phi();
        double phi2=rj1->phi();
	double dp=phi1-phi2;
	if(std::abs(dp)> M_PI ){
		dp=(2*M_PI )-std::abs(dp);// for crcting phi angle crction

	}
	double x = sqrt(pow((eta1-eta2),2)+pow((dp),2));
	return x;       
}
//////////Function Definition: GetDressedElectron///////////////////////
/*const reco::GenParticle* GetDressedElectron( const reco::GenParticle* const BORN_ELE, const reco::GenParticle* const NAKED_ELE, 
const double MX_DELTA_R){
        ///Defining Electrons mass in GeV
        const double ELE_MASS=5.11*pow(10,-4);
        ///Defining a four vector in PtEtaPhiM vector notation. Note headerfile <Math/Vector4D.h> 
        //has to be included
        ROOT::Math::PtEtaPhiMVector DRESSED_ELE_4VEC(NAKED_ELE->pt(), NAKED_ELE->eta(), NAKED_ELE->phi(),ELE_MASS);
        //Asign a pointer to the Born electron.
        const reco::GenParticle* temp_ele = BORN_ELE;
        //Making sure we have don't set the Born electron as the Naked electron
        while(temp_ele!=NAKED_ELE){
                //Returning null pointer if it does not undergo any decay (interaction with environment) at all
                if(temp_ele->numberOfDaughters()==0){
                        return nullptr;
                }
                //If we are not so unfortunate lets go daughter by daughter
                const reco::GenParticle* swap_electron = nullptr;
                //Itertating over the daughters of the BORN Electron
                for (size_t i = 0; i< temp_ele->numberOfDaughters(); ++i){
                //Set the daughters pointer and keep adding the photon four vector to the naked electrons four vector to dress it.      
                        const reco::Candidate* daughter_particle=temp_ele->daughter(i);
                        if (fabs(daughter_particle->pdgId())==11){
                                swap_electron = dynamic_cast<const reco::GenParticle*>(daughter_particle);
                        }
                        else if (fabs(daughter_particle->pdgId())==22){
                                //Checking if deltaR lies with MX_DELTA_R
                                if(dR(daughter_particle, NAKED_ELE)<=MX_DELTA_R){
                                        DRESSED_ELE_4VEC += ROOT::Math::PtEtaPhiMVector(daughter_particle->pt(), daughter_particle->eta(), daughter_particle->phi(), ELE_MASS);
                                }

                                
                        }  
                }
                if (swap_electron){
                        temp_ele=swap_electron;
                }

        }
        //cout<<"Energy before:  "<<DRESSED_ELE_4VEC.E()<<endl;
        // Make a GenParticle from the vector and return it
        reco::Particle::LorentzVector dress_4_vec(DRESSED_ELE_4VEC.Pt(), DRESSED_ELE_4VEC.Eta(), DRESSED_ELE_4VEC.Phi(), DRESSED_ELE_4VEC.E());
        reco::GenParticle* dressed_e = new reco::GenParticle(NAKED_ELE->charge(), dress_4_vec, NAKED_ELE->vertex(), NAKED_ELE->pdgId(), NAKED_ELE->status(), 1);
        //cout<<"Energy after:  "<<dressed_e->energy()<<endl;
        return dressed_e;
}*/
////////////////////////////////////////////////////////////////////////////////////

//////////////Trial/////////////////////////////////////////////////////////////
const reco::GenParticle* GetDressedLepton(
            const reco::GenParticle * const BORN_LEP,
            const reco::GenParticle * const NAKED_LEP,
            const double MAX_DELTA_R,
            const int PDGID
            ) {
        double MASS;
        if (fabs(PDGID) == 11){
        MASS = 5.109989e-4;
        }
        else if (fabs(PDGID) == 13){
        MASS =  1.06e-1;
        }
        else if (fabs(PDGID) == 15){
        MASS = 1.78;
        }
        else {
                cout<<"Error Non-Leptonic PDGId entered. GetDressedLepton will return a nullptr"<<endl;
                return nullptr;
        }
        // Make a 4 vector for the dressed electron
        math::PtEtaPhiMLorentzVector dressed_p4(NAKED_LEP->pt(), NAKED_LEP->eta(), NAKED_LEP->phi(), MASS);

        // Dive down the decay tree from the born electron until we hit the
        // naked electron, saving all the photons and summing them if they are
        // within DeltaR of 0.1 of the naked electon.	
        const reco::GenParticle* tmp_lep = BORN_LEP;
        while (tmp_lep != NAKED_LEP) {
            // For some reason there are no daughters, but the particle is
            // "unstable". Abort and return nullptr.
            if (tmp_lep->numberOfDaughters() == 0) {
                return nullptr;
            }
            // Otherwise look through the daughters and find an electron
            const reco::GenParticle* swap_lep = nullptr;
            for (size_t i = 0; i < tmp_lep->numberOfDaughters(); ++i) {
                const reco::Candidate* test_particle = tmp_lep->daughter(i);
                // If we find electron, we save it as the next item to recurse over
                if (fabs(test_particle->pdgId()) == fabs( PDGID)) {
                    swap_lep = dynamic_cast<const reco::GenParticle*>(test_particle);
                } 
                // If we find a photon, add its 4 vector if it is within some
                // distance of the naked electron
                else if (fabs(test_particle->pdgId()) == 22) {
                    const double DELTA_R = deltaR(test_particle->eta(), test_particle->phi(), NAKED_LEP->eta(), NAKED_LEP->phi());
                    //cout<<"Status of photon: "<<test_particle->status()<<endl;
                    if (DELTA_R < MAX_DELTA_R) {
                        dressed_p4 += math::PtEtaPhiMLorentzVector(
                                test_particle->pt(),
                                test_particle->eta(),
                                test_particle->phi(),
                                MASS
                                );
                    }
                }
            }
            // Now that we done searching this level of the decay tree, move to
            // the next
            if (swap_lep) {
                tmp_lep = swap_lep;
            }
        }

        // Make a GenParticle from the vector and return it
        reco::GenParticle* dressed_lep = new reco::GenParticle(
                NAKED_LEP->charge(),
                dressed_p4,
                NAKED_LEP->vertex(),
                NAKED_LEP->pdgId(),
                NAKED_LEP->status(),
                1
            );

        return dressed_lep;
    }


////Function to store the dR values of all photons emitted in the born -> naked process: Using a modified version of the GeDressed Function////
vector<double> GetdR( const reco::GenParticle * const BORN_LEP, const reco::GenParticle * const NAKED_LEP, const int PDGID){
        // Make a 4 vector for the dressed electron
        //math::PtEtaPhiMLorentzVector dressed_p4(NAKED_ELECTRON->pt(), NAKED_ELECTRON->eta(), NAKED_ELECTRON->phi(), ELECTRON_MASS);

        // Dive down the decay tree from the born electron until we hit the
        // naked electron, saving all the photons and summing them if they are
        // within DeltaR of 0.1 of the naked electon.
        
        vector<double> delta_R;
        const reco::GenParticle* tmp_lep = BORN_LEP;
        while (tmp_lep != NAKED_LEP) {
            // For some reason there are no daughters, but the particle is
            // "unstable". Abort and return nullptr.
            //
            //cout<<tmp_lep->numberOfDaughters()<<endl;
            if (tmp_lep->numberOfDaughters() != 0) {
                //return vector<double>();
            //
            // Otherwise look through the daughters and find an electron
            const reco::GenParticle* swap_lep = nullptr;
            for (size_t i = 0; i < tmp_lep->numberOfDaughters(); ++i) {
                const reco::Candidate* test_particle = tmp_lep->daughter(i);
                // If we find electron, we save it as the next item to recurse over
                //
                if (fabs(test_particle->pdgId()) == fabs(PDGID)) {
                    //
                    swap_lep = dynamic_cast<const reco::GenParticle*>(test_particle);
                    //
                }
                // If we find a photon, add its 4 vector if it is within some
                // distance of the naked electron
                else if (fabs(test_particle->pdgId()) == 22) {
                    //
			
                    const double DELTA_R = deltaR(test_particle->eta(), test_particle->phi(), NAKED_LEP->eta(), NAKED_LEP->phi());
                    //
                    delta_R.push_back(DELTA_R);
                }
            }
            // Now that we done searching this level of the decay tree, move to
            // the next
            if (swap_lep) {
                //
                tmp_lep = swap_lep;
                //
            }
        }
        else{
                break;
        }

        // Make a GenParticle from the vector and return it
        /*reco::GenParticle* dressed_e = new reco::GenParticle(
                NAKED_ELECTRON->charge(),
                dressed_p4,
                NAKED_ELECTRON->vertex(),
                NAKED_ELECTRON->pdgId(),
                NAKED_ELECTRON->status(),
                1
            );*/
        }
        return delta_R;
        
        
}


////////////////////////////////



////////////Function Definition: GetNakedLepton///////////////////////////////////
const reco::GenParticle* GetNakedLepton(const reco::GenParticle* const BORN_LEP, const int PDGID){
        //Initialising the pointer
        const reco::GenParticle* naked_lep=BORN_LEP;
        //Now iterate over the daughters untill find a stable electron
        while (naked_lep ->status()!=1){
                //if it does not have any daughters then return nullptr
                if(naked_lep->numberOfDaughters()==0){
                        return nullptr;
                }
                //if not then we continue our search
                for(size_t i=0; i< naked_lep->numberOfDaughters(); i++){
                        const reco::Candidate* test_part = naked_lep->daughter(i);
                        if (fabs(test_part->pdgId())==fabs(PDGID)){
                                naked_lep=dynamic_cast<const reco::GenParticle*>(test_part);
                                break;
                        }

                }
                
        }
        return naked_lep;

}


bool descending(double i, double j) { return i > j; }

vector<int> index_list(vector<float> v){
	vector<float> l=v;
	vector<int> index;
	sort(l.begin(),l.end(), descending);
        for (size_t j=0; j<l.size(); j++){
        	for (size_t i=0; i<v.size(); i++){
                if (l[j]==v[i]){
                	index.push_back(i);
                        break;
                }
                }
        }
        return index;
}




class ZEE_RecHit_NTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZEE_RecHit_NTuplizer(const edm::ParameterSet&);
      ~ZEE_RecHit_NTuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

     std::vector<float> Hit_ES_Eta[2];
     std::vector<float> Hit_ES_Phi[2];
     std::vector<float> Hit_ES_X[2];
     std::vector<float> Hit_ES_Y[2];
     std::vector<float> Hit_ES_Z[2];
     std::vector<float> ES_RecHitEn[2];


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


//   cluster tools
      EcalClusterLazyTools *clustertools;
      noZS::EcalClusterLazyTools *clustertools_NoZS;

//   Identify if the SC lies in EB OR EE based on its seed
     bool isEB = 0;
     bool isEE = 0; // !isEB not sufficient since later will try to include the preshower as well

// Get the hits from the ES
//     std::vector<GlobalPoint> GetESPlaneRecHits(const reco::SuperCluster& sc, unsigned int planeIndex) const;
     void GetESPlaneRecHits(const reco::SuperCluster& sc, const CaloGeometry* &geo, unsigned int elenum, unsigned int planeIndex);

//   clear the vectors 
     void ClearTreeVectors();
// Initialisation
     void initialisation();
// ----------member data ---------------------------
     TTree* T;
     //My variables
     double qS;
     double alpQCD;
     double alpQED;
     int idFirst;
     int idSecond;
     double xFirst;
     double xSecond;
     double scalePDF;

        int event_num=0;

     int nMEP;
     int nMEPF;
     unsigned int sPID;
     double wt;
     double wtP;
     bool hPF;
     bool hDRV;

     // Variables for Run info.
     int run;
     int event;
     int lumi;

     double number_of_dressed_ele;
     int num_gen_jets;
     int num_reco_jets;
///////////My Electron Variables//////////////////////////////
     //int numElectrons_
     // Electron variables
     int nElectrons_;
     Float_t rho;
     bool isMC_;
     /*std::vector<float> iEta[2];
     std::vector<float> iPhi[2];
     std::vector<float> Hit_Eta[2];
     std::vector<float> Hit_Phi[2];
     std::vector<float> Hit_X[2];
     std::vector<float> Hit_Y[2];
     std::vector<float> Hit_Z[2];

        
     std::vector<float> RecHitFrac[2];
     std::vector<float> RecHitEn[2];

     std::vector<float> Ele_pt_;
     std::vector<float> Ele_eta_;
     std::vector<float> Ele_phi_;
     std::vector<float> Ele_energy_;
     std::vector<float> Ele_ecal_energy_;
     std::vector<float> Ele_ecal_mustache_energy_;

     std::vector<float> Ele_R9;
     std::vector<float> Ele_S4;
     std::vector<float> Ele_SigIEIE;
     std::vector<float> Ele_SigIPhiIPhi;
     std::vector<float> Ele_SCEtaW;
     std::vector<float> Ele_SCPhiW;
     std::vector<float> Ele_CovIEtaIEta;
     std::vector<float> Ele_CovIEtaIPhi;
     std::vector<float> Ele_ESSigRR;
     std::vector<float> Ele_SCRawE;
     std::vector<float> Ele_SC_ESEnByRawE;
     std::vector<float> Ele_HadOverEm;

     std::vector<float> HFEMClust_pt_;
     std::vector<float> HFEMClust_eta_;
     std::vector<float> HFEMClust_phi_;
     std::vector<float> HFEMClust_energy_;
    //energy in long or short fibers various cluster sizes
     std::vector<float> HFEMClust_eLong1x1_;
     std::vector<float> HFEMClust_eShort1x1_;
     std::vector<float> HFEMClust_eLong3x3_;
     std::vector<float> HFEMClust_eShort3x3_;
     std::vector<float> HFEMClust_eLong5x5_;
     std::vector<float> HFEMClust_eShort5x5_;
    //total energy in various clusters
     std::vector<float> HFEMClust_e1x1_;
     std::vector<float> HFEMClust_e3x3_;
     std::vector<float> HFEMClust_e5x5_;
    //Identification Variables
     std::vector<float> HFEMClust_eSeL_; //Longitudinal variable: E(3x3,short fibers)/E(3x3,long fibers)
     std::vector<float> HFEMClust_eCOREe9_; // Transverse Variable: E(Core of cluster)/E(3x3)
     std::vector<float> HFEMClust_e9e25_; // Shower Exclusion Variable: E(3x3)/E(5x5)
     std::vector<float> HFEMClust_eCore_; // energy in central highest energy cells (at least 50% energy of previous total energy startign with seed cell)


*/	
	

     std::vector<float> Ele_Gen_Pt;
     std::vector<float> Ele_Gen_Eta;
     std::vector<float> Ele_Gen_Phi;
     std::vector<float> Ele_Gen_E;
     /*
     std::vector<float> eta_lead_vec;
     std::vector<float> eta_sublead_vec;
     std::vector<float> eta_lead_brl;
     std::vector<float> eta_lead_ecp;
     std::vector<float> eta_sublead_brl;
     std::vector<float> eta_sublead_ecp;
     */

/*
     std::vector<int> passLooseId_;
     std::vector<int> passMediumId_;
     std::vector<int> passTightId_;
     std::vector<int> passMVAMediumId_;

     std::vector<int> isTrue_;
*/
////////////////////Defining variables to hold the properties of dressed electrons///////////////////////////
     std::vector<float> Dressed_Ele_Pt;
     std::vector<float> Dressed_Ele_Eta;
     std::vector<float> Dressed_Ele_Phi;
     std::vector<float> Dressed_Ele_E;
     std::vector<double> dR_values;
     std::vector<double> dR_0;
     std::vector<double> dR_1;

     std::vector<float> Born_Ele_Pt;
     std::vector<float> Born_Ele_Eta;
     std::vector<float> Born_Ele_Phi;
     std::vector<float> Born_Ele_E;

     std::vector<float> Naked_Ele_Pt;
     std::vector<float> Naked_Ele_Eta;
     std::vector<float> Naked_Ele_Phi;
     std::vector<float> Naked_Ele_E;

     std::vector<float> Z_Truth_Pt;
     std::vector<float> Z_Truth_Eta;
     std::vector<float> Z_Truth_Phi;
     std::vector<float> Z_Truth_E;

     std::vector<float> recon_ele_Pt_vec;
     std::vector<float> recon_ele_Eta_vec;
     std::vector<float> recon_ele_Phi_vec;
     std::vector<float> recon_ele_E_vec;

     std::vector<float> gsf_electrons_Pt_vec;
     std::vector<float> gsf_electrons_Eta_vec;
     std::vector<float> gsf_electrons_Phi_vec;
     std::vector<float> gsf_electrons_E_vec;

////////GenJet vectors//////////////////////////////

     std::vector<double> Gen_Jet_Pt_vec;
     std::vector<double> Gen_Jet_Eta_vec;
     std::vector<double> Gen_Jet_Phi_vec;
     std::vector<double> Gen_Jet_E_vec;
     std::vector<int> Gen_Jet_Mother_vec;
     //std::vector<double>dR_values_jets;
//     std::vector<int> Gen_Jet_Num_vec;
     //std::vector<float> Hadron_Flavour_vec;

////////RecoJet vectors////////////////////////////
     std::vector<double> Reco_Jet_Pt_vec;
     std::vector<double> Reco_Jet_Eta_vec;
     std::vector<double> Reco_Jet_Phi_vec;
     std::vector<double> Reco_Jet_E_vec;
     std::vector<int> Reco_Jet_Mother_vec;
//     std::vector<int> Reco_Jet_Num_vec;



      // -----------------Handles--------------------------
 /*     edm::Handle<EcalRecHitCollection> EBRechitsHandle;
      edm::Handle<EcalRecHitCollection> EERechitsHandle;
      edm::Handle<EcalRecHitCollection> ESRechitsHandle;
      edm::Handle<reco::RecoEcalCandidateCollection> HFElectrons;
      edm::Handle<reco::SuperClusterCollection> HFSC;
      edm::Handle<reco::HFEMClusterShapeAssociationCollection> HFClusterCollection;*/
      edm::Handle<edm::View<pat::Electron>> electrons;  
      edm::Handle<GenEventInfoProduct> xhandle;
      edm::Handle<edm::View<reco::GenParticle> >gphandle;
      //edm::Handle<edm::View<pat::Jet>> jets_handle;
      
/////////////Z peak///////////////////////////////////////
      edm::Handle<reco::GenParticleCollection> gpchandle;
/////////////////Jets///////////////////////////////
      edm::Handle<reco::GenJetCollection> gjchandle;
      edm::Handle<std::vector<pat::Jet>> recjhandle;

//      edm::Handle<edm::View<reco::GenParticle> > genParticles;
//      edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
//      edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
      //---------------- Input Tags-----------------------
 /*     edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionEBToken_;
      edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionEEToken_;
      edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionESToken_;
      edm::EDGetTokenT<reco::RecoEcalCandidateCollection> HFElectronsToken_;
      edm::EDGetTokenT<reco::SuperClusterCollection> HFSCToken_;
      edm::EDGetTokenT<reco::HFEMClusterShapeAssociationCollection> HFClusterToken_;*/
      edm::EDGetToken electronsToken_;
      edm::EDGetTokenT<GenEventInfoProduct> xtoken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > gptoken_;
      edm::EDGetTokenT<reco::GenParticleCollection> gpctoken_;
      edm::EDGetTokenT<reco::GenJetCollection> gjctoken_;
      edm::EDGetTokenT<std::vector<pat::Jet>> recjtoken_;
     // edm::EDGetTokenT<edm::View<pat::Jet>> jets_;
//      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
//      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
//      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;


};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
ZEE_RecHit_NTuplizer::ZEE_RecHit_NTuplizer(const edm::ParameterSet& iConfig):
/*   recHitCollectionESToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsES"))),
   HFElectronsToken_(consumes<reco::RecoEcalCandidateCollection>(edm::InputTag("hfRecoEcalCandidate"))),
   HFSCToken_(consumes<reco::SuperClusterCollection>(edm::InputTag("hfEMClusters"))),
   HFClusterToken_(consumes<reco::HFEMClusterShapeAssociationCollection>(edm::InputTag("hfEMClusters"))),*/
   xtoken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
	gpctoken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
        gjctoken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJetsAK8"))),
        recjtoken_(consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJetsAK8")))
        //jets_(consumes<pat::Jet>(edm::InputTag("slimmedJetsAK8")))
{   
   //now do what ever initialization is needed
   isMC_= iConfig.getParameter<bool>("isMC");

   electronsToken_ = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
//   genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
   gptoken_ = mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
   usesResource("TFileService");
}


ZEE_RecHit_NTuplizer::~ZEE_RecHit_NTuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)



}


//
// member functions
//

// ------------ method called for each event  ------------
void ZEE_RecHit_NTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;
   event_num=event_num+1;
   /*if(event_num%100==0){
        cout<<"Event number: "<<event_num<<endl;
   }*/
   
/*   iEvent.getByToken(HFElectronsToken_, HFElectrons);
   iEvent.getByToken(HFSCToken_, HFSC);
   iEvent.getByToken(HFClusterToken_, HFClusterCollection);*/
   iEvent.getByToken(electronsToken_, electrons);
   iEvent.getByToken(recjtoken_, recjhandle);
   //iEvent.getByToken(jets_, jets_handle);
   if(isMC_){
   	//std::cout<<"This is simulated data. Finding the Generator level information"<<std::endl;
   	initialisation();
   	iEvent.getByToken(xtoken_, xhandle);
   	iEvent.getByToken(gptoken_, gphandle);
	iEvent.getByToken(gpctoken_,gpchandle);
        iEvent.getByToken(gjctoken_, gjchandle);
        
   /*std::cout<<"qScale:"<<xhandle->qScale()<<std::endl;
   std::cout<<"alpha qcd:"<<xhandle->alphaQCD()<<std::endl;
   std::cout<<"alpha qed:"<<xhandle->alphaQED()<<std::endl;*/
        //
   	qS=xhandle->qScale();
   	alpQCD=xhandle->alphaQCD();
   	alpQED=xhandle->alphaQED();
   	nMEP=xhandle->nMEPartons();
   	nMEPF=xhandle->nMEPartons();
   	sPID=xhandle->signalProcessID();
   	wt=xhandle->weight();
   	wtP=xhandle->weightProduct();
   	hPF=xhandle->hasPDF();
   	hDRV=xhandle->hasDJRValues();

        const gen::PdfInfo* mypdf = xhandle->pdf();
        idFirst = mypdf->id.first;
        idSecond = mypdf->id.second;
        xFirst = mypdf->x.first;
        xSecond = mypdf->x.second;
        scalePDF = mypdf->scalePDF;
        //cout<<"xhandle fine"<<endl;
   }
        
  // else{
//	std::cout<<"This is real data. Generator level data will be skipped"<<std::endl;
  // }
//Clear all vectors to be written to the tree
   ClearTreeVectors();
   run=0;  event=0;  lumi=0;
  // cout<<"electrons: "<<electrons->size()<<endl;
   //cout<<"reco-jet handle: "<<recjhandle->size()<<endl;
   
///////////////////////Codeblock: Defining the pointers to born, naked and dressed electrons/////////////////////
  
const reco::GenParticle* bornEle_0 = nullptr;
const reco::GenParticle* bornEle_1 = nullptr;
const reco::GenParticle* nakedEle_0 = nullptr;
const reco::GenParticle* nakedEle_1 = nullptr;
const reco::GenParticle* dressedEle_0 = nullptr;
const reco::GenParticle* dressedEle_1 = nullptr;
const reco::GenParticle* z_boson=nullptr;

///////////////////////////Fill Electron/Photon related stuff/////////////////////////////////////////////////////
   nElectrons_ = 0;
   number_of_dressed_ele =0;
   /*for (size_t i = 0; i < electrons->size(); ++i){
	if(electrons->size() > 2) break;
	const auto ele = electrons->ptrAt(i);
        if( ele->pt() < 1. ) continue;
	if(!(ele->ecalDrivenSeed())) continue;
	if(ele->parentSuperCluster().isNull()) continue;

        nElectrons_++;
        Ele_pt_.push_back( ele->pt() );
        Ele_eta_.push_back( ele->superCluster()->eta() );
        Ele_phi_.push_back( ele->superCluster()->phi() );
        Ele_energy_.push_back( ele->energy() );
	Ele_ecal_energy_.push_back( ele->correctedEcalEnergy() );
        Ele_R9.push_back(ele->full5x5_r9());
        Ele_SigIEIE.push_back(ele->full5x5_sigmaIetaIeta());
        Ele_SigIPhiIPhi.push_back(ele->full5x5_sigmaIphiIphi());
        Ele_SCEtaW.push_back(ele->superCluster()->etaWidth());
        Ele_SCPhiW.push_back(ele->superCluster()->phiWidth());
	Ele_HadOverEm.push_back(ele->hadronicOverEm());
        const CaloClusterPtr seed_clu = ele->superCluster()->seed();
	Ele_SCRawE.push_back(ele->superCluster()->rawEnergy());
        Ele_SC_ESEnByRawE.push_back( (ele->superCluster()->preshowerEnergy())/(ele->superCluster()->rawEnergy()) );
	}   */

        
/////////////////////My analysis related to electron////////////////////////
/*	numElectron_=0;
        for (i=0; i < gphandle->size(); i++){
	if(gphandle */

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////HFCluster variables/////////////////////////////////////////////

//   std::cout<<std::endl<<" HFClusters = "<<HFClusterCollection->size()<<" HFElectrons = "<<HFSC->size()<<std::endl;
  /* for (unsigned int iii = 0; iii < HFSC->size(); ++iii){
//	std::cout<<std::endl<<" Iteration number = "<<iii++<<std::endl;
	const SuperCluster& supClus = (*HFSC)[iii];	
	reco::SuperClusterRef theClusRef = edm::Ref<reco::SuperClusterCollection>(HFSC, iii);	
	const HFEMClusterShapeRef clusShapeRef = HFClusterCollection->find(theClusRef)->val;
	const HFEMClusterShape& clusShape=*clusShapeRef;

        HFEMClust_pt_.push_back(supClus.energy()/cosh(supClus.eta()) );
        HFEMClust_eta_.push_back(supClus.eta());
        HFEMClust_phi_.push_back(supClus.phi());
        HFEMClust_energy_.push_back(supClus.energy());
        HFEMClust_eLong1x1_.push_back(clusShape.eLong1x1());
        HFEMClust_eShort1x1_.push_back(clusShape.eShort1x1());
        HFEMClust_eLong3x3_.push_back(clusShape.eLong3x3());
        HFEMClust_eShort3x3_.push_back(clusShape.eShort3x3());
        HFEMClust_eLong5x5_.push_back(clusShape.eLong5x5());
        HFEMClust_eShort5x5_.push_back(clusShape.eShort5x5());
        HFEMClust_e1x1_.push_back(clusShape.e1x1());
        HFEMClust_e3x3_.push_back(clusShape.e3x3());
        HFEMClust_e5x5_.push_back(clusShape.e5x5());
        HFEMClust_eSeL_.push_back(clusShape.eSeL());
        HFEMClust_eCOREe9_.push_back(clusShape.eCOREe9());
        HFEMClust_e9e25_.push_back(clusShape.e9e25());
        HFEMClust_eCore_.push_back(clusShape.eCore());

   }*/

///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////
   if (isMC_){
	for(edm::View<GenParticle>::const_iterator part = gphandle->begin(); part!=gphandle->end(); ++part){
		if(part ->status()==1 && abs(part ->pdgId())==11){
			Ele_Gen_Pt.push_back(part->pt());
			Ele_Gen_Eta.push_back(part->eta());
			Ele_Gen_Phi.push_back(part->phi());
			Ele_Gen_E.push_back(part->energy());
			}
		}
	
/////////////////////////////Debug: Ele_Gen_Eta//////////////////////////////////////////////////////////////////////////
    /*    //cout<<"1"<<endl;
        size_t e1=-9999;
        size_t e2=-9999;
        size_t e3=-9999;
        size_t e4=-9999;
        const reco::GenParticle* leadele = nullptr;
        const reco::GenParticle* subleadele = nullptr;




        if(Ele_Gen_E.size()>=2){
                e1=index_list(Ele_Gen_E)[0];
                e2=index_list(Ele_Gen_E)[1];

                for(size_t j =0; j< gphandle->size(); j=j+1){
                        const reco::GenParticle* part1=&gphandle->at(j);
                        if(part1 ->status()==1 && abs(part1 ->pdgId())==11 ){
                                //cout<<"Ele_Gen_E: "<<Ele_Gen_E[e1]<<endl;
                                //cout<<"part energy: "<<part1->energy()<<endl;
                                float energy1=part1->energy();
                                float energy2=Ele_Gen_E[e1];
                                float energy3=Ele_Gen_E[e2];
                                if(energy1==energy2){
                                        //cout<<"j: "<<j<<endl;
                                        //cout<<"part energy: "<<part1->energy()<<endl;
                                        e3=j;
                                        //cout<<"e3: "<<e3<<endl;
                                }
                                if(energy1==energy3){
                                        e4=j;
                                }
                        }
                }
                //cout<<"End"<<endl;

                //cout<<"e1: "<<e1<<endl;
                //cout<<"e2: "<<e2<<endl;
                //cout<<"e3: "<<e3<<endl;
                //cout<<"e4: "<<e4<<endl;
                /*ROOT::Math::PtEtaPhiEVector ele_1((Ele_Gen_Pt)[e1], (Ele_Gen_Eta)[e1],(Ele_Gen_Phi)[e1], (Ele_Gen_E)[e1]);
                ROOT::Math::PtEtaPhiEVector ele_2((Ele_Gen_Pt)[e2], (Ele_Gen_Eta)[e2],(Ele_Gen_Phi)[e2], (Ele_Gen_E)[e2]);
                assert(ele_1.E()>ele_2.E());
                if (ele_1.Pt()>=10 && ele_2.Pt()>=10 && fabs(ele_1.Eta())<2.4 && fabs(ele_2.Eta()<2.4)){
                        eta_lead_vec.push_back(ele_1.Eta());
                }*/
               //cout<<"2"<<endl;
                //int iterator=0;
                //cout<<"size: "<<gphandle->size()<<endl;
              /*  for(size_t i =0; i< gphandle->size(); i=i+1){
                        //cout<<"i: "<<i<<endl;
                        //cout<<"e1: "<<e1<<endl;
                        //cout<<"e2: "<<e2<<endl;
                        //iterator=iterator+1;
                        const reco::GenParticle* part=&gphandle->at(i);
                        //cout<<"3"<<endl;
                        /*if(part->status()==1){
                                cout<<"Ok1"<<endl;
                        }
                        if(abs(part->pdgId())==11){
                                cout<<"ok2"<<endl;
                        }
                        /*if(fabs(part->eta())<2.4){
                                cout<<"ok3"<<endl;
                        }*/
                        //cout<<"status: "<<part->status()<<endl;
                        //cout<<"pid: "<<part->pdgId()<<endl;
                        //cout<<"eta: "<<part->eta()<<endl;
                  /*      if(part->status()==1 && abs(part->pdgId())==11 ){
                                //cout<<"4"<<endl;
                                //cout<<"i: "<<i<<endl;
                                //cout<<"e1: "<<e1<<endl;
                                //cout<<"e2: "<<e2<<endl;
                                if (i==e3){
                                        //cout<<"5"<<endl;
                                        auto max_it = std::max_element(Ele_Gen_E.begin(), Ele_Gen_E.end());
                                        //cout<<"part pt: "<<part->pt()<<endl;
                                        //cout<<"max element: "<<*max_it<<endl;
                                        float part_energy=part->energy();
                                        if(part_energy!=*max_it){
                                                cout<<"Problem"<<endl;
                                                cout<<"part energy: "<<part->energy()<<endl;
                                                cout<<"max element: "<<*max_it<<endl;
                                        }
                                        leadele=part;
                                }
                                if (i==e4){
                                        //cout<<"6"<<endl;
                                        subleadele=part;
                                }

                        }
                }
                if (leadele != nullptr && subleadele !=nullptr){
                        ROOT::Math::PtEtaPhiEVector ele_1(leadele->pt(), leadele->eta(),leadele->phi(), leadele->energy());
                        ROOT::Math::PtEtaPhiEVector ele_2(subleadele->pt(), subleadele->eta(),subleadele->phi(), subleadele->energy());
                        //cout<<"7"<<endl;
                        if(ele_1.Pt()>10. && ele_2.Pt()>10.){
                                //cout<<ele_1.Eta();
                                eta_lead_vec.push_back(ele_1.Eta());
                                eta_sublead_vec.push_back(ele_2.Eta());
                                if(fabs(ele_1.Eta())<1.4){
                                        eta_lead_brl.push_back(ele_1.Eta());
                                }
                                if(fabs(ele_2.Eta())<1.4){
                                        eta_sublead_brl.push_back(ele_2.Eta());
                                }
                                if(fabs(ele_1.Eta())>=1.4){
                                        eta_lead_ecp.push_back(ele_1.Eta());
                                }
                                if(fabs(ele_2.Eta())>=1.4){
                                        eta_sublead_ecp.push_back(ele_2.Eta());
                                }
                        }     
                }
        
        }
        

        */





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        /*for(edm::View<GenParticle>::const_iterator g_part = gphandle->begin(); g_part!=gphandle->end(); ++g_part){
                if(g_part->status()!=1 && abs(g_part->pdgId())==11){
                        for(size_t j = 0; j<g_part->numberOfMothers(); ++j){
                                if(g_part->mother(j)->pdgId()==23){
                                        cout<<"We found a born electron"<<endl;
                                        for(size_t k=0; k<g_part->numberOfDaughters(); ++k){
                                                const edm::View<GenParticle>::const_iterator tp= g_part->daughter(k);
                                                cout<<"Daughter Number "<<k<<"has PDGID"<<tp->pdgId()<<endl;
                                        }
                
                                        //bornEle_0=g_part;
                                        //nakedEle_0=GetNakedElectron(g_part);
                                        //if(nakedEle_0 !=nullptr){
                                        //        cout<<"Got the naked electron too."<<endl;
                                        }
                                }
                        }
                }

        }
        /*for (unsigned int i = 0; i <gpchandle->size(); ++i){
        //Grabbing pointer of one particle of the collection
	const reco::GenParticle* gpart = &gpchandle->at(i);
        //If the grabbed particle is electron and we have atleast still one unfound dressed electron
        cout<<"Status:  "<<gpart->status()<<"PDGID:  "<<gpart->pdgId()<<endl;
   }*/
////////////////Code Block: Finding Z Mass by Bottom Up method/////////////////
////////////Lets grab a particle////////////////////////
        /*if(gpchandle.isValid()){
                //cout<<"Valid"<<endl;
                if(gpchandle->empty()){
                        //cout<<"Empty"<<endl;
                }
        }*/
        
        for (unsigned int i = 0; i <gpchandle->size(); ++i){
              //cout<<"We are inside the loop"<<endl;
              const reco::GenParticle* gpart = &gpchandle->at(i);
              //cout<<"gpchandle fine"<<endl;
              //cout<<"gpart pdgId: "<<gpart->pdgId()<<endl;
              /////If we have not found atleast one of the born electrons yet and the particle at hand is an electron or positron/////
              
              if(fabs(gpart->pdgId())==11 && (dressedEle_0 ==nullptr || dressedEle_1 == nullptr)){
                //cout<<"Number of mothers:  "<<gpart->numberOfMothers()<<endl;
                /////Loop over their mothers to see if its mother is a Z and the particle itself is not a final state particle 
                
                for(size_t j = 0; j<gpart->numberOfMothers(); ++j){
                       // cout<<"mother pdgid: "<<gpart->mother(j)->pdgId()<<" status: "<<gpart->status()<<endl;
                       /////////////////////////////////////
                       /*if(gpart->mother(j)->pdgId()==23){
                        if(gpart->status()==3){
                                cout<<"We got one"<<endl;
                        }
                        else{
                                cout<<"Nope"<<endl;

                       }
                       }*/
                      
                       /////////////////////////////////////////////
                        if (gpart->mother(j)->pdgId()==23 && gpart->status() != 1){
                               
                                //cout<<"We got a born electron"<<endl;
                                ///If we have not found the first born electron////
                                if(bornEle_0==nullptr){
                                        bornEle_0=gpart;
                                        
                                        //cout<<"Status of born0: "<<bornEle_0->status()<<endl;
                               //       cout<<"Kudos we are good to move ahead"<<endl;
                               ///Getting the corresponding naked electron
                                        nakedEle_0=GetNakedLepton(bornEle_0, 11);
                                        
                                        //cout<<"Its pt is: "<<(nakedEle_0->pt())<<endl;
                                //If we have born and naked electron we can find the dressed electron
                                        if(bornEle_0==nullptr){
                                                cout<<"Empty born"<<endl;
                                        }
                                        if(nakedEle_0==nullptr){
                                                cout<<"Emptry naked"<<endl;
                                        }
                                        if (bornEle_0 && nakedEle_0){
                                                
                                                dressedEle_0=GetDressedLepton(bornEle_0, nakedEle_0, 0.1, 11);
                                                
                                                if(bornEle_0==nullptr){
                                                cout<<"Empty born"<<endl;
                                                }
                                                if(nakedEle_0==nullptr){
                                                cout<<"Emptry naked"<<endl;
                                                }
                                                dR_0=GetdR(bornEle_0, nakedEle_0, 11);
                                                
                                                for (size_t l = 0; l < dR_0.size(); l++){
                                                        //cout<<"l: "<<l<<"dR: "<<dR_0[l]<<endl;
                                                        dR_values.push_back(dR_0[l]);
                                                        
                                                }
                                                /*for (size_t l=0; dR_0.size(); l++){
                                                        //cout<<dR_0[l];
                                                        dR_values.push_back(dR_0[l]);
                                                }*/
                                                //dR_values.insert(dR_values.end(), dR_0.begin(), dR_0.end());
                                                /*cout<<"Size:  "<<dR_values.size();
                                                for (size_t l = 0; l<dR_values.size(); l++){
                                                        cout<<"l: "<<l<<" values: "<<dR_values[l]<<endl;
                                                }*/
                                               
                                        }
                             //           cout<<"The current energy for de0 is "<<(dressedEle_0->energy())<<endl;
                                        
                                        
                                }
                                /////If we already found the first born electron then lets find and do the samething for the second born electron
                                else{
                                        bornEle_1=gpart;
                                        //cout<<"status born1: "<<bornEle_1->status()<<endl;
                                        nakedEle_1=GetNakedLepton(bornEle_1, 11);
                                        
                                        if (bornEle_1 && nakedEle_1){
                                                
                                                dressedEle_1=GetDressedLepton(bornEle_1, nakedEle_1, 0.1, 11);
                                                
                                                dR_1=GetdR(bornEle_1, nakedEle_1, 11);
                                                for (size_t b=0; b< dR_1.size(); b++){
                                                        //cout<<dR_1[b]<<endl;
                                                        dR_values.push_back(dR_1[b]);
                                                        
                                                }
                                        }
                           //             cout<<"The current energy for de1 is "<<(dressedEle_1->energy())<<endl;
                                                               
                
                                }
                        }
                }
              }
                
                ///We only expect one Z and hence one electron and positron pair, so once we have it we can break out of the loop 
                if(dressedEle_0 !=nullptr && dressedEle_1 !=nullptr){
                        //cout<<"Found "<<"for 0: "<<dressedEle_0->energy()<<"and for 1 "<<dressedEle_1->energy()<<endl;
                        break;
                }
        }
        /*if(dressedEle_0 !=nullptr && dressedEle_1!=nullptr){
        cout<<"Dressedele0 energy: "<<(dressedEle_0->energy())<<endl;
        cout<<"Dressedele1 energy: "<<(dressedEle_1->energy())<<endl;
        }*/
        //////Ensuring that the 0 is the most energetic electron among the two/////
        
        if(dressedEle_0 !=nullptr && dressedEle_1 !=nullptr){
                if (dressedEle_0->pt() < dressedEle_1->pt()){
                        //cout<<"Swapping active"<<endl;
                        std::swap(dressedEle_0, dressedEle_1);
                        std::swap(nakedEle_0, nakedEle_1);
                        std::swap(bornEle_0, bornEle_1);
                                        
                }
        
        
                Dressed_Ele_Pt.push_back(dressedEle_0->pt());
                Dressed_Ele_Eta.push_back(dressedEle_0->eta());
                Dressed_Ele_Phi.push_back(dressedEle_0->phi());
                Dressed_Ele_E.push_back(dressedEle_0->energy());

                Dressed_Ele_Pt.push_back(dressedEle_1->pt());
                Dressed_Ele_Eta.push_back(dressedEle_1->eta());
                Dressed_Ele_Phi.push_back(dressedEle_1->phi());
                Dressed_Ele_E.push_back(dressedEle_1->energy());


                Born_Ele_Pt.push_back(bornEle_0->pt());
                Born_Ele_Eta.push_back(bornEle_0->eta());
                Born_Ele_Phi.push_back(bornEle_0->phi());
                Born_Ele_E.push_back(bornEle_0->energy());

                Born_Ele_Pt.push_back(bornEle_1->pt());
                Born_Ele_Eta.push_back(bornEle_1->eta());
                Born_Ele_Phi.push_back(bornEle_1->phi());
                Born_Ele_E.push_back(bornEle_1->energy());
                
	
                Naked_Ele_Pt.push_back(nakedEle_0->pt());
                Naked_Ele_Eta.push_back(nakedEle_0->eta());
                Naked_Ele_Phi.push_back(nakedEle_0->phi());
                Naked_Ele_E.push_back(nakedEle_0->energy());

                Naked_Ele_Pt.push_back(nakedEle_1->pt());
                Naked_Ele_Eta.push_back(nakedEle_1->eta());
                Naked_Ele_Phi.push_back(nakedEle_1->phi());
                Naked_Ele_E.push_back(nakedEle_1->energy());

                for (size_t l = 0; l<electrons->size(); ++l){
                        const auto ele1 = electrons->ptrAt(l);
                        //cout<<"electrons fine"<<endl;
                        if (deltaR(dressedEle_0->eta(), dressedEle_0->phi(), ele1->eta(), ele1->phi()) <0.1 || deltaR(dressedEle_1->eta(), dressedEle_1->phi(), ele1->eta(), ele1->phi()) < 0.1){
                                recon_ele_Pt_vec.push_back(ele1->pt());
                                recon_ele_Eta_vec.push_back(ele1->eta());
                                recon_ele_Phi_vec.push_back(ele1->phi());
                                recon_ele_E_vec.push_back(ele1->energy());
                        }
                }
        }
        //cout<<"1"<<endl;
   }
///////////Code block: To store the gdgsfelectrons////////////////////////////////////////////////////////
        
        for(size_t l=0; l<electrons->size(); ++l){
                const auto reco_electron1=electrons->ptrAt(l);
                gsf_electrons_Pt_vec.push_back(reco_electron1->pt());
                gsf_electrons_Eta_vec.push_back(reco_electron1->eta());
                gsf_electrons_Phi_vec.push_back(reco_electron1->phi());
                gsf_electrons_E_vec.push_back(reco_electron1->energy());
        }
        //cout<<"2"<<endl;
////Code block: To find the True Z////////////////
        //cout<<"Going inside the loop"<<endl;
if(isMC_){
        for (unsigned int k = 0; k <gpchandle->size(); ++k){
       // cout<<"We are inside the loop"<<endl;
                const reco::GenParticle* gpart1 = &gpchandle->at(k);
                //cout<<"gpart pdgId: "<<gpart->pdgId()<<endl;
                if (gpart1->pdgId()==23 && z_boson == nullptr){
                        z_boson = gpart1;
                        Z_Truth_Pt.push_back(z_boson->pt());
                        Z_Truth_Eta.push_back(z_boson->eta());
                        Z_Truth_Phi.push_back(z_boson->phi());
                        Z_Truth_E.push_back(z_boson->energy());
                        //cout<<"Found z"<<endl;
                }
        
        }
}

        //cout<<"3"<<endl;
////////Code block: GenJets/////////////////////////////////
        //cout<<"We are here"<<endl;
        //cout<<"GenJetCollection size: "<<gjchandle->size()<<endl;
        //
if(isMC_){
        for(unsigned int g = 0; g < gjchandle->size(); ++g){
                const reco::GenJet* gj = &gjchandle->at(g);
                //cout<<"gjchandle fine"<<endl;
        //////////////Pt cut on GenJet/////////////////////////////////

                //cout<<"GenJet Pt: "<<gj->pt()<<endl;
                num_gen_jets=num_gen_jets+1.0;
                //cout<<"Num Gen Jet current: "<<num_gen_jets<<endl;
                Gen_Jet_Pt_vec.push_back(gj->pt());
                Gen_Jet_Eta_vec.push_back(gj->eta());
                Gen_Jet_Phi_vec.push_back(gj->phi());
                Gen_Jet_E_vec.push_back(gj->energy());
                //Hadron_Flavour_vec.push_back(gj->hadronFlavour());
                //int mother = -9999;
                //cout<<"Number of mothers: "<<gj->size()<<endl;
                //cout<<"3.5"<<endl;
                /////////////////////////////////
                //cout<<__LINE__<<endl;
                const auto& gjconstituents = gj->getJetConstituents();
                //cout<<__LINE__<<endl;
                if(!gjconstituents.empty()){
                        for(const auto& gjc : gjconstituents){
                                //cout<<"1"<<endl;
                                if(gjc.isNonnull()){
                                        if(gjc->numberOfMothers()!=0){
                                        //cout<<__LINE__<<endl;
                                        const reco::Candidate* mthr = gjc ->mother();
                                        Gen_Jet_Mother_vec.push_back(mthr->pdgId());
                                        }
                                }
                                

                                

                        }

                }
                
               //////////////////////////////////////
                //cout<<"Fine"<<endl;

        }
}
        //cout<<"Num GenJets: "<<num_gen_jets<<"  "<<"Size of Pt vector: "<<Gen_Jet_Pt_vec.size()<<endl;
        //cout<<"GenJet Number: "<<num_gen_jets<<endl;
        //cout<<"GenJet pt size: "<<Gen_Jet_Pt_vec.size()<<endl;
        //cout<<"4"<<endl;
/////////////Code block: RecoJets/////////////////////
        //
        for(unsigned int q = 0; q< recjhandle->size(); ++q){
                const pat::Jet* pfj = &recjhandle->at(q);
                //cout<<"recjhandle fine"<<endl;
                num_reco_jets=num_reco_jets+1.0;
                //cout<<"Number of Reco Jets: "<<num_reco_jets<<endl;
                Reco_Jet_Pt_vec.push_back(pfj->pt());
                Reco_Jet_Eta_vec.push_back(pfj->eta());
                Reco_Jet_Phi_vec.push_back(pfj->phi());
                Reco_Jet_E_vec.push_back(pfj->energy());

                //int reco_mother=-9999;

                const auto& rjconstituents = pfj->getJetConstituents();
                for(const auto& rjc : rjconstituents){
                        const reco::Candidate* rec_mthr = rjc -> mother();
                        if (rec_mthr){
                        cout<<"Mother found: "<<rec_mthr->pdgId()<<endl;
                        }
                       // Reco_Jet_Mother_vec.push_back(rec_mthr->pdgId());
                       else{
                        cout<<"No mother found"<<endl;
                       }
                }
                /*const auto& rjconst= pfj->getPFConstituents();
                for(const auto& rjp : rjconst){
                        const reco::Candidate* rec_mthr = rjp->mother();
                        if(rec_mthr){
                                Reco_Jet_Mother_vec.push_back(rec_mthr->pdgId());
                                cout<<"Found Reco Mother"<<endl;
                        }
                }*/
                

        }
        //cout<<"5"<<endl;

       ////////////////////////dR_values GenJet and RecoJet///////////////////////////////
        /*for(unsigned int g = 0; g < gjchandle->size(); ++g){
                const reco::GenJet* gj = &gjchandle->at(g);
        //////////////Pt cut on GenJet/////////////////////////////////

                if(gj->pt()>=20){
                        for(unsigned int q = 0; q< recjhandle->size(); ++q){
                                const reco::PFJet* pfj = &recjhandle->at(q);
                                if(pfj->pt()>=20){
                                        dR_values_jets.push_back(dR_gen_reco_jets(gj, pfj));
                                }

                        }
                }
        }*/
        //cout<<"RecoJet Num:"<< num_reco_jets<<endl;
        //cout<<"RecoJet Pt size"<<Reco_Jet_Pt_vec.size()<<endl;

//////Codeblock: dR matching with Gsfelectron/////////////////////
        /*for (size_t l = 0; l<electrons->size(); ++l){
                const auto ele1 = electrons->ptrAt(l);
                cout<<"Pt: "<<Dressed_Ele_Pt[0]<<" Eta "<<Dressed_Ele_Eta[0]<<" Phi: "<<Dressed_Ele_Phi[0]<<" Energy: "<<Dressed_Ele_E[0]<<" Size: "<< Dressed_Ele_E.size()<<endl;
                if (deltaR(dressedEle_0->eta(), dressedEle_0->phi(), ele1->eta(), ele1->phi()) <0.1 || deltaR(dressedEle_1->eta(), dressedEle_1->phi(), ele1->eta(), ele1->phi()) < 0.1){
                        recon_ele_Pt_vec.push_back(ele1->pt());
                        recon_ele_Eta_vec.push_back(ele1->eta());
                        recon_ele_Phi_vec.push_back(ele1->phi());
                        recon_ele_E_vec.push_back(ele1->energy());
                }
        }*/

///////////////////////Codeblock: Z Peak from Mother Daughter matching using Bottom Up Approach////////////////////////////////////
   /*for (unsigned int i = 0; i <gpchandle->size(); ++i){
        //Grabbing pointer of one particle of the collection
	const reco::GenParticle* gpart = &gpchandle->at(i);
        //If the grabbed particle is electron and we have atleast still one unfound dressed electron
	if (fabs(gpart->pdgId())==11 && (dressedEle_0 == nullptr || dressedEle_1 == nullptr)){
		for (size_t j=0; j< gpart->numberOfMothers(); ++j){

			if(gpart->mother(j)->pdgId()==23 && gpart->status()!=1){
			
				if (bornEle_0 ==nullptr){
					bornEle_0=gpart;

                                        //Initialising the pointer

                                        const reco::GenParticle* naked_ele_0 = bornEle_0;
                                        //Now iterate over the daughters untill find a stable electron
                                        while (naked_ele_0 ->status()!=1){
                                                //if it does not have any daughters then return nullptr
                                                if(naked_ele->numberOfDaughters()==0){
                                                        return nullptr;
                                                }
                                        //if not then we continue our search
                                        for(size_t i=0; i< naked_ele_0->numberOfDaughters(); i++){
                                                const reco::Candidate* test_part = naked_ele_0->daughter(i);
                                                if (fabs(test_part->pdgId())==11){
                                                        naked_ele_0=dynamic_cast<const reco::GenParticle*>(test_part);
                                                        break;
                                                }

                                        }
                                        
                                }

                        }


					cout<<"naked_ele_0"<<typeid(*naked_ele).name()<<endl;
					dressedEle_0=GetDressedElectron(bornEle_0, nakedEle_0, 0.3);
				}
				else{
					bornEle_1=gpart;
					nakedEle_1=GetNakedElectron(bornEle_1);
					dressedEle_1=GetDressedElectron(bornEle_1, nakedEle_1, 0.3);
				}
			}
		}
	}
        //We expect one Z boson in one event so as soon as we find both the Dressed Electron we can come out of the loop
   	if(dressedEle_0 != nullptr && dressedEle_1 !=nullptr){
		break;
	}
   }
   //Making sure more energetic electron is at 0th index
   if (dressedEle_0->energy() < dressedEle_1->energy()){
	std::swap(dressedEle_0, dressedEle_1);
	std::swap(nakedEle_0, nakedEle_1);
	std::swap(bornEle_0, bornEle_1);

   }

  if (dressedEle_0 != nullptr && dressedEle_1 !=nullptr){
	std::cout<<"Got it"<<"Dressed_Ele_0_Pt:"<<dressedEle_0->pt()<<endl;
        std::cout<<"Got it"<<"Dressed_Ele_1_Pt:"<<dressedEle_1->pt()<<endl;
  }*/
/////////////////////////Run, event, lumi//////////////////////////////////
   run=iEvent.id().run();
   event=iEvent.id().event();
   lumi=iEvent.luminosityBlock();
///////////////////////////////////////////////////////////////////////////



   T->Fill(); // Write out the events


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void ZEE_RecHit_NTuplizer::beginJob()
{
        edm::Service<TFileService> fs;
        T=fs->make<TTree>("T","MyTuple");
/*	T->Branch("iEtaEle1"  ,  &(iEta[0]));
        T->Branch("iPhiEle1"  ,  &(iPhi[0]));
	T->Branch("Hit_ES_Eta_Ele1"  ,  &(Hit_ES_Eta[0]));
        T->Branch("Hit_ES_Phi_Ele1"  ,  &(Hit_ES_Phi[0]));
        T->Branch("Hit_ES_X_Ele1"  ,  &(Hit_ES_X[0]));
        T->Branch("Hit_ES_Y_Ele1"  ,  &(Hit_ES_Y[0]));
        T->Branch("Hit_ES_Z_Ele1"  ,  &(Hit_ES_Z[0]));
        T->Branch("ES_RecHitEnEle1"  ,  &(ES_RecHitEn[0]));

	T->Branch("Hit_Eta_Ele1"  ,  &(Hit_Eta[0]));
        T->Branch("Hit_Phi_Ele1"  ,  &(Hit_Phi[0]));
	T->Branch("Hit_X_Ele1"  ,  &(Hit_X[0]));
        T->Branch("Hit_Y_Ele1"  ,  &(Hit_Y[0]));
	T->Branch("Hit_Z_Ele1"  ,  &(Hit_Z[0]));
        T->Branch("RecHitEnEle1"  ,  &(RecHitEn[0]));
	T->Branch("RecHitFracEle1"  ,  &(RecHitFrac[0]));
	T->Branch("iEtaEle2"  ,  &(iEta[1]));
        T->Branch("iPhiEle2"  ,  &(iPhi[1]));
	T->Branch("Hit_ES_Eta_Ele2"  ,  &(Hit_ES_Eta[1]));
        T->Branch("Hit_ES_Phi_Ele2"  ,  &(Hit_ES_Phi[1]));
        T->Branch("Hit_ES_X_Ele2"  ,  &(Hit_ES_X[1]));
        T->Branch("Hit_ES_Y_Ele2"  ,  &(Hit_ES_Y[1]));
        T->Branch("Hit_ES_Z_Ele2"  ,  &(Hit_ES_Z[1]));
	T->Branch("ES_RecHitEnEle2"  ,  &(ES_RecHitEn[1]));

	T->Branch("Hit_Eta_Ele2"  ,  &(Hit_Eta[1]));
        T->Branch("Hit_Phi_Ele2"  ,  &(Hit_Phi[1]));
	T->Branch("Hit_X_Ele2"  ,  &(Hit_X[1]));
        T->Branch("Hit_Y_Ele2"  ,  &(Hit_Y[1]));
        T->Branch("Hit_Z_Ele2"  ,  &(Hit_Z[1]));
        T->Branch("RecHitEnEle2"  ,  &(RecHitEn[1]));
        T->Branch("RecHitFracEle2"  ,  &(RecHitFrac[1]));

        T->Branch("nElectrons",  &nElectrons_ , "nEle/I");
        T->Branch("pt"  ,  &Ele_pt_);
        T->Branch("eta" ,  &Ele_eta_ );
        T->Branch("phi" ,  &Ele_phi_ );
        T->Branch("energy", &Ele_energy_);
	T->Branch("energy_ecal", &Ele_ecal_energy_);
	T->Branch("energy_ecal_mustache", &Ele_ecal_mustache_energy_);

        T->Branch("passMediumId" ,  &passMediumId_ );
        T->Branch("passTightId"  ,  &passTightId_ );
        T->Branch("passMVAMediumId", &passMVAMediumId_);

        T->Branch("Ele_R9"  ,  &Ele_R9);
        T->Branch("Ele_S4"  ,  &Ele_S4);
        T->Branch("Ele_SigIEIE"  ,  &Ele_SigIEIE);
        T->Branch("Ele_SigIPhiIPhi" , &Ele_SigIPhiIPhi);
        T->Branch("Ele_SCEtaW"  ,  &Ele_SCEtaW);
        T->Branch("Ele_SCPhiW"  ,  &Ele_SCPhiW);
        T->Branch("Ele_CovIEtaIEta"  ,  &Ele_CovIEtaIEta);
        T->Branch("Ele_CovIEtaIPhi"  ,  &Ele_CovIEtaIPhi);
        T->Branch("Ele_ESSigRR"  ,  &Ele_ESSigRR);
        T->Branch("Ele_SCRawE"  ,  &Ele_SCRawE);
        T->Branch("Ele_SC_ESEnByRawE"  ,  &Ele_SC_ESEnByRawE);
	T->Branch("Ele_HadOverEm"  ,  &Ele_HadOverEm);
	
	T->Branch("HFEMClust_pt", &HFEMClust_pt_);
	T->Branch("HFEMClust_eta", &HFEMClust_eta_);
        T->Branch("HFEMClust_phi", &HFEMClust_phi_);
        T->Branch("HFEMClust_energy", &HFEMClust_energy_);
        T->Branch("HFEMClust_eLong1x1", &HFEMClust_eLong1x1_);
        T->Branch("HFEMClust_eShort1x1", &HFEMClust_eShort1x1_);
        T->Branch("HFEMClust_eLong3x3", &HFEMClust_eLong3x3_);
        T->Branch("HFEMClust_eShort3x3", &HFEMClust_eShort3x3_);
        T->Branch("HFEMClust_eLong5x5", &HFEMClust_eLong5x5_);
        T->Branch("HFEMClust_eShort5x5", &HFEMClust_eShort5x5_);
        T->Branch("HFEMClust_e1x1", &HFEMClust_e1x1_);
        T->Branch("HFEMClust_e3x3", &HFEMClust_e3x3_);
        T->Branch("HFEMClust_e5x5", &HFEMClust_e5x5_);
        T->Branch("HFEMClust_eSeL", &HFEMClust_eSeL_);
        T->Branch("HFEMClust_eCOREe9", &HFEMClust_eCOREe9_);
        T->Branch("HFEMClust_e9e25", &HFEMClust_e9e25_);
        T->Branch("HFEMClust_eCore", &HFEMClust_eCore_);
        T->Branch("rho", &rho, "rho/F");*/

        /*T->Branch("run",&run,"run/I");
        T->Branch("event",&event,"event/I");
        T->Branch("lumi",&lumi,"lumi/I");*/
/////////////////////////////////////////////////////////
	if(isMC_){
		T->Branch("qScale",&qS,"qS/D");
		T->Branch("alphaQCD",&alpQCD,"alpQCD/D");
		T->Branch("alphaQED",&alpQED,"alpQED/D");
                T->Branch("idFirst", &idFirst, "idFirst/I");
                T->Branch("idSecond", &idSecond, "idSecond/I");
                T->Branch("xFirst", &xFirst, "xFirst/D");
                T->Branch("xSecond", &xSecond, "xSecond/D");
                T->Branch("scalePDF", &scalePDF, "scalePDF/D");


		T->Branch("nMEPartons",&nMEP,"nMEP/I");
		T->Branch("nMEPartonsFiltered",&nMEPF,"nMEPF/I");
		T->Branch("signalProcessID",&sPID,"sPID/i");
		T->Branch("weight",&wt,"wt/D");
		T->Branch("weightProduct",&wtP,"wtP/D");
		T->Branch("hasPDF",&hPF,"hPF/O");
		T->Branch("hasDJRValues",&hDRV,"hDRV/O");
                
		T->Branch("Ele_Gen_Pt" , &Ele_Gen_Pt);
                T->Branch("Ele_Gen_Eta" , &Ele_Gen_Eta);
                T->Branch("Ele_Gen_Phi" , &Ele_Gen_Phi);
                T->Branch("Ele_Gen_E" , &Ele_Gen_E);

	
                T->Branch("DressedEle_Pt", &Dressed_Ele_Pt);
                T->Branch("DressedEle_Eta", &Dressed_Ele_Eta);
                T->Branch("DressedEle_Phi", &Dressed_Ele_Phi);
                T->Branch("DressedEle_E", &Dressed_Ele_E);
                
                T->Branch("NakedEle_Pt", &Naked_Ele_Pt);
                T->Branch("NakedEle_Eta", &Naked_Ele_Eta);
                T->Branch("NakedEle_Phi", &Naked_Ele_Phi);
                T->Branch("NakedEle_E", &Naked_Ele_E);

                T->Branch("BornEle_Pt", &Born_Ele_Pt);
                T->Branch("BornEle_Eta", &Born_Ele_Eta);
                T->Branch("BornEle_Phi", &Born_Ele_Phi);
                T->Branch("BornEle_E", &Born_Ele_E);

                T->Branch("ZTruth_Pt", &Z_Truth_Pt);
                T->Branch("ZTruth_Eta", &Z_Truth_Eta);
                T->Branch("ZTruth_Phi", &Z_Truth_Phi);
                T->Branch("ZTruth_E", &Z_Truth_E);

                T->Branch("Recon_Ele_Pt", &recon_ele_Pt_vec);
                T->Branch("Recon_Ele_E", &recon_ele_E_vec);
                T->Branch("Recon_Ele_Phi", &recon_ele_Phi_vec);
                T->Branch("Recon_Ele_Eta", &recon_ele_Eta_vec);
                T->Branch("DR_values", &dR_values);


                ////////////////GenJet////////////////////////
                T->Branch("Number_GenJets", &num_gen_jets,"num_gen_jets/I");
                T->Branch("GenJet_Pt", &Gen_Jet_Pt_vec);
                T->Branch("GenJet_Eta", &Gen_Jet_Eta_vec);
                T->Branch("GenJet_Phi", &Gen_Jet_Phi_vec);
                T->Branch("GenJet_Energy", &Gen_Jet_E_vec);
                //T->Branch("GenJet_Hadron_Flavour", &Hadron_Flavour_vec);
                T->Branch("GenJet_Mother", &Gen_Jet_Mother_vec);
                //T->Branch("dR_Jets", &dR_values_jets);
        }

                ////////////GsfElectrons///////////////////
                T->Branch("Gsf_Electrons_Pt", &gsf_electrons_Pt_vec);
                T->Branch("Gsf_Electrons_Eta", &gsf_electrons_Eta_vec);
                T->Branch("Gsf_Electrons_Phi", &gsf_electrons_Phi_vec);
                T->Branch("Gsf_Electrons_Energy", &gsf_electrons_E_vec);
                /////////////RecJet/////////////////////////
                T->Branch("Number_RecoJets", &num_reco_jets, "num_reco_jets/I");
                T->Branch("RecoJet_Pt", &Reco_Jet_Pt_vec);
                T->Branch("RecoJet_Eta", &Reco_Jet_Eta_vec);
                T->Branch("RecoJet_Phi", &Reco_Jet_Phi_vec);
                T->Branch("RecoJet_Energy", &Reco_Jet_E_vec);
                T->Branch("RecoJet_Mother", &Reco_Jet_Mother_vec);

               /* T->Branch("eta_lead", &eta_lead_vec);
                T->Branch("eta_sublead", &eta_sublead_vec);
                T->Branch("eta_lead_brl", &eta_lead_brl);
                T->Branch("eta_lead_ecp", &eta_lead_ecp);
                T->Branch("eta_sublead_brl", &eta_sublead_brl);
                T->Branch("eta_sublead_ecp", &eta_sublead_ecp);
                */
        
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZEE_RecHit_NTuplizer::endJob()
{
}


// Extract the rechits of the SC from the ES layers
/*void ZEE_RecHit_NTuplizer::GetESPlaneRecHits(const reco::SuperCluster& sc, const CaloGeometry* &geo, unsigned int elenum, unsigned int planeIndex) {
        double RawenergyPlane = 0.;
        double pfRawenergyPlane = 0.;
//      if(!_ESRechitsHandle.isValid())
//              return RawenergyPlane;
              
//        if (!sc.preshowerClusters().isAvailable()) //protection for miniAOD
//                break;

	int NumHits = 0;

	const CaloSubdetectorGeometry* ecalESGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalPreshower));


        for(auto iES = sc.preshowerClustersBegin(); iES != sc.preshowerClustersEnd(); ++iES) {//loop over preshower clusters
                const std::vector< std::pair<DetId, float> > hits = (*iES)->hitsAndFractions();
                for(std::vector<std::pair<DetId, float> >::const_iterator rh = hits.begin(); rh != hits.end(); ++rh) { // loop over recHits of the cluster
                        //      std::cout << "print = " << (*iES)->printHitAndFraction(iCount);
                        //      ++iCount;
                        for(ESRecHitCollection::const_iterator esItr = ESRechitsHandle->begin(); esItr != ESRechitsHandle->end(); ++esItr) {//loop over ES rechits to find the one in the cluster
                                ESDetId rhid = ESDetId(esItr->id());
                                if(rhid == (*rh).first) { // found ESrechit
//                                        std::cout << " ES energy = " << esItr->energy() << " pf energy = " << (*rh).second << std::endl;
                                        if((int) rhid.plane() == (int) planeIndex) {
						std::shared_ptr<const CaloCellGeometry> geom = ecalESGeom->getGeometry(rhid);
						Hit_ES_Eta[elenum].push_back( geom->etaPos() );
                                                Hit_ES_Phi[elenum].push_back( geom->phiPos() );
						Hit_ES_X[elenum].push_back( geom->getPosition().x() );
						Hit_ES_Y[elenum].push_back( geom->getPosition().y() );
						Hit_ES_Z[elenum].push_back( geom->getPosition().z() ) ;
						ES_RecHitEn[elenum].push_back(esItr->energy());
//						std::cout << "Preshower" <<std::setprecision(4) << " Eta = " <<geom->etaPos() << " : " <<" Phi = "<< geom->phiPos() << " 3D position" << geom->getPosition().z() << std::endl;
                                                RawenergyPlane += esItr->energy();
                                                pfRawenergyPlane += rh->second;
						NumHits++;
                                        }
                                        break;
                                }
                        }
                }

//		std::cout<<std::endl<<" Number of ES hits in plane 1 = "<<NumHits<<std::endl;
        }

 //       return RawenergyPlane;
}*/
//Initialisation/////
void ZEE_RecHit_NTuplizer::initialisation()
{
        qS =-9999;
        alpQCD =-9999;
        alpQED =-9999;
        idFirst=-9999;
        idSecond=-9999;
        xFirst=-9999.0;
        xSecond=-9999.0;
        scalePDF=-9999.0;



        nMEP=-9999;
        nMEPF=-9999;
        sPID=-9999;
        wt=-9999;
        wtP=-9999;
        hPF=true;
        hDRV=true;
        
        num_gen_jets=0;
        num_reco_jets=0;
}


//Clear tree vectors each time analyze method is called
void ZEE_RecHit_NTuplizer::ClearTreeVectors()
{
	nElectrons_ = 0;
/*	iEta[0].clear();
	iPhi[0].clear();

	
	Hit_ES_Eta[0].clear();
        Hit_ES_Phi[0].clear();
        Hit_ES_X[0].clear();
        Hit_ES_Y[0].clear();
        Hit_ES_Z[0].clear();
	ES_RecHitEn[0].clear();


        Hit_Eta[0].clear();
	Hit_Phi[0].clear();
	Hit_X[0].clear();
	Hit_Y[0].clear();
	Hit_Z[0].clear();
	RecHitEn[0].clear();
	RecHitFrac[0].clear();
	iEta[1].clear();
        iPhi[1].clear();

	Hit_ES_Eta[1].clear();
        Hit_ES_Phi[1].clear();
        Hit_ES_X[1].clear();
        Hit_ES_Y[1].clear();
        Hit_ES_Z[1].clear();
	ES_RecHitEn[1].clear();

	Hit_Eta[1].clear();
        Hit_Phi[1].clear();
        Hit_X[1].clear();
        Hit_Y[1].clear();
        Hit_Z[1].clear();
        RecHitEn[1].clear();
        RecHitFrac[1].clear();
	Ele_pt_.clear();
	Ele_eta_.clear();
	Ele_phi_.clear();
	Ele_energy_.clear();
	Ele_ecal_energy_.clear();
	Ele_ecal_mustache_energy_.clear();		
	Ele_R9.clear();
	Ele_S4.clear();
	Ele_SigIEIE.clear();
	Ele_SigIPhiIPhi.clear();
	Ele_SCEtaW.clear();
	Ele_SCPhiW.clear();
	Ele_CovIEtaIEta.clear();
	Ele_CovIEtaIPhi.clear();
	Ele_ESSigRR.clear();
	Ele_SCRawE.clear();
	Ele_SC_ESEnByRawE.clear();
	Ele_HadOverEm.clear();

	HFEMClust_pt_.clear();
	HFEMClust_eta_.clear();
	HFEMClust_phi_.clear();
	HFEMClust_energy_.clear();
	HFEMClust_eLong1x1_.clear();
	HFEMClust_eShort1x1_.clear();
	HFEMClust_eLong3x3_.clear();
	HFEMClust_eShort3x3_.clear();
	HFEMClust_eLong5x5_.clear();
	HFEMClust_eShort5x5_.clear();
	HFEMClust_e1x1_.clear();
	HFEMClust_e3x3_.clear();
	HFEMClust_e5x5_.clear();
	HFEMClust_eSeL_.clear();
	HFEMClust_eCOREe9_.clear();
	HFEMClust_e9e25_.clear(); 
	HFEMClust_eCore_.clear();
*/
	if(isMC_){

	Ele_Gen_Pt.clear();
	Ele_Gen_Eta.clear();
	Ele_Gen_Phi.clear();
	Ele_Gen_E.clear();
      /*  eta_lead_vec.clear();
        eta_sublead_vec.clear();
        eta_lead_brl.clear();
        eta_lead_ecp.clear();
        eta_sublead_brl.clear();
        eta_sublead_ecp.clear();
        */




		
	Dressed_Ele_Pt.clear();
	Dressed_Ele_Eta.clear();
	Dressed_Ele_Phi.clear();
	Dressed_Ele_E.clear();
        
        dR_0.clear();
        dR_1.clear();
        dR_values.clear();

        Naked_Ele_Pt.clear();
        Naked_Ele_Eta.clear();
        Naked_Ele_Phi.clear();
        Naked_Ele_E.clear();

        Born_Ele_Pt.clear();
        Born_Ele_Eta.clear();
        Born_Ele_Phi.clear();
        Born_Ele_E.clear();

        Z_Truth_Pt.clear();
        Z_Truth_Eta.clear();
        Z_Truth_Phi.clear();
        Z_Truth_E.clear();

        recon_ele_E_vec.clear();
        recon_ele_Pt_vec.clear();
        recon_ele_Phi_vec.clear();
        recon_ele_Eta_vec.clear();

        Gen_Jet_Pt_vec.clear();
        Gen_Jet_Eta_vec.clear();
        Gen_Jet_Phi_vec.clear();
        Gen_Jet_E_vec.clear();
        Gen_Jet_Mother_vec.clear();
        //dR_values_jets.clear();

        //Hadron_Flavour_vec.clear();

        }
        Reco_Jet_Pt_vec.clear();
        Reco_Jet_Eta_vec.clear();
        Reco_Jet_Phi_vec.clear();
        Reco_Jet_E_vec.clear();
        Reco_Jet_Mother_vec.clear();

        gsf_electrons_Pt_vec.clear();
        gsf_electrons_Eta_vec.clear();
        gsf_electrons_Phi_vec.clear();
        gsf_electrons_E_vec.clear();

	/*passMediumId_.clear();
	passTightId_ .clear();
	passMVAMediumId_.clear();

	isTrue_.clear();*/
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZEE_RecHit_NTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZEE_RecHit_NTuplizer);
