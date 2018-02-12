#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "RConfig.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h" 
#include "TMultiGraph.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TSystem.h"
#include <iostream>
#include "TLorentzVector.h"
#include "math.h"
#include "ExRootClasses.h"
#include "ExRootTreeReader.h"
#include "TLeaf.h"
#include "classes/DelphesClasses.h"
#include <cmath>


#include "/home/salv/eigen/Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::ofstream;
using std::cout;
using std::endl;
using std::string;
using std::cin;

int get_nevents(string fileName);
MatrixXd get_eventdisplay(string fileName, int event);
void ladebalken(int i, int max);
void speichere(std::string name, MatrixXd data);
void speichere(std::string name, VectorXd data);
VectorXd get_event_MET(string fileName, int event);
TLorentzVector get_event_MET_neutrino(string fileName, int event);
MatrixXd get_eventdisplay_particle(string fileName, int event, int PID);
int get_numberOfPhotons(string fileName);
int get_numberOfJets(string fileName);
int get_numberOfbJets(string fileName);
TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[], double mass);
TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[]);
int get_nTruth(string fileName, int PID);
VectorXd vertex_match(TClonesArray* TCP, int iEvent);
double get_weight(TBranch* branchE, int iEvent);



int error=0;
bool debug=true;

int main(int argc, char const *argv[])
{

	double mmu = 0.105658;
	double me  = 0.000511;
	double mw  = 80.385;

	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/libExRootAnalysis.so");
	gSystem->Load("/home/salv/Dokumente/Masterarbeit/Delphes/libDelphes.so");

	
	string fileName = "samples/tua_LH_decay_LO_wlep_PY8_DELATL_50.root";
	int nEvents = get_nevents(fileName);
	MatrixXd ev = get_eventdisplay(fileName,1);
	speichere("einevent",ev);


	string fileNames[3];
	fileNames[0] = "samples/tua_LH_decay_LO_wlep_PY8_DELATL_50.root";
	fileNames[1] = "samples/tua_LH_interference_wlep_tja_ta_PY8_DELATL_50.root";
	fileNames[2] = "samples/tua_LH_production_wlep_tja_ta_PY8_DELATL_50.root";

	string filePPTNames[3];filePPTNames[0] = "data/dec_Photon_PT";filePPTNames[1] = "data/int_Photon_PT";filePPTNames[2] = "data/pro_Photon_PT";
	string filePEtaNames[3];filePEtaNames[0] = "data/dec_Photon_Eta";filePEtaNames[1] = "data/int_Photon_Eta";filePEtaNames[2] = "data/pro_Photon_Eta";
	string filePPhiNames[3];filePPhiNames[0] = "data/dec_Photon_Phi";filePPhiNames[1] = "data/int_Photon_Phi";filePPhiNames[2] = "data/pro_Photon_Phi";
	

	string fileWPTNames[3];fileWPTNames[0] = "data/dec_WBoson_PT";fileWPTNames[1] = "data/int_WBoson_PT";fileWPTNames[2] = "data/pro_WBoson_PT";
	string fileWEtaNames[3];fileWEtaNames[0] = "data/dec_WBoson_Eta";fileWEtaNames[1] = "data/int_WBoson_Eta";fileWEtaNames[2] = "data/pro_WBoson_Eta";
	string fileWPhiNames[3];fileWPhiNames[0] = "data/dec_WBoson_Phi";fileWPhiNames[1] = "data/int_WBoson_Phi";fileWPhiNames[2] = "data/pro_WBoson_Phi";
	string fileWMNames[3];fileWMNames[0] = "data/dec_WBoson_M";fileWMNames[1] = "data/int_WBoson_M";fileWMNames[2] = "data/pro_WBoson_M";
	

	string filetPTNames[3];filetPTNames[0] = "data/dec_TopQuark_PT";filetPTNames[1] = "data/int_TopQuark_PT";filetPTNames[2] = "data/pro_TopQuark_PT";
	string filetEtaNames[3];filetEtaNames[0] = "data/dec_TopQuark_Eta";filetEtaNames[1] = "data/int_TopQuark_Eta";filetEtaNames[2] = "data/pro_TopQuark_Eta";
	string filetPhiNames[3];filetPhiNames[0] = "data/dec_TopQuark_Phi";filetPhiNames[1] = "data/int_TopQuark_Phi";filetPhiNames[2] = "data/pro_TopQuark_Phi";
	string filetMNames[3];filetMNames[0] = "data/dec_TopQuark_M";filetMNames[1] = "data/int_TopQuark_M";filetMNames[2] = "data/pro_TopQuark_M";

	string filebPTNames[3];filebPTNames[0] = "data/dec_BQuark_PT";filebPTNames[1] = "data/int_BQuark_PT";filebPTNames[2] = "data/pro_BQuark_PT";
	string filebEtaNames[3];filebEtaNames[0] = "data/dec_BQuark_Eta";filebEtaNames[1] = "data/int_BQuark_Eta";filebEtaNames[2] = "data/pro_BQuark_Eta";
	string filebPhiNames[3];filebPhiNames[0] = "data/dec_BQuark_Phi";filebPhiNames[1] = "data/int_BQuark_Phi";filebPhiNames[2] = "data/pro_BQuark_Phi";
	string filebMNames[3];filebMNames[0] = "data/dec_BQuark_M";filebMNames[1] = "data/int_BQuark_M";filebMNames[2] = "data/pro_BQuark_M";

	string fileuPTNames[3];fileuPTNames[0] = "data/dec_UQuark_PT";fileuPTNames[1] = "data/int_UQuark_PT";fileuPTNames[2] = "data/pro_UQuark_PT";
	string fileuEtaNames[3];fileuEtaNames[0] = "data/dec_UQuark_Eta";fileuEtaNames[1] = "data/int_UQuark_Eta";fileuEtaNames[2] = "data/pro_UQuark_Eta";
	string fileuPhiNames[3];fileuPhiNames[0] = "data/dec_UQuark_Phi";fileuPhiNames[1] = "data/int_UQuark_Phi";fileuPhiNames[2] = "data/pro_UQuark_Phi";
	string fileuMNames[3];fileuMNames[0] = "data/dec_UQuark_M";fileuMNames[1] = "data/int_UQuark_M";fileuMNames[2] = "data/pro_UQuark_M";

	string filebJPTNames[3];filebJPTNames[0] = "data/dec_bJet_PT";filebJPTNames[1] = "data/int_bJet_PT";filebJPTNames[2] = "data/pro_bJet_PT";
	string filebJEtaNames[3];filebJEtaNames[0] = "data/dec_bJet_Eta";filebJEtaNames[1] = "data/int_bJet_Eta";filebJEtaNames[2] = "data/pro_bJet_Eta";
	string filebJPhiNames[3];filebJPhiNames[0] = "data/dec_bJet_Phi";filebJPhiNames[1] = "data/int_bJet_Phi";filebJPhiNames[2] = "data/pro_bJet_Phi";
	string filebJMNames[3];filebJMNames[0] = "data/dec_bJet_M";filebJMNames[1] = "data/int_bJet_M";filebJMNames[2] = "data/pro_bJet_M";

	string fileJPTNames[3];fileJPTNames[0] = "data/dec_Jet_PT";fileJPTNames[1] = "data/int_Jet_PT";fileJPTNames[2] = "data/pro_Jet_PT";
	string fileJEtaNames[3];fileJEtaNames[0] = "data/dec_Jet_Eta";fileJEtaNames[1] = "data/int_Jet_Eta";fileJEtaNames[2] = "data/pro_Jet_Eta";
	string fileJPhiNames[3];fileJPhiNames[0] = "data/dec_Jet_Phi";fileJPhiNames[1] = "data/int_Jet_Phi";fileJPhiNames[2] = "data/pro_Jet_Phi";
	string fileJMNames[3];fileJMNames[0] = "data/dec_Jet_M";fileJMNames[1] = "data/int_Jet_M";fileJMNames[2] = "data/pro_Jet_M";
	
	string fileLPTNames[3];fileLPTNames[0] = "data/dec_Lepton_PT";fileLPTNames[1] = "data/int_Lepton_PT";fileLPTNames[2] = "data/pro_Lepton_PT";
	string fileLEtaNames[3];fileLEtaNames[0] = "data/dec_Lepton_Eta";fileLEtaNames[1] = "data/int_Lepton_Eta";fileLEtaNames[2] = "data/pro_Lepton_Eta";
	string fileLPhiNames[3];fileLPhiNames[0] = "data/dec_Lepton_Phi";fileLPhiNames[1] = "data/int_Lepton_Phi";fileLPhiNames[2] = "data/pro_Lepton_Phi";
	string fileLMNames[3];fileLMNames[0] = "data/dec_Lepton_M";fileLMNames[1] = "data/int_Lepton_M";fileLMNames[2] = "data/pro_Lepton_M";

	string fileMETNames[3];fileMETNames[0] = "data/dec_MET";fileMETNames[1] = "data/int_MET";fileMETNames[2] = "data/pro_MET";
	string fileMETEtaNames[3];fileMETEtaNames[0] = "data/dec_MET_Eta";fileMETEtaNames[1] = "data/int_MET_Eta";fileMETEtaNames[2] = "data/pro_MET_Eta";
	string fileMETPhiNames[3];fileMETPhiNames[0] = "data/dec_MET_Phi";fileMETPhiNames[1] = "data/int_MET_Phi";fileMETPhiNames[2] = "data/pro_MET_Phi";
	
	string fileWeightNames[3];fileWeightNames[0] = "data/dec_weight";fileWeightNames[1] = "data/int_weight";fileWeightNames[2] = "data/pro_weight";

	string filetPRNames[3];filetPRNames[0] = "data/dec_TopQuark_Photon_R";filetPRNames[1] = "data/int_TopQuark_Photon_R";filetPRNames[2] = "data/pro_TopQuark_Photon_R";
	string filetbRNames[3];filetbRNames[0] = "data/dec_TopQuark_BQuark_R";filetbRNames[1] = "data/int_TopQuark_BQuark_R";filetbRNames[2] = "data/pro_TopQuark_BQuark_R";
	string filetWRNames[3];filetWRNames[0] = "data/dec_TopQuark_WBoson_R";filetWRNames[1] = "data/int_TopQuark_WBoson_R";filetWRNames[2] = "data/pro_TopQuark_WBoson_R";
	string filetjRNames[3];filetjRNames[0] = "data/dec_TopQuark_LeadingJet_R";filetjRNames[1] = "data/int_TopQuark_LeadingJet_R";filetjRNames[2] = "data/pro_TopQuark_LeadingJet_R";

	string filetPMNames[3];filetPMNames[0] = "data/dec_TopQuark_Photon_M";filetPMNames[1] = "data/int_TopQuark_Photon_M";filetPMNames[2] = "data/pro_TopQuark_Photon_M";
	string filetbMNames[3];filetbMNames[0] = "data/dec_TopQuark_BQuark_M";filetbMNames[1] = "data/int_TopQuark_BQuark_M";filetbMNames[2] = "data/pro_TopQuark_BQuark_M";
	string filetWMNames[3];filetWMNames[0] = "data/dec_TopQuark_WBoson_M";filetWMNames[1] = "data/int_TopQuark_WBoson_M";filetWMNames[2] = "data/pro_TopQuark_WBoson_M";
	string filetjMNames[3];filetjMNames[0] = "data/dec_TopQuark_LeadingJet_M";filetjMNames[1] = "data/int_TopQuark_LeadingJet_M";filetjMNames[2] = "data/pro_TopQuark_LeadingJet_M";

	string filePtRNames[3];filePtRNames[0] = "data/dec_Photon_TopQuark_R";filePtRNames[1] = "data/int_Photon_TopQuark_R";filePtRNames[2] = "data/pro_Photon_TopQuark_R";
	string filePbRNames[3];filePbRNames[0] = "data/dec_Photon_BQuark_R";filePbRNames[1] = "data/int_Photon_BQuark_R";filePbRNames[2] = "data/pro_Photon_BQuark_R";
	string filePWRNames[3];filePWRNames[0] = "data/dec_Photon_WBoson_R";filePWRNames[1] = "data/int_Photon_WBoson_R";filePWRNames[2] = "data/pro_Photon_WBoson_R";
	string filePjRNames[3];filePjRNames[0] = "data/dec_Photon_LeadingJet_R";filePjRNames[1] = "data/int_Photon_LeadingJet_R";filePjRNames[2] = "data/pro_Photon_LeadingJet_R";

	string filePtMNames[3];filePtMNames[0] = "data/dec_Photon_TopQuark_M";filePtMNames[1] = "data/int_Photon_TopQuark_M";filePtMNames[2] = "data/pro_Photon_TopQuark_M";
	string filePbMNames[3];filePbMNames[0] = "data/dec_Photon_BQuark_M";filePbMNames[1] = "data/int_Photon_BQuark_M";filePbMNames[2] = "data/pro_Photon_BQuark_M";
	string filePWMNames[3];filePWMNames[0] = "data/dec_Photon_WBoson_M";filePWMNames[1] = "data/int_Photon_WBoson_M";filePWMNames[2] = "data/pro_Photon_WBoson_M";
	string filePjMNames[3];filePjMNames[0] = "data/dec_Photon_LeadingJet_M";filePjMNames[1] = "data/int_Photon_LeadingJet_M";filePjMNames[2] = "data/pro_Photon_LeadingJet_M";


	TFile* files[3];

	//#pragma omp parallel for
	for (int iFile = 0; iFile < 3; ++iFile)
	{
		
		//TFile* file = new TFile(fileNames[iFile].c_str(),"READ");
		//TTree* tree = (TTree*)file->Get("Delphes");
		files[iFile] 		= new TFile(fileNames[iFile].c_str(),"READ");
		TTree* 		tree 	= (TTree*)files[iFile]->Get("Delphes");

		
		TBranch *bE 		= tree->GetBranch("Event");
		TBranch *bP 		= tree->GetBranch("Particle");
		TBranch *bPhoton	= tree->GetBranch("Photon");
		TBranch *bJet   	= tree->GetBranch("Jet");
		TBranch *bElectron	= tree->GetBranch("Electron");
		TBranch *bMuon		= tree->GetBranch("Muon");
		TBranch *bMET		= tree->GetBranch("MissingET");


		TClonesArray *TCP 		= 0;
		TClonesArray *TCPhoton 	= 0;
		TClonesArray *TCJet 	= 0;
		TClonesArray *TCElectron= 0;
		TClonesArray *TCMuon 	= 0;
		TClonesArray *TCMET 	= 0;
		TClonesArray *TCE 		= 0;


		bP->SetAddress(&TCP);
		bPhoton->SetAddress(&TCPhoton);
		bJet->SetAddress(&TCJet);
		bElectron->SetAddress(&TCElectron);
		bMuon->SetAddress(&TCMuon);
		bMET->SetAddress(&TCMET);
		bE->SetAddress(&TCE);


		int nEvents = get_nevents(fileNames[iFile].c_str());

		VectorXd VPhotonPT = VectorXd::Zero(nEvents);VectorXd VPhotonEta = VectorXd::Zero(nEvents);VectorXd VPhotonPhi = VectorXd::Zero(nEvents);
		VectorXd VTopQuarkPT = VectorXd::Zero(nEvents);VectorXd VTopQuarkEta = VectorXd::Zero(nEvents);VectorXd VTopQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTopQuarkM = VectorXd::Zero(nEvents);
		VectorXd VWBosonPT = VectorXd::Zero(nEvents);VectorXd VWBosonEta = VectorXd::Zero(nEvents);VectorXd VWBosonPhi = VectorXd::Zero(nEvents);VectorXd VWBosonM = VectorXd::Zero(nEvents);
		VectorXd VBQuarkPT = VectorXd::Zero(nEvents);VectorXd VBQuarkEta = VectorXd::Zero(nEvents);VectorXd VBQuarkPhi = VectorXd::Zero(nEvents);VectorXd VBQuarkM = VectorXd::Zero(nEvents);
		VectorXd VUQuarkPT = VectorXd::Zero(nEvents);VectorXd VUQuarkEta = VectorXd::Zero(nEvents);VectorXd VUQuarkPhi = VectorXd::Zero(nEvents);VectorXd VUQuarkM = VectorXd::Zero(nEvents);
		
		VectorXd VbJetPT = VectorXd::Zero(nEvents);VectorXd VbJetEta = VectorXd::Zero(nEvents);VectorXd VbJetPhi = VectorXd::Zero(nEvents);VectorXd VbJetM = VectorXd::Zero(nEvents);
		VectorXd VJetPT = VectorXd::Zero(nEvents);VectorXd VJetEta = VectorXd::Zero(nEvents);VectorXd VJetPhi = VectorXd::Zero(nEvents);VectorXd VJetM = VectorXd::Zero(nEvents);
		VectorXd VLeptonPT = VectorXd::Zero(nEvents);VectorXd VLeptonEta = VectorXd::Zero(nEvents);VectorXd VLeptonPhi = VectorXd::Zero(nEvents);VectorXd VLeptonM = VectorXd::Zero(nEvents);
		VectorXd VMET = VectorXd::Zero(nEvents);VectorXd VMETEta = VectorXd::Zero(nEvents);VectorXd VMETPhi = VectorXd::Zero(nEvents);

		VectorXd VTopQuark_Photon_R = VectorXd::Zero(nEvents);VectorXd VTopQuark_bQuark_R = VectorXd::Zero(nEvents);VectorXd VTopQuark_WBoson_R = VectorXd::Zero(nEvents);VectorXd VTopQuark_LedJet_R = VectorXd::Zero(nEvents);
		VectorXd VTopQuark_Photon_M = VectorXd::Zero(nEvents);VectorXd VTopQuark_bQuark_M = VectorXd::Zero(nEvents);VectorXd VTopQuark_WBoson_M = VectorXd::Zero(nEvents);VectorXd VTopQuark_LedJet_M = VectorXd::Zero(nEvents);

		VectorXd VPhoton_TopQuark_R = VectorXd::Zero(nEvents);VectorXd VPhoton_bQuark_R = VectorXd::Zero(nEvents);VectorXd VPhoton_WBoson_R = VectorXd::Zero(nEvents);VectorXd VPhoton_LedJet_R = VectorXd::Zero(nEvents);
		VectorXd VPhoton_TopQuark_M = VectorXd::Zero(nEvents);VectorXd VPhoton_bQuark_M = VectorXd::Zero(nEvents);VectorXd VPhoton_WBoson_M = VectorXd::Zero(nEvents);VectorXd VPhoton_LedJet_M = VectorXd::Zero(nEvents);

		VectorXd VWeight = VectorXd::Zero(nEvents);
		
		//int NPhotons = 0;
		//int NbJets = 0;
		//int numberOfPhotons = get_numberOfPhotons(fileNames[iFile].c_str());
		//int numberOfbJets = get_numberOfbJets(fileNames[iFile].c_str());

		int iE 	 = 0;
		int iELe = 0;
		int iEta = 0;
		int iMas = 0;
		int iWei = 0;
		for (int iEvent = 0; iEvent < nEvents; ++iEvent)
		{
			//string a; cin >> a;
			bP->GetEntry(iEvent);
			bPhoton->GetEntry(iEvent);
			bJet->GetEntry(iEvent);
			bElectron->GetEntry(iEvent);
			bMuon->GetEntry(iEvent);
			bMET->GetEntry(iEvent);
			bE->GetEntry(iEvent);

			//####################################################
			//####################Weight##########################
			//####################################################

			HepMCEvent *E_Event = (HepMCEvent*)TCE->At(0);
			VWeight(iEvent) = E_Event->Weight;

			//####################################################
			//##############Photon Zeugs##########################
			//####################################################

			int nPhotons = TCPhoton->GetEntries();
			TLorentzVector TV_Photon;
			for (int iPhoton = 0; iPhoton < nPhotons; ++iPhoton)
			{
				Photon *P_Photon = (Photon*)TCPhoton->At(iPhoton);
				if (P_Photon->PT > TV_Photon.Pt())
				{
					TV_Photon.SetPtEtaPhiE(P_Photon->PT,P_Photon->Eta,P_Photon->Phi,P_Photon->E);
				}
			}


			//####################################################
			//##############LJet Zeugs############################
			//####################################################
			//##############BJet Zeugs############################
			//####################################################

			TLorentzVector TV_bJet;
			TLorentzVector TV_Jet;
			int nbtags = 0;
			int nJets = TCJet->GetEntries();
			for (int iJets = 0; iJets < nJets; ++iJets)
			{
			 	Jet	*P_Jet 	= (Jet*)TCJet->At(iJets);
			 	if(P_Jet->BTag == 1){
			 		nbtags++;
			 		if(TV_bJet.Pt() < P_Jet->PT){
						TV_bJet.SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
	 				}
			 	}	else if(TV_Jet.Pt() < P_Jet->PT){
					TV_Jet.SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
			 	}
			}

			//####################################################
			//#########################Elektronen#################
			//####################################################

			int nElectrons = TCElectron->GetEntries();
			int nMuons = TCMuon->GetEntries();
			TLorentzVector TV_Lepton;

			for (int iElectron = 0; iElectron < nElectrons; ++iElectron)
			{
				Electron* P_Electron = (Electron*)TCElectron->At(iElectron);
				if (TV_Lepton.Pt() < P_Electron->PT)
				{
					TV_Lepton.SetPtEtaPhiM((P_Electron->PT),(P_Electron->Eta),(P_Electron->Phi),me);
				}
			}
			for (int iMuon = 0; iMuon < nMuons; ++iMuon)
			{
				Muon* P_Muon = (Muon*)TCMuon->At(iMuon);
				if (TV_Lepton.Pt() < P_Muon->PT)
				{
					TV_Lepton.SetPtEtaPhiM((P_Muon->PT),(P_Muon->Eta),(P_Muon->Phi),mmu);
				}
			}

			//####################################################
			//#########################MET und W##################
			//####################################################
			TLorentzVector TV_METp;
			TLorentzVector TV_METm;
			TLorentzVector TV_WBoson;
			MissingET* P_MET = (MissingET*)TCMET->At(0);

			double arg = ((mw*mw)/(2*(P_MET->MET)*(TV_Lepton.Pt())))+cos((P_MET->Phi)-(TV_Lepton.Phi()));
			if (arg < 1){arg = 1;}
			double METEtap = acosh(arg)+VLeptonEta(iEvent);
			double METEtam = -acosh(arg)+VLeptonEta(iEvent);
			
			//cout << METEtap << " " << METEtam << " " << VMETEta(iEvent) << endl;
			TV_METp.SetPtEtaPhiE(P_MET->MET,METEtap,P_MET->Phi,P_MET->MET);
			TV_METm.SetPtEtaPhiE(P_MET->MET,METEtam,P_MET->Phi,P_MET->MET);

			if(fabs(mw-(TV_Lepton+TV_METm).M())>fabs(mw-(TV_Lepton+TV_METp).M()))
			{
				TV_WBoson = TV_Lepton+TV_METp;
			} else {
				TV_WBoson = TV_Lepton+TV_METm;
			}


			//####################################################
			//##################Top Reco Zeugs E##################
			//####################################################

			TLorentzVector TV_Top;
			TV_Top = TV_WBoson+TV_bJet;

			iE++;
			if(nElectrons+nMuons==0){continue;}
			iELe++;
			if(fabs(TV_Jet.Eta())>2.5 || fabs(TV_bJet.Eta())>2.5){continue;}
			iEta++;
			if(TV_WBoson.M() <= 0 || TV_Top.M() <= 0){continue;}
			iMas++;
			if(VWeight(iEvent)==0){continue;}
			iWei++;

			//####################################################
			//#######################Saving#######################
			//####################################################

			VPhotonPT(iEvent) 	= 	TV_Photon.Pt();
			VPhotonEta(iEvent) 	= 	TV_Photon.Eta();
			VPhotonPhi(iEvent) 	= 	TV_Photon.Phi();					

			VTopQuarkM(iEvent)	=	TV_Top.M();
			VTopQuarkPT(iEvent)	=	TV_Top.Pt();
			VTopQuarkEta(iEvent)=	TV_Top.Eta();
			VTopQuarkPhi(iEvent)=	TV_Top.Phi();
			
			VMET(iEvent) 		= 	P_MET->MET;
			VMETPhi(iEvent) 	= 	P_MET->Phi;
			VMETEta(iEvent) 	= 	P_MET->Eta;

			VWBosonPT(iEvent)	=	TV_WBoson.Pt();
			VWBosonEta(iEvent)	=	TV_WBoson.Eta();
			VWBosonPhi(iEvent)	=	TV_WBoson.Phi();
			VWBosonM(iEvent)	=	TV_WBoson.M();

			VLeptonPT(iEvent) 	=	TV_Lepton.Pt();
			VLeptonEta(iEvent)	=	TV_Lepton.Eta();
			VLeptonPhi(iEvent)	=	TV_Lepton.Phi();
			VLeptonM(iEvent)	=	TV_Lepton.M();

			VbJetPT(iEvent)		=	TV_bJet.Pt();
			VbJetEta(iEvent)	=	TV_bJet.Eta();
			VbJetPhi(iEvent)	=	TV_bJet.Phi();
			VbJetM(iEvent)		=	TV_bJet.M();

			VJetPT(iEvent)		=	TV_Jet.Pt();
			VJetEta(iEvent)		=	TV_Jet.Eta();
			VJetPhi(iEvent)		=	TV_Jet.Phi();
			VJetM(iEvent)		=	TV_Jet.M();

			VTopQuark_Photon_R(iEvent)	= TV_Top.DeltaR(TV_Photon);
			VTopQuark_bQuark_R(iEvent)	= TV_Top.DeltaR(TV_bJet);
			VTopQuark_WBoson_R(iEvent)	= TV_Top.DeltaR(TV_WBoson);
			VTopQuark_LedJet_R(iEvent)	= TV_Top.DeltaR(TV_Jet);
			VTopQuark_Photon_M(iEvent)	= (TV_Top+TV_Photon).M();
			VTopQuark_bQuark_M(iEvent)	= (TV_Top+TV_bJet).M();
			VTopQuark_WBoson_M(iEvent)	= (TV_Top+TV_WBoson).M();
			VTopQuark_LedJet_M(iEvent)	= (TV_Top+TV_Jet).M();
			VPhoton_TopQuark_R(iEvent)	= TV_Photon.DeltaR(TV_Top);
			VPhoton_bQuark_R(iEvent)	= TV_Photon.DeltaR(TV_bJet);
			VPhoton_WBoson_R(iEvent)	= TV_Photon.DeltaR(TV_WBoson);
			VPhoton_LedJet_R(iEvent)	= TV_Photon.DeltaR(TV_Jet);
			VPhoton_TopQuark_M(iEvent)	= (TV_Photon+TV_Top).M();
			VPhoton_bQuark_M(iEvent)	= (TV_Photon+TV_bJet).M();
			VPhoton_WBoson_M(iEvent)	= (TV_Photon+TV_WBoson).M();
			VPhoton_LedJet_M(iEvent)	= (TV_Photon+TV_Jet).M();
		}
		speichere(filePPTNames[iFile], VPhotonPT);
		speichere(filePEtaNames[iFile],VPhotonEta);
		speichere(filePPhiNames[iFile],VPhotonPhi);
		speichere(filetPTNames[iFile], VTopQuarkPT);
		speichere(filetEtaNames[iFile],VTopQuarkEta);
		speichere(filetPhiNames[iFile],VTopQuarkPhi);
		speichere(filetMNames[iFile],  VTopQuarkM);
		speichere(filebJPTNames[iFile], VbJetPT);
		speichere(filebJEtaNames[iFile],VbJetEta);
		speichere(filebJPhiNames[iFile],VbJetPhi);
		speichere(filebJMNames[iFile],  VbJetM);
		speichere(fileJPTNames[iFile], VJetPT);
		speichere(fileJEtaNames[iFile],VJetEta);
		speichere(fileJPhiNames[iFile],VJetPhi);
		speichere(fileJMNames[iFile],  VJetM);
		speichere(fileWPTNames[iFile], VWBosonPT);
		speichere(fileWEtaNames[iFile],VWBosonEta);
		speichere(fileWPhiNames[iFile],VWBosonPhi);
		speichere(fileWMNames[iFile],  VWBosonM);
		speichere(filetPRNames[iFile], VTopQuark_Photon_R);
		speichere(filetbRNames[iFile], VTopQuark_bQuark_R);
		speichere(filetWRNames[iFile], VTopQuark_WBoson_R);
		speichere(filetjRNames[iFile], VTopQuark_LedJet_R);
		speichere(filetPMNames[iFile], VTopQuark_Photon_M);
		speichere(filetbMNames[iFile], VTopQuark_bQuark_M);
		speichere(filetWMNames[iFile], VTopQuark_WBoson_M);
		speichere(filetjMNames[iFile], VTopQuark_LedJet_M);
		speichere(filePtRNames[iFile], VPhoton_TopQuark_R);
		speichere(filePbRNames[iFile], VPhoton_bQuark_R);
		speichere(filePWRNames[iFile], VPhoton_WBoson_R);
		speichere(filePjRNames[iFile], VPhoton_LedJet_R);
		speichere(filePtMNames[iFile], VPhoton_TopQuark_M);
		speichere(filePbMNames[iFile], VPhoton_bQuark_M);
		speichere(filePWMNames[iFile], VPhoton_WBoson_M);
		speichere(filePjMNames[iFile], VPhoton_LedJet_M);
		speichere(fileWeightNames[iFile], VWeight);


		ofstream dat((fileNames[iFile]+"_Cut_Effs.txt").c_str());
		dat.is_open();
		dat << "Events: " << iE << "\n";
		dat << "Eff Lepton Cut: " << ((float(iELe)/float(iE))*100) << "%" << "\n";
		dat << "Eff Eta Cut: " << ((float(iEta)/float(iE))*100) << "%" << "\n";
		dat << "Eff Mass Cut: " << ((float(iMas)/float(iE))*100) << "%" << "\n";
		dat << "Eff Weight Cut: " << ((float(iWei)/float(iE))*100) << "%" << "\n";
		dat.close();

		cout << iE << endl;

		files[iFile]->Close();
		//file->Close();
	}

	return 0;

}



VectorXd vertex_match(TClonesArray* TCP,int iFile){
//	int npart_sel = 50;
	int nPart = TCP->GetEntries();
	MatrixXd particle_infos(nPart,7);

	int nt = 0;int i_nt = 0;
	int nu = 0;int i_nu = 0;
	int ng = 0;int i_ng = 0;
	int ngam = 0;int i_ngam = 0;

	for (int iPart = 0; iPart < nPart; ++iPart)
	{

		GenParticle* P_Particle = (GenParticle*)TCP->At(iPart);
		GenParticle* P_ParticleHelp = (GenParticle*)TCP->At(iPart);

		particle_infos(iPart,0) = P_Particle->PID;
		if(abs(P_Particle->M1) == 1){particle_infos(iPart,1) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->M1));
			particle_infos(iPart,1) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->M2) == 1){particle_infos(iPart,2) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->M2));
			particle_infos(iPart,2) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->D1) == 1){particle_infos(iPart,3) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->D1));
			particle_infos(iPart,3) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->D2) == 1){particle_infos(iPart,4) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->D2));
			particle_infos(iPart,4) = P_ParticleHelp->PID;
		}
		particle_infos(iPart,5) = iPart;
		particle_infos(iPart,6) = abs(P_Particle->PT);
		if(abs(particle_infos(iPart,0)) == 6)nt++;
		if(abs(particle_infos(iPart,0)) == 2 || abs(particle_infos(iPart,0)) == 4)nu++;
		if(abs(particle_infos(iPart,0)) == 21 || abs(particle_infos(iPart,0)) == 9)ng++;
		if(abs(particle_infos(iPart,0)) == 22)ngam++;
	}


	MatrixXd t_infos = MatrixXd::Zero(nt,7);
	MatrixXd u_infos = MatrixXd::Zero(nu,7);
	MatrixXd g_infos = MatrixXd::Zero(ng,7);
	MatrixXd gam_infos = MatrixXd::Zero(ngam,7);

	for (int iPart = 0; iPart < nPart; ++iPart)
	{
		if(abs(particle_infos(iPart,0)) == 6){t_infos.row(i_nt) = particle_infos.row(iPart);i_nt++;}
		if(abs(particle_infos(iPart,0)) == 22){gam_infos.row(i_ngam) = particle_infos.row(iPart);i_ngam++;}
		if(abs(particle_infos(iPart,0)) == 2 || abs(particle_infos(iPart,0)) == 4){u_infos.row(i_nu) = particle_infos.row(iPart);i_nu++;}
		if(abs(particle_infos(iPart,0)) == 21 || abs(particle_infos(iPart,0)) == 9){g_infos.row(i_ng) = particle_infos.row(iPart);i_ng++;}
	}


	VectorXd i_vertex_matched = VectorXd::Zero(10);
	bool matched_t = false; 
	bool matched_u = false;   //i_nu=i_nu/2;
	bool matched_b = false;   //i_ngam=i_ngam/2;
	bool matched_gam = false; //i_ng=i_ng/2; 
	GenParticle* P_ParticleHelp = (GenParticle*)TCP->At(0);
	for (int it = 0; it < i_nt && (!matched_t || !matched_b); ++it)
	{
		if ((abs(t_infos(it,1))==2 || abs(t_infos(it,2))==2 || abs(t_infos(it,3))==2 || abs(t_infos(it,4))==2) && !matched_t)
		{
			i_vertex_matched(1) = t_infos(it,5);
			matched_t = true;
		}
/*		if ((abs(t_infos(it,1))==22 || abs(t_infos(it,2))==22 || abs(t_infos(it,3))==22 || abs(t_infos(it,4))==22) && !matched_t)
		{
			i_vertex_matched(1) = t_infos(it,5);
			matched_t = true;
		}*/
		if (((abs(t_infos(it,1))==21 || abs(t_infos(it,2))==21) && !matched_t) && iFile > 0)
		{
			i_vertex_matched(1) = t_infos(it,5);
			matched_t = true;
		}
		if ( (abs(t_infos(it,4))==24) && (abs(t_infos(it,3))==5))
		{
			i_vertex_matched(3) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(3));
			i_vertex_matched(4) = P_ParticleHelp->D1;
			i_vertex_matched(5) = P_ParticleHelp->D2;
			matched_b = true;
		}
		if ( (abs(t_infos(it,3))==24) && (abs(t_infos(it,4))==5))
		{
			i_vertex_matched(3) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(3));
			i_vertex_matched(4) = P_ParticleHelp->D2;
			i_vertex_matched(5) = P_ParticleHelp->D1;
			matched_b = true;
		}
	}

	for (int iu = 0; iu < i_nu && !matched_u; ++iu)
	{
		if(abs(u_infos(iu,1))==6 || abs(u_infos(iu,2))==6 || abs(u_infos(iu,3))==6 || abs(u_infos(iu,4))==6)
		{
			i_vertex_matched(2) = u_infos(iu,5);
			matched_u = true;
		}
	}
	
	for (int igam = 0; igam < i_ngam && !matched_gam; ++igam)
	{
		if ( (abs(gam_infos(igam,1))==21  ||  abs(gam_infos(igam,1))==9 ) || (abs(gam_infos(igam,2))==21  ||  abs(gam_infos(igam,2))==9 ) )   
		{
			i_vertex_matched(0) = gam_infos(igam,5);
			matched_gam = true;
		}
		if (abs(gam_infos(igam,1))==2 || abs(gam_infos(igam,2))==2 || abs(gam_infos(igam,3))==2 || abs(gam_infos(igam,4))==2)
		{
			i_vertex_matched(0) = gam_infos(igam,5);
			matched_gam = true;
		}
		if (abs(gam_infos(igam,1))==6 || abs(gam_infos(igam,2))==6 || abs(gam_infos(igam,3))==6 || abs(gam_infos(igam,4))==6)
		{
			i_vertex_matched(0) = gam_infos(igam,5);
			matched_gam = true;
		}
	}
	i_vertex_matched(6) = int(matched_gam);	i_vertex_matched(7) = int(matched_t);	i_vertex_matched(8) = int(matched_u);	i_vertex_matched(9) = int(matched_b);
	
	//cout << i_vertex_matched(6) << " Photon" << endl;
	//cout << i_vertex_matched(7) << " Top" << endl;
	//cout << i_vertex_matched(8) << " Up" << endl;
	//cout << i_vertex_matched(9) << " B" << endl;
	//if(debug)cout << particle_infos << endl;//.block(0,0,100,7) << endl;
	//if(debug && matched_u)cout << "Up gemached" << endl;
	//if(debug && matched_t)cout << "Top gemached" << endl;
	//if(debug && matched_gam)cout << "Gamma gemached" << endl;
	//if(debug && matched_b)cout << "b+W gemached" << endl;
	//if(debug){cout << "Der Vertex Vektor:"<< endl;cout << i_vertex_matched << endl;}


	return i_vertex_matched;
}



MatrixXd get_eventdisplay(string fileName, int event){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");

	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();
	MatrixXd display = MatrixXd::Zero(numberOfParticles,10);


    for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{

		GenParticle *P = (GenParticle*)TCP->At(ipart);
		display(ipart,0)=P->PID;
		display(ipart,1)=P->Px;
		display(ipart,2)=P->Py;
		display(ipart,3)=P->Pz;
		display(ipart,4)=P->PT;
		display(ipart,5)=P->Mass;
		display(ipart,6)=P->E;
		display(ipart,7)=P->Eta;
		display(ipart,8)=P->Phi;
		display(ipart,9)=ipart;
	}

	file->Close();
	return display;
}

MatrixXd get_eventdisplay_particle(string fileName, int event, int PID){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	int nPIDs = 0;
	for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{
	   	GenParticle *P = (GenParticle*)TCP->At(ipart);
		if (P->PID == PID)nPIDs++;
	}
	if(nPIDs == 0){cout << "Kein Teilchen gefunden!" << endl; return MatrixXd::Zero(1,1);}
	MatrixXd display = MatrixXd::Zero(nPIDs,12);
	nPIDs = 0;
    for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{
	   	GenParticle *P = (GenParticle*)TCP->At(ipart);
		if (P->PID == PID){
		display(nPIDs,0)=P->PID;
		display(nPIDs,1)=P->Px;
		display(nPIDs,2)=P->Py;
		display(nPIDs,3)=P->Pz;
		display(nPIDs,4)=P->PT;
		display(nPIDs,5)=P->Mass;
		display(nPIDs,6)=P->E;
		display(nPIDs,7)=P->Eta;
		display(nPIDs,8)=P->Phi;
		display(nPIDs,9)=ipart;
		display(nPIDs,10)=P->D1;
		display(nPIDs,11)=P->D2;
		nPIDs++;
		}
	}

	file->Close();
	return display;
}

int get_nTruth(string fileName, int PID){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(0);  
	int nEvents = branchP->GetEntries();

	int nPIDs = 0;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		branchP->GetEntry(iEvent);  
		int nParticles = TCP->GetEntries();
		for (int ipart = 0; ipart < nParticles; ++ipart)
		{
		   	GenParticle *P = (GenParticle*)TCP->At(ipart);
			if(abs(P->PID) == PID)nPIDs++;
		}
	}

	file->Close();
	return nPIDs;
}


int get_nevents(string fileName){
	TFile* h_file = new TFile(fileName.c_str(),"READ");
	TTree* h_tree = (TTree*)h_file->Get("Delphes");
	int numberOfEntries = h_tree->GetEntries();
	h_file->Close();
	return numberOfEntries;

}

double get_weight(TBranch* branchE, int iEvent){
	TClonesArray *TCE = 0;
	branchE->SetAddress(&TCE);
	branchE->GetEntry(iEvent);  
	HepMCEvent *E = (HepMCEvent*)TCE->At(0);
	return E->Weight;
}


VectorXd get_event_MET(string fileName, int event){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	VectorXd MET = VectorXd::Zero(3);

	for (int iparticle = 0; iparticle < numberOfParticles; ++iparticle)
	{
		GenParticle *P = (GenParticle*)TCP->At(iparticle);
		MET(0)=MET(0)-P->Px;MET(1)=MET(1)-P->Py;
	}
	MET(2)=sqrt(pow(MET(0),2)+pow(MET(1),2));

	file->Close();
	return MET;
}


TLorentzVector get_event_MET_neutrino(string fileName, int event){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	TLorentzVector MET;

	for (int iparticle = 0; iparticle < numberOfParticles; ++iparticle)
	{
		GenParticle *P = (GenParticle*)TCP->At(iparticle);
		if (P->PID == 12 || P->PID == -12 || P->PID == 14 || P->PID == -14 || P->PID == 18 || P->PID == -18)
		{
			TLorentzVector nu(P->Px,P->Py,P->Pz,P->E);
			MET = MET+nu;
		}
	}
	MET.SetE(MET.Vect().Mag());

	file->Close();
	return MET;
}



int get_numberOfPhotons(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bPhoton = tree->GetBranch("Photon");
	TClonesArray* TCPhoton = 0;
	bPhoton->SetAddress(&TCPhoton);
	int nEvents = bPhoton->GetEntries();

	int nPhotons = 0;

	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bPhoton->GetEntry(iEvent);
		nPhotons += TCPhoton->GetEntries();
	}
	file->Close();
	return nPhotons;
}

int get_numberOfJets(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bJet = tree->GetBranch("Jet");
	TClonesArray* TCJet = 0;
	bJet->SetAddress(&TCJet);
	int nEvents = bJet->GetEntries();

	int nJet = 0;

	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bJet->GetEntry(iEvent);
		nJet += TCJet->GetEntries();
	}
	file->Close();
	return nJet;
}

int get_numberOfbJets(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bJet = tree->GetBranch("Jet");
	TClonesArray* TCJet = 0;
	bJet->SetAddress(&TCJet);
	int nEvents = bJet->GetEntries();
	int nbJet = 0;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bJet->GetEntry(iEvent);
		int nJet = TCJet->GetEntries();
		for (int iJet = 0; iJet < nJet; ++iJet)
		{
			Jet* P_Jet = (Jet*)TCJet->At(iJet);
			if(P_Jet->BTag == 1)nbJet++;
		}
	}
	file->Close();
	return nbJet;
}


void ladebalken(int i, int max){
	double progress = (1.*i)/(1.*max)*100;
	#pragma omp critical
	std::cout << "\rSchritt " << i << " von " << max << " geschafft! " << "Also " << progress << " %";
	if(i == max-1 || i==max){std::cout << "\rSchritt " << max << " von " << max << " geschafft. Fertig!" << std::endl;}
	
	return;
}


void speichere(std::string name, MatrixXd data){
	cout << ("\nSpeichere Datei " +name+ ".txt ab:").c_str() << endl;
	ofstream dat((name+".txt").c_str());
	dat.is_open();
	dat << data << "\n";
	dat.close();
	cout << "Datei " +name+ ".txt abgespeichert!\n" << endl;
	return;
}



void speichere(std::string name, VectorXd data){
	cout << ("\nSpeichere Datei " +name+ ".txt ab:").c_str() << endl;
	ofstream dat((name+".txt").c_str());
	dat.is_open();
	dat << data << "\n";
	dat.close();
	cout << "Datei " +name+ ".txt abgespeichert!\n" << endl;
	return;
}




TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[], double mass)
{
	if(nVec == 1)
	{
		return vecs[0];
	}
	else if(nVec == 2)
	{
		return vecs[0]+vecs[1];
	}
	else if(nVec == 3)
	{
		TLorentzVector perms[3];
		double perm_m[3];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();
		perms[2] = vecs[1]+vecs[2];perm_m[2]=perms[2].M();
		if (fabs(perm_m[0]-mass) > fabs(perm_m[1]))
		{
			if (fabs(perm_m[1]-mass) > fabs(perm_m[2]))
			{
				return perms[2];
			}else
			{
				return perms[1];
			}
		}
		if (fabs(perm_m[0]-mass) > fabs(perm_m[2]))
		{
			return perms[2];
		}
		return perms[0];
	}
	else if(nVec == 4)
	{
		TLorentzVector perms[6];
		double perm_m[6];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();perm_m[0]=fabs(perm_m[0]-mass);
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();perm_m[1]=fabs(perm_m[1]-mass);
		perms[2] = vecs[0]+vecs[3];perm_m[2]=perms[2].M();perm_m[2]=fabs(perm_m[2]-mass);
		perms[3] = vecs[1]+vecs[2];perm_m[3]=perms[3].M();perm_m[3]=fabs(perm_m[3]-mass);
		perms[4] = vecs[1]+vecs[3];perm_m[4]=perms[4].M();perm_m[4]=fabs(perm_m[4]-mass);
		perms[5] = vecs[2]+vecs[3];perm_m[5]=perms[5].M();perm_m[5]=fabs(perm_m[5]-mass);
		for (int i = 0; i < 6; ++i)
		{
			if((perm_m[i] < perm_m[(i+1)%6]) && (perm_m[i] < perm_m[(i+2)%6]) && (perm_m[i] < perm_m[(i+3)%6]) && (perm_m[i] < perm_m[(i+4)%6]) && (perm_m[i] < perm_m[(i+5)%6]))
			{
				return perms[i];
			}
		}
	}
	else if(nVec == 5)
	{
		TLorentzVector perms[10];
		double perm_m[10];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();perm_m[0]=fabs(perm_m[0]-mass);
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();perm_m[1]=fabs(perm_m[1]-mass);
		perms[2] = vecs[0]+vecs[3];perm_m[2]=perms[2].M();perm_m[2]=fabs(perm_m[2]-mass);
		perms[3] = vecs[0]+vecs[4];perm_m[3]=perms[3].M();perm_m[3]=fabs(perm_m[3]-mass);
		perms[4] = vecs[1]+vecs[2];perm_m[4]=perms[4].M();perm_m[4]=fabs(perm_m[4]-mass);
		perms[5] = vecs[1]+vecs[3];perm_m[5]=perms[5].M();perm_m[5]=fabs(perm_m[5]-mass);
		perms[6] = vecs[1]+vecs[4];perm_m[6]=perms[6].M();perm_m[6]=fabs(perm_m[6]-mass);
		perms[7] = vecs[2]+vecs[3];perm_m[7]=perms[7].M();perm_m[7]=fabs(perm_m[7]-mass);
		perms[8] = vecs[2]+vecs[4];perm_m[8]=perms[8].M();perm_m[8]=fabs(perm_m[8]-mass);
		perms[9] = vecs[3]+vecs[4];perm_m[9]=perms[9].M();perm_m[9]=fabs(perm_m[9]-mass);
		for (int i = 0; i < 10; ++i)
		{
			if((perm_m[i] < perm_m[(i+1)%10]) && (perm_m[i] < perm_m[(i+2)%10]) && (perm_m[i] < perm_m[(i+3)%10]) && (perm_m[i] < perm_m[(i+4)%10]) && (perm_m[i] < perm_m[(i+5)%10]) && (perm_m[i] < perm_m[(i+6)%10]) && (perm_m[i] < perm_m[(i+7)%10]) && (perm_m[i] < perm_m[(i+8)%10]) && (perm_m[i] < perm_m[(i+9)%10]))
			{
				return perms[i];
			}
		}
	}
	else{
		cout << "ERROR: Not supported number of Jets" << endl;
		return vecs[0];
	}
}


TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[])
{
	if(error < 2){error++;
	cout << "INFO: No Mass was given choosing highest PT" << endl;}
	if (nVec < 2)
	{	
		if(error < 2){error++;
		cout << "INFO: Only one Vector given -> Outout = Input" << endl;}
		return vecs[0];
	}

	TLorentzVector reco[2];
	TLorentzVector swap;
	reco[0] = vecs[0];
	reco[1] = vecs[1];
	for (int i = 2; i < nVec; ++i)
	{
		if (reco[0].Pt() < reco[1].Pt())
		{
			swap = reco[1]; reco[1]=reco[0]; reco[0]=swap;
		}
		if (vecs[i].Pt() > reco[1].Pt())
		{
			reco[1] = vecs[i];
		}
	}
	swap = reco[0]+reco[1];
	return swap;
}


