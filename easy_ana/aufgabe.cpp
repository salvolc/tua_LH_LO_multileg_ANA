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


#include "/home/salv/eigen/Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::ofstream;
using std::cout;
using std::endl;
using std::string;

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
VectorXd vertex_match(TClonesArray* TCP);


int error=0;

int main(int argc, char const *argv[])
{


	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MadGraph/MG5_aMC_v2_6_0/ExRootAnalysis/libExRootAnalysis.so");
	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MadGraph/Delphes-3.4.1/libDelphes.so");

	
	string fileName = "../samples/tca_RH_ttbar_PY8_ATLASDELPHES_50000.root";

	int nEvents = get_nevents(fileName);
	MatrixXd ev = get_eventdisplay(fileName,nEvents-1);
	speichere("einevent",ev);


	string fileNames[4];
	fileNames[0] = "../samples/tua_LH_ttbar_PY8_ATLASDELPHES_50000.root";
	fileNames[1] = "../samples/tua_RH_ttbar_PY8_ATLASDELPHES_50000.root";
	fileNames[2] = "../samples/tca_LH_ttbar_PY8_ATLASDELPHES_50000.root";
	fileNames[3] = "../samples/tca_RH_ttbar_PY8_ATLASDELPHES_50000.root";

	// string filePPTNames[4];filePPTNames[0] = "data/tua_LH_Photon_PT";filePPTNames[1] = "data/tua_RH_Photon_PT";filePPTNames[2] = "data/tca_LH_Photon_PT";filePPTNames[3] = "data/tca_RH_Photon_PT";
	// string filebJetPTNames[4];	filebJetPTNames[0] = "data/tua_LH_bJet_PT";filebJetPTNames[1] = "data/tua_RH_bJet_PT";filebJetPTNames[2] = "data/tca_LH_bJet_PT";filebJetPTNames[3] = "data/tca_RH_bJet_PT";
	// string fileWPTENames[4];	fileWPTENames[0] = "data/tua_LH_WE_PT";fileWPTENames[1] = "data/tua_RH_WE_PT";fileWPTENames[2] = "data/tca_LH_WE_PT";fileWPTENames[3] = "data/tca_RH_WE_PT";
	// string fileWPTNames[4];		fileWPTNames[0] = "data/tua_LH_W_PT";fileWPTNames[1] = "data/tua_RH_W_PT";fileWPTNames[2] = "data/tca_LH_W_PT";fileWPTNames[3] = "data/tca_RH_W_PT";
	// string fileWMNames[4];		fileWMNames[0] = "data/tua_LH_W_M";fileWMNames[1] = "data/tua_RH_W_M";fileWMNames[2] = "data/tca_LH_W_M";fileWMNames[3] = "data/tca_RH_W_M";
	// string filetPTNames[4];		filetPTNames[0] = "data/tua_LH_t_PT";filetPTNames[1] = "data/tua_RH_t_PT";filetPTNames[2] = "data/tca_LH_t_PT";filetPTNames[3] = "data/tca_RH_t_PT";
	// string filetMNames[4];		filetMNames[0] = "data/tua_LH_t_M";filetMNames[1] = "data/tua_RH_t_M";filetMNames[2] = "data/tca_LH_t_M";filetMNames[3] = "data/tca_RH_t_M";
	

	string fileTPPTNames[4];fileTPPTNames[0] = "data/tua_LH_Photon_PT_truth";fileTPPTNames[1] = "data/tua_RH_Photon_PT_truth";fileTPPTNames[2] = "data/tca_LH_Photon_PT_truth";fileTPPTNames[3] = "data/tca_RH_Photon_PT_truth";
	string fileTPEtaNames[4];fileTPEtaNames[0] = "data/tua_LH_Photon_Eta_truth";fileTPEtaNames[1] = "data/tua_RH_Photon_Eta_truth";fileTPEtaNames[2] = "data/tca_LH_Photon_Eta_truth";fileTPEtaNames[3] = "data/tca_RH_Photon_Eta_truth";
	string fileTPPhiNames[4];fileTPPhiNames[0] = "data/tua_LH_Photon_Phi_truth";fileTPPhiNames[1] = "data/tua_RH_Photon_Phi_truth";fileTPPhiNames[2] = "data/tca_LH_Photon_Phi_truth";fileTPPhiNames[3] = "data/tca_RH_Photon_Phi_truth";
	

	//string fileTWPTNames[4];fileTWPTNames[0] = "data/tua_LH_WBoson_PT_truth";fileTWPTNames[1] = "data/tua_RH_WBoson_PT_truth";fileTWPTNames[2] = "data/tca_LH_WBoson_PT_truth";fileTWPTNames[3] = "data/tca_RH_WBoson_PT_truth";
	//string fileTWEtaNames[4];fileTWEtaNames[0] = "data/tua_LH_WBoson_Eta_truth";fileTWEtaNames[1] = "data/tua_RH_WBoson_Eta_truth";fileTWEtaNames[2] = "data/tca_LH_WBoson_Eta_truth";fileTWEtaNames[3] = "data/tca_RH_WBoson_Eta_truth";
	//string fileTWPhiNames[4];fileTWPhiNames[0] = "data/tua_LH_WBoson_Phi_truth";fileTWPhiNames[1] = "data/tua_RH_WBoson_Phi_truth";fileTWPhiNames[2] = "data/tca_LH_WBoson_Phi_truth";fileTWPhiNames[3] = "data/tca_RH_WBoson_Phi_truth";
	//string fileTWMNames[4];fileTWMNames[0] = "data/tua_LH_WBoson_M_truth";fileTWMNames[1] = "data/tua_RH_WBoson_M_truth";fileTWMNames[2] = "data/tca_LH_WBoson_M_truth";fileTWMNames[3] = "data/tca_RH_WBoson_M_truth";
	

	string fileTtPTNames[4];fileTtPTNames[0] = "data/tua_LH_TopQuark_PT_truth";fileTtPTNames[1] = "data/tua_RH_TopQuark_PT_truth";fileTtPTNames[2] = "data/tca_LH_TopQuark_PT_truth";fileTtPTNames[3] = "data/tca_RH_TopQuark_PT_truth";
	string fileTtEtaNames[4];fileTtEtaNames[0] = "data/tua_LH_TopQuark_Eta_truth";fileTtEtaNames[1] = "data/tua_RH_TopQuark_Eta_truth";fileTtEtaNames[2] = "data/tca_LH_TopQuark_Eta_truth";fileTtEtaNames[3] = "data/tca_RH_TopQuark_Eta_truth";
	string fileTtPhiNames[4];fileTtPhiNames[0] = "data/tua_LH_TopQuark_Phi_truth";fileTtPhiNames[1] = "data/tua_RH_TopQuark_Phi_truth";fileTtPhiNames[2] = "data/tca_LH_TopQuark_Phi_truth";fileTtPhiNames[3] = "data/tca_RH_TopQuark_Phi_truth";
	string fileTtMNames[4];fileTtMNames[0] = "data/tua_LH_TopQuark_M_truth";fileTtMNames[1] = "data/tua_RH_TopQuark_M_truth";fileTtMNames[2] = "data/tca_LH_TopQuark_M_truth";fileTtMNames[3] = "data/tca_RH_TopQuark_M_truth";

	for (int iFile = 3; iFile < 4; ++iFile)
	{
		TFile* file = new TFile(fileNames[iFile].c_str(),"READ");
		TTree* tree = (TTree*)file->Get("Delphes");

		TBranch *bP 		= tree->GetBranch("Particle");
		TClonesArray *TCP 	= 0;
		bP->SetAddress(&TCP);
		bP->GetEntry(0);

		int nPhotons = 0;
		int nMatched = 0;
		for (int iEvent = nEvents-1; iEvent < nEvents; ++iEvent)
		{
			int nParticles = TCP->GetEntries();
			for (int iParticle = 0; iParticle < nParticles; ++iParticle)
			{
				GenParticle * P_Particle = (GenParticle*)TCP->At(iParticle);
				if (P_Particle->PID == 22) 
				{
					nPhotons++;
				}
			}
			VectorXd match = vertex_match(TCP);
			if (match(3)==1){nMatched++;}

		}
		cout << nPhotons << endl;
		cout << nMatched << endl;
		file->Close();
	}



	// for (int iFile = 0; iFile < 4; ++iFile)
	// {
		
	// 	TFile* file = new TFile(fileNames[iFile].c_str(),"READ");
	// 	TTree* tree = (TTree*)file->Get("Delphes");


	// 	TBranch *bP 		= tree->GetBranch("Particle");
	// 	TBranch *bPhoton	= tree->GetBranch("Photon");
	// 	TBranch *bJet   	= tree->GetBranch("Jet");
	// 	TBranch *bElectron	= tree->GetBranch("Electron");
	// 	TBranch *bMuon		= tree->GetBranch("Muon");
	// 	TBranch *bMET		= tree->GetBranch("MissingET");



	// 	TClonesArray *TCP 		= 0;
	// 	TClonesArray *TCPhoton 	= 0;
	// 	TClonesArray *TCJet 	= 0;
	// 	TClonesArray *TCElectron= 0;
	// 	TClonesArray *TCMuon 	= 0;
	// 	TClonesArray *TCMET 	= 0;

	// 	bP->SetAddress(&TCP);
	// 	bPhoton->SetAddress(&TCPhoton);
	// 	bJet->SetAddress(&TCJet);
	// 	bElectron->SetAddress(&TCElectron);
	// 	bMuon->SetAddress(&TCMuon);
	// 	bMET->SetAddress(&TCMET);


	// 	bP->GetEntry(0);
	// 	bPhoton->GetEntry(0);
	// 	bJet->GetEntry(0);
	// 	bElectron->GetEntry(0);
	// 	bMuon->GetEntry(0);
	// 	bMET->GetEntry(0);

	// 	bool truth = true;


	// 	//int nTElectrons	= get_nTruth(fileNames[iFile].c_str(),11); int iTElectrons	= 0;
	// 	//int nTWBoson	= get_nTruth(fileNames[iFile].c_str(),24); int iTWBoson		= 0;
	// 	//cout << iTWBoson << endl;
	// 	VectorXd VTPhotonPT = VectorXd::Zero(nEvents);VectorXd VTPhotonEta = VectorXd::Zero(nEvents);VectorXd VTPhotonPhi = VectorXd::Zero(nEvents);
	// 	//VectorXd VTWBosonPT = VectorXd::Zero(nTWBoson);VectorXd VTWBosonEta = VectorXd::Zero(nTWBoson);VectorXd VTWBosonPhi = VectorXd::Zero(nTWBoson);VectorXd VTWBosonM = VectorXd::Zero(nTWBoson);
	// 	VectorXd VTTopQuarkPT = VectorXd::Zero(nEvents);VectorXd VTTopQuarkEta = VectorXd::Zero(nEvents);VectorXd VTTopQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTTopQuarkM = VectorXd::Zero(nEvents);



	// 	for (int iEvent = 0; iEvent < nEvents && truth; ++iEvent)
	// 	{
	// 		bP->GetEntry(iEvent);

	// 		int nParticles = TCP->GetEntries();

	// 		VectorXd i_matched = vertex_match(TCP);
	// 		if(i_matched(3)==0){continue;}

	// 		GenParticle* P_Particle = (GenParticle*)TCP->At(i_matched(0));


	// 		VTPhotonPT(iEvent) = P_Particle->PT;
	// 		VTPhotonEta(iEvent) = P_Particle->Eta;
	// 		VTPhotonPhi(iEvent) = P_Particle->Phi;

	// 		P_Particle = (GenParticle*)TCP->At(i_matched(1));
	// 		VTTopQuarkPT(iEvent) = P_Particle->PT;
	// 		VTTopQuarkEta(iEvent) = P_Particle->Eta;
	// 		VTTopQuarkPhi(iEvent) = P_Particle->Phi;
	// 		VTTopQuarkM(iEvent) = P_Particle->Mass;
	// 		// for (int iparticle = 0; iparticle < nParticles; ++iparticle)
	// 		// {
	// 		// 	if (abs(P_Particle->PID) == 24){
	// 		// 		GenParticle* P_Particle = (GenParticle*)TCP->At(iparticle);
	// 		// 		VTWBosonPT(iTWBoson) = P_Particle->PT;
	// 		// 		VTWBosonEta(iTWBoson) = P_Particle->Eta;
	// 		// 		VTWBosonPhi(iTWBoson) = P_Particle->Phi;
	// 		// 		VTWBosonM(iTWBoson) = P_Particle->Mass;
	// 		// 		iTWBoson++;}
	// 		// }
	// 		if (iEvent%10000==0)
	// 		{
	// 		cout << iEvent << endl;
	// 		}
	// 	}


	// 	speichere(fileTPPTNames[iFile],VTPhotonPT);
	// 	speichere(fileTPEtaNames[iFile],VTPhotonEta);
	// 	speichere(fileTPPhiNames[iFile],VTPhotonPhi);

	// 	speichere(fileTWPTNames[iFile],VTWBosonPT);
	// 	speichere(fileTWEtaNames[iFile],VTWBosonEta);
	// 	speichere(fileTWPhiNames[iFile],VTWBosonPhi);
	// 	speichere(fileTWMNames[iFile],VTWBosonM);

	// 	speichere(fileTtPTNames[iFile],VTTopQuarkPT);
	// 	speichere(fileTtEtaNames[iFile],VTTopQuarkEta);
	// 	speichere(fileTtPhiNames[iFile],VTTopQuarkPhi);
	// 	speichere(fileTtMNames[iFile],VTTopQuarkM);


	// 	file->Close();
	// }

	return 0;

}




MatrixXd get_eventdisplay(string fileName, int event){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");

	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();
	MatrixXd display = MatrixXd::Zero(numberOfParticles,9);


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
	MatrixXd display = MatrixXd::Zero(nPIDs,9);
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


VectorXd vertex_match(TClonesArray* TCP){
	int nPart = TCP->GetEntries();

	bool matched1=false;
	bool matched2=false;

	MatrixXd particle_infos(nPart,7);

	int nt = 0;int i_nt = 0; 

	int nb = 0;int i_nb = 0; 


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
		if(abs(particle_infos(iPart,0)) == 5)nb++;
	}


	MatrixXd t_infos = MatrixXd::Zero(nt,7);
	MatrixXd b_infos = MatrixXd::Zero(nb,7);

	for (int iPart = 0; iPart < nPart; ++iPart)
	{
		if(abs(particle_infos(iPart,0)) == 6){t_infos.row(i_nt) = particle_infos.row(iPart);i_nt++;}
		if(abs(particle_infos(iPart,0)) == 5){b_infos.row(i_nb) = particle_infos.row(iPart);i_nb++;}
	}

	cout << t_infos << endl;

	VectorXd i_vertex_matched = VectorXd::Zero(4);
	bool matched_t = false;
	GenParticle* P_ParticleHelp = (GenParticle*)TCP->At(0);

	for (int it = 0; it < i_nt; ++it)
	{
		if ( (t_infos(it,3)==22) && (abs(t_infos(it,4))==4  ||  abs(t_infos(it,4))==2 ) )   
		{
			i_vertex_matched(1) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(1));
			i_vertex_matched(0) = P_ParticleHelp->D1;
			i_vertex_matched(2) = P_ParticleHelp->D2;
			matched_t = true;
		}
		if ( (t_infos(it,4)==22) && (abs(t_infos(it,3))==4  ||  abs(t_infos(it,3))==2 ) )
		{
			i_vertex_matched(1) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(1));
			i_vertex_matched(0) = P_ParticleHelp->D2;
			i_vertex_matched(2) = P_ParticleHelp->D1;
			matched_t = true;
		}
		
	}



	if(matched_t){i_vertex_matched(3)=1;}else{i_vertex_matched(3)=0;}
	return i_vertex_matched;
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


