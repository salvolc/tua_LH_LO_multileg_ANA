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
VectorXd vertex_match(TClonesArray* TCP);


int error=0;

int main(int argc, char const *argv[])
{


	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/libExRootAnalysis.so");
	gSystem->Load("/home/salv/Dokumente/Masterarbeit/Delphes/libDelphes.so");

	
	string fileName = "samples/tua_LH_decaymode_commmadecay_LO.root";
	int nEvents = get_nevents(fileName);
	MatrixXd ev = get_eventdisplay(fileName,1);
	speichere("einevent",ev);


	string fileNames[3];
	fileNames[0] = "samples/tua_LH_decaymode_commmadecay_LO.root";
	fileNames[1] = "samples/tua_LH_interference_LO_multileg_tua.root";
	fileNames[2] = "samples/tua_LH_production_mulileg_tua_LO.root";

	string fileTPPTNames[3];fileTPPTNames[0] = "data/dec_Photon_PT_truth";fileTPPTNames[1] = "data/int_Photon_PT_truth";fileTPPTNames[2] = "data/pro_Photon_PT_truth";
	string fileTPEtaNames[3];fileTPEtaNames[0] = "data/dec_Photon_Eta_truth";fileTPEtaNames[1] = "data/int_Photon_Eta_truth";fileTPEtaNames[2] = "data/pro_Photon_Eta_truth";
	string fileTPPhiNames[3];fileTPPhiNames[0] = "data/dec_Photon_Phi_truth";fileTPPhiNames[1] = "data/int_Photon_Phi_truth";fileTPPhiNames[2] = "data/pro_Photon_Phi_truth";
	

	string fileTWPTNames[3];fileTWPTNames[0] = "data/dec_WBoson_PT_truth";fileTWPTNames[1] = "data/int_WBoson_PT_truth";fileTWPTNames[2] = "data/pro_WBoson_PT_truth";
	string fileTWEtaNames[3];fileTWEtaNames[0] = "data/dec_WBoson_Eta_truth";fileTWEtaNames[1] = "data/int_WBoson_Eta_truth";fileTWEtaNames[2] = "data/pro_WBoson_Eta_truth";
	string fileTWPhiNames[3];fileTWPhiNames[0] = "data/dec_WBoson_Phi_truth";fileTWPhiNames[1] = "data/int_WBoson_Phi_truth";fileTWPhiNames[2] = "data/pro_WBoson_Phi_truth";
	string fileTWMNames[3];fileTWMNames[0] = "data/dec_WBoson_M_truth";fileTWMNames[1] = "data/int_WBoson_M_truth";fileTWMNames[2] = "data/pro_WBoson_M_truth";
	

	string fileTtPTNames[3];fileTtPTNames[0] = "data/dec_TopQuark_PT_truth";fileTtPTNames[1] = "data/int_TopQuark_PT_truth";fileTtPTNames[2] = "data/pro_TopQuark_PT_truth";
	string fileTtEtaNames[3];fileTtEtaNames[0] = "data/dec_TopQuark_Eta_truth";fileTtEtaNames[1] = "data/int_TopQuark_Eta_truth";fileTtEtaNames[2] = "data/pro_TopQuark_Eta_truth";
	string fileTtPhiNames[3];fileTtPhiNames[0] = "data/dec_TopQuark_Phi_truth";fileTtPhiNames[1] = "data/int_TopQuark_Phi_truth";fileTtPhiNames[2] = "data/pro_TopQuark_Phi_truth";
	string fileTtMNames[3];fileTtMNames[0] = "data/dec_TopQuark_M_truth";fileTtMNames[1] = "data/int_TopQuark_M_truth";fileTtMNames[2] = "data/pro_TopQuark_M_truth";

	string fileTbPTNames[3];fileTbPTNames[0] = "data/dec_BQuark_PT_truth";fileTbPTNames[1] = "data/int_BQuark_PT_truth";fileTbPTNames[2] = "data/pro_BQuark_PT_truth";
	string fileTbEtaNames[3];fileTbEtaNames[0] = "data/dec_BQuark_Eta_truth";fileTbEtaNames[1] = "data/int_BQuark_Eta_truth";fileTbEtaNames[2] = "data/pro_BQuark_Eta_truth";
	string fileTbPhiNames[3];fileTbPhiNames[0] = "data/dec_BQuark_Phi_truth";fileTbPhiNames[1] = "data/int_BQuark_Phi_truth";fileTbPhiNames[2] = "data/pro_BQuark_Phi_truth";
	string fileTbMNames[3];fileTbMNames[0] = "data/dec_BQuark_M_truth";fileTbMNames[1] = "data/int_BQuark_M_truth";fileTbMNames[2] = "data/pro_BQuark_M_truth";

	string fileTuPTNames[3];fileTuPTNames[0] = "data/dec_UQuark_PT_truth";fileTuPTNames[1] = "data/int_UQuark_PT_truth";fileTuPTNames[2] = "data/pro_UQuark_PT_truth";
	string fileTuEtaNames[3];fileTuEtaNames[0] = "data/dec_UQuark_Eta_truth";fileTuEtaNames[1] = "data/int_UQuark_Eta_truth";fileTuEtaNames[2] = "data/pro_UQuark_Eta_truth";
	string fileTuPhiNames[3];fileTuPhiNames[0] = "data/dec_UQuark_Phi_truth";fileTuPhiNames[1] = "data/int_UQuark_Phi_truth";fileTuPhiNames[2] = "data/pro_UQuark_Phi_truth";
	string fileTuMNames[3];fileTuMNames[0] = "data/dec_UQuark_M_truth";fileTuMNames[1] = "data/int_UQuark_M_truth";fileTuMNames[2] = "data/pro_UQuark_M_truth";

	string fileWeightNames[3];fileWeightNames[0] = "data/dec_weight_truth";fileWeightNames[1] = "data/int_weight_truth";fileWeightNames[2] = "data/pro_weight_truth";


	for (int iFile = 0; iFile < 3; ++iFile)
	{
		
		TFile* file = new TFile(fileNames[iFile].c_str(),"READ");
		TTree* tree = (TTree*)file->Get("Delphes");


		TBranch *bP 		= tree->GetBranch("Particle");
		TClonesArray *TCP 	= 0;
		bP->SetAddress(&TCP);
		bP->GetEntry(0);

		int nEvents = get_nevents(fileNames[iFile].c_str());
		MatrixXd ev = get_eventdisplay(fileNames[iFile].c_str(),1);

		bool truth = true;

		VectorXd VTPhotonPT = VectorXd::Zero(nEvents);VectorXd VTPhotonEta = VectorXd::Zero(nEvents);VectorXd VTPhotonPhi = VectorXd::Zero(nEvents);
		VectorXd VTTopQuarkPT = VectorXd::Zero(nEvents);VectorXd VTTopQuarkEta = VectorXd::Zero(nEvents);VectorXd VTTopQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTTopQuarkM = VectorXd::Zero(nEvents);
		VectorXd VTWBosonPT = VectorXd::Zero(nEvents);VectorXd VTWBosonEta = VectorXd::Zero(nEvents);VectorXd VTWBosonPhi = VectorXd::Zero(nEvents);VectorXd VTWBosonM = VectorXd::Zero(nEvents);
		VectorXd VTBQuarkPT = VectorXd::Zero(nEvents);VectorXd VTBQuarkEta = VectorXd::Zero(nEvents);VectorXd VTBQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTBQuarkM = VectorXd::Zero(nEvents);
		VectorXd VTUQuarkPT = VectorXd::Zero(nEvents);VectorXd VTUQuarkEta = VectorXd::Zero(nEvents);VectorXd VTUQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTUQuarkM = VectorXd::Zero(nEvents);


		for (int iEvent = 0; iEvent < nEvents && truth; ++iEvent)
		{
			bP->GetEntry(iEvent);

			int nParticles = TCP->GetEntries();

/*			string a;
			cin >> a;*/
			VectorXd i_matched = VectorXd::Zero(4); i_matched = vertex_match(TCP);
			if(i_matched(6)==0){continue;}

			GenParticle* P_Particle = (GenParticle*)TCP->At(i_matched(0));
			VTPhotonPT(iEvent) = P_Particle->PT;
			VTPhotonEta(iEvent) = P_Particle->Eta;
			VTPhotonPhi(iEvent) = P_Particle->Phi;

			P_Particle = (GenParticle*)TCP->At(i_matched(1));
			VTTopQuarkPT(iEvent) = P_Particle->PT;
			VTTopQuarkEta(iEvent) = P_Particle->Eta;
			VTTopQuarkPhi(iEvent) = P_Particle->Phi;
			VTTopQuarkM(iEvent) = P_Particle->Mass;

			P_Particle = (GenParticle*)TCP->At(i_matched(4));
			VTBQuarkPT(iEvent) = P_Particle->PT;
			VTBQuarkEta(iEvent) = P_Particle->Eta;
			VTBQuarkPhi(iEvent) = P_Particle->Phi;
			VTBQuarkM(iEvent) = P_Particle->Mass;

			P_Particle = (GenParticle*)TCP->At(i_matched(2));
			VTUQuarkPT(iEvent) = P_Particle->PT;
			VTUQuarkEta(iEvent) = P_Particle->Eta;
			VTUQuarkPhi(iEvent) = P_Particle->Phi;
			VTUQuarkM(iEvent) = P_Particle->Mass;

			P_Particle = (GenParticle*)TCP->At(i_matched(5));
			VTWBosonPT(iEvent) = P_Particle->PT;
			VTWBosonEta(iEvent) = P_Particle->Eta;
			VTWBosonPhi(iEvent) = P_Particle->Phi;
			VTWBosonM(iEvent) = P_Particle->Mass;
			if(iEvent%2000==0){ladebalken(iEvent*(iFile+1),nEvents*3);}
			//cout << "Event " << iEvent << " Photon PT " << VTPhotonPT(iEvent) << " TopQuarkPT " << VTTopQuarkPT(iEvent) << endl;
		}


		speichere(fileTPPTNames[iFile],VTPhotonPT);
		speichere(fileTPEtaNames[iFile],VTPhotonEta);
		speichere(fileTPPhiNames[iFile],VTPhotonPhi);

		speichere(fileTtPTNames[iFile],VTTopQuarkPT);
		speichere(fileTtEtaNames[iFile],VTTopQuarkEta);
		speichere(fileTtPhiNames[iFile],VTTopQuarkPhi);
		speichere(fileTtMNames[iFile],VTTopQuarkM);

		speichere(fileTbPTNames[iFile],VTBQuarkPT);
		speichere(fileTbEtaNames[iFile],VTBQuarkEta);
		speichere(fileTbPhiNames[iFile],VTBQuarkPhi);
		speichere(fileTbMNames[iFile],VTBQuarkM);

		speichere(fileTuPTNames[iFile],VTBQuarkPT);
		speichere(fileTuEtaNames[iFile],VTBQuarkEta);
		speichere(fileTuPhiNames[iFile],VTBQuarkPhi);
		speichere(fileTuMNames[iFile],VTBQuarkM);

		speichere(fileTWPTNames[iFile],VTWBosonPT);
		speichere(fileTWEtaNames[iFile],VTWBosonEta);
		speichere(fileTWPhiNames[iFile],VTWBosonPhi);
		speichere(fileTWMNames[iFile],VTWBosonM);


		file->Close();
	}

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

	//cout << particle_infos.block(0,0,50,7) << endl;

	VectorXd i_vertex_matched = VectorXd::Zero(8);
	bool matched_t = false; 
	bool matched_u = false; i_nu=i_nu/2;
	bool matched_b = false; i_ngam=i_ngam/2;
	bool matched_g = false; i_ng=i_ng/2; 
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
		if ( (t_infos(it,4)==24) && abs(t_infos(it,3))==5)
		{
			i_vertex_matched(3) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(3));
			i_vertex_matched(4) = P_ParticleHelp->D1;
			i_vertex_matched(5) = P_ParticleHelp->D2;
			matched_b = true;
		}
		if ( (t_infos(it,3)==24) && abs(t_infos(it,4))==5)
		{
			i_vertex_matched(3) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(3));
			i_vertex_matched(4) = P_ParticleHelp->D2;
			i_vertex_matched(5) = P_ParticleHelp->D1;
			matched_b = true;
		}
	}

	if(!matched_t){//Falls t nicht so gefunden wird haben wir wohl ein up was nach t gamma geht
	for (int iu = 0; iu < i_nu; ++iu)
	{
		if ( (u_infos(iu,3)==22) && (abs(u_infos(iu,4))==6) )   
		{
			i_vertex_matched(2) = u_infos(iu,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(2));
			i_vertex_matched(0) = P_ParticleHelp->D1;
			i_vertex_matched(1) = P_ParticleHelp->D2;
			matched_u = true;
		}
		if ( (u_infos(iu,4)==22) && (abs(u_infos(iu,3))==6) )
		{
			i_vertex_matched(2) = u_infos(iu,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(2));
			i_vertex_matched(0) = P_ParticleHelp->D2;
			i_vertex_matched(1) = P_ParticleHelp->D1;
			matched_u = true;
		}
	}
	}
	if(!matched_t && !matched_u){
	for (int ig = 0; ig < i_ng; ++ig)
	{
/*		if((abs(g_infos(ig,3)==6))){cout << "hallo1" << endl; cout << abs(g_infos(ig,3)) << endl;}
		if((abs(g_infos(ig,4))==4) || (abs(g_infos(ig,4))==2)){cout << "hallo2" << endl; cout << abs(g_infos(ig,4)) << endl;}
*/		if ( (abs(g_infos(ig,3))==6) && ( (abs(g_infos(ig,4))==4) || (abs(g_infos(ig,4))==2) ) )   
		{
			i_vertex_matched(7) = g_infos(ig,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(7));
			i_vertex_matched(1) = P_ParticleHelp->D1;
			i_vertex_matched(2) = P_ParticleHelp->D2;
			matched_g = true;
		}
		if ( abs((g_infos(ig,4))==6) && ( (abs(g_infos(ig,3))==4) || (abs(g_infos(ig,3))==2) ) )
		{
			i_vertex_matched(7) = g_infos(ig,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(7));
			i_vertex_matched(1) = P_ParticleHelp->D2;
			i_vertex_matched(2) = P_ParticleHelp->D1;
			matched_g = true;
		}
	}
	for (int igam = 0; igam < i_ngam; ++igam)
	{
		if ( (abs(gam_infos(igam,1))==21  ||  abs(gam_infos(igam,1))==9 ) || (abs(gam_infos(igam,2))==21  ||  abs(gam_infos(igam,2))==9 ) )   
		{
			i_vertex_matched(0) = gam_infos(igam,5);
		}
	}
	}
	if(matched_u && error != 0)cout << "Up gemached" << endl;
	if(matched_t && error != 0)cout << "Top gemached" << endl;
	if(matched_g && error != 0)cout << "Gluon gemached" << endl;
	if(matched_t || matched_u || matched_g){i_vertex_matched(6)=1;}else{i_vertex_matched(6)=0;}
	if(matched_b){i_vertex_matched(7)=1;}else{i_vertex_matched(7)=0;}
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


