#include "TROOT.h"  // For ROOT types like Double_t
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include <iostream>
#include <cmath>
#include "TH2D.h"
#include <vector>
#include <stdexcept>

void DM_manual(){

	
//////////////////////////////////////////////////////////////////////////////////////////////////
// Binding term values and ranges at density ratio of 5, recommended resolutions		//
// 1MeV: 0 to 2200 hexaquarks, res of 1 							//
// 10keV: 2.3*10^6 to  3.0*10^6 res of 100							//				       
//////////////////////////////////////////////////////////////////////////////////////////////////



  
          TLorentzVector d_m1, d_n1;
          TLorentzVector *d_m2 = new TLorentzVector();
          TLorentzVector *d_n2 = new TLorentzVector();

        if (d_m2 == nullptr) {
                d_m2 = new TLorentzVector();  // Allocate memory if nullptr
        }

        if (d_n2 == nullptr) {
                d_n2 = new TLorentzVector();  // Allocate memory if nullptr
        }


Double_t density_ratio = 3.0;

Double_t min_hex = 0;//15.4*TMath::Power(10,6);//1.945084*TMath::Power(10,12);
Double_t max_hex = 2400;//15.8*TMath::Power(10,6);//1.945094*TMath::Power(10,12);

	//2500*TMath::Power(10,6);//1.945094*TMath::Power(10,12);

Double_t res = 10.0;//*TMath::Power(10,4);//2.0*TMath::Power(10,3);//0.1*(max_hex_val / min_hex_val); //depending on range can increase

int bin = int((max_hex - min_hex) / res);
if ((bin <= 0)|| bin >=100) bin = 200; // Default to a reasonable value if calculation fails

Double_t volume_term = 0.001; //Volume term in GeV
 TString title = Form("Mass-Velocity Phase Space, Binding Strength (%.1e GeV), Density ratio (%.2f)",volume_term, density_ratio);
 ////

//Histograms

 TH1D *beta_dist = new TH1D("beta_dist","Beta Dist", bin, 0.001, 0.06);//0.0001, 0.01);
 
//TF1 * bind = new TF1("bind"   , " 0.001*(x+1)*(x)-0.001*0.73*TMath::Power((3),1/3)*(x+1)*(x)*TMath::Power((x+1),-1/3)", 1000, 10000);
 
//Angular Dependence
//TH2D* h_d_m2=new TH2D("h_d_m2","Angular dependence of Incoming D*_m",200,0,1,200,0,180);
//TH2D* h_d_n2=new TH2D("h_d_n2","Angular dependence of Scattered D*_n",200,0,1,200,0,180);

//Energy/velocity distribution 
// TH1F *v_d_n2 =new TH1F("v_d_n2","Energy Velocity phase space",200,0,1);
// TH1D *v_d_m2 =new TH1D("v_d_m2","Energy Velocity phase space for Scattered DM particle ",300,2.5,17.5);
 TH2D *bind_1D =new TH2D("bind_1D","1D Plot", bin, min_hex, max_hex, 100, -3,3);

// Mass velocity phase space 
TH2D *M_N = new TH2D("M_N", title, bin, min_hex,max_hex, bin, min_hex,max_hex);
//50,1000,TMath::Power(10,19),50,3,19
//set min 1*1-10 and max 1 for beta after fill before plot  
//16 ,min_hex_val,max_hex_val, 16, min_hex_val,max_hex_val);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Creating particles and terms
//loop over m and n for varying masses
for (double k = min_hex; k<max_hex; k = k+res){
	//double_t n=TMath::Power(10,(k));
	double n = round(k);
	for (double j = min_hex; j<max_hex; j = j+res){
	//	double_t m=TMath::Power(10,(j));
		double m = round(j);

		//Binding energy formula for n1 and m1 incoming, n2 and m2 outgoing
          Double_t bind_e_m1 = volume_term*(m)*(m-1) - 0.001*0.73*TMath::Power((density_ratio),1.0/3)*(m)*(m-1)*TMath::Power((m),-1.0/3);
          Double_t bind_e_n1 = volume_term*(n)*(n-1) - 0.001*0.73*TMath::Power((density_ratio),1.0/3)*(n)*(n-1)*TMath::Power((n),-1.0/3);
          Double_t bind_e_m2 = volume_term*(m+1)*(m) - 0.001*0.73*TMath::Power((density_ratio),1.0/3)*(m+1)*(m)*TMath::Power((m+1),-1.0/3);
          Double_t bind_e_n2 = volume_term*(n-1)*(n-2) - 0.001*0.73*TMath::Power((density_ratio),1.0/3)*(n-1)*(n-2)*TMath::Power((n-1),-1.0/3);
//            std::cout<<m<<"m" "    " "n"<<n<<"    "  "Energy n " <<2.38 - bind_e_n1/n<<" energy m" <<2.38 -  bind_e_m1/m <<endl;

	  //SetPxyPyPzM
	  // Initially 0 unless table condensate
	  
	  d_m2->SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
          d_n2->SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
 
          d_n1.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
          d_m1.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

	  //total energy formula 
	  Double_t e_m1 = (m*2.38) - bind_e_m1;
          Double_t e_m2 = ((m+1)*2.38) - bind_e_m2;
          Double_t e_n2 = ((n-1)*2.38) - bind_e_n2;
          Double_t e_n1 = (n*2.38) - bind_e_n1;
	  Double_t Masses_1[2] = {0,0};
          Double_t Masses_2[2] = {0,0}; 



	  //conidtionals for stab;le condensate creation between maximum mass energy and minimum decay 
	  if ((e_m1 > m*0.5) && (e_m1 < m*2.38 )){

	          if ((e_m2 > (m+1)*0.5) && (e_m2 < (m+1)*2.38 )){

			  if ((e_n1 > (n)*0.5) && (e_n1 < (n)*2.38 )){
        
				  if ((e_n2 > (n-1)*0.5) && (e_n2 < (n-1)*2.38 )){

				    //Enforced energy conservation
				    double tmpdiff = e_m1 + e_n1- e_n2-e_m2;
				  if ((tmpdiff)> 0.0){

				    //at rest defined masses reactants and products 
				   d_m1.SetPxPyPzE(0.0, 0.0, 0.0, e_m1);
				   d_n2->SetPxPyPzE(0.0, 0.0, 0.0, e_n2);
				   d_m2->SetPxPyPzE(0.0, 0.0, 0.0, e_m2);
				   d_n1.SetPxPyPzE(0.0, 0.0, 0.0, e_n1);
				   Double_t Masses_1[2] = {(e_m1),(e_n1)}; // DM_BEC with m and  Hexaquarks            
				   Double_t Masses_2[2] = {(e_m2),(e_n2)}; // DM_BEC with m+1, n-1 after

				   TLorentzVector V1 = d_m1 + d_n1;


				   // Creating decay vertices
				   TGenPhaseSpace Vertex_1;

				   TLorentzVector V2 =  *d_m2 + *d_n2;
				   Vertex_1.SetDecay(V1,2,Masses_2);
				   Double_t weight = Vertex_1.Generate();
				   
				   d_m2 = Vertex_1.GetDecay(0);
				   d_n2 = Vertex_1.GetDecay(1);
				   //std::cout<< d_m<<endl;

				   if (!d_m2 || !d_n2) {
				     std::cerr << "Decay products are uninitialized!" << std::endl;
				   }


				   if (d_m2->IsZombie()) {
				     std::cerr << "Decay product is a zombie object!" << std::endl;
				   }

				   // beta definition

				   double beta_m =d_m2->Beta();
        

				   double beta_n =d_n2->Beta();

				   ///////////////////////////////////////////////////////////////////////
				   //Histgrams fill and plot
    
				   if (((beta_m<1)&&(beta_m > 0))){
				     M_N->Fill(m, n, beta_m);
				     bind_1D-> Fill(m, bind_e_m1/m);
				     beta_dist -> Fill( beta_m);
				   }

				 }
				 }
			}
		}
	}



} 
}
TCanvas* c7 = new TCanvas("c7","1D",5,5,1000,1000);


//bind_1D->Draw();
beta_dist -> Draw();
beta_dist -> GetXaxis()->SetTitle("#beta, |v/c|");
beta_dist -> GetYaxis()->SetTitle("Frequency");
beta_dist -> SetTitle("Outgoing Velocity distribution");
gStyle->SetLabelFont(42, "XYZ"); // Helvetica (like a journal)
gStyle->SetTitleFont(42, "XYZ");
gStyle->SetTitleSize(0.05, "XYZ");
gStyle->SetLabelSize(0.04, "XYZ");
gStyle->SetPadTickX(1); // Ticks on top and bottom
gStyle->SetPadTickY(1); // Ticks on left and right

c7->SetGridx();
c7->SetGridy();


beta_dist->GetXaxis()->CenterTitle();  // X-axis
beta_dist->GetYaxis()->CenterTitle();  // Y-axis
beta_dist->GetZaxis()->CenterTitle();  // Z-axis (if using COLZ)

gStyle->SetOptStat(0); // Turn off statistics box
gStyle->SetNdivisions(515, "XY"); // More control over tick marks


TCanvas* c6 = new TCanvas("c6","distribution m and n",5,5,1000,1000);

 M_N->GetXaxis()->SetNdivisions(5); // X-axis: 5 numbered divisions
M_N->GetZaxis()->SetNdivisions(5); // Z-axis: 5 numbered divisions
M_N->GetYaxis()->SetNdivisions(5); // Y-axis: 5 numbered divisions
M_N->GetXaxis()->CenterTitle();  // X-axis
M_N->GetYaxis()->CenterTitle();  // Y-axis
M_N->GetZaxis()->CenterTitle();  // Z-axis (if using COLZ)

c6->Update();

gStyle->SetPaintTextFormat("1.0e");  // 3 significant figures
    TGaxis::SetExponentOffset(0.0, 0.04, "Z");
    TGaxis::SetExponentOffset(-0.06, -0.04, "X");

c6->SetGridx();
c6->SetGridy();


c6->Modified();

M_N->SetMarkerStyle(21);
M_N->GetXaxis()->SetTitle("m, number of hexaquarks");
M_N->GetYaxis()->SetTitle("n, number of hexaquarks");
M_N->GetZaxis()->SetTitle("#beta (v/c)");
M_N->GetZaxis()->SetRange(TMath::Power(10,-3),TMath::Power(10,-2));

gStyle->SetTextFont(42);           // Default text
M_N->GetXaxis()->SetLabelFont(42);
M_N->GetXaxis()->SetTitleFont(42);
M_N->GetYaxis()->SetLabelFont(42);
M_N->GetYaxis()->SetTitleFont(42);


M_N->Draw("colz");
//c6->SetLogz();
c6->Update();
	

} 
