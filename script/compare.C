#include <sstream>

void compare( TString fn1, TString l1, TString fn2, TString l2, TString fn3 = "", TString l3 = ""){
	
	TCanvas* c = new TCanvas( "c1", "c1", 600, 800 );
	TFile * f1 = new TFile( fn1 );
	TFile * f2 = new TFile( fn2 );
	TFile * f3 = NULL;
	if ( fn3 != "" )
		f3 = new TFile( fn3 );

	c->Print( "compare.pdf[" );

		gStyle->SetOptStat( 0 );

	stringstream sstr;
	for ( int i = 1; i < 39; i++ ){
		TLegend * leg = new TLegend( 0.6, 0.75, 0.9, 0.9 );
		sstr.str("");    	sstr << "corno" << i;

	
	    TH1D* tmp1 = (TH1D*) f1->Get( sstr.str().c_str() );
	 
	    TH1D* tmp2 = (TH1D*) f2->Get( sstr.str().c_str() );

	    

	    tmp2->SetLineColor( 3 );
	    tmp2->SetLineWidth( 2 );

	    leg->AddEntry( tmp1, l1, "l" );
	    leg->AddEntry( tmp2, l2, "l" );

	    tmp1->SetTitle( "Slewing Correction");
	    tmp1->GetYaxis()->SetRangeUser( -5, 2);
	    tmp1->SetLineWidth( 2 );
	    tmp1->Draw();
	    tmp2->SetTitle( "Slewing Correction");
	    tmp2->Draw("same");

	    if ( f3 ){
	    	TH1D* tmp3 = (TH1D*) f3->Get( sstr.str().c_str() );
	    	tmp3->SetTitle( "Slewing Correction");
	    	leg->AddEntry( tmp3, l3, "l" );
	    	tmp3->SetLineColor( 2 );
	    	tmp3->SetLineWidth( 2 );
			tmp3->Draw("same");	    	

		}

	    leg->Draw();

	    c->Print( "compare.pdf" );
	    delete leg;
	}
	
	
	



	c->Print( "compare.pdf]" );

}