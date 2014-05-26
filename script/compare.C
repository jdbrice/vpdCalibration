#include <sstream>

void compare( TString fn1, TString fn2){
	
	TCanvas* c = new TCanvas( "c1", "c1", 600, 800 );
	TFile * f1 = new TFile( fn1 );
	TFile * f2 = new TFile( fn2 );

	c->Print( "cor.pdf[" );

	stringstream sstr;
	for ( int i = 1; i < 39; i++ ){
		TLegend * leg = new TLegend( 0.6, 0.6, 0.9, 0.9 );
		sstr.str("");    	sstr << "cor" << i;

		//f1->cd();
	    TH1D* tmp1 = (TH1D*) f1->Get( sstr.str().c_str() );
	    //f2->cd();
	    TH1D* tmp2 = (TH1D*) f2->Get( sstr.str().c_str() );

	    tmp2->SetLineColor( 3 );

	    leg->AddEntry( tmp1, fn1, "l" );
	    leg->AddEntry( tmp2, fn2, "l" );

	    tmp1->GetYaxis()->SetRangeUser( -15, 15);
	    tmp1->Draw();
	    tmp2->Draw("same");

	    leg->Draw();

	    c->Print( "cor.pdf" );
	    delete leg;
	}
	
	
	



	c->Print( "cor.pdf]" );

}