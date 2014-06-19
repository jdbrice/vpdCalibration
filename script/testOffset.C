


void testOffset( ){
	
	TFile * f = new TFile( "qa.root" );

	f->cd( "initialOffset");
	TH1D* tdc = (TH1D*)f->Get( "tdc");

	//tdc->Draw( "colz" );

}