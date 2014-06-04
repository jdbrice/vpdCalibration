
Double_t detectorResolution(Double_t *x, Double_t *par){
	Double_t resval = 0.0;
	if (x[0]>0){
		resval	= TMath::Sqrt( x[0] / ( 1.0 + x[0] ) );
		return ( par[0] / resval );
	} else {
		return 0.0;
	}
}


void detRes(){

	TFile * f = new TFile( "qa.root" );

	f->cd( "channel0" );

	TF1 * g = new TF1( "g", "gaus", -1.0, 1.0 );

	it1cutAvgNB->FitSlicesY(g, 0, -1, 10);

	TH1D* sig = (TH1D*)it1cutAvgNB_2->Clone("sigma");
	
	TF1 * fr = new TF1( "fr", detectorResolution, 0, 19, 1);

	sig->Fit( "fr", "NQR" );

	sig->Draw("P");





}