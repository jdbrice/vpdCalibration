


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

	TH1D* sigma = new TH1D( "sigma", "sigma 0", 38, 0.5, 39.5 );

	for ( int i = 0; i < 38; i++){

		stringstream sstr;
		sstr << "channel" << i;
		f->cd( sstr.str().c_str() );

		TF1 * g = new TF1( "g", "gaus", -1.0, 1.0 );

		it5cutAvgN->FitSlicesY(g, 0, -1, 10);


		TH1D* sig = (TH1D*)it5cutAvgN_2->Clone("sigma");
		
		TF1 * fr = new TF1( "fr", detectorResolution, 0, 19, 1);

		sig->Fit( "fr", "QR" );

		sig->Draw("P");

		cout << "res: " << fr->GetParameter( 0 ) << endl;
		sigma->SetBinContent( i+1, fr->GetParameter( 0 ) );
	}

	sigma ->Draw("HP");





}