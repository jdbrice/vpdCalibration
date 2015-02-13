

void tofQAMaker( string file, int xRange = 10 ){

	TFile * f = new TFile( file.c_str(), "READ" );


	TH1D * tpc = (TH1D*)f->Get( "zvertex" );
	TH1D * vpd = (TH1D*)f->Get( "zvertexVPD" );
	TH1D * delta = (TH1D*)f->Get( "zvertexDelta" );

	Reporter rp( "vertexQA.pdf" );

	TF1* gaus = new TF1( "g", "gaus" );
	gaus->SetRange( -xRange, xRange );

	rp.newPage();
	tpc->Draw();
	rp.savePage();
	rp.newPage();
	vpd->Draw();
	rp.savePage();
	rp.newPage();
	delta->GetXaxis()->SetRangeUser( -25, 25 );
	gStyle->SetOptFit( 111 );
	delta->Draw();
	delta->Fit( gaus, "R" );
	rp.savePage();

}