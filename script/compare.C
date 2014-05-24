

void compare( TString fn1, TString fn2){
	
	TCanvas* c = new TCanvas( "c1", "c1", 600, 800 );
	TFile * f1 = new TFile( fn1 );
	TFile * f2 = new TFile( fn2 );

	c->Print( "cor.pdf[" );
	TLegend * leg = new TLegend( 0.6, 0.6, 0.9, 0.9 );
	f1->cd();
	cor1->Draw();
	cor1->GetYaxis()->SetRangeUser( -15, 15);
	leg->AddEntry( cor1, fn1, "l" );
	f2->cd();
	cor1->SetLineColor( 3 );
	cor1->Draw("same");
	leg->AddEntry( cor1, fn2, "l" );
	leg->Draw();


	c->Print( "cor.pdf" );

	f1->cd();
	cor2->Draw();
	cor2->GetYaxis()->SetRangeUser( -15, 15);
	f2->cd();
	cor2->SetLineColor( 3 );
	cor2->Draw("same");
	c->Print( "cor.pdf" );

	f1->cd();
	cor3->Draw();
	cor3->GetYaxis()->SetRangeUser( -15, 15);
	f2->cd();
	cor3->SetLineColor( 3 );
	cor3->Draw("same");
	c->Print( "cor.pdf" );

	f1->cd();
	cor4->Draw();
	cor4->GetYaxis()->SetRangeUser( -15, 15);
	f2->cd();
	cor4->SetLineColor( 3 );
	cor4->Draw("same");
	c->Print( "cor.pdf" );

	f1->cd();
	cor5->Draw();
	cor5->GetYaxis()->SetRangeUser( -15, 15);
	f2->cd();
	cor5->SetLineColor( 3 );
	cor5->Draw("same");
	c->Print( "cor.pdf" );

	f1->cd();
	cor6->Draw();
	cor6->GetYaxis()->SetRangeUser( -15, 15);
	f2->cd();
	cor6->SetLineColor( 3 );
	cor6->Draw("same");
	c->Print( "cor.pdf" );

	f1->cd();
	cor7->Draw();
	cor7->GetYaxis()->SetRangeUser( -15, 15);
	f2->cd();
	cor7->SetLineColor( 3 );
	cor7->Draw("same");
	c->Print( "cor.pdf" );

	f1->cd();
	cor8->Draw();
	cor8->GetYaxis()->SetRangeUser( -15, 15);
	f2->cd();
	cor8->SetLineColor( 3 );
	cor8->Draw("same");
	c->Print( "cor.pdf" );

	c->Print( "cor.pdf]" );

}