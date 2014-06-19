

#include "random.h"



Double_t detectorResolution(Double_t *x, Double_t *par){
  Double_t resval = 0.0;
  if (x[0]>0){
    resval  = TMath::Sqrt( x[0] / ( 1.0 + x[0] ) );
    return ( par[0] / resval );
  } else {
    return 0.0;
  }
}

void exploreStat(){
	


  	TFile * f = new TFile( "stat.root", "RECREATE" );

  	TRandom2 * r1 = new TRandom2( clock() );

    TH2D* stat = new TH2D( "stat", "stat", 200, 0.5, 200.5, 500, -10, 10 );

    for ( int N = 1; N < 201; N++ ){

      for ( int j = 0; j < 5000; j++ ){

        double sample = r1->Gaus( 0, 2 );

        double total = 0;
        for (int i=0; i<N; ++i) {
          double number = r1->Gaus( 0, r1->Gaus( 0, 5 ) );
          total += number;
        }
        double avg = total / (double)N;

        book->fill( "stat", N, sample - avg );

      }

    }

  	gROOT->SetOptStat( 1111 );
  	stat->Draw( "colz" );

    TF1 * fr = new TF1( "fr", detectorResolution, 0, 190, 1);
    book->get2D("stat" )->FitSlicesY();
    TH1D* fsy = (TH1D*) gDirectory->Get( "stat_2" );
    book->add( "statSig", fsy );
    fsy->Fit( "fr", "R" );


  	f->Write();

}

