#include <TFile.h>
#include <TNtuple.h>
#include <TRandom1.h>
#include <math.h>

int main()
{
  TRandom *r = new TRandom1();

  // Test 1: test multiple files
  // Produces a flat background 0-10 MeV and 0-6000mm radius
  TFile *f1 = new TFile("data1_1.root","RECREATE");
  TNtuple *n1 = new TNtuple("ntuple","data from ascii file","energy:radius");
  for (int i=0;i<10000;i++){
    n1->Fill(r->Rndm()*10,pow(r->Rndm(),0.3333)*6000);
  }
  f1->Write();
  f1->Close();

  TFile *f2 = new TFile("data1_2.root","RECREATE");
  TNtuple *n2 = new TNtuple("ntuple","data from ascii file","energy:radius");
  for (int i=0;i<10000;i++){
    n2->Fill(r->Rndm()*10,pow(r->Rndm(),0.3333)*6000);
  }
  f2->Write();
  f2->Close();

  // Test 2: test multiple contributions
  // Produces a flat background 0-5 MeV and a spike at 8 MeV near the edge
  TFile *f3 = new TFile("data2_1.root","RECREATE");
  TNtuple *n3 = new TNtuple("ntuple","data from ascii file","energy:radius");
  for (int i=0;i<10000;i++){
    n3->Fill(r->Rndm()*5,pow(r->Rndm(),0.3333)*6000);
  }
  f3->Write();
  f3->Close();

  TFile *f4 = new TFile("data2_2.root","RECREATE");
  TNtuple *n4 = new TNtuple("ntuple","data from ascii file","energy:radius");
  for (int i=0;i<10000;i++){
    n4->Fill(8,pow(r->Rndm(),0.1)*6000);
  }
  f4->Write();
  f4->Close();

  // Test 3: test signal with zero expected
  // Produces a flat background from 0-1 MeV
  
  TFile *f5 = new TFile("data3_1.root","RECREATE");
  TNtuple *n5 = new TNtuple("ntuple","data from ascii file","energy:radius");
  for (int i=0;i<10000;i++){
    n5->Fill(r->Rndm(),pow(r->Rndm(),0.33333)*6000);
  }
  f5->Write();
  f5->Close();

  // Test 4: test signal with expected but all cut
  // Produces a flat background 0-5 MeV all at the edge

  TFile *f6 = new TFile("data4_1.root","RECREATE");
  TNtuple *n6 = new TNtuple("ntuple","data from ascii file","energy:radius");
  for (int i=0;i<10000;i++){
    n6->Fill(r->Rndm()*5,6000);
  }
  f6->Write();
  f6->Close();
  return 0;
}
