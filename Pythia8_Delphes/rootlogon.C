{
  // Please replace with location of your libDelphes.so
  gSystem->Load("/home/hamzeh-khanpour/Delphes-3.5.0/libDelphes.so");
  // Please replace with location of your Delphes folder
  gROOT->ProcessLine(".include /home/hamzeh-khanpour/Delphes-3.5.0");
  // Please replace with location of your Delphes/external folder
  gROOT->ProcessLine(".include /home/hamzeh-khanpour/Delphes-3.5.0/external");
}
