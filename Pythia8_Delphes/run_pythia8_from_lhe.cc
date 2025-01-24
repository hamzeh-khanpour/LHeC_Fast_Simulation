//
//  g++ run_pythia8_from_lhe.cc -o run_pythia8_from_lhe -I../include -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -L../lib -Wl,-rpath,../lib -lpythia8 -ldl
//
//  ./run_pythia8_from_lhe

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {
    // Initialize Pythia instance
    Pythia pythia;

    // Enable reading LHEF
    pythia.readString("Beams:frameType = 4"); // Specifies LHE input
    pythia.readString("Beams:LHEF = /home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_SM/Events/run_01/aa_ww_semi_leptonic_SM.lhe"); // Path to your LHE file

    // Initialize Pythia
    pythia.init();

    // Loop over events
    for (int iEvent = 0; iEvent < 1000; ++iEvent) {
        if (!pythia.next()) continue; // Generate the next event

        // Print basic event information
        pythia.event.list();
    }

    // Print statistics
    pythia.stat();

    return 0;
}
