import pythia8

# Initialize Pythia
pythia = pythia8.Pythia()

# Configure Pythia to read the LHE file
pythia.readString("Beams:frameType = 4")  # Frame type for LHE input
pythia.readString("Beams:LHEF = /home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_SM/Events/run_01/aa_ww_semi_leptonic_SM.lhe")  # Path to your LHE file

# Enable Parton Showering and Hadronization
pythia.readString("HadronLevel:all = on")  # Enable hadronization
pythia.readString("PartonLevel:FSR = on")  # Enable final state radiation
pythia.readString("PartonLevel:ISR = on")  # Enable initial state radiation
pythia.readString("PartonLevel:MPI = on")  # Enable multiple parton interactions

# Initialize Pythia
pythia.init()

# Loop over events
for iEvent in range(10):  # Adjust the number of events as needed
    if not pythia.next():
        continue  # Skip event if it fails

    print(f"Event {iEvent + 1}")
    # Loop over particles in the event
    for i in range(pythia.event.size()):
        particle = pythia.event[i]
        if particle.isFinal():  # Only final state particles
            print(f"  Particle {i}: ID={particle.id()}, "
                  f"pT={particle.pT():.2f} GeV, eta={particle.eta():.2f}, phi={particle.phi():.2f}")

# Print statistics
pythia.stat()
