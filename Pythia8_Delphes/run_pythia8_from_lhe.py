import pythia8

# Initialize Pythia
pythia = pythia8.Pythia()

# Set up to read LHE file
pythia.readString("Beams:frameType = 4")  # Specifies LHE input
pythia.readString("Beams:LHEF = /home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_SM/Events/run_01/aa_ww_semi_leptonic_SM.lhe")  # Path to your LHE file

# Initialize Pythia
pythia.init()

# Process events
for iEvent in range(1000):  # Generate 1000 events
    if not pythia.next():
        continue  # Skip invalid events

    # Access particles in the event
    for particle in pythia.event:
        print(f"ID: {particle.id()}, px: {particle.px()}, py: {particle.py()}, pz: {particle.pz()}")

# Print statistics
pythia.stat()
