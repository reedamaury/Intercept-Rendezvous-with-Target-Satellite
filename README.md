# Intercept-Rendezvous-with-Target-Satellite

In a joint Blue Origin/Stanford capstone project, sponsored by Blue Origin, I worked with a small team to develop CubeSat mission aimed at identifying, targeting, attaching to, and deorbiting an inactive, uncooperative satellite in Amazon's Kuiper Constellation. In the development of this mission, it was necessary to select a launch provider.

To do so, we established requirements to be met by our selected launch vehicle, namely the ​ΔV required to inject into the target orbit and rendezvous with the target inactive satellite.  Our approach involved employing a Lambert solver, which solves a two-point boundary value problem. Given the launch site coordinates at the initial time (initial position), the time of flight (TOF), and the satellite's coordinates at the final time (final position), the Lambert solver outputs a short-way, conic section trajectory and the associated impulsive ​ΔV burn. A second impulsive burn was required for rendezvous with the target satellite.

To find the minimum ​ΔV , we employed a sweeping algorithm, reminiscent of the one used to generate porkchop plots for interplanetary trajectories. For the specified launch date, we looked at a 6 hour launch window, with launch times in intervals of 1 minute. For each launch time, we swept through 120 minutes of TOF's (incrementing the TOF by 2 minutes in each iteration), and calculated the ​ΔV for each launch time and TOF combination. Our algorithm also included a conditionals, which excluded unrealistic parabolic trajectories and trajectories that impacted the Earth before reaching the target satellite. For a launch from Cape Canaveral on April 2nd, 2024, targeting a hypothetical inactive Kuiper satellite, we obtained the following results: 


Launch Date: 4/2/2024

Launch Time: 04:24:00 UTC

Time of Flight:  20 mins

Arrival Time: 04:44:00 UTC

ΔV requirement: 8.617 km/s

Launch Vehicle: SpaceX Falcon 9


A 3D visualization of the scenario is included in the LambertRendezvous3dPlot script. 

While this approach is a good first step in approximating the ΔV requirement, it is still very ideal. We are currently working on a more advanced simulation, where we inject into an orbit that slightly off from the desired orbit and then employ the Clohessy-Wiltshire equations for rendezvous with target satellite. 
