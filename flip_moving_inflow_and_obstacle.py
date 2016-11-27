#
# This FLIP example is an adaption of flip06_obstacle.py
# Instead of static flow and obstacle objects, it is using a moving objects
# (for flow and obstacle) and continous inflow
#
from manta import *

dim	= 3
res	= 64
gs	 = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

narrowBand     = 3
minParticles   = pow(2,dim)
frames         = 500
particleNumber = 2
randomness     = 0.1
bWidth         = 1
inVel          = True

# Adaptive time stepping
s.frameLength = 0.8				 # length of one frame (in "world time")
s.cfl         = 3.0				 # maximal velocity per cell and timestep, 3 is fairly strict
s.timestep    = s.frameLength
s.timestepMin = s.frameLength / 4.  # time step range
s.timestepMax = s.frameLength * 4.

# prepare grids, particles and mesh
flags     = s.create(FlagGrid)
phi       = s.create(LevelsetGrid)
phiIn     = s.create(LevelsetGrid)
phiParts  = s.create(LevelsetGrid)
phiObs    = s.create(LevelsetGrid)
phiObsIn  = s.create(LevelsetGrid)
vel       = s.create(MACGrid)
velOld    = s.create(MACGrid)
velParts  = s.create(MACGrid)
pressure  = s.create(RealGrid)
fractions = s.create(MACGrid)
tmpVec3   = s.create(VecGrid)
pp        = s.create(BasicParticleSystem)
pVel      = pp.create(PdataVec3)
mesh      = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi = s.create(IntGrid)

phi.initFromFlags(flags)

# optional inflow velocity
if inVel:
    fluidVel = Cylinder(parent=s, center=gs*vec3(0.99,0.8,0.5), radius=res*0.05, z=gs*vec3(0, 0.01, 0))
    fluidSetVel = vec3(-5,-5,0)

lastFrame = -1
if 1 and (GUI):
	gui = Gui()
	gui.show()
	#gui.pause()

#main loop
while s.frame < frames:

	## BEGIN INFLOW
	if lastFrame<50:
		
		# Set initial velocity for fluid flow (optional)
		if inVel:
			fluidVel.applyToGrid( grid=vel , value=fluidSetVel )
			mapGridToPartsVec3(source=vel, parts=pp, target=pVel )
		
		# Create a moving inflow object and a moving obstacle object
		# moving inflow object
		cylinder = Cylinder(parent=s, center=gs*vec3(1-(0.01*lastFrame),0.8,0.5), radius=res*0.05, z=gs*vec3(0, 0.01, 0))
		phiIn = cylinder.computeLevelset()
		
		# moving obstacle object
		sphere = Sphere(parent=s, center=gs*vec3(0.01*lastFrame,0.3,0.5), radius=res*0.2)
		phiObsIn = sphere.computeLevelset()
		
		# Because of moving objects we have to reset flags for every step. Also resetting phiObs
		flags.initDomain(boundaryWidth=bWidth, phiWalls=phiObs)
		
		# Create actual grid grid from both obstacle and fluid inflow - order of join/subtract operations important!
		phiObs.join(phiObsIn)
		phiIn.subtract(phiObs)
		phi.join(phiIn)

		# Now that phiObs is up-to-date for this step, we can also update the flags
		updateFractions(flags=flags, phiObs=phiObs, fractions=fractions, boundaryWidth=bWidth)
		setObstacleFlags(flags=flags, phiObs=phiObs, fractions=fractions)

		# Sample new particles, but only into parts of phiIn that do not contain any particles anymore
		sampleLevelsetWithParticles(phi=phiIn, flags=flags, parts=pp, discretization=particleNumber, randomness=randomness, refillEmpty=True)
		
		# Set flow cells to fluid (1), and other non-obstacle cells to empty (4)
		# Important! This has to be called after the sampling. Otherwise the 'refillEmpty' option does not work (function would break too early out of sampling loop)
		flags.updateFromLevelset(phi)

	## END INFLOW

	# From here, same step routine as in flip06_obstacle.py
	maxVel = vel.getMaxValue()
	s.adaptTimestep( maxVel )
	mantaMsg('\nFrame %i, time-step size %f' % (s.frame, s.timestep))
	
	# FLIP 
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False, stopInObstacle=False )
	pushOutofObs( parts=pp, flags=flags, phiObs=phiObs )

	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1) # first order is usually enough
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

	# create level set of particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phiParts )

	# combine level set of particles with grid level set
	phi.addConst(1.); # shrink slightly
	phi.join( phiParts );
	extrapolateLsSimple(phi=phi, distance=narrowBand+2, inside=True ) 
	extrapolateLsSimple(phi=phi, distance=3 )
	phi.setBoundNeumann(1) # make sure no particles are placed at outer boundary
	flags.updateFromLevelset(phi)

	# combine particles velocities with advected grid velocities
	mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3)
	extrapolateMACFromWeight( vel=velParts , distance=2, weight=tmpVec3 )
	combineGridVel(vel=velParts, weight=tmpVec3 , combineVel=vel, phi=phi, narrowBand=(narrowBand-1), thresh=0)
	velOld.copyFrom(vel)

	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=(0,-0.001,0))

	extrapolateMACSimple( flags=flags, vel=vel , distance=2, intoObs=True )
	setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)	

	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi, fractions=fractions )

	extrapolateMACSimple( flags=flags, vel=vel , distance=4, intoObs=True )
	setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)

	if (dim==3):
		# mis-use phiParts as temp grid to close the mesh
		phiParts.copyFrom(phi)
		phiParts.setBound(0.5,0)
		phiParts.createMesh(mesh)

	# set source grids for resampling, used in adjustNumber!
	pVel.setSource( vel, isMAC=True )
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, exclude=phiObs, narrowBand=narrowBand ) 
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )
	
	s.step()
	lastFrame = s.frame

	# cleanup: objects cylinder and sphere are created in loop but never freed (probably not ideal)
	if 'sphere'   in globals() : del sphere
	if 'cylinder' in globals() : del cylinder



