#
# This FLIP example combines narrow band flip, 2nd order wall boundary conditions,
# adaptive time stepping, and moving obstacles
# 
from manta import *

# simulate with obstacle fraction?
withFractions = False
# sphere or box obstacle?
withSphereObs = False

dim    = 3
res    = 40
gs     = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

narrowBand    = 3
minParticles  = pow(2,dim)
frames        = 200

# Adaptive time stepping
s.frameLength = 1.0                 # length of one frame (in "world time")
s.cfl         = 4.0                 # maximal velocity per cell and timestep, 3 is fairly strict
s.timestep    = s.frameLength 
s.timestepMin = s.frameLength / 2.  # time step range
s.timestepMax = s.frameLength * 2.

# prepare grids and particles
flags     = s.create(FlagGrid)
phi       = s.create(LevelsetGrid)
phiParts  = s.create(LevelsetGrid)
phiObs    = s.create(LevelsetGrid)

vel       = s.create(MACGrid)
velOld    = s.create(MACGrid)
velParts  = s.create(MACGrid)
obsVel  = s.create(MACGrid)

pressure  = s.create(RealGrid)
tmpVec3   = s.create(VecGrid)
fractions = s.create(MACGrid) if withFractions else None

pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 
mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# scene setup
bWidth=1
flags.initDomain(boundaryWidth=bWidth, phiWalls=phiObs if withFractions else None)
phi.setConst(999.)

# standing dam
fluidbox1 = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.3,1)) 
phi.join( fluidbox1.computeLevelset() )

# init obstacle properties
endVelFrame = 20 # when to stop moving the obstacle
obsSize = 0.2
obsPosSphere = vec3(0.2,0.2,0.5)
obsPosBoxP0 = vec3(0.15-obsSize*0.5,0.1,0.4)
obsPosBoxP1 = vec3(0.15+obsSize*0.5,0.5,0.6)
obsVelVec = vec3(0.6,0.0,0.0) * (1./Real(endVelFrame)) * float(res) # velocity in grid units for 100 steps
obsVel.setConst(obsVelVec)
obsVel.setBound(value=Vec3(0.), boundaryWidth=bWidth+1) # make sure walls are static
#obs = "dummy"; phiObs = "dummy2"

# init obstacle
if withSphereObs:
	obs = Sphere( parent=s, center=gs*obsPosSphere, radius=res*obsSize)
else:
	obs = Box( parent=s, p0=gs*obsPosBoxP0, \
						 p1=gs*obsPosBoxP1)
phiObs.join( obs.computeLevelset() )

flags.updateFromLevelset(phi)
phi.subtract( phiObs );
sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.05 )

# also sets boundary flags for phiObs
if withFractions:
	updateFractions( flags=flags, phiObs=phiObs, fractions=fractions, boundaryWidth=bWidth )
setObstacleFlags(flags=flags, phiObs=phiObs, fractions=fractions if withFractions else None )

lastFrame = -1
if 1 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

#main loop
while s.frame < frames:
	t = s.timeTotal
	maxVel = vel.getMax()
	s.adaptTimestep( maxVel )
	
	# FLIP 
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False, stopInObstacle=False )
	pushOutofObs( parts=pp, flags=flags, phiObs=phiObs )

	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1) # first order is usually enough
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

	if t<=endVelFrame:
		phiObs.clear()
		flags.initDomain(boundaryWidth=bWidth, phiWalls=phiObs if withFractions else None )

		del obs
		if withSphereObs:
			obs = Sphere( parent=s, center=gs*obsPosSphere + float(t) * obsVelVec, radius=res*obsSize)
		else:
			obs = Box( parent=s, p0=gs*obsPosBoxP0 + float(t) * obsVelVec, \
			                     p1=gs*obsPosBoxP1 + float(t) * obsVelVec)
		phiObs.join(obs.computeLevelset())

		if withFractions:
			updateFractions( flags=flags, phiObs=phiObs, fractions=fractions, boundaryWidth=bWidth )
		setObstacleFlags(flags=flags, phiObs=phiObs, fractions=fractions if withFractions else None)
	elif t==endVelFrame+1:
		# stop moving
		obsVel.setConst(Vec3(0.))
		
	# create level set of particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phiParts )

	# combine level set of particles with grid level set
	phi.addConst(1.); # shrink slightly
	phi.join( phiParts );
	extrapolateLsSimple(phi=phi, distance=narrowBand+2, inside=True ) 
	extrapolateLsSimple(phi=phi, distance=3 )
	phi.setBoundNeumann(0) # make sure no particles are placed at outer boundary, warning - larger values can delete thin sheets at outer walls...
	flags.updateFromLevelset(phi)

	# combine particles velocities with advected grid velocities
	mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3)
	extrapolateMACFromWeight( vel=velParts , distance=2, weight=tmpVec3 )
	combineGridVel(vel=velParts, weight=tmpVec3 , combineVel=vel, phi=phi, narrowBand=(narrowBand-1), thresh=0)
	velOld.copyFrom(vel)

	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=(0,-0.001,0))

	extrapolateMACSimple( flags=flags, vel=vel , distance=2, intoObs=True )
	setWallBcs(flags=flags, vel=vel, fractions=fractions if withFractions else None, phiObs=phiObs, obvel=obsVel)

	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi, fractions=fractions if withFractions else None )

	extrapolateMACSimple( flags=flags, vel=vel , distance=4, intoObs=True )
	setWallBcs(flags=flags, vel=vel, fractions=fractions if withFractions else None, phiObs=phiObs, obvel=obsVel)
	
	# make sure we have proper velocities
	extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel*1.5)+2) )

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

	#s.printMemInfo()
	lastFrame = s.frame;


