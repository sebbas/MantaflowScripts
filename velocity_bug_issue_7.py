#
# Velocity extrapolation issue from Blender as a Mantaflow script
# https://github.com/sebbas/BlenderMantaflow/issues/7
#

from manta import *
import os, shutil, math, sys, gc

dim_s1     = 3
res_s1     = 64
gravity_s1 = vec3(0, 0, -1)
gs_s1      = vec3(64, 23, 49)

if dim_s1 == 2:
    gs_s1.z    = 1
    gravity_s1 = vec3(0,-1,0)

# ISSUE: Compare doOpen_s1=False and boundConditions_s1='' with doOpen_s1=True and boundConditions_s1='xXyY'
doOpen_s1              = False
boundConditions_s1     = ''
boundaryWidth_s1       = 1

using_highres_s1   = False
using_adaptTime_s1 = True

s1 = Solver(name='solver_s1', gridSize=gs_s1, dim=dim_s1)

dt_default_s1  = 0.1
dt_factor_s1   = 1
fps_s1         = 24
dt0_s1         = dt_default_s1 * (25.0 / fps_s1) * dt_factor_s1
s1.frameLength = dt0_s1
s1.timestepMin = dt0_s1 / 10
s1.timestepMax = dt0_s1
s1.cfl         = 2
s1.timestep    = (s1.timestepMax+s1.timestepMin)*0.5

flags_s1      = s1.create(FlagGrid)
phiParts_s1   = s1.create(LevelsetGrid)
phi_s1        = s1.create(LevelsetGrid)
phiIn_s1      = s1.create(LevelsetGrid)
pressure_s1   = s1.create(RealGrid)

phiObs_s1     = s1.create(LevelsetGrid)
phiObsIn_s1   = s1.create(LevelsetGrid)
fractions_s1  = s1.create(MACGrid)

vel_s1        = s1.create(MACGrid)
x_vel_s1      = s1.create(RealGrid)
y_vel_s1      = s1.create(RealGrid)
z_vel_s1      = s1.create(RealGrid)

velOld_s1     = s1.create(MACGrid)
velParts_s1   = s1.create(MACGrid)
mapWeights_s1 = s1.create(MACGrid)

pp_s1         = s1.create(BasicParticleSystem)
pVel_pp1      = pp_s1.create(PdataVec3)
mesh_s1       = s1.create(Mesh)

# Acceleration data for particle nbs
pindex_s1     = s1.create(ParticleIndexSystem)
gpi_s1        = s1.create(IntGrid)

forces_s1     = s1.create(MACGrid)
x_force_s1    = s1.create(RealGrid)
y_force_s1    = s1.create(RealGrid)
z_force_s1    = s1.create(RealGrid)

phi_s1.initFromFlags(flags_s1)
phiIn_s1.initFromFlags(flags_s1)

narrowBandWidth_s1         = 3
combineBandWidth_s1        = narrowBandWidth_s1 - 1
adjustedNarrowBandWidth_s1 = 4 # only used in adjustNumber to control band width

particleNumber_s1 = 2
minParticles_s1   = pow(particleNumber_s1, dim_s1)
radiusFactor_s1   = 2
randomness_s1     = 0.1
maxVel_s1         = 1 # just declared here, do not set

using_drops_s1   = False
using_bubbles_s1 = False
using_floats_s1  = False
using_tracers_s1 = False

# inflow
obsBox = Box( parent=s1, p0=gs_s1*vec3(0,0,0), p1=gs_s1*vec3(0.2,1,0.75))
phiObsIn_s1.join( obsBox.computeLevelset() )

fluidBox = Box( parent=s1, p0=gs_s1*vec3(0.1,0,0.85), p1=gs_s1*vec3(0.15,1.0,0.9))
phiIn_s1.join( fluidBox.computeLevelset() )

# cleanup
def liquid_post_step_low_1():
    forces_s1.clear()
    phiObs_s1.setConst(9999)

def liquid_step_1():
    pp_s1.advectInGrid(flags=flags_s1, vel=vel_s1, integrationMode=IntRK4, deleteInObstacle=False, stopInObstacle=False)
	
    pushOutofObs(parts=pp_s1, flags=flags_s1, phiObs=phiObs_s1)
	
    advectSemiLagrange(flags=flags_s1, vel=vel_s1, grid=phi_s1, order=1) # first order is usually enough
    advectSemiLagrange(flags=flags_s1, vel=vel_s1, grid=vel_s1, order=2, openBounds=doOpen_s1, boundaryWidth=boundaryWidth_s1)
    
    # create level set of particles
    gridParticleIndex(parts=pp_s1, flags=flags_s1, indexSys=pindex_s1, index=gpi_s1)
    unionParticleLevelset(pp_s1, pindex_s1, flags_s1, gpi_s1, phiParts_s1)
    
    # combine level set of particles with grid level set
    phi_s1.addConst(1.) # shrink slightly
    phi_s1.join(phiParts_s1)
    extrapolateLsSimple(phi=phi_s1, distance=narrowBandWidth_s1+2, inside=True)
    extrapolateLsSimple(phi=phi_s1, distance=3)
    phi_s1.setBoundNeumann(boundaryWidth_s1) # make sure no particles are placed at outer boundary
    
    if doOpen_s1:
        resetOutflow(flags=flags_s1, phi=phi_s1, parts=pp_s1, index=gpi_s1, indexSys=pindex_s1)
    flags_s1.updateFromLevelset(phi_s1)
    
    # combine particles velocities with advected grid velocities
    mapPartsToMAC(vel=velParts_s1, flags=flags_s1, velOld=velOld_s1, parts=pp_s1, partVel=pVel_pp1, weight=mapWeights_s1)
    extrapolateMACFromWeight(vel=velParts_s1, distance=2, weight=mapWeights_s1)
    combineGridVel(vel=velParts_s1, weight=mapWeights_s1, combineVel=vel_s1, phi=phi_s1, narrowBand=combineBandWidth_s1, thresh=0)
    velOld_s1.copyFrom(vel_s1)
    
    # forces & pressure solve
    addGravity(flags=flags_s1, vel=vel_s1, gravity=gravity_s1)
    addForceField(flags=flags_s1, vel=vel_s1, force=forces_s1)
	
    extrapolateMACSimple(flags=flags_s1, vel=vel_s1, distance=2, phiObs=phiObs_s1, intoObs=True)
    setWallBcs(flags=flags_s1, vel=vel_s1, fractions=fractions_s1, phiObs=phiObs_s1)
    
    solvePressure(flags=flags_s1, vel=vel_s1, pressure=pressure_s1, phi=phi_s1, fractions=fractions_s1)
    
    extrapolateMACSimple(flags=flags_s1, vel=vel_s1, distance=4, phiObs=phiObs_s1, intoObs=True)
    setWallBcs(flags=flags_s1, vel=vel_s1, fractions=fractions_s1, phiObs=phiObs_s1)
    
    if (dim_s1==3):
        # mis-use phiParts as temp grid to close the mesh
        phiParts_s1.copyFrom(phi_s1)
        phiParts_s1.setBound(0.5,0)
        phiParts_s1.createMesh(mesh_s1)
    
    # Create interpolated version of original phi grid for later use in (optional) high-res step
    if using_highres_s1:
        interpolateGrid(target=phi_xl1, source=phiParts_s1)
    
    # set source grids for resampling, used in adjustNumber!
    pVel_pp1.setSource(vel_s1, isMAC=True)
    adjustNumber(parts=pp_s1, vel=vel_s1, flags=flags_s1, minParticles=1*minParticles_s1, maxParticles=2*minParticles_s1, phi=phi_s1, exclude=phiObs_s1, radiusFactor=radiusFactor_s1, narrowBand=adjustedNarrowBandWidth_s1)
    flipVelocityUpdate(vel=vel_s1, velOld=velOld_s1, flags=flags_s1, parts=pp_s1, partVel=pVel_pp1, flipRatio=0.97)

def manta_step_1(framenr):
    s1.frame = framenr
    s1.timeTotal = s1.frame * dt0_s1
    last_frame_s1 = s1.frame

    while s1.frame == last_frame_s1:
        
        flags_s1.initDomain(boundaryWidth=boundaryWidth_s1, phiWalls=phiObs_s1, outflow=boundConditions_s1)
        
        phiObs_s1.join(phiObsIn_s1)
        phi_s1.join(phiIn_s1)
        phi_s1.subtract(phiObsIn_s1)

        updateFractions(flags=flags_s1, phiObs=phiObs_s1, fractions=fractions_s1, boundaryWidth=boundaryWidth_s1)
        setObstacleFlags(flags=flags_s1, phiObs=phiObs_s1, fractions=fractions_s1)
        
        sampleLevelsetWithParticles(phi=phiIn_s1, flags=flags_s1, parts=pp_s1, discretization=particleNumber_s1, randomness=randomness_s1, refillEmpty=True)
        flags_s1.updateFromLevelset(phi_s1)
        
        if using_adaptTime_s1:
            maxVel_s1 = vel_s1.getMaxValue()
            s1.adaptTimestep(maxVel_s1)
        
        liquid_step_1()
        
        if using_highres_s1:
            xl1.timestep = s1.timestep
            liquid_step_high_1()
        s1.step()
    
    liquid_post_step_low_1()

if (GUI):
    gui=Gui()
    gui.show()
    gui.pause()

start_frame = 1
end_frame = 1000

# All low and high res steps
while start_frame <= end_frame:
    manta_step_1(start_frame)
    start_frame += 1
