import lm
from pyLM.units import *
import numpy as np
import scipy as sp
from scipy import ndimage
import math
import itertools

def createERCellularAutomaton(size, 
                              dimensions=2, 
                              survivalRate = 0.5, 
                              deathLimitFxn = 3,
                              birthNumberFxn = 4,
                              steps = 10,
                              limits = []):
    '''
    size - dimension in lattice sites
    dimensions - dimensionality of the lattice (1, 2, or 3)
    survivalRate - the rate of survival; a function that takes the step and returns a number (2D: 0-9, 3D: 0-27)
    limits - A list of functions; a function that takes the step number and returns a birth number (2D: 0-9, 3D: 0-27)
    steps - number of iterations of the "game"
    '''
    # Generate lattices objects
    print("Creating grids")
    shape = dimensions*[size]
    lattice0 = np.zeros(shape=shape, dtype=int)
    lattice1 = np.zeros(shape=shape, dtype=int)
    
    # Randomly set up first lattice
    print("Creating random lattice")
    lattice0 = np.random.choice([0,1], 
                                shape, 
                                p=[1.0-survivalRate, survivalRate])
    
    # Carve out initial value
    print("Carving initial value")
    x = np.zeros(shape)
    ind = np.indices(x.shape) - 140
    xdiff = ind[0]
    ydiff = ind[1]
    zdiff = ind[2]
    dist = xdiff**2 + ydiff**2 + zdiff**2
    boolLattice = np.ones(shape, dtype=int)
    for lFxn in limits:
        boolLattice[lFxn(dist)] = 0
    lattice0 = lattice0 & boolLattice
    totalCytoplasmic = np.sum(boolLattice)
    print("Total cytoplasmic sites:",totalCytoplasmic)
    #return lattice0
    
    def countAliveNeightbors(lat, idx):
        count = np.sum(lat[(idx[0]-1):(idx[0]+2), 
                           (idx[1]-1):(idx[1]+2), 
                           (idx[2]-1):(idx[2]+2)])
        return count
        
        count = 0
        for idxOffset in itertools.product([-1,0,1], 
                                           repeat=dimensions):
            try:
                newIdx = np.array(idx)+np.array(idxOffset)
                if lat[tuple(newIdx)]:
                    count += 1
            except Exception as e:
                pass
        return count
    
    for s in range(steps):
        # Compute rates for this step
        deathLimit = deathLimitFxn(s)
        birthNumber = birthNumberFxn(s)
        birthNumberMat = birthNumber*((np.sqrt(dist)+62.0)/140)**(1.0/3.0)
        deathLimitMat = deathLimit*(1.0 + (np.sqrt(dist)-73)/1000)
        
        print("Computing Iteration",s)
        print("Birth number:",birthNumber,"Death limit:",deathLimit)

        
        changed = 0
        born = 0
        for idx in itertools.product(np.arange(size), 
                                     repeat=dimensions):
            if not boolLattice[idx]:
                continue

            # Compute Life game
            count = countAliveNeightbors(lattice0, idx)
            if lattice0[idx] == 1:
                if count < deathLimit:
                    changed -=1 
                    lattice1[idx] = 0
                else:
                    lattice1[idx] = 1
            else:
                if count > birthNumberMat[idx]:
                    changed +=1 
                    born += 1
                    lattice1[idx] = 1
                else:
                    lattice1[idx] = 0

        # Swap lattices
        print("Changed:",changed, "Born:", born)
        print("Alive0:",np.sum(lattice0), "(%0.1f%%)"%(100.0*np.sum(lattice0)/totalCytoplasmic))
        print("Alive1:",np.sum(lattice1), "(%0.1f%%)"%(100.0*np.sum(lattice1)/totalCytoplasmic))
        tmpLattice = lattice0
        lattice0 = lattice1
        lattice1 = tmpLattice
        
    return lattice0



def geometryModel(sim,
                  cellCenter = lm.point(*micron(9,9,9)),
                  cellRadius = micron(8.9),
                  membraneThickness= micron(0.128),
                  nuclSize= micron(4.15),
                  speckleRadius=micron(0.35),
                  cajalRadius=micron(0.5),
                  poreRadius=micron(0.083),
                  n_cajals=4,
                  n_speckles=20,
                  n_NPCs=1515,
                  n_mito=2000,
                  fb = 0.9,
                  fd = 0.8,
                  steps=11,
                  lsize=288,
                  buildER=True,
                  limits=[lambda x: x<= 65.0**2,lambda x: x> 139.0**2]
    ):
    '''
    all length units are in micrometers
    Definition of the parameters:
    cellCenter - the coordinates of the center of the cell in simulation box
    membraneThickness - the thickness of the plasma membrane
    nuclSize - the radius of the spherical nucleus
    speckleRadius - the radius of the individual speckles 
    cajalRadius - the radius of the individual Cajal bodies
    poreRadius - The radius of the individual nuclear pore complexes (NPCs)
    n_cajals - number of Cajal bodies implemented as spheres
    n_speckles - number of nuclear speckles implemented as spheres
    n_NPCs - number of NPCs on the nuclear envelope
    n_mito - number of mitochondria implemented as spherocylinders
    fb - birth threshold in the cellular automata algorithm for the ER generation
    fd - death threshold in the cellular automata algorithm for the ER generation
    steps - number of steps in the cellular automata algorithm for the ER generation
    lsize - size of the total simulation lattice in the cellular automata algorithm for the ER generation
    buildER - specify whether or not generate a new ER lattice or use an existing one
    limits - the upper and lower boundaries determining where ER should fit
    '''

    ################ Create cell components ##########################
    CellWall  = lm.Sphere(cellCenter, cellRadius, sim.siteTypes['CellWall'])
    Cytoplasm = lm.Sphere(cellCenter, cellRadius-membraneThickness, sim.siteTypes['Cytoplasm'])             
    NucEnv    = lm.Sphere(cellCenter, nuclSize, sim.siteTypes['NucEnv'])
    Nucleus   = lm.Sphere(cellCenter, nuclSize-membraneThickness, sim.siteTypes['Nucleus'])
    CellWall.thisown = 0
    Cytoplasm.thisown = 0
    NucEnv.thisown = 0
    Nucleus.thisown = 0

    ################ Cajal body: "n_cajal" with radius "cajalRadius" #############################
    dist = [] ; listdist = []
    Cajalpre = []
    accepted = []
    np_accepted = []
    cc = 1 ; point = 0
    while cc < n_cajals + 1:
        x,y,z = np.random.uniform(-0.5,0.5,3)
        Rp = nuclSize - cajalRadius - micron(0.1)  
        xpp = x*Rp + cellCenter.x
        ypp = y*Rp + cellCenter.y  
        zpp = z*Rp + cellCenter.z 
        caj = [xpp,ypp,zpp]
        listdist = []
        point += 1
        if (point == 1 ):
            accepted.append(caj)
            np_accepted = np.array(accepted)
        for i in range (0, cc):
            dist =  np.sqrt(np.sum((np.array(caj)-np_accepted[i])**2))
            listdist.append(dist) 
        np_dist = np.array(listdist)
        if(np.all(np_dist > 2*cajalRadius + micron(0.05))):
            accepted.append(caj)
            np_accepted = np.array(accepted)
            listdist = []
            cc += 1
    NCajal = np.array(accepted)

    cc = [] 
    Cajal = lm.UnionSet(sim.siteTypes['Cajal'])
    Cajal.thisown = 0
    for k in range(0, n_cajals):
        w = lm.Sphere(lm.point(NCajal[k][0],
                               NCajal[k][1],
                               NCajal[k][2]),
                      cajalRadius,
                      sim.siteTypes['Cajal'])
        w.thisown=0
        cc.append(w)
        Cajal.addShape(w)
    print("Cajal done!")

    ################ Nuclear Speckles: "n_speckles" with radius "speckleRadius" ##########################
    counter = 0
    taken = 0 
    listdist = [] ; np_acceptspeck = [] ; acceptspeck = [] 
    while counter < n_speckles + 1:
        u = np.random.uniform(-1, 1)
        thet = np.random.uniform(0, 2.0*math.pi)
        rad = np.random.uniform(speckleRadius*1e6 , nuclSize*1e6 - (speckleRadius*1e6 + 0.1))
        xm = micron(rad*math.sqrt(1-u**2)*np.cos(thet)) + cellCenter.x
        ym = micron(rad*math.sqrt(1-u**2)*np.sin(thet)) + cellCenter.y
        zm = micron(rad*u) + cellCenter.z
        point = [xm, ym, zm]
        listdist = [] ; dist = 0 ; np_dist = []
        taken += 1
        if ( taken == 1 ):
            acceptspeck.append(point)
            np_acceptspeck = np.array(acceptspeck)
        for i in range (0, counter):
            dist =  np.sqrt(np.sum((point-np_acceptspeck[i])**2))
            listdist.append(dist)
        np_dist = np.array(listdist)
        tt = 0
        if(np.all(np_dist > 2*speckleRadius + micron(0.1))):
            for k in range (0, n_cajals):
                if (np.sqrt(np.sum((np.array(point)-NCajal[k])**2)) > speckleRadius + cajalRadius + micron(0.1)):
                    tt += 1
            if (tt == 4):
                acceptspeck.append(point)
                np_acceptspeck = np.array(acceptspeck)
                counter += 1
                listdist = []
    NSpeck = np.array(acceptspeck)

    collection = []
    j = 0
    Speckle = lm.UnionSet(sim.siteTypes['Speckle'])
    Speckle.thisown = 0
    for k in range(0, n_speckles):
        q = lm.Sphere(lm.point(NSpeck[k][j],
                               NSpeck[k][j+1],
                               NSpeck[k][j+2]),
                      speckleRadius,
                      sim.siteTypes['Speckle'])
        q.thisown=0
        collection.append(q)
        Speckle.addShape(q)
    print("Speckles done!")

    ################# NPC: with the number of "n_NPCs" #########################
    NPCpre = []
    taken = 0

    ldist = [] ; npcdist = 0 ; np_npcdist = []
    acceptnpc = []
    npc = 0

    while npc < n_NPCs + 1:   
        x,y,z = np.random.uniform(-0.5,0.5,3)
        r = math.sqrt(x**2 + y**2 + z**2) 
        Rp = (nuclSize - sim.latticeSpacing)/r
        xpp = x*Rp + cellCenter.x
        ypp = y*Rp + cellCenter.y
        zpp = z*Rp + cellCenter.z
        pores = [xpp,ypp,zpp]
        NPCpre.append(pores)
        ldist = [] ; np_npcdist = []
        taken += 1
        if (taken == 1 ):
            acceptnpc.append(pores)
            np_acceptnpc = np.array(acceptnpc)
        for i in range (0, npc):
            npcdist =  np.sqrt(np.sum((np.array(pores)-np_acceptnpc[i])**2))
            ldist.append(npcdist)
        np_npcdist = np.array(ldist)
        if(np.all(np_npcdist > 0.2e-6)):
            acceptnpc.append(pores)
            np_acceptnpc = np.array(acceptnpc)
            ldist = []
            npc += 1

    NP = np.array(NPCpre)

    collect = []
    j = 0
    NPC = lm.UnionSet(sim.siteTypes['NPC'])
    NPC.thisown = 0
    for k in range(0, n_NPCs):
        s = lm.Sphere(lm.point(NP[k][j],
                               NP[k][j+1],
                               NP[k][j+2]),
                      poreRadius,
                      sim.siteTypes['NPC'])
        s.thisown=0
        collect.append(s)
        NPC.addShape(s)

    print("NPCs done!")

    ####################### Golgi Apparatus  #################################################
    m1 = lm.Sphere(lm.point(*micron(15,15,15)), micron(3.5), sim.siteTypes['Golgi'])             
    m11 = lm.Sphere(lm.point(*micron(15,15,15)), micron(3.308), sim.siteTypes['Golgi'])             
    m12 = lm.Sphere(lm.point(*micron(15,15,15)), micron(3.116), sim.siteTypes['Golgi'])             
    m13 = lm.Sphere(lm.point(*micron(15,15,15)), micron(2.924), sim.siteTypes['Golgi'])             
    m14 = lm.Sphere(lm.point(*micron(15,15,15)), micron(2.732), sim.siteTypes['Golgi'])             

    m2 = lm.Cone(cellCenter, micron(2), micron(6.5), sim.siteTypes['Golgi'], lm.vector(1.0,1.0,1.0))
    m22 = lm.Cone(cellCenter, micron(2.5), micron(8), sim.siteTypes['Golgi'], lm.vector(1.0,1.0,1.0))
    m23 = lm.Cone(cellCenter, micron(3), micron(9.5), sim.siteTypes['Golgi'], lm.vector(1.0,1.0,1.0))
    m24 = lm.Cone(cellCenter, micron(3.5), micron(10.5), sim.siteTypes['Golgi'], lm.vector(1.0,1.0,1.0))


    Gol1= lm.Difference(m1, m11, sim.siteTypes['Golgi'])
    Gol2= lm.Difference(m11, m12, sim.siteTypes['Golgi'])
    Gol3= lm.Difference(m12, m13, sim.siteTypes['Golgi'])
    Gol4= lm.Difference(m13, m14, sim.siteTypes['Golgi'])

    G1 = lm.Intersection(Gol1, m2, sim.siteTypes['Golgi'])
    G2 = lm.Intersection(Gol2, m22, sim.siteTypes['Golgi'])
    G3 = lm.Intersection(Gol3, m23, sim.siteTypes['Golgi'])
    G4 = lm.Intersection(Gol4, m24, sim.siteTypes['Golgi'])

    G11 = lm.Union(G1, G2, sim.siteTypes['Golgi'])
    G22 = lm.Union(G3, G4, sim.siteTypes['Golgi'])
    Golgi = lm.Union(G11, G22, sim.siteTypes['Golgi'])
    Golgi.thisown = 0
    print("Golgi done!")

    ######################## Mitochondria with number of "n_mito" ###############################3
    mit = []
    pp1 = []
    pp2 = []
    Mito = lm.UnionSet(sim.siteTypes['Mito'])
    Mito.thisown=0
    for i in range(1, n_mito+1):
        u = np.random.uniform(-1, 1)
        thet = np.random.uniform(0, 2.0*math.pi)
        rad = np.random.uniform(4.706, 8.5)
        xm = micron(rad*math.sqrt(1-u**2)*np.cos(thet)) + micron(9) 
        ym = micron(rad*math.sqrt(1-u**2)*np.sin(thet)) + micron(9) 
        zm = micron(rad*u) + micron(9)
        centers = [xm, ym, zm]
        if i < 666 :
            point1 = [xm - micron(0.125), ym-micron(0.125) , zm]
            point2 = [xm + micron(0.125), ym+micron(0.125) , zm]
        elif  i >= 666 and i < 1332 :
            point1 = [xm, ym - micron(0.125), zm- micron(0.125)]
            point2 = [xm, ym + micron(0.125), zm+ micron(0.125)]
        elif  i >= 1332 :
            point1 = [xm - micron(0.125), ym , zm - micron(0.125)]
            point2 = [xm + micron(0.125), ym , zm + micron(0.125)]
            
        pp1.append(point1)
        pp2.append(point2)
    mp1 = np.array(pp1)
    mp2 = np.array(pp2)

    for h in range(0, n_mito):
        m = lm.Capsule(lm.point(mp1[h][0], mp1[h][1], mp1[h][2]),
                       lm.point(mp2[h][0], mp2[h][1], mp2[h][2]),
                       micron(0.249),
                       sim.siteTypes['Mito'])
        m.thisown = 0
        mit.append(m)
        Mito.addShape(m)
    print("Mito done!")



    # Add all geometries to the simulation box
    sim.lm_builder.addRegion(CellWall)
    sim.lm_builder.addRegion(Cytoplasm)
    sim.lm_builder.addRegion(NucEnv)
    sim.lm_builder.addRegion(Nucleus)
    sim.lm_builder.addRegion(NPC)
    sim.lm_builder.addRegion(Cajal)
    sim.lm_builder.addRegion(Speckle)
    sim.lm_builder.addRegion(Mito)
    sim.lm_builder.addRegion(Golgi)

    ######################### ER with parameters described above #########################################
    print("Building ER")

    if buildER:
        deathLimitFxn  = lambda s: int(fd*27*(steps-0.25*s)**2/steps**2)
        birthNumberFxn = lambda s: int(fb*27) 
        im = createERCellularAutomaton(lsize, 
                                        dimensions=3,
                                       survivalRate=0.85, 
                                       deathLimitFxn  = deathLimitFxn,
                                       birthNumberFxn = birthNumberFxn, 
                                       steps=steps,
                                       limits = limits)
        
        # Save the grid to a 3D lattice as a numpy array; this can be loaded with np.load()
        np.save("ER.npy", np.array(im, dtype="float"))
        print("Done building ER")
     
    print("Loading ER")
    ER = np.load ("ER.npy")
    print("Discretizing")
    lattice = sim.getLattice() 
    print("Setting ER")
    for x in range(ER.shape[0]):
        for y in range(ER.shape[1]):
            for z in range(ER.shape[2]):
                idx = (x,y,z)
                if ER[x,y,z]:
                    sim.setLatticeSite(idx, "ER")
 
    print("ER Done")
 
