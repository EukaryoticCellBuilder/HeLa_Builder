

def diffusionModel(sim,
                   d_premRNA = 6.1e-13,
                   d_S = 2.6e-13,
                   d_SpremRNA = 2.6e-13,
                   d_mRNA = 6.1e-13
                   ):
    
    ##########################
    # Get handles to regions #
    ##########################
    nucleus   = sim.modifyRegion('Nucleus')
    speckle   = sim.modifyRegion('Speckle')
    npc       = sim.modifyRegion('NPC')
    cytoplasm = sim.modifyRegion('Cytoplasm')
    cajal     = sim.modifyRegion('Cajal')
    
    #######################
    # Set diffusion rates #
    #######################
    nucleus.setDiffusionRate(species='premRNA', rate= d_premRNA) 
    nucleus.setDiffusionRate(species='S', rate= d_S) 
    speckle.setDiffusionRate(species='premRNA', rate= d_premRNA) 
    speckle.setDiffusionRate(species='S', rate= d_S) 

    nucleus.setDiffusionRate(species='SpremRNA', rate= d_SpremRNA)   
    speckle.setDiffusionRate(species='SpremRNA', rate= d_SpremRNA)   

    nucleus.setDiffusionRate(species='mRNA', rate= d_mRNA)   
    speckle.setDiffusionRate(species='mRNA', rate= d_mRNA)   

    #setTwoWayTransitionRate

    sim.setTwoWayTransitionRate(species='premRNA', one='Nucleus', two='Speckle', rate= d_premRNA)
    sim.setTwoWayTransitionRate(species='premRNA', one='Speckle', two='Nucleus', rate= 0.01*d_premRNA)

    sim.setTwoWayTransitionRate(species='S', one='Nucleus', two='Speckle', rate= d_S)
    sim.setTwoWayTransitionRate(species='S', one='Speckle', two='Nucleus', rate= 0.01*d_S)

    sim.setTwoWayTransitionRate(species='SpremRNA', one='Nucleus', two='Speckle', rate= d_SpremRNA)
    sim.setTwoWayTransitionRate(species='SpremRNA', one='Speckle', two='Nucleus', rate= 0.01*d_SpremRNA)

    sim.setTwoWayTransitionRate(species='mRNA', one='Nucleus', two='Speckle', rate= d_mRNA)
    sim.setTwoWayTransitionRate(species='mRNA', one='Speckle', two='Nucleus', rate= d_mRNA)
    print ("Diffusions are set!")	
	
