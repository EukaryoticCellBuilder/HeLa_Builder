

def reactionModel(sim,
                  # Reaction rates 
                  Ktr=0.0034,
                  Kon=4.66e7,
                  Koff=0.062,
                  Ksplice=0.05
                  ):
    
    ###################################################################
    # Adjust rates by volume #
    # The reaction rates have to be in the unit of 1/s. 
    # Therefore, in most cases the experimental rates have to be scaled.
    ####################################################################
    scale = 157863.12  # = 6.022e23 * 64^3 e-24

    Kon = Kon/scale
    
    ##########################
    # Get handles to regions #
    ##########################
    nucleus   = sim.modifyRegion('Nucleus')
    cytoplasm = sim.modifyRegion('Cytoplasm')
    npc       = sim.modifyRegion('NPC')
    speckle   = sim.modifyRegion('Speckle')

    #####################
    # Add all Reactions #
    #####################
    nucleus.addReaction(reactant=('gene'), product=('gene', 'premRNA'), rate= Ktr) 
    nucleus.addReaction(reactant=('premRNA', 'S'), product='SpremRNA', rate= Kon)
    nucleus.addReaction(reactant=('SpremRNA'), product=('premRNA', 'S'), rate= Koff)
    nucleus.addReaction(reactant=('SpremRNA'), product=('S', 'mRNA'), rate= Ksplice) 
    print("Reactions are set!")
    
def particleModel(sim,
                  n_gene=200,
                  n_premRNA=0,
                  n_S=1000,
                  n_mRNA=0,
                  n_SpremRNA=0
                  ):

    ###################################################
    # Add all Species with their corresponding counts #
    ###################################################

    sim.addParticles(species='gene', region='Nucleus', count= n_gene)  
    sim.addParticles(species='premRNA', region='Nucleus', count= n_premRNA)
    sim.addParticles(species='S', region='Nucleus', count= n_S)
    sim.addParticles(species='SpremRNA', region='Nucleus', count=n_SpremRNA)
    sim.addParticles(species='mRNA', region='Nucleus', count=n_mRNA)
    
    print("Particles were added!")
