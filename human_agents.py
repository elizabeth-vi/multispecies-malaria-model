from index_names import Compartments, Species, Transitions, Treatments
import random

class Pathogen(object):
    """
    the `state' of this Disease, wrt
    :param.current: current compartment.current (S, I, or R)
    :param.time: time of next event
    :param next: compartment transition to at time.time
    """
    __slots__ = ('current', 'time', 'next')

    def __init__(self, transition_time):
        """
        :param transition_time: time of transition from default allocation to susceptible class
        """
        self.current = Compartments.S
        self.time = transition_time  # int(time_end + 1)  # effectively infinity
        self.next = Compartments.S  # no default change from this compartment type


class Agent(object):
    """
    individuals in the human population
    """
    __slots__ = 'state', 'memory', 'sex', 'G6PD_level', 'transition_overrides', 'G_treatment', 'param_overrides'

    def __init__(self, transition_time):
        """
        initialise Agent as per the Pathogen object
        :param transition_time: time of transition from default allocation to susceptible class
        """
        self.state = (Pathogen(transition_time=transition_time), Pathogen(transition_time=transition_time))  # todo: not hardcode this for 2 pathogens
        self.memory = ([], [])  # memory by species

        #Probabilities of [Severe, Intermediate, Normal] G6PD enzyme levels for XX and XY chromosomes
        p_XY = [0.137, 0, 1-0.137]
        p_XX = [0.05, 0.158, (1-0.05-0.158)]
        
        #Assign sex chromosomes
        r = random.random()
        self.sex = random.randint(0,1) # XY=0, XX=1

        #Assign true G6PD status

        #Ineligible
        if r<0.15: #Use a different r
            G6PD_level = None 

        #male XY
        if self.sex == 0: #XY
            if r < p_XY[0]:   #Hemizygous male
                G6PD_level = 0
            else:               #Normal male
                G6PD_level = 1
        
        # female XX
        if self.sex == 1: #XX
            if r < p_XX[0]:         #Homozygous female
                G6PD_level = 0
            elif r < sum(p_XX[0:2]):#Heterozygous female
                G6PD_level = 0.5
            else:                       #Normal female
                G6PD_level = 1

        self.G6PD_level = G6PD_level

        #Treatment is a subclass of "G"
        self.G_treatment = None


        #Override population default transition if a particular value is set here
        # Ie transition_overrides[Transitions.L_whatever] = 0.4
        self.transition_overrides = {}
        self.param_overrides = {}

    #For radical cure treatment specifically
    def start_G_treatment(self, params, treatment):

        G_params = ["pTfP", "pP", "psi","pG","pN","c"] #List of potential parameters that will change with varying radical cure treatment

        self.G_treatment = treatment

        #Rate out of G for specified treatment
        self.transition_overrides[Transitions.G_done] = params[treatment].psi[Species.vivax]
        
        #Add params to override defaults
        for param in G_params:
            self.param_overrides[param] = getattr(params[treatment],param)


    def finish_G_treatment (self):

        G_params = ["pTfP", "pP", "psi","pG","pN","c"] #List of potential parameters that will change with varying radical cure treatment

        del self.transition_overrides[Transitions.G_done]
        for param in G_params:
            del self.param_overrides[param]
        self.G_treatment = None
        

    #EDIT TO GIVE PROPER DISTRIBUTION
    def G6PD_test(self):
        
        #sex -add
        sensitivity = [0.99, 0.44] #Sensitivity for [Severe, Intermediate] True Pos 
        specificity = [0.99, 0.97] #Specificity for [Males, Females] True neg

        #ineligible
        if self.G6PD_level == None:
            return None

        r = random.random()
        if self.G6PD_level < 0.3:
            if r < sensitivity[0]:
                return 0
        elif self.G6PD_level < 0.7:
            if r < sensitivity[1]:
                return 0.5
        else:
            return 1
        return 1
    
    def assign_G_treatment(self, params, policy):
        """
        Assign and commence radical cure treatment, based on current treatment policy and G6PD status
        """

        

        treatment = self.choose_treatment(policy)
        self.start_G_treatment(params, treatment)


    def choose_treatment(self,policy):
        
        G6PD_test = self.G6PD_test()
        if policy == Treatments.Baseline:
            return Treatments.Baseline
        
        elif policy == Treatments.PLD:
            if G6PD_test <= 0.3:
                return Treatments.PG6PD
            else:
                return Treatments.PLD
            
        elif policy == Treatments.PHD:
            if G6PD_test <= 0.7:
                return Treatments.PG6PD
            else:
                return Treatments.PHD
            
        elif policy == Treatments.Taf:
            if G6PD_test <= 0.7:
                return Treatments.PG6PD
            else:
                return Treatments.Taf
        
        print("Error: Policy is not listed")
        quit()
                         
                         
                         
                         
                         
    # def give_treatment(self, treatment):
    #     pass

    #Get a specific transition probability given the global table and this agent's particular overrides
    #transition = transition ID from enum
    # def get_transition_rate(self, transition, global_transitions):
    #     # For each transition category
    #     if (self.transition_overrides.has_key(transition)):
    #         return self.transition_overrides[transition]
        
    #     return global_transitions[transition]



    
    #Not yet used
# class Treatment():
#     def __init__(self, params = {}, transition_overrides = {}, load_dict = x):
#         self.params = params
#         self.transition_overrides = {}