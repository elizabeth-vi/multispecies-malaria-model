from index_names import Compartments, Species, Transitions, Treatments
import random
import numpy as np


class Hypnozoite(object):
    """
    The state of the latent-stage hypnozoites in the liver
    :param number: number of hypnozoites currently in the liver
    :param time_next: time of next natural hypnozoite event (death, activation)
    :param alpha: rate of hypnozoite activation
    :param mu: rate of hypnozoite death
    """

    __slots__ = ('number', 'time_next', 'alpha', 'mu')

    def __init__(self, hyp_init, transition_time, hyp_params, species):
        self.number = hyp_init #integer
        self.time_next = transition_time #number
        self.alpha = hyp_params.alpha_hyp[species] #number
        self.mu = hyp_params.mu_hyp[species] #number

    def infect(self, mean, current_time):
        """
        Add latent hypnozoites upon a new primary infection event
        """

        # ni = self.number
        self.number += np.random.geometric(p=1/(mean+1)) - 1
        # nj = self.number
        # if ni>0:
        #     print("time = "+str(current_time)+str([ni, nj]))
        self.update_transition_time(current_time)

    def treat(self, p_death, current_time):
        # print("Treating. N_hyp start = "+str(self.number))
        self.number = np.random.binomial(self.number,1-p_death)
        self.update_transition_time(current_time)
        # print("N_hyp ending "+str(self.number)+", next time = ",self.time_next)

    def update_transition_time(self, current_time):
        rate = self.number*(self.alpha+self.mu)

        if rate != 0: #avoid divide by zero error
            time = current_time + int(max(1, round(random.expovariate(rate), ndigits=0)))
        else:
            time = 1e10 #arbitary long time
        self.time_next = time
        # print(time)

    
    def update(self, person, species, current_time, p_clinical): 
        """
        Update hypnozoite status after natural event (death / activation)
        Updates human compartment and hypnozoite transition time
        """
        assert self.number >= 1, "cannot have hypnozoite birth/death with no hypnozoites"
        # print(self.number)

        self.number -= 1 #decrement by one
        self.update_transition_time(current_time)

        if random.random() < self.alpha/(self.alpha + self.mu): #activation
            activation = True

            current_compartment = person.state[species].current
            if current_compartment in [Compartments.S, Compartments.A, Compartments.R]:
                if random.random() < p_clinical: #clinical infection
                    person.state[species].current = Compartments.I #transition occurs
                else: #asymptomatic
                    person.state[species].current = Compartments.A #transition occurs
            
        else: #hypnozoite death
            activation = False
            
        return activation
        
class Pathogen(object):
    """
    the `state' of this Disease, wrt
    :param.current: current compartment.current (S, I, or R)
    :param.time: time of next event
    :param next: compartment transition to at time.time
    """
    __slots__ = ('current', 'time', 'next', 'hypnozoites')

    def __init__(self, transition_time, hyp_params, species):
        """
        :param transition_time: time of transition from default allocation to susceptible class
        """
        self.current = Compartments.S
        self.time = transition_time  # int(time_end + 1)  # effectively infinity
        self.next = Compartments.S  # no default change from this compartment type
        self.hypnozoites = Hypnozoite(hyp_init = 0, transition_time=transition_time, hyp_params = hyp_params, species=species)

            



class Agent(object):
    """
    individuals in the human population
    """
    __slots__ = 'state', 'memory', 'sex', 'G6PD_level', 'transition_overrides', 'G_treatment', 'param_overrides', "recent_treatment"

    def __init__(self, transition_time, params):
        """
        initialise Agent as per the Pathogen object
        :param transition_time: time of transition from default allocation to susceptible class
        """
        self.state = (Pathogen(transition_time=transition_time, hyp_params=params, species = Species.falciparum), Pathogen(transition_time=transition_time, hyp_params=params, species=Species.vivax))  # todo: not hardcode this for 2 pathogens
        self.memory = ([], [])  # memory by species
        self.recent_treatment = - (params.relapse_monitor+1) #time of most recent treatment / health service visit, starts with no recent visit
        # hyp_init = [0,0] # hardcoded for 2 pathogens
        # self.hypnozoites = Hypnozoite(hyp_init = hyp_init, transition_time=transition_time, alpha = params.alpha_hyp, mu = params.mu_hyp)

        #Probabilities of [Severe, Intermediate, Normal] G6PD enzyme levels for XX and XY chromosomes
        p_XX = [0.05, 0.158, (1-0.05-0.158)]
        p_XY = [0.137, 0, 1-0.137]
        p_inel = [0.08, 0.01] #[XX, XY]
        
        
        #Assign sex chromosomes
        r1 = random.random()
        r2 = random.random()
        self.sex = random.randint(0,1) # XX=0, XY=1

        #Assign true G6PD status

        # female XX
        if self.sex == 0: #XX
            if r2 < p_XX[0]:                            #Homozygous female
                G6PD_level = random.uniform(0, 0.3)
            elif r2 < sum(p_XX[0:2]):                   #Heterozygous female
                G6PD_level = random.uniform(0.3, 0.7)
            else:                                       #Normal female
                G6PD_level = random.uniform(0.7, 1)

        #male XY
        elif self.sex == 1: #XY
            if r2 < p_XY[0]:                            #Hemizygous male
                G6PD_level = random.uniform(0, 0.3)
            else:                                       #Normal male
                G6PD_level = random.uniform(0.7, 1)

        #Ineligible - overrides above result if assigned
        if r1 < p_inel[self.sex]:
            G6PD_level = None 
        
        self.G6PD_level = G6PD_level

        #Treatment is a subclass of "G"
        self.G_treatment = None


        #Override population default transition if a particular value is set here
        # Ie transition_overrides[Transitions.L_whatever] = 0.4
        self.transition_overrides = {}
        self.param_overrides = {}

    #For radical cure treatment specifically
    def start_G_treatment(self, species, params, treatment, current_time):

        assert species == Species.vivax

        if treatment == Treatments.Baseline:
            return
        # G_params = ["pTfP", "pP", "psi","pG","pN","c","p_rad"] #List of potential parameters that will change with varying radical cure treatment
        # G_params = ["pTfP", "pP", "psi","pG","p_rad"] #List of potential parameters that will change with varying radical cure treatment

        self.G_treatment = treatment
        treatment_params = params.treatment_params[treatment]

        #Rate out of G for specified treatment
        # self.transition_overrides[Transitions.G_done] = params[treatment].psi[species]


        self.transition_overrides[Transitions.G_done] = treatment_params["psi"][species]
        
        #Add params to override defaults
        # for param in G_params:
        #     self.param_overrides[param] = getattr(params[treatment],param)
        for param in treatment_params:
            self.param_overrides[param] = treatment_params[param]

            if param == "pG":
                if self.G6PD_level < treatment_params["haem_threshold"]:
                    self.param_overrides[param] = np.add(treatment_params[param],treatment_params["p_haem"]*params.p_haem_death).tolist()

        #Treat hypnozoites
        # self.state[species].hypnozoites.treat(params[treatment].p_rad, current_time)
        self.state[species].hypnozoites.treat(treatment_params["p_rad"], current_time)


    def finish_G_treatment (self):

        # G_params = ["pTfP", "pP", "psi","pG","pN","c"] #List of potential parameters that will change with varying radical cure treatment

        # del self.transition_overrides[Transitions.G_done]
        # for param in G_params:
        #     del self.param_overrides[param]
        self.transition_overrides = {}
        self.param_overrides = {}
        self.G_treatment = None
        

    def G6PD_test(self):
        """
        Return a G6PD test result based on patient's actual G6PD levels
        """
        #Currently works off brackets but is adaptable to vary based on results within a bracket, if evidence supports
        
        #probability of getting a relevant test result
        testing_prob = [[1.00, 0, 0], #For severe G6PD deficinecy, probability of [severe, intermediate, normal] test result
                        [0.48, 0.42, 0.10], #For intermediate G6PD deficiency
                        [0.01, 0.04, 0.95]] #For normal G6PD levels

        #ineligible
        assert self.G6PD_level is not None, "giving G6PD test to ineligible person"

        r = random.random()

        #Determine severity bracket
        if self.G6PD_level < 0.3: #severe deficiency
            severity = 0
        elif self.G6PD_level < 0.7: #intermediate deficiency
            severity = 1
        else: #Normal levels
            severity = 2
        
        #Determine test result based on bracket
        if r < testing_prob[severity][0]: #severe test result
            test_res = 0
        elif r < testing_prob[severity][0] + testing_prob[severity][1]: #intermediate
            test_res = 0.5
        else: #normal
            test_res = 1

        return test_res
        

    
    def assign_G_treatment(self, params, policy, current_time):
        """
        Assign and commence radical cure treatment, based on current treatment policy and G6PD status
        """

        treatment = self.choose_treatment(policy)
        self.start_G_treatment(species=Species.vivax, params=params, treatment=treatment, current_time=current_time)


    def choose_treatment(self,policy):

        G6PD_test = self.G6PD_test()
        for i in range(len(policy["G6PD_maxes"])):
            if G6PD_test <= policy["G6PD_maxes"][i]:
                return policy["treatments"][i]
        
        print("Error: Policy "+str(policy)+" is not listed")
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