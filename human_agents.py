from index_names import Compartments
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
    __slots__ = 'state', 'memory', 'sex', 'G6PD_level'

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
    
    def give_treatment(self, treatment):
        pass


    # def G6PD_dist():
        
    #     #Probabilities for male and female, [Severe, Intermediate, Normal]
    #     p_male = [0.137, 0, 1-0.137]
    #     p_female = [0.05, 0.158, (1-0.05-0.158)]
        

    #     r = random.random()
    #     sex = random.randint(0,1) # Male=0, Female=1

    #     #male XY
    #     if sex == 0: 
    #         if r < p_male[0]:   #Hemizygous male
    #             G6PD_level = 0
    #         else:               #Normal male
    #             G6PD_level = 1
        
    #     # female XX
    #     if sex == 1: 
    #         if r < p_female[0]:         #Homozygous female
    #             G6PD_level = 0
    #         elif r < sum(p_female[0:2]):#Heterozygous female
    #             G6PD_level = 0.5
    #         else:                       #Normal female
    #             G6PD_level = 1
        
    #     return G6PD_level, sex
    
    # def set_G6PD(agent):
    #     agent.G6PD_level, sex = G6PD_dist()
    #     agent.G6PD_test = G6PD_test(agent.G6PD_level, sex)

    #     return

    
