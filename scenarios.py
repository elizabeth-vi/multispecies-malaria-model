from index_names import Treatments

scenario 0: 
policy_G6PD_max = [0.7, 1.0]
policy_Treatment = [Treatments.PG6PD, Treatments.Taf]

class Scenarios(object):
    def __init__(self,index):
        if index == 0: #baseline
            self.policy_G6PD_max = [1.0]
            self.policy_treatment = [Treatments.Baseline]
        elif index == 1: #rec PLD
            self.policy.G6PD_max = [0.3, 1.0]
            self.policy.treatment = [Treatments.PG6PD, Treatments.PLD]
        elif index == 2: #Rec PHD
            self.policy_G6PD_max = [0.7, 1.0]
            self.policy_treatment = [Treatments.PG6PD, Treatments.PHD]
        elif index == 3: #Rec Taf
            self.policy_G6PD_max = [0.7, 1.0]
            self.policy_treatment = [Treatments.PG6PD, Treatments.Taf]
