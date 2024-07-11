from index_names import Compartments, Transitions, Species, Treatments
import random
import math
from copy import copy

class Disease(object):
    """
    Keep track of human population level aggregates and logic for updating an Agent
    """
    __slots__ = ('pop_counts', 'total_humans', 'species', 'T_deaths', 'G_deaths', 'I_deaths', 'new_infections',
                 'relapses', 'num_relapses', 'relapse_recorded', 'new_inf_clinical', 'new_T', 'new_G')

    def __init__(self, malaria_species, initial_population_counts, time_end):
        """
        Each of the possible Disease states have been pre-allocated as slots, population level aggegrates recorded for use here
        :param malaria_species: which malaria species this instantiation is for
        :param initial_population_counts: initial human population in each compartment
        """
        self.pop_counts = initial_population_counts.copy()  # last one is for `just_died`
        self.species = malaria_species

        # things tracking for later analysis. i.e. `outputs of interest`
        self.T_deaths = []
        self.G_deaths = []
        self.I_deaths = []
        self.new_infections = []
        self.new_inf_clinical = []
        self.relapses = []
        self.num_relapses = [0]*time_end
        self.relapse_recorded = [0]*time_end
        self.new_T = []
        self.new_G = []

    def calculate_time(self, current_time, rate):
        if rate==0:
            return 99999999
        return current_time + int(max(1, round(random.expovariate(rate), ndigits=0)))
    
    def transition_table(self, person, current_time, event_rates, params, treatment_policy):
        """
        Place to specify transition table: Note the fact that a transition is/has just occurred was determined in self.update() and infection events are determined in self.infect_me()
        :param person: current agent
        :param current_time: present simulation time step (integer)
        :param event_rates: transition rates pre-calculated from stochastic_sir_event_rates() using population level stats from previous time step
        :param params: the model parameters
        :return: nothing explicit, but next agent compartment and.time updated
        """
        current_status = person.state[self.species].current  # just transitioned to this
        # params = params_list[treatment_policy]
        
        #Update transition table for person-specific override value (e.g. different radical cure treatment durations)
        for transition in person.transition_overrides:
            event_rates[transition] = person.transition_overrides[transition]

        #Update params for person-specific override value (e.g. different radical cure parameters)
        if person.param_overrides:
            params = copy(params)
            event_rates = copy(event_rates)

            #Update transition table for person-specific override value (e.g. different radical cure treatment durations)
            for transition in person.transition_overrides:
                event_rates[transition] = person.transition_overrides[transition]

            #Update params for person-specific override value (e.g. different radical cure parameters)
            for param in person.param_overrides:
                setattr(params, param, person.param_overrides[param])


        # determine what next transition is, and when (excepting infection)

        if current_status == Compartments.S:
            _next = Compartments.S
            _time = params.time_end + 1  # stay susceptible unless infected

        elif current_status == Compartments.I:
            # I->death/A or I->T/G happens first
            if random.random() < params.pTreat[self.species]: #seek treatment
                time_treat = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.I_treat])
                time_next = 3 * params.time_end #i.e. never
            else: #naturally transition
                time_next = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.I_next])
                time_treat = 3 * params.time_end #i.e. never

            # MDA should affect the I compartment too
            if event_rates[Transitions.A_treat] == 0:  # avoiding a divide by zero error
                time_mda = 3 * params.time_end  # i.e. never, but making sure likely to be bigger than time_recover
            else:
                time_mda = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.A_treat])

            # if transition before mda or treatment progress
            if time_next < time_treat and time_next < time_mda:
                _time = time_next
                if random.random() < params.pI[self.species]:  # malaria death
                    self.I_deaths.append(_time)
                    _next = Compartments.just_died
                else:
                    _next = Compartments.A
            elif time_treat<time_mda:
                _time = time_treat

                # MASKING LOGIC
                if (person.state[Species.falciparum].current in [Compartments.A, Compartments.I]) and (person.state[Species.vivax].current in [Compartments.A, Compartments.I]): #if mixed infection masking is possible
                    
                    if random.random() < params.mask_prob: #treat with prob treat with T given masking may occur
                        _next = Compartments.T
                    else: # if masking doesn't happen treat as usual
                        _next = Compartments.G

                else:
                    #Check if prescribed radical cure treatment (if eligible)
                    if random.random() < 1 - params.pN[self.species]:

                        if person.G6PD_level is None: #ineligible for G6PD treatment, no test administered
                            _next = Compartments.T
                        else: 
                            _next = Compartments.G
                    else:
                        _next = Compartments.T

            else:
                _time = time_mda

                # MASKING LOGIC
                if (person.state[Species.falciparum].current in [Compartments.I, Compartments.A]) and (person.state[Species.vivax].current in [Compartments.I, Compartments.A]): #if mixed infection masking is possible
                    
                    if random.random() < params.mask_prob_mda: # treat with prob treat with T given masking may occur
                        _next = Compartments.T
                    else: # if masking doesn't happen treat as usual
                        _next = Compartments.G

                else:
                    if random.random() < params.pN_mda[self.species]:
                        _next = Compartments.T
                    else:
                        _next = Compartments.G

        elif current_status == Compartments.A:
            # A-> L or R, or T or G
            time_recover = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.A_recover])
            if event_rates[Transitions.A_treat] == 0:  # avoiding a divide by zero error
                time_mda = 3 * params.time_end  # i.e. never, but making sure likely to be bigger than time_recover
            else:
                time_mda = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.A_treat])

            if time_recover < time_mda:
                _time = time_recover
                # if random.random() < params.ph[self.species]:
                #     _next = Compartments.L
                #     assert self.species == Species.vivax, "falciparum being allocated to L in Disease.transitions_table from A"
                # else:
                #     _next = Compartments.R
                _next = Compartments.R
            else:
                # note that this implies that all MDA events are ACT # put in masking
                _time = time_mda
                if random.random() < params.pN_mda[self.species]:
                    _next = Compartments.T
                else:
                    _next = Compartments.G

        elif current_status == Compartments.R:
            # unless infected before this happens, immunity wanes
            _next = Compartments.S
            _time = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.R_S])

        # elif current_status == Compartments.L:
        #     assert self.species == Species.vivax, "falciparum was just assigned to `L`"
        #     # either relapse (I or A) or hypnozoite death (S) (unless infected in the meantime)
        #     time_relapse = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.L_relapse])
        #     time_recover = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.L_S])
        #     if time_relapse < time_recover:
        #         _time = time_relapse
        #         if random.random() < params.pL[self.species]:
        #             _next = Compartments.I
        #         else:
        #             _next = Compartments.A
        #     else:
        #         _time = time_recover
        #         _next = Compartments.S

        elif current_status == Compartments.T:
            self.new_T.append(current_time)  # just got to treatment
            _time = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.T_done])
            if random.random() < params.pT[self.species]:  # malaria death
                self.T_deaths.append(_time)
                _next = Compartments.just_died  # a flag in run_me to update all pathogen compartments and human population size as recoreded in Mozzies
            else:  # 1 - pT
                if random.random() < params.pTfA[self.species]:  # treatment failure
                    _next = Compartments.A
                else:  # 1 - pTfA
                    # if random.random() < params.pA[self.species]:  # recover with hypnozoites
                    #     _next = Compartments.L
                    #     assert self.species == Species.vivax, "shouldn't be assigning compartment L to falciparum in Disease.transition_table, from T"
                    # else:  # 1 - pA :: recover with no hypnozoites
                    #     _next = Compartments.R
                    _next = Compartments.R

        elif current_status == Compartments.G:
            self.new_G.append(current_time)  # record now getting `G` treatment 

            _time = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.G_done])
            
            if random.random() < params.pG[self.species]:  # malaria death
                self.G_deaths.append(_time)
                _next = Compartments.just_died  # remove from compartments do anything with
            else:  # 1 - pG
                if random.random() < params.pTfP[self.species]:  # treatment failure
                    _next = Compartments.A
                else:  # 1 - pTfP
                    # if random.random() < params.pP[self.species]:  # recover with hypnozoites
                    #     _next = Compartments.L
                    #     assert self.species == Species.vivax, "shouldn't be assigning compartment L to falciparum in Disease.transition_table, from compartment G"
                    # else:  # 1 - pP :: recover without hypnozoites
                    #     _next = Compartments.R
                    _next = Compartments.R

        elif current_status == Compartments.just_died:  # person just died
            _next = Compartments.dead
            _time = current_time  # using this to catch if not moved to `dead` and hence agent and disease updated for all pathogens in `run_me`
        else:
            raise ValueError("A compartment that doesn't exist has been assigned and sent to Disease.transmission_table()")

        # make the specified updates
        person.state[self.species].next = _next
        person.state[self.species].time = _time


    def triggering(self, person, current_time, event_rates, params):
        """
        Re-doing transition from current L to elsewhere due to recent falciparum infection, with increased nu_hat = Zf * nu
        :param person:
        :param current_time:
        :param event_rates:
        :param params:
        :return:
        """

        #Update transition table for person-specific override value (e.g. different radical cure treatment durations)
        for transition in person.transition_overrides:
            event_rates[transition] = person.transition_overrides[transition]

        # either relapse (I or A) or hypnozoite death (S) (unless infected in the meantime)
        time_relapse = self.calculate_time(current_time=current_time, rate=params.zf * event_rates[Transitions.L_relapse])
        time_recover = self.calculate_time(current_time=current_time, rate=event_rates[Transitions.L_S])
        if time_relapse < time_recover:
            _time = time_relapse
            if random.random() < params.pL[self.species]:
                _next = Compartments.I
            else:
                _next = Compartments.A
        else:
            _time = time_recover
            _next = Compartments.S

        person.state[Species.vivax].next = _next
        person.state[Species.vivax].time = _time

    def update(self, person, current_time, event_rates, params, policy, event_rates_other_sp = None):
        """
        check if time to next event is now, or infection occurs
        :return: updates state
        """
        assert person.state[self.species].time >= current_time
        current_compartment = person.state[self.species].current
        # params = params_list[policy]

        if person.param_overrides:
            params = copy(params)
            event_rates = copy(event_rates)

            #Update transition table for person-specific override value (e.g. different radical cure treatment durations)
            for transition in person.transition_overrides:
                event_rates[transition] = person.transition_overrides[transition]

            #Update params for person-specific override value (e.g. different radical cure parameters)
            for param in person.param_overrides:
                setattr(params, param, person.param_overrides[param])

        # STANDARD TRANSITIONS
        if person.state[self.species].time == current_time:  # an event has been scheduled for this time
            assert current_compartment not in [Compartments.dead, Compartments.just_died], "Transitioning dead person"
            # decrement population count for current state

            self.pop_counts[current_compartment] -= 1
            # increment population count for next state

            # self.pop_counts[person.state[self.species].next] += 1

            # if current_compartment == Compartments.L and (person.state[self.species].next in [Compartments.I, Compartments.A]):
                # self.relapses.append(current_time)  # record relapse event -- here and not in self.update as a new infection may occur in the meantime

            #start / assign relevant treatment and
            #record treatment time - to monitor recorded relapses
            if person.state[self.species].next in [Compartments.T, Compartments.G]:                
                person.assign_treatment(params, policy)
                person.recent_treatment = current_time

                #Check if person moved from G to T due to G6PD deficiency (treatment assigned above)
                if (person.treatment == params.policy["blood_stage"]) and (person.state[self.species].next == Compartments.G):
                    person.state[self.species].next = Compartments.T
                
                #start treatment and assign override params
                person.start_treatment(self.species, params, current_time)

            #If ending G treatment, undo treatment param changes
            if person.state[self.species].current in [Compartments.T, Compartments.G]:
                person.finish_treatment()


            # add triggering logic (triggering occurs after a pf recovery)
            #UPDATE TRIGGERING WITH HYPNOZOITES IN FUTURE
            if params.flag_triggering_by_pf and self.species == Species.falciparum and person.state[Species.vivax].current == Compartments.L and person.state[self.species].next == Compartments.R:
                self.triggering(person=person, current_time=current_time, event_rates=event_rates_other_sp, params=params)

            # continue as usual
            person.state[self.species].current = person.state[self.species].next  # transition occurs
            self.pop_counts[person.state[self.species].next] += 1
            self.transition_table(person=person, current_time=current_time, event_rates=event_rates, params=params, treatment_policy=policy)  # identify next transition to occur

        else:   # check for infection event
            if current_compartment == Compartments.S: 
                self.infect_me(person=person, rate_infection=event_rates[Transitions.S_inf], prob_I=params.pc[self.species], params=params, policy=policy, current_time=current_time, event_rates=event_rates)
            elif current_compartment == Compartments.R: #or Compartments.A
                self.infect_me(person=person, rate_infection=event_rates[Transitions.R_inf], prob_I=params.pR[self.species], params=params, policy=policy, current_time=current_time, event_rates=event_rates)
            # elif current_compartment == Compartments.L:
            #     assert self.species == Species.vivax, "falciparum person is in `L`"
            #     self.infect_me(person=person, rate_infection=event_rates[Transitions.L_inf], prob_I=params.pL[self.species], params_list=params_list, policy=policy, current_time=current_time, event_rates=event_rates)

        # HYNPOZOITE TRANSITIONS
        if person.state[self.species].hypnozoites.time_next == current_time: #hypnozoite death/activation scheduled for this time

            current_compartment_updated = person.state[self.species].current #update current compartment in case transition occured previously
            # Only changes compartment if in [S,A,R], no effect on dead zombies
            # assert current_compartment not in [Compartments.dead, Compartments.just_died], "Current compartment is "+str(current_compartment.name)
            activation = person.state[self.species].hypnozoites.update(person, self.species, current_time, p_clinical=params.pL[self.species]) #update number of hypnozoites and transition the human compartment
            
            if activation:
                new_compartment = person.state[self.species].current
                if new_compartment != current_compartment_updated:
                    #Check not moving out of dead compartments
                    assert current_compartment_updated not in [Compartments.dead, Compartments.just_died], "Current compartment is "+str(current_compartment.name)

                    self.pop_counts[current_compartment_updated] -= 1 # decrement population count for current state
                    self.pop_counts[new_compartment] += 1 # increment population count for next state
                    self.relapses.append(current_time)  # record relapse event
                    self.num_relapses[current_time] += 1
                    self.new_infections.append(current_time)  # just moved to I or A *now*
                    if new_compartment == Compartments.I: 
                        self.new_inf_clinical.append(current_time)  # just moved to I *now*

                    # continue as usual
                    self.transition_table(person=person, current_time=current_time, event_rates=event_rates, params=params, treatment_policy=policy)  # identify next transition to occur


        #Record relapse for new clinical infections if recently treated
        new_compartment = person.state[self.species].current
        if new_compartment == Compartments.I and current_compartment != Compartments.I:
            if (current_time - person.recent_treatment) < params.relapse_monitor:
                #Append relapse
                self.relapse_recorded[current_time] += 1
                # print("Relapse recorded, from compartment "+str(current_compartment.name)+" to "+(new_compartment.name)+" after "+str(current_time - person.recent_treatment)+ " days")

    def infect_me(self, person, rate_infection, prob_I, params, policy, current_time, event_rates):

        # params = params_list[policy]

        if random.random() < (1 - math.exp(-rate_infection)):
            #Update transition table for person-specific override value (e.g. different radical cure treatment durations)
            for transition in person.transition_overrides:
                event_rates[transition] = person.transition_overrides[transition]

            self.new_infections.append(current_time)  # just moved to I or A *now*
            self.pop_counts[person.state[self.species].current] -= 1  # decrement count for current state
            if random.random() <= prob_I:
                person.state[self.species].current = Compartments.I
                self.new_inf_clinical.append(current_time)  # just moved to I *now*
            else:
                person.state[self.species].current = Compartments.A

            if params.nu_hyp[self.species] > 0: #if avg number of hypnozoites injected from new primary infection, infect with hypnozoites
                assert self.species != Species.falciparum, "Falciparum shouldn't have hypnozoites"
                person.state[self.species].hypnozoites.infect(mean=params.nu_hyp[self.species], current_time=current_time) #inject hypnozoites upon new primary infection
            self.pop_counts[person.state[self.species].current] += 1  # increment new current state
            self.transition_table(person=person, current_time=current_time, event_rates=event_rates, params=params, treatment_policy=policy)  # determine next event and time it occurs

    def stochastic_sir_event_rates(self, params, mozzie, time):
        """
        Calculate the rates events happen for humans in this time step
        Moved to within the disease class so can have species specific parameter values and events
        :param params: model parameters
        :param params: Mozzie class holds the population level counts needed
        :return: event rates for infection and recovery
        """
        event_rate = [None] * params.number_events  # preallocate

        if type(params.b) == list:
            if self.species== Species.falciparum:
                lamx = params.b[self.species] * params.epsilonM[self.species] * (mozzie.mozzie_I_count[self.species] + (1-params.epsilonM[Species.vivax])*mozzie.mozzie_I_count[Species.mixed]) / params.human_population
            else:
                lamx = params.b[self.species] * params.epsilonM[self.species] * (mozzie.mozzie_I_count[self.species] + (1 - params.epsilonM[Species.falciparum]) * mozzie.mozzie_I_count[Species.mixed]) / params.human_population
        else:
            if self.species == Species.falciparum:
                lamx = params.b * params.epsilonM[self.species] * (mozzie.mozzie_I_count[self.species] + (1 - params.epsilonM[Species.vivax]) * mozzie.mozzie_I_count[Species.mixed]) / params.human_population
            else:
                lamx = params.b * params.epsilonM[self.species] * (mozzie.mozzie_I_count[self.species] + (1 - params.epsilonM[Species.falciparum]) *mozzie.mozzie_I_count[Species.mixed]) / params.human_population

        # 0: agent can be infected by mozzies with single or mixed infection
        event_rate[Transitions.S_inf] = lamx

        # 1: I -> A or death, natural recovery/death
        event_rate[Transitions.I_next] = params.sigma[self.species] #preceded by a check if they receive treatment

        # 2: I -> T or G
        event_rate[Transitions.I_treat] = params.tau[self.species] #preceded by a check if they receive treatment

        # 3: A -> R or L
        event_rate[Transitions.A_recover] = params.alpha[self.species]

        # 4: A -> T or G
        event_rate[Transitions.A_treat] = 0
        

        #With MDA. If used in future, adjust rates and transition table

        # if time >= params.mda_t1 and time < params.mda_t2:
        #     event_rate[Transitions.I_treat] = params.c[self.species] * params.tau[self.species] + params.eta[self.species] + params.etaFSAT[self.species] + params.etaMDA[self.species]

        #     event_rate[Transitions.A_treat] = params.eta[self.species] + params.etaFSAT[self.species] + params.etaMDA[self.species]
        # else:
        #     if time>params.mda_t2:
        #         #make sure MDA every 6 months
        #         params.mda_t1 = params.mda_t1 + (365.25/2)
        #         params.mda_t2 = params.mda_t2 + (365.25/2)
        #     event_rate[Transitions.I_treat] = params.c[self.species] * params.tau[self.species] + params.eta[self.species] + params.etaFSAT[self.species]
        #     event_rate[Transitions.A_treat] = params.eta[self.species] + params.etaFSAT[self.species]

        # 5: R -> I or A
        event_rate[Transitions.R_inf] = params.r[self.species] * lamx

        # 6: R -> S
        event_rate[Transitions.R_S] = params.omega[self.species]

        # # 7: L -> I or A
        # event_rate[Transitions.L_inf] = params.r[self.species] * lamx

        # # 8: L -> I or A
        # event_rate[Transitions.L_relapse] = params.nu[self.species]

        # # 9: L -> S
        # event_rate[Transitions.L_S] = params.kappa[self.species]

        # 10: T -> ...
        event_rate[Transitions.T_done] = params.rho[self.species]

        # 11: G -> ...
        # Amended to be person-specific: override if receiving treatment with different duration
        event_rate[Transitions.G_done] = params.psi[self.species]

        #event_rate[Transitions.mixed_inf] =
        return event_rate

