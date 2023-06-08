"""A class of Multiple Tier-based Strictly Local Grammars. Copyright (C) 2019
Alena Aksenova.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.
"""

from copy import deepcopy
from random import shuffle
from mtsl_class import *
from fsm_family import *

class MITSL(MTSL):
    """A class for Multiple Input-local Tier-based Strictly Local grammars and languages.

    Attributes:
        alphabet (list): alphabet used in the language;
        symbols (list): the set of possible m-length strings 
        grammar (list): the list of substructures;
        k (int): locality window across tier symbols;
        m (int): locality window across letters (character length of tier symbols)
        data (list): input data;
        edges (list): start- and end-symbols for the grammar;
        polar ("p" or "n"): polarity of the grammar;
        fsm (FSMFamily): a list of finite state machines that
            corresponds to the grammar;
        tier (list): list of tuples, where every tuple lists elements
            of some tier.
    Learning for k > 2 is not implemented: requires more theoretical work.
    """

    def __init__(
        self, alphabet=None, grammar=None, k=2, data=None, edges=[">", "<"], polar="p"
    ):
        """Initializes the TSL object."""
        super().__init__(alphabet, grammar, k, data, edges, polar)
        self.symbols = set()
        self.fsm = FSMFamily()
        self.m = 2 
        if self.k != 2:
            raise NotImplementedError(
                "The learner for k-MTSL languages is " "still being designed."
            )
        self.tier = None

    def annotate_string(self, string, asData = False):
        """Annotates the string with the start and end symbols.

        Arguments:
            string (str): a string that needs to be annotated.
            asData (boolean) (optional): whether these are being annotated for input data (True) or as a string to be scanned (False). 
        Returns:
            str: annotated version of the string.
        """
        if asData:
            return ">" * (self.m*self.k-1) + string.strip() + "<" * (self.m*self.k-1)#used for annotating data - must attest k*m-1 length sequences of edge characters followed by one alphabet character
        return ">" * self.m + string.strip() + "<" * self.m#used for annotating strings to scan
    
    def extract_alphabet(self):
        """Extracts alphabet from the given data or grammar and saves it into
        the 'alphabet' attribute.
        Also generates all symbols with extract_mgram_symbols()
        CAUTION: if not all symbols were used in the data or grammar,
                the result is not correct: update manually.
        """
        L.extract_alphabet(self)
        self.extract_mgram_symbols()

    def extract_mgram_symbols(self):
        """Generates all m-length symbols from the data and saves in into the 'symbols' attribute """
        self.symbols = list({"".join(gram) for gram in self.generate_all_ngrams(self.alphabet, self.m)}.union({edge * self.m for edge in self.edges}))

    def learn(self, restrictions_remove = [], symbols_remove = []):
        """
        Learns 2-local MITSL grammar for a given sample. The algorithm 
        currently works only for k=2 and is based on MITSL designed 
        by De Santo and Askenova (2021). This implementation 
        is slightly modified to add m-factors to an empty tier, rather
        than removing m-factors from a complete tier.
        Arguments:
            restrictions_remove (list of k-length tuples of m-length symbols representing restrictions, optional): the restrictions, from whose tiers, the symbols_remove symbols should be removed 
            symbols_remove (list of m-length str symbols): the symbols, which should be removed from the given tier
        Results:
            self.grammar is updated with a grammar of the following shape:
            {(tier_1):[bigrams_for_tier_1],
                ...
             (tier_n):[bigrams_for_tier_n]}
        """
        if not self.data:
            raise ValueError("Data needs to be provided.")
        if not self.alphabet:
            raise ValueError(
                "The alphabet is empty. Provide data or "
                "run `grammar.extract_alphabet`."
            )
        if not self.symbols:
            self.extract_mgram_symbols()

        self.data = list(set(self.data))#eliminate duplicates

        possible = set(self.generate_all_ngrams(self.symbols, self.k, addEdges= False))

        attested = set()
        for d in self.data:
            d = self.annotate_string(d, asData = True)
            grams = [d[i:i+self.k*self.m] for i in range(len(d)-self.k*self.m+1)]
            grams = [(gram[:self.m], gram[self.m:]) for gram in grams]
            attested.update(set(grams))

        unattested = list(possible.difference(attested))


        paths = self.all_paths(self.data)

        grammar = []

        b = list(self.symbols)
        for p1, p2 in unattested:
            c = {p1, p2}
            c.update([edge * self.m for edge in self.edges])
            for s in b:
                if s == p1 or s == p2:
                    continue
                relevant_paths = [path for path in paths if path[0] == p1 and s in path[1] and path[2] == p2]
                paths_without_s = [[p[0], p[1].difference({s}), p[2]] for p in relevant_paths]
                are_paths = [p in paths for p in paths_without_s]
                _ = [print('add:', p) for p in paths_without_s if p not in paths and ((p1,p2) in restrictions_remove and s in symbols_remove)]#This is a grammar engineering tool: it shows what paths would need to be added to remove particular Symbols from the tier containing particular Restrictions
                if not all(are_paths):
                    c.add(s)

            grammar.append((c, (p1, p2)))
        gathered = self.gather_grammars(grammar)
        self.grammar = gathered
        self.tier = [i for i in self.grammar]

        if self.check_polarity() == "p":
            self.grammar = self.opposite_polarity()

    def scan(self, string, verbose = False):
        """Scan string with respect to a given MTSL grammar.

        Arguments:
            string (str): a string that needs to be scanned.
            verbose(bool): True to optionally print details about which tier is violated
        Returns:
            bool: well-formedness of the string.
        """
        tier_evals = []
            
        string = self.annotate_string(string)

        bigrams = [(string[i:i+self.m], i) for i in range(len(string)-self.m + 1)]
        for tier in self.grammar:
            t = tier

            restrictions = self.grammar[tier]

            projection = [kfactor for kfactor in bigrams if kfactor[0] in tier]

            this_tier = [(projection[i][0], projection[i+1][0]) for i in range(len(projection) - 1) if (projection[i+1][1] - projection[i][1] > self.m-1)]
            this_tier = [pair in restrictions for pair in this_tier]

            valid = False
            if self.check_polarity() == "p":
                valid = all(this_tier)
            else:
                valid = not any(this_tier)
            if verbose and not valid:
                print(string, tier, projection, sorted(restrictions))
            tier_evals.append(valid)

        return all(tier_evals)

    def gather_grammars(self, grammar):
        """Gathers grammars with the same tier together.

        Arguments:
            grammar (list): a representation of the learned grammar
                where there is a one-to-one mapping between tiers
                and bigrams.
        Returns:
            dict: a dictionary where keys are tiers and values are
                the restrictions imposed on those tiers.
        """
        G = {}
        for i in grammar:
            if tuple(i[0]) in G:
                G[tuple(i[0])] += [i[1]]
            else:
                G[tuple(i[0])] = [i[1]]
        return G

    def path(self, string):
        """Collects a list of paths from a string.

        A path is a
        triplet <a, X, b>, where `a` is a symbol, `b` is a symbol
        that follows `a` in `string`, and `X` is a set of symbols
        in-between `a` and `b`.
        Arguments:
            string (str): a string paths of which need to be found.
        Returns:
            list: list of paths of `string`.
        """
        string = self.annotate_string(string)
        paths = set()

        data = [string[i:i+self.m] for i in range(len(string) - self.m + 1)]
        n = len(data)
        for start in range(0, n-2):
            for end in range(start+2, n):
                first = data[start]
                middle = tuple(sorted(set(data[start+1:end])))
                last = data[end]
                path = (first, middle, last)                
                paths.add(path)

        return paths

    def all_paths(self, dataset):
        """Finds all paths that are present in a list of strings.

        Arguments:
            dataset (list): a list of strings.
        Returns:
            list: a list of paths present in `dataset`.
        """
        paths = set()
        for item in dataset:
            paths.update(self.path(item))

        return [[first, set(middle), last] for first, middle, last in paths]

    def opposite_polarity(self):
        """Generates a grammar of the opposite polarity.

        Returns:
            dict: a dictionary containing the opposite ngram lists
                for every tier of the grammar.
        """
        if not self.grammar:
            raise ValueError(
                "Grammar needs to be provided. It can also "
                "be learned using `grammar.learn()`."
            )
        opposite = {}
        for i in self.grammar:
            possib = self.generate_all_ngrams(list(i), self.k, addEdges=False)#check as tuple
            opposite[i] = [j for j in possib if j not in self.grammar[i]]

        return opposite

    def switch_polarity(self):
        """Changes polarity of the grammar, and rewrites grammar to the
        opposite one."""
        self.grammar = self.opposite_polarity()
        self.change_polarity()

    def map_restrictions_to_fsms(self):
        """Builds FSM family corresponding to the given grammar"""
        restr_to_fsm = list()

        grammar = self.grammar if self.check_polarity() == "p" else self.opposite_polarity()
        for tier, ngrams in grammar.items():
            g = ngrams
            fsm = FSM(self.edges[0] * self.m, self.edges[1] * self.m)
            fsm.sl_to_fsm(g)
            restr_to_fsm.append([tier, g, fsm])
        return restr_to_fsm


    def fsmize(self):
        """Builds FSM family corresponding to the given grammar and saves in it
        the fsm attribute."""
        restr_to_fsm = self.map_restrictions_to_fsms()
        self.fsm.family = [i[2] for i in restr_to_fsm]

    def generate_sample(self, n=10, repeat=True, safe=True):
        """Generates a data sample of the required size, with or without
        repetitions depending on `repeat` value.

        Arguments:
            n (int): the number of examples to be generated;
            repeat (bool): allows (rep=True) or prohibits (rep=False)
               repetitions within the list of generated items;
            safe (bool): automatically breaks out of infinite loops,
                for example, when the grammar cannot generate the
                required number of data items, and the repetitions
                are set to False.
        Returns:
            list: generated data sample.
        """
        if not self.alphabet:
            raise ValueError("Alphabet cannot be empty.")
        if not self.fsm.family:
            self.fsmize()

        tier_smap = self.tier_state_maps()
        if not any([len(tier_smap[x]) for x in tier_smap]):
            raise (
                ValueError(
                    "There are ngrams in the grammar that are"
                    " not leading anywhere. Clean the grammar "
                    " or run `grammar.clean_grammar()`."
                )
            )

        main_smap = self.general_state_map(tier_smap)
        data = [self.generate_item(tier_smap, main_smap) for i in range(n)]

        if not repeat:
            data = set(data)
            useless_loops = 0
            prev_len = len(data)

            while len(data) < n:
                data.add(self.generate_item(tier_smap))

                if prev_len == len(data):
                    useless_loops += 1
                else:
                    useless_loops = 0

                if safe and useless_loops > 500:
                    print(
                        "The grammar cannot produce the requested "
                        "number of strings. Check the grammar, "
                        "reduce the number, or allow repetitions."
                    )
                    break

        return list(data)

    def tier_image(self, string):
        """
        Creates tier images of a string with respect to the different
        tiers listed in the grammar.
        Returns:
            dict: a dictionary of the following shape:
                { (tier_1):"string_image_given_tier_1",
                    ...,
                  (tier_n):"string_image_given_tier_n"
                }
        """

        tiers = {}
        for i in self.grammar:
            curr_tier = list()
            for index, s in enumerate(string):
                if s in self.edges or s in i:
                    curr_tier += [(s, index)]
            tiers[i] = curr_tier
        return tiers

    def merge_symbols(self, word):
        """Merges tuples of m-length symbols into one string by deleting overlapping letters between symbols
        Arguments:
            word (tuple of m-length strings): The word, represented as a sequence of m-length symbols
        Returns:
            (str): The word, with the overlapping parts of each symbol removed
        """
        return "".join((s[-1] for s in word[1:-self.m]))

    def generate_item(self, tier_smap, main_smap = None):
        """Generates a well-formed string with respect to the given grammar.
        Arguments:
            tier_smap (dict): The dictionary of transitions within the FSMs that correspond to the tier grammars
            main_smap (dict) (optional): The dictionary of transitions within all FSMs of the FSM family, will be computed if not provided
        Returns:
            str: a well-formed string.
        """
        word = (self.edges[0]*2,) * (self.k - 1) 
        main_smap = main_smap if main_smap is not None else self.general_state_map(tier_smap)#I am saving time by not recalculating this for each sample
        tier_images = self.tier_image(word)

        while word[-1] != self.edges[1] * self.m:
            symbols = list(main_smap[word[-(self.k - 1) :][0]])
            symbols = [symbol for symbol in symbols if symbol[:(self.m-1)] == word[-1][-(self.m-1):]]#This restriction ensures that adjacent symbols can be merged together at the end
            shuffle(symbols)#just generating a list once and then going through it in a shuffled order is more efficient than picking a different random one each time
            exhausted = True
            for maybe in symbols:
                good = True
                for tier in tier_smap:
                    if maybe in tier:
                        old_image = [oldSymbol[0] for oldSymbol in tier_images[tier] if oldSymbol[1] < len(word) - (self.m - 1)]#this ignores any previous symbols that would overlap with the next symbol
                        while len(old_image) < self.k - 1:
                            old_image = [self.edges[0]*2] + old_image
                        if maybe not in tier_smap[tier][tuple(old_image[-(self.k - 1) :])]:
                            good = False
                            break
                if good:
                    word += (maybe,)
                    tier_images = self.tier_image(word)
                    exhausted = False#we found a symbol that works here, so we did not exhaust our options
                    break
            if exhausted:
                print('We have generated ourselves into a corner with the word', str(word), 'please check your grammar.')
                return ""

        #print(word)

        newword = self.merge_symbols(word)

        if self.scan(newword):
            return newword
        else:
            return self.generate_item(tier_smap)

    def tier_state_maps(self):
        """
        Generates a dictionary of transitions within the FSMs
        that correspond to the tier grammars.
        Returns:
            dict: the dictionary of the form
                {
                 (tier_1):{"keys":[list of next symbols]},
                 (tier_2):{"keys":[list of next symbols]},
                   ...
                 (tier_n):{"keys":[list of next symbols]},
                }, where keys are (k-1)-long tier representations.
        Warning: the list of next symbols is tier-specific,
            so this estimates the rough options: refer to
            generate_item for the filtering of wrongly
            generated items.
        """
        restr_to_fsm = self.map_restrictions_to_fsms()
        tier_smaps = {}

        for curr_tier in restr_to_fsm:
            sl = SL()
            sl.change_polarity(self.check_polarity())
            sl.edges = [edge * self.m for edge in self.edges]
            sl.k = self.k
            sl.alphabet = list(curr_tier[0])
            sl.grammar = curr_tier[1]
            sl.fsm = curr_tier[2]
            tier_smaps[tuple(sl.alphabet)] = sl.state_map()

        return tier_smaps

    def general_state_map(self, smaps):
        """
        Generates a dictionary of transitions within all
        FSMs of the FSM family.
        Returns:
            dict: the dictionary of the form
                {"keys":[list of next symbols]}, where 
                keys are (k-1)-long strings.
        Warning: the list of next symbols is tier-specific,
            so this estimates the rough options: refer to
            generate_item for the filtering of wrongly
            generated items.
        """
        local_smaps = deepcopy(smaps)

        for tier in local_smaps:
            non_tier = [i for i in self.symbols if i not in tier]
            for entry in local_smaps[tier]:
                local_smaps[tier][entry].extend(non_tier)

        local_smaps = list(local_smaps.values())

        main_smap = deepcopy(local_smaps[0])

        for other in local_smaps[1:]:
            for entry in other:

                if entry not in main_smap:
                    main_smap[entry] = other[entry]
                else:
                    inter = [i for i in main_smap[entry] if i in other[entry]]
                    main_smap[entry] = inter

        free_ones = []
        for i in self.symbols:
            for j in self.grammar:
                if i in j:
                    break
            free_ones.append(i)

        ext_alphabet = deepcopy(self.symbols)
        for x in free_ones:
            main_smap[x] = ext_alphabet

        return main_smap

    def clean_grammar(self):
        """Removes useless ngrams from the grammar.

        If negative, it just removes duplicates. If positive, it detects
        ngrams to which one cannot get     from the initial symbol and
        from which one cannot get     to the final symbol, and removes
        them.
        """
        for tier in self.grammar:
            sl = SL()
            sl.change_polarity(self.check_polarity())
            sl.edges = [edge * self.m for edge in self.edges]
            sl.alphabet = self.symbols
            sl.k = self.k
            sl.grammar = self.grammar[tier]
            sl.fsmize()
            sl.clean_grammar()
            self.grammar[tier] = deepcopy(sl.grammar)





