import os, re
import count_freqs

class Tagger(object):
    # emissions is dictionary with structure:
    #    emissions[tag] = {word1: count}
    #    ie, emissions[I-GENE] = {'hydrolase': 2, 'opsin': 1}
    trained_tag_counts = {}
    trained_word_counts = {}
    emission_counts = {}
    ngrams = {}

    # pi_cache used to improve performance of Viterbi algorithm
    pi_cache = {}

    # define the specialty categories to group infrequent words
    category_keywords = ['_NUMERIC_', '_ALLCAPS_', '_LASTCAP_', '_RARE_']
    rare_cnt_threshold = 5

    def get_rare_keyword(self, word):
        ''' Return the most appropriate RARE category keyword based on properties
            of the given word
        '''
        if re.search('[0-9]', word):
            # 0. Includes a number
            return self.category_keywords[0]
        elif re.search('^[A-Z]+$', word):
            # 1. All Capital Letters
            return self.category_keywords[1]
        elif re.search('^.*[A-Z]$', word):
            # 2. Ends in Capital Letter
            return self.category_keywords[2]
        else:
            # 3. Remaining Catch-all Rare word
            return self.category_keywords[3]

    def get_word_or_keyword(self, word):
        ''' If the given word is rare, return the appropriate RARE keyword,
            else just return the given word
        '''
        if word not in self.trained_word_counts or self.trained_word_counts[word] < self.rare_cnt_threshold:
            # This word was not in training set, so assume it is a RARE keyword
            word = self.get_rare_keyword(word)
        return word

    def get_possible_tags(self, k, N):
        ''' List the possible tags at location `k` in sentence of length `N`
        '''
        if k < 0:
            return ['*']
        elif k >= N:
            return ['STOP']
        else:
            return self.trained_tag_counts.keys()

    def get_sentence_tags(self, sentence):
        ''' Run viterbi algorithm to get arg max tags for the given
            (space-separated) sentence
        '''
        max_prob = 0
        max_tags = []
        N = len(sentence)
        # clear pi_cache
        self.pi_cache = {}
        self.bp_cache = [None]*N
        # For each possible tag at word location N-1
        for v in self.get_possible_tags(N-1, N):
            tags, prob = self.pi(N, v, 'STOP', sentence)
            if prob >= max_prob:
                max_prob = prob
                max_tags = tags
        # Remove initial start tags from max_tags list
        max_tags = max_tags[2:]
        return (max_tags, max_prob)
        

    def pi(self, k, u, v, sentence):
        ''' helper function for Viterbit algorithm
            This method is recursive. It traverses a sentence in reverse
            to calculate the most likely sequence of tags, based on
            likely tag-sequences (from HMM counts) and likely emissions
        '''
        # If checking the STOP point of sentence, no emission involved STOP always yields STOP
        if k==len(sentence):
            emit_prob=1
        # Check for reasons to halt recursive algorithm
        elif k>=0:
            # Halt if tag v never emits word
            word = self.get_word_or_keyword(sentence[k])
            emit_prob = self.get_emission_prob(word, v)
            if emit_prob == 0:
                return ([], 0)

        # Check for cached pi(k,u,v) value
        if (k,u,v) in self.pi_cache:
            return self.pi_cache[(k,u,v)]
        # Check for base case
        elif k==-1:
            if u=='*' and v=='*':
                max_prob = 1
            else:
                max_prob = 0
            max_tags = ['*', '*']
        # Return most likely tag for k-2nd word, iterate over all tags
        else:
            max_prob = 0
            max_tags = []
            for w in self.get_possible_tags(k-2, len(sentence)):
                tags, prob = self.pi(k-1, w, u, sentence)
                if prob != 0:
                    prob *= self.get_trigram_prob(v, w, u)
                    prob *= emit_prob
                    if prob > max_prob:
                        max_prob = prob
                        max_tags = tags+[v]
        if (k,u,v) not in self.pi_cache:
            self.pi_cache[(k,u,v)] = (max_tags, max_prob)
        return (max_tags, max_prob)

    def get_word_tag(self, word):
        ''' Calculate the most likely tag for the given input word and return tag
            * calculation based solely on emission; arg(y) max p(y=>word)
            @return tuple (tag, probability) example: ("I-GENE", 0.00045)
        '''
        word = self.get_word_or_keyword(word)
        max_tag = None
        max_prob = 0
        for tag in self.emission_counts:
            prob = self.get_emission_prob(word, tag)
            if prob > max_prob:
                max_prob = prob
                max_tag = tag
        return (max_tag, max_prob)

    def get_trigram_prob(self, tag, prior2, prior1):
        ''' Return the probability of the given tag following the given 2 prior tags
            q(tag|prior2, prior1) = count(prior2 prior1 tag) / count(prior2 prior1)
        '''
        uvs = ' '.join([prior2, prior1, tag])
        if uvs in self.ngrams:
            uvs_cnt = self.ngrams[uvs]
        else:
            return 0
        uv = ' '.join([prior2, prior1])
        uv_cnt = self.ngrams[uv]
        return float(uvs_cnt)/uv_cnt


    def get_emission_prob(self, word, tag):
        ''' Return the probability of the given tag emitting the given word
            e(x|s) = count(s=>x) / count(s)
        '''
        if tag not in self.trained_tag_counts or word not in self.emission_counts[tag]:
            return 0
        else:
            emit_cnt = self.emission_counts[tag][word]
            tag_cnt = self.trained_tag_counts[tag]
            return float(emit_cnt)/tag_cnt

    def flag_rare_words(self):
        ''' For any word that appears infrequently in training set
            (with ANY tag), replace it with a RARE keyword in emissions data.
        '''
        for word in self.trained_word_counts.copy():
            count = self.trained_word_counts[word]
            # If word is RARE, update RARE word and emissions counts
            if count < self.rare_cnt_threshold:
                rare_keyword = self.get_rare_keyword(word)
                # Update the RARE keyword is in word_counts
                if rare_keyword not in self.trained_word_counts:
                    self.trained_word_counts[rare_keyword] = 0
                self.trained_word_counts[rare_keyword] += count
                
                # Update the RARE keyword is in emissions
                for tag in self.emission_counts.copy():
                    if word in self.emission_counts[tag]:
                        emit_cnt = self.emission_counts[tag][word]
                        if rare_keyword not in self.emission_counts[tag]:
                            self.emission_counts[tag][rare_keyword] = 0
                        self.emission_counts[tag][rare_keyword] += emit_cnt

    def read_tag_count_file(self, filename):
        ''' Read the given tag counts file and store results locally to Tagger instance
        '''
        try:
            file = open(filename, 'r')
        except:
            raise Exception('Cannot open file: %s' % filename)

        for line in file.readlines():
            row = line.split(' ')
            if row[1] == 'WORDTAG':
                self.process_wordtag(row)
            elif 'GRAM' in row[1]:
                self.process_ngram(row)

    def process_wordtag(self, row):
        ''' Store the counts data for a WORDTAG row (ie, "4 WORDTAG I-GENE obsin")
        '''
        count = int(row[0].strip())
        tag = row[2].strip()
        word = row[3].strip()
        if tag not in self.trained_tag_counts:
            self.trained_tag_counts[tag] = 0
            self.emission_counts[tag] = {}
        if word not in self.emission_counts[tag]:
            self.emission_counts[tag][word] = 0
        self.trained_tag_counts[tag] += count
        self.emission_counts[tag][word] += count
        # Make sure word is included in tracking word_counts
        if word not in self.trained_word_counts:
            self.trained_word_counts[word] = 0
        self.trained_word_counts[word] += count

    def process_ngram(self, row):
        ''' Store the counts data for an N-GRAM row (ie, "15 3-GRAM I-GENE I-GENE O")
        '''
        count = int(row[0].strip())
        sequence = ' '.join(row[2:]).strip()
        if sequence not in self.ngrams:
            self.ngrams[sequence] = 0
        self.ngrams[sequence] += count

    def get_sentences(self, file_lines):
        ''' Parse the formatted file to break it down into sentences
            ie, array of words
        '''
        sentences = []
        sentence = []
        for row in file_lines:
            word = row.strip()
            if word=='':
                sentences.append(sentence)
                sentence = []
            else:
                sentence.append(word)
        return sentences

    def tag_file(self, input_filename, output_filename):
        ''' For each word, in each sentence in input_filename, find the 
            most likely tag and output results to output_filename.
            @param string input_filename. (File format ["This", "Gene", "myosin"])
            @param string output_filename. (File format ["This O", "Gene O", "myosin I-GENE"])
        '''
        # Open input file for reading
        try:
            ifile = open(input_filename, 'r')
        except:
            raise Exception('Cannot open file: %s' % input_filename)
        # Open output file for reading
        try:
            ofile = open(output_filename, 'w')
        except:
            raise Exception('Cannot open file: %s' % output_filename)

        # Read each sentence in input_file and write with proper tags to output
        sentences = self.get_sentences(ifile.readlines())
        for i,s in enumerate(sentences):
            print 'Tagging sentence', str(i), '(# words= ', len(s), ')'
            tags, prob = self.get_sentence_tags(s)
            for i in range(len(s)):
                ofile.write(' '.join([s[i], tags[i]]) + '\n')
            ofile.write('\n')
        ifile.close()
        ofile.close()
        return
