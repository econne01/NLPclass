import os
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
    category_keywords = ['_RARE_']
    rare_cnt_threshold = 5

    def get_rare_keyword(self, word):
        ''' Return the most appropriate RARE category keyword based on properties
            of the given word
        '''
        return self.category_keywords[0]

    def get_word_or_keyword(self, word):
        ''' If the given word is rare, return the appropriate RARE keyword,
            else just return the given word
        '''
        if word not in self.trained_word_counts or self.trained_word_counts[word] < self.rare_cnt_threshold:
            # This word was not in training set, so assume it is a RARE keyword
            word = self.get_rare_keyword(word)
        return word

    def get_sentence_tags(self, sentence):
        ''' Run viterbi algorithm to get arg max tags for the given
            (space-separated) sentence
        '''
        max_prob = 0
        max_tags = ['Default']
        n = len(sentence)
        # clear pi_cache
        self.pi_cache = {}
        # Find max ending tag-pair using (recursive) Viterbi algorithm
        for u in self.trained_tag_counts:
            word = self.get_word_or_keyword(sentence[-2])
            # Check for reasons to skip whole Viterbi algorithm
            # Case 1, u tag never emits required word
            if self.get_emission_prob(word, u) == 0:
                continue
            for v in self.trained_tag_counts:
                # Case 2, no tag-sequence ends in u,v,STOP
                if self.get_trigram_prob('STOP', u, v) == 0:
                    continue
                # Case 3
                word = self.get_word_or_keyword(sentence[-1])
                if self.get_emission_prob(word, v) == 0:
                    continue

                tags, prob = self.pi(n, u, v, sentence)
                if prob * self.get_trigram_prob('STOP', u, v) > max_prob:
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
        word = self.get_word_or_keyword(sentence[k-1])

        if k>0:
            # Check for reasons to halt recursive algorithm
            if (self.get_emission_prob(word, v) == 0 #Case 1, u tag never emits required word
                    or ' '.join([u,v]) not in self.ngrams == 0 #Case 2, no tag-sequence of (u,v) in training set
                    or self.get_emission_prob(word, v) == 0 #Case 3, tag v never emits word
                    ):
                return (['Break case'], 0)

        # Check for cached pi(k,u,v) value
        if (k,u,v) in self.pi_cache:
            return self.pi_cache[(k,u,v)]
        # Check for base case
        elif k==0:
            if u=='*' and v=='*':
                max_prob = 1
                max_tags = ['*', '*']
            else:
                max_prob = 0
                max_tags = ['*', '*']
        # If k is 2nd word or less, then we need to set w='*'
        elif k in [1,2]:
            w = '*'
            max_tags, max_prob = self.pi(k-1, w, u, sentence)
            if max_prob != 0:
                max_prob *= self.get_trigram_prob(v, w, u)
                max_prob *= self.get_emission_prob(word, v)
                max_tags = max_tags+[v]
        # Return most likely tag for k-2nd word, iterate over all tags
        else:
            max_prob = 0
            max_tags = ['Initial Pi']
            for w in self.trained_tag_counts:
                tags, prob = self.pi(k-1, w, u, sentence)
                if prob != 0:
                    prob *= self.get_trigram_prob(v, w, u)
                    prob *= self.get_emission_prob(word, v)
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
            uvs_cnt = 0
        uv = ' '.join([prior2, prior1])
        if uv in self.ngrams:
            uv_cnt = self.ngrams[uv]
        else:
            uv_cnt = 0
            # Cannot divide by zero, so just return 0 here instead of error-ing out
            return 0
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
        sentence = []
        sent_cnt = 1
        for i,row in enumerate(ifile.readlines()):
            if i%500==0:
                print 'Processing row %s...' %i

            word = row.strip()
            if word=='':
                print 'Tagging sentence', sent_cnt, '(# words= ', len(sentence), ')'
                tags, prob = self.get_sentence_tags(sentence)
                print '    prob=', prob
                for i in range(len(sentence)):
                    ofile.write(' '.join([sentence[i], tags[i]]) + '\n')
                ofile.write('\n')
                sentence = []
                sent_cnt += 1
            else:
                sentence.append(word)
        ifile.close()
        ofile.close()
        return
