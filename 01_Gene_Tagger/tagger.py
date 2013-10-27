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

    # define the specialty categories to group infrequent words
    category_keywords = ['_RARE_']
    rare_cnt_threshold = 5

    def get_rare_keyword(self, word):
        ''' Return the most appropriate RARE category keyword based on properties
            of the given word
        '''
        return self.category_keywords[0]

    def get_tag(self, word):
        ''' Calculate the most likely tag for the given input word and return tag
            @return tuple (tag, probability) example: ("I-GENE", 0.00045)
        '''
        if word not in self.trained_word_counts or self.trained_word_counts[word] < self.rare_cnt_threshold:
            # This word was not in training set, so assume it is a RARE keyword
            word = self.get_rare_keyword(word)
        max_tag = None
        max_prob = 0
        for tag in self.emission_counts:
            prob = self.get_emission_prob(word, tag)
            if prob > max_prob:
                max_prob = prob
                max_tag = tag
        return (max_tag, max_prob)

    def get_emission_prob(self, word, tag):
        ''' Return the probability of the given tag emitting the given word
            e(x|s) = count(s=>x) / count(s)
        '''
        if tag not in self.trained_tag_counts or word not in self.emission_counts[tag]:
            return 0
        else:
            emit_cnt = self.emission_counts[tag][word]
            tag_cnt = self.trained_tag_counts[tag]
            #tag_cnt = sum([self.trained_tag_counts[tag][w] for w in self.trained_tag_counts[tag]])
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
        # Read each word in input_file and write the tag to output
        for i,row in enumerate(ifile.readlines()):
            if i%500==0:
                print 'Processing row %s...' %i

            word = row.strip()
            if word=='':
                ofile.write('\n')
            else:
                tag, max_prob = self.get_tag(word)
                ofile.write(' '.join([word, tag]) + '\n')
        ifile.close()
        ofile.close()
        return
