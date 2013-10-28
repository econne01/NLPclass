import os
from tagger import Tagger

if __name__ == '__main__':
    print 'Begin Tagger main method'

    train_filename = os.getcwd() + r'/gene.train'
    count_filename = os.getcwd() + r'/gene.counts'
    dev_filename = os.getcwd() + r'/gene.dev'
    output_filename = os.getcwd() + r'/gene_dev.p2.out'

    tagger = Tagger()
    tagger.read_tag_count_file(count_filename)
    tagger.flag_rare_words()

    sentence = 'We found that RhoA can initiate a linear kinase cascade leading to the activation of ERK6 ( p38 gamma ), a recently identified number of the p38 family of MAPKs .'
    #tagger.get_sentence_tags(sentence.split(' '))
    tagger.tag_file(dev_filename, output_filename)

    print 'Task complete'
