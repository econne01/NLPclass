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

    tagger.tag_file(dev_filename, output_filename)

    print 'Task complete'
