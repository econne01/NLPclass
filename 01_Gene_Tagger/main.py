from tagger import Tagger

if __name__ == '__main__':
    print 'Begin Tagger main method'

    train_filename = r'C:\Users\Eric\projects\academic\NLP\01_Tagging_Chemistry\gene.train'
    count_filename = r'C:\Users\Eric\projects\academic\NLP\01_Tagging_Chemistry\gene.counts'
    dev_filename = r'C:\Users\Eric\projects\academic\NLP\01_Tagging_Chemistry\gene.test'
    output_filename = r'C:\Users\Eric\projects\academic\NLP\01_Tagging_Chemistry\gene_test.p1.out'

    tagger = Tagger()
    tagger.read_tag_count_file(count_filename)
    tagger.flag_rare_words()
    
    tagger.tag_file(dev_filename, output_filename)

    print 'Task complete'