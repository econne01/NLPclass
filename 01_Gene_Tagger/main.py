from tagger import Tagger

if __name__ == '__main__':
    print 'Begin Tagger main method'

    train_filename = r'C:\Users\Eric\projects\academic\NLPclass\01_Gene_Tagger\gene.train'
    count_filename = r'C:\Users\Eric\projects\academic\NLPclass\01_Gene_Tagger\gene.counts'
    dev_filename = r'C:\Users\Eric\projects\academic\NLPclass\01_Gene_Tagger\gene.dev'
    output_filename = r'C:\Users\Eric\projects\academic\NLPclass\01_Gene_Tagger\gene_dev.p1.out'

    tagger = Tagger()
    tagger.read_tag_count_file(count_filename)
    tagger.flag_rare_words()
    
    sentence = 'We found that RhoA can initiate a linear kinase cascade leading to the activation of ERK6 ( p38 gamma ), a recently identified number of the '
    sentence = ['BACKGROUND', ':', 'Ischemic', 'heart', 'disease', 'is', 'the', 'primary', 'cause', 'of', 'morbidity', 'and', 'mortality', 'among', 'diabetics']
    #tagger.get_sentence_tags(sentence)
    tagger.tag_file(dev_filename, output_filename)

    print 'Task complete'