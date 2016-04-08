# wannovar file title, update on 2016-03-20
default_annovar_info_names = ['1000g2014oct_all', '1000g2014oct_eas', '1000g2014oct_sas', '1000g_afr',
                              '1000g_all', '1000g_amr', '1000g_eas', '1000g_eur', '1000g_sas', 'cadd',
                              'esp6500siv2_aa', 'esp6500siv2_all', 'esp6500siv2_ea', 'exac_all', 'popfreqmax',
                              'aachange.refgene', 'cadd_phred', 'cadd_raw', 'cg46', 'clinvar_20150330',
                              'clinvar_db', 'clinvar_dbid',
                              'clinvar_dis', 'clinvar_id', 'clinvar_sig', 'clinvar_status', 'cosmic_dis',
                              'cosmic_id', 'cytoband', 'dbsnp', 'esp6500si_aa', 'esp6500si_all', 'esp6500si_ea',
                              'esp6500siv2_all', 'exac_afr', 'exac_amr', 'exac_eas', 'exac_fin', 'exac_freq',
                              'exac_nfe', 'exac_oth', 'exac_sas', 'exonicfunc.refgene', 'fathmm_pred',
                              'fathmm_score', 'func.refgene', 'gene.refgene', 'genedetail.refgene', 'gerp++_rs',
                              'gwas_beta', 'gwas_dis', 'gwas_or', 'gwas_p', 'gwas_pubmed', 'gwas_snp',
                              'lr_pred', 'lr_score', 'lrt_pred', 'lrt_score', 'mutationassessor_pred',
                              'mutationassessor_score', 'mutationtaster_pred', 'mutationtaster_score',
                              'nci60', 'phylop100way_vertebrate', 'phylop46way_placental', 'polyphen2_hdiv_pred',
                              'polyphen2_hdiv_score', 'polyphen2_hvar_pred', 'polyphen2_hvar_score',
                              'radialsvm_pred', 'radialsvm_score', 'sift_pred', 'sift_score', 'siphy_29way_logodds',
                              'snp138', 'vest3_score']



default_annovar_ids = ['chr', 'start', 'end', 'ref', 'alt']

default_vcf_title = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format']

default_skip = ['otherinfo']


def get_sample_id_from_sample_file(route):
    f = open(route)
    title = True
    idl = []
    for line in f:
        if title:
            title = False
            continue
        id = line.split('\t')[1]
        idl.append(id.lower())
    return idl



def check_annovar_title(route, sample_file):
    f = open(route)
    title = f.readline()
    title_list = title.strip().split('\t')
    dic = {}
    if title_list[:5] == default_annovar_ids \
       and title_list[5:69] == default_annovar_info_names \
       and title_list[72:81] == default_vcf_title:
        dic = {'ids': [i for i in range(5)],
               'infos': [i for i in range(5, 69)],
               'vcfs': [i for i in range(72, 81)],
               'samples': [i for i in range(81, len(title_list))]}
    else:
        sample_ids = get_sample_id_from_sample_file(sample_file)
        dic = {'ids': [], 'infos': [], 'vcfs': [], 'samples': []}
        unknown = []
        for i in title_list:
            ilower = i.lower()
            if ilower == 'ref':
                if title_list[title_list.index(i)-1].lower() == 'end':
                    dic['ids'].append(title_list.index(i))
                elif title_list[title_list.index(i) - 1].lower() == 'id':
                    dic['vcfs'].append(title_list.index(i))
                else:
                    print "Error, the %d column name of your title %s is invalid. Doesn't match vcf nor annovar."\
                          % (title_list.index(i), 'ref')
                    raise
            elif ilower == 'alt':
                if title_list[title_list.index(i) - 2].lower() == 'end':
                    dic['ids'].append(title_list.index(i))
                elif title_list[title_list.index(i) - 2].lower() == 'id':
                    dic['vcfs'].append(title_list.index(i))
                else:
                    print "Error, the %d column name of your title %s is invalid. Doesn't match vcf nor annovar." \
                          % (title_list.index(i), 'ref')
                    raise

            elif ilower in default_annovar_ids:
                dic['ids'].append(title_list.index(i))
            elif ilower in default_annovar_info_names:
                dic['infos'].append(title_list.index(i))
            elif ilower in default_vcf_title:
                dic['vcfs'].append(title_list.index(i))
            elif ilower in default_skip:
                continue
            elif ilower in sample_ids:
                dic['samples'].append(title_list.index(i))
            else:
                unknown.append(i)
        if unknown != []:
            print 'Warming: %d columns are unknown in input file:' % len(unknown)
            print unknown
            user_input = raw_input('Will you ignore these unknown columns? (Y/N)')
            if user_input.upper() != 'Y':
                print 'Exit without filtering.\nPlease revise the column name in your input file or add new names into config.py'
                raise
    return dic


genotype_empty_symble = ['0/0', './.']


gene_name_column_in_wannovar = 'gene.refgene'


info_na_char = ['.', '']

