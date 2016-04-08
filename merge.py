import linecache
import os
import filter
import data


def get_anno(path):
    f = open(path)
    d = {}
    cache_data = linecache.getlines(path)
    title = True
    for line in range(len(cache_data)):
        if title:
            title = False
            t = cache_data[line]
            continue
        ll = cache_data[line].strip().split('\t')
        key = tuple(ll[-8:-3])
        d[key] = ll
    linecache.clearcache()
    return d, t


def count_alle_freq(sample_info, sample_id, cases, ctrls):
    allele_count = {}
    case_count = 0.0
    ctrl_count = 0.0

    for sample in range(len(sample_info)):
        gt = filter.search_words('./.', sample_info[sample])
        sample_id_ = sample_id[sample]
        if sample_id_ in cases:
            type = 'case'
            case_count += 1
        elif sample_id_ in ctrls:
            type = 'ctrl'
            ctrl_count += 1
        else:
            print 'Error: Unknown sample', sample_id_
            raise
        if gt:
            gtl = gt.split('/')
            for i in gtl:
                if i not in ['0', '.']:
                    try:
                        allele_count[i][type] += 1
                    except KeyError:
                        allele_count[i] = {'case': 0, 'ctrl': 0}
                        allele_count[i][type] += 1

    allele_count_l = []
    allele_freq_l = []
    for geno in allele_count:
        allele_count_l.append('case:'+str(allele_count[geno]['case'])+';'+'ctrl:'+str(allele_count[geno]['ctrl']))

        allele_freq_l.append('%.3f' % (allele_count[geno]['case']/case_count/2) +
                             ':'+'%.3f' % (allele_count[geno]['ctrl']/ctrl_count/2))
    allele_count_s = ','.join(allele_count_l)
    allele_freq_s = ','.join(allele_freq_l)
    sample_count_s = 'case:'+str(int(case_count))+';'+'ctrl:'+str(int(ctrl_count))
    if allele_count_s == '':
        allele_count_s = '.'
    if allele_freq_s == '':
        allele_freq_s = '.'
    return [allele_count_s, sample_count_s, allele_freq_s]


def add_genotype_from_vcf_to_annovar(vcf_f, anno_f, new_f, sample_file):
    anno_d, anno_title = get_anno(anno_f)
    nf = open(new_f, 'w')
    print 'get anno info'
    sg = data.SampleGroup(sample_file)
    cases = sg.cases_id
    ctrls = sg.ctrls_id
    if os.path.exists(vcf_f):
        cache_data = linecache.getlines(vcf_f)
        print 'get vcf info'
        for line in range(1, len(cache_data)):
            if cache_data[line][:2] == '##':
                continue
            elif cache_data[line][0] == '#':
                vcf_title = cache_data[line].strip().split('\t')[8:] + ['allele_count', 'sample_count', 'allele_freq']
                ntl = anno_title.strip().split('\t') + vcf_title
                nt = '\t'.join(ntl) + '\n'
                nf.write(nt)
                continue
            ll = cache_data[line].strip().split('\t')
            key = tuple(ll[:5])
            if key in anno_d:
                nl = anno_d[key]+ll[8:]+count_alle_freq(sample_info=ll[9:], sample_id=vcf_title[1:-3], cases=cases, ctrls=ctrls)
                nf.write('\t'.join(nl)+'\n')
        return
    else:
        print('the path [{}] is not exist!'.format(vcf_f))


add_genotype_from_vcf_to_annovar('/Users/rubisco/data/deaf_seq/vcf/mergeall.vcf',
                                 '/Users/rubisco/PycharmProjects/FilterVCF/work/filtered.annovar.txt',
                                 '/Users/rubisco/PycharmProjects/FilterVCF/work/filtered.gt.txt',
                                 '/Users/rubisco/PycharmProjects/FilterVCF/work/SampleInfo.txt')
