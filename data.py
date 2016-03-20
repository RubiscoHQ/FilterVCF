import os
import config


class Vcf():
    def __init__(self, route):
        if os.path.exists(route):
            self.head = []
            self.variants = []

            f = open(route)
            for line in f:
                if line[:2] == '##':
                    self.head.append(line)
                elif line[0] == '#':
                    self.title = line
                    title_list = line[1:].strip().split('\t')
                    self.sample_list = title_list[9:]
                else:
                    info_dict = {}
                    info_list = line.strip().split('\t')
                    for i in range(9):
                        info_dict[title_list[i]] = info_list[i]
                    info_dict['samples'] = {}
                    for sample in self.sample_list:
                        info_dict['samples'][sample] = info_list[title_list.index(sample)]
                    self.variants.append(Variant(info_dict, _format='vcf'))
        else:
            print 'Error: file do not exist', route
            raise IOError
        return

    def return_vcf(self):
        return


def check_annovar_title():
    return

class Annovar():
    def __init__(self, route, sample_columns, info_columns, vcf_columns, id_columns):
        if os.path.exists(route):
            self.variants = []

            f = open(route)
            title = True
            for line in f:
                if title:
                    self.title = line
                    title_list = line.strip().split('\t')
                    self.sample_list = []
                    for i in sample_columns:
                        self.sample_list.append(title_list[i])
                    title = False
                else:
                    info_dict = {}
                    info_list = line.strip().split('\t')

                    info_dict['infos'] = {}
                    for i in info_columns:
                        info_dict['infos'][title_list[i]] = info_list[i]

                    info_dict['samples'] = {}
                    for i in sample_columns:
                        info_dict['samples'][title_list[i]] = info_list[i]

                    info_dict['vcfs'] = {}
                    for i in vcf_columns:
                        info_dict['vcfs'][title_list[i]] = info_list[i]

                    info_dict['ids'] = {}
                    for i in range(5):
                        info_dict['ids'][title_list[i]] = info_list[i]

                    self.variants.append(Variant(info_dict, _format='wannovar'))
        else:
            print 'Error: file do not exist', route
            raise IOError
        return

    def return_annovar(self):
        return

    def return_vcf(self):
        return


class Variant():
    def __init__(self, info_dict, _format='vcf'):
        if _format == 'vcf':
            self.infos = {}
            self.info_format = 'vcf'

            self.infos['chrom'] = info_dict['CHROM']
            self.infos['vcfpos'] = info_dict['POS']
            self.infos['id'] = info_dict['ID']
            self.infos['vcfref'] = info_dict['REF']
            self.infos['vcfalt'] = info_dict['ALT']
            self.infos['qual'] = info_dict['QUAL']
            self.infos['filter'] = info_dict['FILTER']

            info_list = info_dict['INFO'].split(';')
            for i in info_list:
                lt = i.split('=')
                try:
                    self.infos[lt[0]] = lt[1]
                except IndexError:
                    self.infos[lt[0]] = None

            self.vcfformat = info_dict['FORMAT']
            format_list = self.vcfformat.split(':')

            self.samples = {}
            for sample in info_dict['samples']:
                sample_info = info_dict['samples'][sample].split(':')
                self.samples[sample] = {}
                for i in range(len(format_list)):
                    self.samples[sample][format_list[i]] = sample_info[i]

        elif _format == 'wannovar':
            self.info_format = 'annovar'
            self.infos = info_dict['infos']

            self.infos['chrom'] = info_dict['ids']['Chr']
            self.infos['pos_start'] = info_dict['ids']['Start']
            self.infos['pos_end'] = info_dict['ids']['End']
            self.infos['vcfref'] = info_dict['ids']['Ref']
            self.infos['vcfalt'] = info_dict['ids']['Alt']
            self.infos['vcfpos'] = info_dict['vcfs']['POS']
            self.infos['id'] = info_dict['vcfs']['ID']
            self.infos['qual'] = info_dict['vcfs']['QUAL']
            self.infos['filter'] = info_dict['vcfs']['FILTER']

            info_list = info_dict['vcfs']['INFO'].split(';')
            for i in info_list:
                lt = i.split('=')
                try:
                    self.infos[lt[0]] = lt[1]
                except IndexError:
                    self.infos[lt[0]] = None

            self.vcfformat = info_dict['vcfs']['FORMAT']
            format_list = self.vcfformat.split(':')

            self.samples = {}
            for sample in info_dict['samples']:
                sample_info = info_dict['samples'][sample].split(':')
                self.samples[sample] = {}
                for i in range(len(format_list)):
                    self.samples[sample][format_list[i]] = sample_info[i]
        else:
            print "Error: Unknow variant type:", _format
            raise
        return


class Sample():
    def __init__(self, health_type, infos):
        self.type = health_type

    pass


class Gene():
    pass














