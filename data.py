import os
import config


class Vcf():
    def __init__(self, route):
        if os.path.exists(route):
            self.head = []
            self.variants = []
            self.route = route

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
        f = open(os.path.splitext(self.route)[0]+'.return.vcf', 'w')
        for line in self.head:
            f.write(line)

        title_list = self.title.strip().split('\t')[:9]+sorted(self.sample_list)
        f.write('\t'.join(title_list)+'\n')
        for var in self.variants:
            f.write(var.return_vcf_line())
        return


class Annovar():
    def __init__(self, route, sample_columns, info_columns, vcf_columns, id_columns):
        if os.path.exists(route):
            self.variants = []
            self.route = route
            self.id_columns = id_columns
            self.info_columns = info_columns
            self.vcf_columns = vcf_columns
            self.sample_columns = sample_columns

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
        f = open(os.path.splitext(self.route)[0]+'.return.annovar.txt', 'w')

        title_list = self.title.strip().split('\t')
        new_title_list = [title_list[i] for i in self.id_columns] + \
                         [title_list[i] for i in self.info_columns] + \
                         [title_list[i] for i in self.vcf_columns] + \
                         [title_list[i] for i in self.sample_columns]
        f.write('\t'.join(new_title_list)+'\n')

        for var in self.variants:
            f.write(var.return_annovar_line(title_list, self.info_columns, self.vcf_columns, self.sample_columns))
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
                    self.infos[lt[0]] = ''

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
            self.infos['annovarref'] = info_dict['ids']['Ref']
            self.infos['annovaralt'] = info_dict['ids']['Alt']
            self.infos['vcfref'] = info_dict['vcfs']['REF']
            self.infos['vcfalt'] = info_dict['vcfs']['ALT']
            self.infos['vcfpos'] = info_dict['vcfs']['POS']
            self.infos['id'] = info_dict['vcfs']['ID']
            self.infos['qual'] = info_dict['vcfs']['QUAL']
            self.infos['filter'] = info_dict['vcfs']['FILTER']

            info_list = info_dict['vcfs']['INFO'].split(';')
            self.vcf_infos_names = []
            # save vcf info column item names in order to recover when do return_annoavr_line
            for i in info_list:
                lt = i.split('=')
                self.vcf_infos_names.append(lt[0])
                try:
                    self.infos[lt[0]] = lt[1]
                except IndexError:
                    self.infos[lt[0]] = ''

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

    def return_vcf_line(self):
        infos_list = []
        for item in self.infos:
            infos_list.append(item+'='+self.infos[item])
        infos_str = ';'.join(infos_list)

        line_list = [self.infos['chrom'], self.infos['vcfpos'], self.infos['id'],
                     self.infos['vcfref'], self.infos['vcfalt'], self.infos['qual'],
                     self.infos['filter'], infos_str, self.vcfformat]

        sample_list= sorted(self.samples.iteritems(), key=lambda d: d[0])
        for sample in sample_list:
            sample_info = []
            for item in self.vcfformat.split(':'):
                sample_info.append(sample[1][item])
            line_list.append(':'.join(sample_info))

        line = '\t'.join(line_list) + '\n'
        return line

    def return_annovar_line(self, title, info_columns, vcf_columns, sample_columns):
        line_list = [self.infos['chrom'], self.infos['pos_start'], self.infos['pos_end'],
                     self.infos['annovarref'], self.infos['annovaralt']]
        # id columns

        for i in info_columns:
            line_list.append(self.infos[title[i]])
        # info columns

        line_list += [self.infos['chrom'], self.infos['vcfpos'], self.infos['id'],
                      self.infos['vcfref'], self.infos['vcfalt'], self.infos['qual'],
                      self.infos['filter']]
        vcf_infos = []
        for i in self.vcf_infos_names:
            vcf_infos.append(i+'='+self.infos[i])
        line_list.append(';'.join(vcf_infos))
        line_list.append(self.vcfformat)
        # vcf columns

        for i in sample_columns:
            sample = title[i]
            sample_info = []
            for item in self.vcfformat.split(':'):
                sample_info.append(self.samples[sample][item])
            line_list.append(':'.join(sample_info))
        # sample_columns

        line = '\t'.join(line_list) + '\n'
        return line


class Sample():
    def __init__(self, health_type, infos):
        self.type = health_type

    pass


class Gene():
    pass













