import sys
from enum import Enum


class RankTemp:
    def __init__(self, name: str, index: int, gtdb_prefix: str) -> object:
        self.name = name
        self.index = index
        self.gtdb_prefix = gtdb_prefix


class Rank(Enum):
    DOMAIN = RankTemp('Domain', 0, 'd__')
    PHYLUM = RankTemp('Phylum', 1, 'd__')
    CLASS = RankTemp('Class', 2, 'c__')
    ORDER = RankTemp('Order', 3, 'o__')
    FAMILY = RankTemp('Family', 4, 'f__')
    GENUS = RankTemp('Genus', 5, 'a__')
    SPECIES = RankTemp('Species', 6, 's__')
    STRAIN = RankTemp('Strain', 7, 't__')


class Lineage:
    def __init__(self):
        self.list = [""] * len(Rank)

    def Get(self, rank: Rank):
        return self.list[rank.value.index]

    def Set(self, rank: Rank, value: str):
        self.list[rank.value.index] = value

    def ToString(self):
        return '\t'.join(self.list)


class ProfileEntry:
    def __init__(self, name: str, lineage: Lineage, abundance: float):
        self.name = name
        self.lineage = lineage
        self.abundance = abundance

    def ToString(self):
        return "{}\t{}\t{}".format(self.name, self.lineage.ToString(), self.abundance)

    def GetName(self):
        return self.name

    def GetLineage(self):
        return self.lineage

    def GetAbundance(self):
        return self.abundance


class Profile:
    def __init__(self):
        self.name = ''
        self.name_to_abundance = dict()

    def Add(self, entry: ProfileEntry):
        self.name_to_abundance[entry.name] = entry

    def Get(self, name: str) -> ProfileEntry:
        if name not in self.name_to_abundance:
            print("{} not present".format(name))
            return None
        return self.name_to_abundance[name]

    def Rows(self):
        return self.name_to_abundance.values()

    def ToString(self):
        outstr = "Name {}".format(self.name)
        for name,row in self.name_to_abundance.items():
            outstr += row.ToString()
            outstr += '\n'
        return outstr



class ProfileFactory:
    default_position_to_rank = {
        0: Rank.DOMAIN,
        1: Rank.PHYLUM,
        2: Rank.CLASS,
        3: Rank.ORDER,
        4: Rank.FAMILY,
        5: Rank.GENUS,
        6: Rank.SPECIES,
        7: Rank.STRAIN
    }

    def __init__(self, pos_to_rank=default_position_to_rank):
        self.position_to_rank = pos_to_rank
        self.separator_list = [';', '|']

    def DynamicSplit(self, lineage_str: str):
        split_list = []
        for sep in self.separator_list:
            splitlen = len(lineage_str.split(sep))

            if splitlen <= len(Rank):
                split_list.append((sep, splitlen))

        split_list.sort(reverse=True, key=lambda x: x[1])
        # print(split_list[0], split_list[0][0], lineage_str, lineage_str.split(split_list[0][0]))
        return lineage_str.split(split_list[0][0])

    def GetLineage(self, lineage_str: str):
        lineage = Lineage()
        tokens = lineage_str.split(';')
        tokens = self.DynamicSplit(lineage_str)
        for i in range(len(tokens)):
            lineage.Set(self.position_to_rank[i], tokens[i])
        return lineage

    def GetProfileEntry(self, name: str, lineage_str: str, abundance: float):
        lineage = self.GetLineage(lineage_str)
        return ProfileEntry(name, lineage, abundance)

    def __call__(self, file_path: str, name_column: int = 0, lineage_column: int = 1, abundance_column: int = 2,
                 rank_column: int = -1):
        profile = Profile()
        print("File: {}".format(file_path))
        with open(file_path, 'r') as file:
            line_num = 0
            for line in file:
                line = line.rstrip()
                tokens = line.split('\t')
                if len(tokens) < 3:
                    print("SHORTER: {}".format(tokens))

                if name_column >= len(tokens):
                    print("skip.. {} >= {} {}".format(name_column, len(tokens), line))
                    continue
                if lineage_column >= len(tokens):
                    print("skip..  {} >= {} {}".format(lineage_column, len(tokens), line))
                    continue
                if abundance_column >= len(tokens):
                    print("skip..  {} >= {} {}".format(abundance_column, len(tokens), line))
                    continue
                if rank_column > -1 and rank_column >= len(tokens):
                    print("skip..  {} >= {} {}".format(rank_column, len(tokens), line))
                    continue

                name = str(line_num) if name_column == -1 else tokens[name_column]
                lineage = tokens[lineage_column]
                abundance = 0

                try:
                    abundance = float(tokens[abundance_column])
                except ValueError:
                    print("Warning: {} is not convertible to float".format(abundance))
                    continue
                if abundance == 0:
                    continue

                if rank_column > -1:
                    rank = tokens[rank_column]
                    if str.lower(rank) != "species":
                        continue

                profile.Add(self.GetProfileEntry(
                    name, lineage, abundance)
                )

                pentry = profile.Get(name)
                # print("line {} {} {}".format(pentry.GetName(), pentry.GetLineage().Get(Rank.CLASS),
                #                              pentry.GetAbundance()))

                line_num += 1

        abundance_sum = sum(row.GetAbundance() for row in profile.Rows())

        if abundance_sum > 99 and abundance_sum < 101:
            for row in profile.Rows():
                row.abundance /= 100
            abundance_sum = sum(row.GetAbundance() for row in profile.Rows())
            print(abundance_sum)

        # if abundance_sum > 50:
        #     print(abundance_sum)
        #     exit(0)

        return profile
