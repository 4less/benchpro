import os.path
import sys
from enum import Enum
from loguru import logger


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


default_ranks = [Rank.DOMAIN, Rank.PHYLUM, Rank.CLASS, Rank.ORDER, Rank.FAMILY, Rank.GENUS,
                 Rank.SPECIES]

rank_map = dict([(b.value.name, b) for a, b in enumerate(Rank)])


class Lineage:
    def __init__(self):
        self.list = [""] * len(Rank)

    def Get(self, rank: Rank):
        if self.list[rank.value.index] == "":
            logger.warning(f"{self.list} is empty at rank index {rank.value.index} for rank {rank}")
        return self.list[rank.value.index]

    def Set(self, rank: Rank, value: str):
        self.list[rank.value.index] = value

    def ToString(self, sep='\t'):
        return sep.join(self.list)

    def of_rank(self, rank):
        return self.list[rank.value.index] and (rank.value.index + 1 == len(self.list) or not self.list[rank.value.index + 1])

    def lowest_rank(self):
        index = self.list.index("") - 1
        return default_ranks[index]


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

    def AddAbundance(self, num):
        self.abundance += num

    def of_rank(self, rank):
        return self.lineage.of_rank(rank)

    def lowest_rank(self):
        return self.lineage.lowest_rank()


class Profile:
    def __init__(self, name):
        self.name = name
        self.name_to_abundance = dict()

    def Add(self, entry: ProfileEntry):
        self.name_to_abundance[entry.name] = entry

    def Get(self, name: str) -> ProfileEntry:
        if name not in self.name_to_abundance:
            logger.error("{} not present".format(name))
            return None
        return self.name_to_abundance[name]

    def Rows(self, rank=None):
        # print("Rows {}".format(rank))
        # print("length: {}".format(len(self.name_to_abundance.values())))
        # if rank:
        #     for v in self.name_to_abundance.values():
        #         print("Value: {}".format(v))
        #         print("{} of rank {}? {}".format(v.ToString(), rank, v.of_rank(rank)))
        return self.name_to_abundance.values() if not rank else {value for value in self.name_to_abundance.values() if value.of_rank(rank)}

    def infer_missing_abundances(self, rank):
        rank_rows = self.Rows(rank)
        target_ranks = [r for r in default_ranks if r.value.index < rank.value.index]

        for trank in target_ranks:
            for row in rank_rows:
                name = row.GetLineage().Get(trank)
                if name not in self.name_to_abundance:
                    lineage = Lineage()
                    for i in range(0, trank.value.index + 1):
                        lineage.Set(default_ranks[i], row.GetLineage().Get(default_ranks[i]))
                    self.name_to_abundance[name] = ProfileEntry(name, lineage, 0)
                self.name_to_abundance[name].AddAbundance(row.GetAbundance())
                logger.info(f"{name}: {row.GetAbundance()}")




    def count_ranks(self):
        rd = dict()
        for row in self.Rows():
            lowest_rank = row.lowest_rank()
            if lowest_rank not in rd:
                rd[lowest_rank] = 0
            rd[lowest_rank] += 1
        return rd

    def ToString(self):
        outstr = "Name {}".format(self.name)
        for name, row in self.name_to_abundance.items():
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
        tokens = self.DynamicSplit(lineage_str)
        for i in range(len(tokens)):
            lineage.Set(self.position_to_rank[i], tokens[i])
        return lineage

    def GetProfileEntry(self, name: str, lineage_str: str, abundance: float):
        lineage = self.GetLineage(lineage_str)

        return ProfileEntry(name, lineage, abundance)

    def __call__(self, file_path: str, profile_name: str, name_column: int = 0, lineage_column: int = 1, abundance_column: int = 2,
                 rank_column: int = -1):
        profile = Profile(profile_name)

        # print(f"Profile Name: {file_path}")
        # input()
        print(file_path)
        print("File Exists : {}".format(os.path.exists(file_path)))
        with open(file_path, 'r') as file:
            line_num = 0
            for line in file:
                if len(line) == 0:
                    logger.warning("line is empty {}".format(line))
                    continue

                line = line.rstrip()
                tokens = line.split('\t')
                if len(tokens) < 3:
                    logger.warning("SHORTER: {}".format(tokens))

                if name_column >= len(tokens):
                    logger.warning("skip.. {} >= {} {}".format(name_column, len(tokens), line))
                    continue
                if lineage_column >= len(tokens):
                    logger.warning("skip..  {} >= {} {}".format(lineage_column, len(tokens), line))
                    continue
                if abundance_column >= len(tokens):
                    logger.warning("skip..  {} >= {} {}".format(abundance_column, len(tokens), line))
                    continue
                if rank_column > -1 and rank_column >= len(tokens):
                    logger.warning("skip..  {} >= {} {}".format(rank_column, len(tokens), line))
                    continue

                name = str(line_num) if name_column == -1 else tokens[name_column]
                lineage = tokens[lineage_column]
                abundance = 0

                try:
                    abundance = float(tokens[abundance_column])
                except ValueError:
                    logger.warning("Warning: {} is not convertible to float".format(tokens[abundance_column]))
                    continue
                if abundance == 0:
                    continue

                # if rank_column > -1:
                #     rank = tokens[rank_column]
                #     if str.lower(rank) != "species":
                #         logger.warning("skip no species {} because {}".format(line, rank_column))
                #         continue

                profile.Add(self.GetProfileEntry(
                    name, lineage, abundance)
                )

                line_num += 1

        rank_counts = profile.count_ranks()

        if len(rank_counts) == 1:
            rank = [rank for rank in rank_counts][0]
            profile.infer_missing_abundances(rank)


        for rank in default_ranks:
            abundance_sum = sum(row.GetAbundance() for row in profile.Rows(rank))
            print("Abundance sum {} -> {}".format(profile.name, abundance_sum))
            if abundance_sum > 50:
                for row in profile.Rows(rank):
                    row.abundance /= 100



        return profile
