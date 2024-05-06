import os.path
import sys
from enum import Enum
from loguru import logger
from ete3 import NCBITaxa


class RankTemp:
    def __init__(self, name: str, index: int, gtdb_prefix: str) -> object:
        self.name = name
        self.index = index
        self.gtdb_prefix = gtdb_prefix


class Rank(Enum):
    __order__ = 'DOMAIN PHYLUM CLASS ORDER FAMILY GENUS SPECIES STRAIN'
    DOMAIN = RankTemp('Domain', 0, 'd__')
    PHYLUM = RankTemp('Phylum', 1, 'p__')
    CLASS = RankTemp('Class', 2, 'c__')
    ORDER = RankTemp('Order', 3, 'o__')
    FAMILY = RankTemp('Family', 4, 'f__')
    GENUS = RankTemp('Genus', 5, 'g__')
    SPECIES = RankTemp('Species', 6, 's__')
    STRAIN = RankTemp('Strain', 7, 't__')


default_ranks = [Rank.DOMAIN, Rank.PHYLUM, Rank.CLASS, Rank.ORDER, Rank.FAMILY, Rank.GENUS,
                 Rank.SPECIES]

rank_map = dict([(b.value.name, b) for a, b in enumerate(Rank)])
rank_map["strain"] = Rank.STRAIN
rank_map["species"] = Rank.SPECIES
rank_map["genus"] = Rank.GENUS
rank_map["family"] = Rank.FAMILY
rank_map["order"] = Rank.ORDER
rank_map["class"] = Rank.CLASS
rank_map["phylum"] = Rank.PHYLUM
rank_map["superkingdom"] = Rank.DOMAIN
rank_map["domain"] = Rank.DOMAIN

idx2rank = [
    Rank.DOMAIN,
    Rank.PHYLUM,
    Rank.CLASS,
    Rank.ORDER,
    Rank.FAMILY,
    Rank.GENUS,
    Rank.SPECIES,
    Rank.STRAIN
]

class Lineage:
    def __init__(self):
        self.list = [None] * len(Rank)

    def Get(self, rank: Rank):
        return self.list[rank.value.index]

    def Set(self, rank: Rank, value: str):
        self.list[rank.value.index] = value

    def GetRank(self):
        none_idx = len(self.list) if None not in self.list else self.list.index(None)
        return idx2rank[none_idx - 1]

    def GetTaxon(self):
        return self.list[self.GetRank().value.index]

    def IsTaxonKnown(self):
        return self.GetTaxon() != "" and not self.GetTaxon().endswith("__")

    def ToString(self, sep='\t'):
        lineage_list = [self.GetRank().name] + [t for t in self.list if t != None]
        return sep.join(lineage_list)

    def of_rank(self, rank) -> bool:
        # print("Rank {}".format(rank))
        # print("1: {}".format( self.list[rank.value.index]))
        # print("2: {}".format(len(self.list)))
        # print("3: {}".format(self.list[rank.value.index + 1]))
        # print(rank.value.index + 1 == len(self.list) or not self.list[rank.value.index + 1])
        # print(len(self.list))
        # print([taxon == '' for taxon in self.list[rank.value.index + 1:]])
        # print(all([taxon == '' for taxon in self.list[rank.value.index + 1:]]))
        # if self.list[rank.value.index] != "" and all([taxon == '' for taxon in self.list[rank.value.index + 1:]]):
        #     print(self.ToString())
        #     input()
        return self.list[rank.value.index] != None and (len(self.list) == rank.value.index + 1 or self.list[rank.value.index + 1] == None)
        # return self.list[rank.value.index] != "" and all([taxon == '' for taxon in self.list[rank.value.index + 1:]])

    def lowest_rank(self):
        if self.list[-1] != "":
            return Rank.STRAIN
        # index = self.list.index("") - 1
        index = 0
        for i,v in enumerate(reversed(self.list)):
            if v != "":
                index = len(self.list) - i - 1
                break
        return default_ranks[index]
    
    def IsGTDB(self) -> bool:
        return any(self.list[rank.value.index].startswith(rank.value.gtdb_prefix) for rank in Rank if self.list[rank.value.index] != None)
    
    def HasRank(self, rank: Rank) -> bool:
        return self.list[rank.value.index] != "" and (self.list[rank.value.index].startswith(rank.value.gtdb_prefix) or not self.IsGTDB())
    
    def HasUnidentifiedRank(self, rank: Rank) -> bool:
        return self.list[rank.value.index] != "" and (not self.list[rank.value.index].startswith(rank.value.gtdb_prefix) and self.IsGTDB())

    def HasUnidentifiedTaxon(self, rank: Rank) -> bool:
        return self.list[rank.value.index] == rank.value.gtdb_prefix and self.IsGTDB()
    

class ProfileEntry:
    def __init__(self, name: str, lineage: Lineage, abundance: float):
        self.name = name
        self.lineage = lineage
        self.abundance = abundance

    def IsUnclassified(self):
        return self.lineage[0] == "Unclassified"

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
    
    def HasRank(self, rank: Rank) -> bool:
        return self.lineage.HasRank(rank)


class Profile:
    def __init__(self, name):
        self.name = name
        self.name_to_abundance = dict()


    def IsEmpty(self):
        return len(self.name_to_abundance) == 0

    def Add(self, entry: ProfileEntry):
        self.name_to_abundance[entry.name] = entry

    def Get(self, name: str) -> ProfileEntry:
        if name not in self.name_to_abundance:
            logger.error("{} not present".format(name))
            return None
        return self.name_to_abundance[name]

    def Rows(self, rank=None, exclude_unknown=False):
        return [row for row in self.name_to_abundance.values() if (row.GetLineage().GetRank() == rank or rank == None) and (not exclude_unknown or row.GetLineage().IsTaxonKnown())]

    def ValidRows(self, rank=None):
        return self.name_to_abundance.values() if not rank else {
                row for row in self.name_to_abundance.values() if row.GetLineage().GetRank() == rank or rank == None
            }

    def infer_missing_abundances(self, rank):
        rank_rows = self.Rows(rank)
        all_rows = self.Rows()
        target_ranks = [r for r in default_ranks if r.value.index < rank.value.index]
    
        rank_sum = sum(v.abundance for v in rank_rows)
        all_sum = sum(v.abundance for v in all_rows)
        # for row in self.Rows():
        #     row.abundance /= rank_sum

        print("Target ranks: {}".format(target_ranks))
        print("{} vs {}".format(rank_sum, all_sum))
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

            new_sum = sum(v.abundance for v in self.Rows(trank))

            print("Rank {} new sum {} vs {}".format(trank, new_sum, rank_sum))
            # input()




    def count_ranks(self):
        rd = dict()
        for row in self.Rows():
            lowest_rank = row.GetLineage().GetRank()
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

        dbfile = "/usr/users/QIB_fr017/fritsche/Projects/benchpro/taxa.sqlite"
        dumpfile = "/usr/users/QIB_fr017/fritsche/Projects/benchpro/taxdump.tar.gz"
        dumpfile = "/usr/users/QIB_fr017/fritsche/Projects/benchpro/taxdump_cami2_toy.tar.gz"
        if os.path.exists(dbfile):
            self.ncbi = NCBITaxa(dbfile=dbfile)
        else:
            if dumpfile is None:
                raise IOError(
                    "The db is empty, you must specify the taxdump.tar.gz file to create the database.")
            self.ncbi = NCBITaxa(dbfile=dbfile, taxdump_file=dumpfile)

    def DynamicSplit(self, lineage_str: str):
        split_list = []
        for sep in self.separator_list:
            splitlen = len(lineage_str.split(sep))

            if splitlen <= len(Rank):
                split_list.append((sep, splitlen))

        split_list.sort(reverse=True, key=lambda x: x[1])
        # print(split_list[0], split_list[0][0], lineage_str, lineage_str.split(split_list[0][0]))
        return lineage_str.split(split_list[0][0])

    def GetLineage(self, lineage_str: str, rank=None):
        lineage = Lineage()
        tokens = self.DynamicSplit(lineage_str)
        n = len(tokens) if rank == None else rank.value.index + 1

        if n > len(tokens):
            ncbi_ranks = self.ncbi.get_rank([t for t in tokens if t.isnumeric()])
            new_tokens = n * ['']
            for t in reversed(tokens):
                if t.isnumeric() and int(t) in ncbi_ranks and ncbi_ranks[int(t)] in rank_map:
                    new_rank = rank_map[ncbi_ranks[int(t)]]
                    new_tokens[new_rank.value.index] = t

            tokens = new_tokens

        for i in range(n):
            lineage.Set(self.position_to_rank[i], tokens[i])
        return lineage

    def GetProfileEntry(self, name: str, lineage_str: str, abundance: float, rank=None):
        lineage = self.GetLineage(lineage_str, rank)

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

                if line.startswith("@") or line.startswith("#"):
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
                    logger.exception("Warning: {} is not convertible to float".format(tokens[abundance_column]))
                    exit
                if abundance == 0:
                    continue

                if rank_column > -1:
                    rank = tokens[rank_column]
                    
                    if rank in rank_map:
                        rank = rank_map[rank]
                    else:
                        rank = None

                else:
                    rank = None

                profile.Add(self.GetProfileEntry(
                    name, lineage, abundance, rank)
                )


                line_num += 1

        rank_counts = profile.count_ranks()



        if len(rank_counts) == 1:
            # Present rank to infer abundances from
            rank = [rank for rank in rank_counts][0]
            print("Infer missing abundances: {}\n{}".format(rank, rank_counts))
            profile.infer_missing_abundances(rank)

        normalize_abundances = False
        for rank in default_ranks:
            abundance_sum = sum(row.GetAbundance() for row in profile.Rows(rank))
            print("Abundance sum {} {} -> {} > 1.5? {}".format(rank, profile.name, abundance_sum, abundance_sum > 1.5))

            if abundance_sum > 1.5:
                normalize_abundances = True

        if normalize_abundances:
            for row in profile.Rows():
                row.abundance /= 100


        return profile
