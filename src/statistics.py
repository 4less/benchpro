from src.profile import Profile, default_ranks
from src.binary_statistics import BinaryStatistics
from src.profile import RankTemp
from src.profile import Rank
import scipy as sc
import sys
import math
from loguru import logger


class AbundanceStatistics:
    def __init__(self):
        self.pearson_correlation_intersection = 0
        self.pearson_correlation_union = 0
        self.bray_curtis = 0
        self.l2 = 0


def get_abundance_dict(profile, rank):
    abundance_dict = dict()
    for row in profile.Rows():
        name = row.GetLineage().Get(rank)
        if name == "":
            continue
        if name not in abundance_dict:
            abundance_dict[name] = 0.0
        abundance_dict[name] += row.GetAbundance()

    return abundance_dict


def bray_curtis(a, b):
    total_a = sum(a)
    total_b = sum(b)
    shared_a = sum(va for va, vb in zip(a, b) if va > 0 and vb > 0)
    shared_b = sum(vb for va, vb in zip(a, b) if va > 0 and vb > 0)
    #
    # print("Bray_Curtis")
    # print("gold: {}".format(a))
    # print("pred: {}".format(b))
    # print("total_gold: {}, total_pred: {}\tshared_gold: {}, pred_gold: {}".format(total_a, total_b, shared_a, shared_b))
    # print("min@ {}".format(min(shared_a, shared_b)))
    # print("(total_a + total_b): {}".format((total_a + total_b)))
    bc = (2*min(shared_a, shared_b) / (total_a + total_b))

    # print("Bray: {}".format(bc))
    return bc

def euclidian(a,b):
    sum = 0
    for a1, a2 in zip(a, b):
        sum += (a1 - a2) ** 2
    return math.sqrt(sum)

def get_rank_statistics(gold_profile: Profile, prediction_profile: Profile,
                        ranks=default_ranks, workbook_logger=None):

    logger.info(f"{get_rank_statistics} for {gold_profile.name} and {prediction_profile.name}")
    rank_dict = dict()
    abundance_stats = dict()

    for rank in ranks:
        logger.info(f"Process rank: {rank}")
        stats = BinaryStatistics()
        abundance = AbundanceStatistics()


        prediction_set = set(row.GetLineage().Get(rank) for row in prediction_profile.Rows(rank))
        gold_set = set(row.GetLineage().Get(rank) for row in gold_profile.Rows(rank))

        shared = prediction_set.intersection(gold_set)
        pred_only = prediction_set.difference(gold_set)
        gold_only = gold_set.difference(prediction_set)
        pred_no_shared = prediction_set.difference(shared)
        gold_no_shared = gold_set.difference(shared)

        if len(shared) == 0:
            logger.warning("No shared taxa between gold profile and prediction profile")
            logger.warning("prediction_set: {}".format(prediction_set))
            logger.warning("gold_set: {}" .format(gold_set))

        gold_set = set(row.GetLineage().Get(rank) for row in gold_profile.Rows())
        stats.tp = len(shared)
        stats.fp = len(pred_no_shared)
        stats.fn = len(gold_no_shared)



        ##########################################################################################
        # Pearson correlation on shared taxa

        gold_dict = get_abundance_dict(gold_profile, rank)
        prediction_dict = get_abundance_dict(prediction_profile, rank)
        gold_vec = [gold_dict[key] for key in shared]
        pred_vec = [prediction_dict[key] for key in shared]

        if workbook_logger:
            workbook_logger.add_sample(gold_profile.name, gold_dict, prediction_dict)

        if len(shared) > 1:
            # print(sc.stats.pearsonr(gold_vec, pred_vec))
            abundance.pearson_correlation_intersection = sc.stats.pearsonr(gold_vec, pred_vec).statistic
        else:
            abundance.pearson_correlation_intersection = 0

        ##########################################################################################
        # Pearson correlation on union taxa

        # gold_vec += [value for key,value in gold_dict]
        # gold_vec += [0.0 for _ in prediction_dict]
        # pred_vec += [0.0 for _ in gold_dict]
        # pred_vec += [value for key, value in prediction_dict.items()]

        if len(gold_vec) > 1:
            # print(sc.stats.pearsonr(gold_vec, pred_vec))
            abundance.pearson_correlation_union = sc.stats.pearsonr(gold_vec, pred_vec).statistic
        else:
            abundance.pearson_correlation_union = 0


        ##########################################################################################
        # Bray-Curtis
        abundance.bray_curtis = bray_curtis(gold_vec, pred_vec)
        ##########################################################################################
        # Euclidean (L2)
        abundance.l2 = euclidian(gold_vec, pred_vec)
        ##########################################################################################

        # print("False negatives: {}".format(gold_set.difference(prediction_set)))
        # print("False positives: {}".format(prediction_set.difference(gold_set)))

        abundance_stats[rank] = abundance
        rank_dict[rank] = stats

        try:
            stats.Validate()
        except:
            logger.exception(f"{gold_profile.name}\t{prediction_profile.name}\t{rank}")
            exit(9)

    return rank_dict, abundance_stats


def print_binary_stats(stats_dict, name: str = '', dataset: str = '', sep: str = '\t', file_handle=sys.stdout):
    prefix = '{}{}'.format(name, sep) if name else ''
    dataset_prefix = '{}{}'.format(dataset, sep) if name else ''

    for rank, stats in stats_dict.items():
        rank = '{}{}'.format(rank.value.name, sep) if rank else ''

        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'TP', stats.tp))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'FP', stats.fp))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'FN', stats.fn))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'Sensitivity', stats.Sensitivity()))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'Precision', stats.Precision()))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'F1', stats.F1()))


def print_abundance_stats(stats_dict, name: str = '', dataset: str = '', sep: str = '\t', file_handle=sys.stdout):
    prefix = '{}{}'.format(name, sep) if name else ''
    dataset_prefix = '{}{}'.format(dataset, sep) if name else ''

    for rank, stats in stats_dict.items():
        rank = '{}{}'.format(rank.value.name, sep) if rank else ''

        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'PearsonCorrelationIntersect',
                                                  stats.pearson_correlation_intersection))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'PearsonCorrelationUnion',
                                                  stats.pearson_correlation_union))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'BrayCurtis',
                                                  stats.bray_curtis))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'L2',
                                                  stats.l2))


def write_statistics(output: str, all_statistics_dict, all_abundance_statistics_dict):
    with open(output, 'w') as out:
        out.write('{}\t{}\t{}\t{}\t{}\n'.format(
            "tool", "sample", "rank", "metric", "value"
        ))
        for tool, statistics_dict_list in all_statistics_dict.items():
            abundance_dict_list = all_abundance_statistics_dict[tool]
            dataset_num = 0
            for statistics_dict, abundance_dict in zip(statistics_dict_list, abundance_dict_list):
                dataset = str(dataset_num)
                print_binary_stats(statistics_dict, tool, dataset, file_handle=out)
                print_abundance_stats(abundance_dict, tool, dataset, file_handle=out)
                dataset_num += 1
                