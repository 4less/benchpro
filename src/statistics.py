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
        self.bray_curtis_intersection = 0
        self.bray_curtis_union = 0
        self.l2_intersection = 0
        self.l2_union = 0


def get_abundance_dict(profile, rank):
    abundance_dict = dict()
    for row in profile.Rows(rank):
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
        # Get abundance dicts for abundance statistics

        gold_dict = get_abundance_dict(gold_profile, rank)
        prediction_dict = get_abundance_dict(prediction_profile, rank)
        gold_vec = [gold_dict[key] for key in shared]
        pred_vec = [prediction_dict[key] for key in shared]

        logger.info("Shared taxa {}".format(len(gold_vec)))
        logger.info("Shared taxa: GoldStd Sum {}, Prediction Sum {}".format(sum(gold_vec), sum(pred_vec)))

        gold_vec_union = gold_vec + [value for key, value in gold_dict.items() if key not in prediction_dict]
        pred_vec_union = pred_vec + [0.0 for key in gold_dict if key not in prediction_dict]
        pred_vec_union += [value for key, value in prediction_dict.items() if key not in gold_dict]
        gold_vec_union += [0.0 for key in prediction_dict if key not in gold_dict]

        logger.info("Union taxa {}".format(len(gold_vec_union)))
        logger.info("Union taxa abundances {}".format(gold_vec_union))
        logger.info("Union taxa: GoldStd Sum {}, Prediction Sum {}".format(sum(gold_vec_union), sum(pred_vec_union)))
        if workbook_logger:
            workbook_logger.add_sample(gold_profile.name, gold_dict, prediction_dict)


        ##########################################################################################
        # Pearson correlation
        abundance.pearson_correlation_intersection = sc.stats.pearsonr(gold_vec, pred_vec).statistic if len(gold_vec) > 1 else 0
        abundance.pearson_correlation_union = sc.stats.pearsonr(gold_vec_union, pred_vec_union).statistic if len(gold_vec) > 1 else 0

        logger.info("Pearson Correlation Shared {} Union {}".format(abundance.pearson_correlation_intersection, abundance.pearson_correlation_union))

        ##########################################################################################
        # Bray-Curtis
        abundance.bray_curtis_intersection = bray_curtis(gold_vec, pred_vec)
        abundance.bray_curtis_union = bray_curtis(gold_vec_union, pred_vec_union)

        logger.info("Bray-Curtis Shared {} Union {}".format(abundance.bray_curtis_intersection, abundance.bray_curtis_union))

        ##########################################################################################
        # Euclidean (L2)
        abundance.l2_intersection = euclidian(gold_vec, pred_vec)
        abundance.l2_union = euclidian(gold_vec_union, pred_vec_union)

        logger.info("L2 Shared {} Union {}".format(euclidian(gold_vec, pred_vec), euclidian(pred_vec_union, gold_vec_union)))
        ##########################################################################################

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
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'BrayCurtisIntersect',
                                                  stats.bray_curtis_intersection))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'BrayCurtisUnion',
                                                  stats.bray_curtis_union))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'L2Intersect',
                                                  stats.l2_intersection))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'L2Union',
                                                  stats.l2_union))


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
                