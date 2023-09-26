import argparse
from src.profile import ProfileFactory, Rank, rank_map, default_ranks
from src.meta import Meta, write_template_meta
from loguru import logger
import sys

class Options:

    def __init__(self):
        self.gold_std_files = []
        self.predictions = dict()
        self.prediction_columns = dict()
        self.gold_columns = dict()
        self.names = []
        self.ranks = []
        self.metafile_path = None
        self.metafile = None

        self.output = ''

    def LoadPrediction(self, positional_argument: str):
        for name, tool in zip(self.names, positional_argument):
            self.predictions[name] = tool.rstrip(',').split(',')

    def ParseGoldColumns(self, gold_columns_str: str):
        tokens=gold_columns_str.split('|')
        self.gold_columns['name'] = int(tokens[0])
        self.gold_columns['lineage'] = int(tokens[1])
        self.gold_columns['abundance'] = int(tokens[2])
        if len(tokens) > 3:
            self.gold_columns['rank'] = int(tokens[3])
        else:
            self.gold_columns['rank'] = -1

    def ParsePredictionColumns(self, prediction_columns_str: str):
        tool_strs = prediction_columns_str.split(',')

        if len(tool_strs) == 1:
            tool_strs = tool_strs * len(self.predictions.keys())
        if len(tool_strs) != len(self.predictions.keys()):
            print("asdsadasdasd")
            exit(9)
        for tool, tool_str in zip(self.predictions.keys(), tool_strs):
            self.prediction_columns[tool] = dict()
            tokens = tool_str.split('|')
            self.prediction_columns[tool]['name'] = -1 if tokens[0] == 'X' else int(tokens[0])
            self.prediction_columns[tool]['lineage'] = int(tokens[1])
            self.prediction_columns[tool]['abundance'] = int(tokens[2])
            if len(tokens) > 3:
                self.prediction_columns[tool]['rank'] = int(tokens[3])
            else:
                self.prediction_columns[tool]['rank'] = -1

def parse_columns(columns_str: str):
    dic = dict()
    tokens = columns_str.split('|')
    dic['name'] = -1 if tokens[0] == 'X' else int(tokens[0])
    dic['lineage'] = int(tokens[1])
    dic['abundance'] = int(tokens[2])
    if len(tokens) > 3:
        dic['rank'] = int(tokens[3])
    else:
        dic['rank'] = -1

    return dic

def get_args():
    if '--meta_template' in sys.argv:
        parser = argparse.ArgumentParser()
        parser.add_argument('--meta_template', dest='meta_template', type=str, default=None,
                            help='Provide path to write template file')
    if '--meta' in sys.argv:
        parser = argparse.ArgumentParser()
        parser.add_argument('--meta', dest='meta_file', type=str, default=None,
                            help='Meta file containing all the information on where to find files etc to run the benchmarks. For more information please run with --meta_help.')
        parser.add_argument('--output', dest='output', type=str,
                            help='Output statistics file.')
        parser.add_argument('--detailed_output', dest='detailed_output', type=str, default=None,
                            help='Detailed output statistics file.')
        parser.add_argument('--ranks', dest='ranks', type=str, default=None,
                            help='target ranks')

    else:
        parser = argparse.ArgumentParser(description='Benchmark the performance of metagenomic profilers. Works also with '
                                                     'GTDB.')
        parser.add_argument('--gold', dest='gold_standard_files', type=str,
                            help='File path to gold standard file')
        parser.add_argument('--names', dest='names', type=str,
                            help='names of individual profiles')
        parser.add_argument('--ranks', dest='ranks', type=str, default=None,
                            help='target ranks')
        parser.add_argument('--output', dest='output', type=str,
                            help='Output statistics file.')
        parser.add_argument('--detailed_output', dest='detailed_output', type=str, default=None,
                            help='Detailed output statistics file.')
        parser.add_argument('--meta', dest='meta_file', type=str, default=None,
                            help='Meta file containing all the information on where to find files etc to run the benchmarks. For more information please run with --meta_help.')
        parser.add_argument('--meta_template', dest='meta_template', type=str, default=None,
                            help='Provide path to write template file')
        parser.add_argument('--meta_help', action='store_true',
                            help='Get help for meta_file')
        parser.add_argument('--gold_cols', dest='gold_columns', type=str, default='0|1|2',
                            help='Give column format for  gold standard profile. name_col|lineage_col|abundance_col(|rank_col :optional) (default: 0|1|2)')
        parser.add_argument('--pred_cols', dest='prediction_columns', type=str, default='0|1|2',
                            help='Give column format for  prediction profiles. name_col|lineage_col|abundance_col(|rank_col :optional) (default: 0|1|2). Separate by tool by deliminating with ,. (example: 0|1|2,1|2|3)')
        parser.add_argument('tool_predictions', nargs='+')
    args = parser.parse_args()
    return args


def get_options():
    args = get_args()

    if hasattr(args, "meta_template"):
        write_template_meta(args.meta_template)
        logger.info(f"Written template of meta file to {args.meta_template}.")
        exit(2)

    options = Options()

    print(args.meta_file)

    if hasattr(args, "meta_file"):
        options.meta = Meta(args.meta_file)
        options.output = args.output
        options.detailed_output = args.detailed_output

        if not args.ranks:
            options.ranks = default_ranks
        else:
            options.ranks = [rank_map[t] for t in args.ranks.split(',')]
        return options

    options.gold_std_files = args.gold_standard_files.rstrip(',').split(',')
    options.names = args.names.rstrip(',').split(',')
    options.output = args.output
    options.detailed_output = args.detailed_output
    options.ranks = args.ranks
    predictions = args.tool_predictions

    options.meta_file = args.meta
    options.meta_help = args.meta_help
    options.meta = None


    if not options.ranks:
        options.ranks = default_ranks
    else:
        options.ranks = [rank_map[t] for t in options.ranks.split(',')]

    if len(options.names) != len(predictions):
        print("length  of names {} does not match number of tools provided {}".format(len(options.names),
                                                                                      len(predictions)))
        exit(2)

    options.LoadPrediction(args.tool_predictions)
    options.ParseGoldColumns(args.gold_columns)
    options.ParsePredictionColumns(args.prediction_columns)


    for tool, files in options.predictions.items():
        if len(options.gold_std_files) != len(files):
            print(
                "length  of files {} for tool {} does not match with number of gold std profiles {}"
                .format(len(files), tool, len(options.gold_std_files)))


    return options
