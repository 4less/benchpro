import sys
import argparse

parser = argparse.ArgumentParser(description='Benchmark the performance of metagenomic profilers. Works also with '
                                             'GTDB.')
parser.add_argument('-1', '--first', dest='first_file', type=str,
                    help='File path to gold standard file')
parser.add_argument('-2', '--second', dest='second_file', type=str,
                    help='File path to gold standard file')
parser.add_argument('-d', '--delim', dest='delim', type=str, default="\t",
                    help='Global delimiter')
parser.add_argument('-d1', '--delim1', dest='delim_first', type=str,
                    help='Delimiter for first')
parser.add_argument('-d2', '--delim2', dest='delim_second', type=str,
                    help='Delimiter for second')
parser.add_argument('-j', '--join', dest='join', type=int, default=0,
                    help='Delimiter for first')
parser.add_argument('-j1', '--join1', dest='join_first', type=int,
                    help='Delimiter for first')
parser.add_argument('-j2', '--join2', dest='join_second', type=int,
                    help='Delimiter for second')
parser.add_argument('-c1', '--columns1', dest='columns_first', type=str, default=None,
                    help='Columns for first')
parser.add_argument('-c2', '--columns2', dest='columns_second', type=str, default=None,
                    help='Columns for second')

args = parser.parse_args()




first_file = args.first_file
second_file = args.second_file

delimiter_first = args.delim
delimiter_second = args.delim

delimiter_first = delimiter_first.replace("\\t", "\t")
delimiter_second = delimiter_second.replace("\\t", "\t")
delimiter_output = delimiter_first

join_first = args.join
join_second = args.join

columns_first = list(map(int, args.columns_first.split(','))) if args.columns_first is not None else None
columns_second = list(map(int, args.columns_second.split(','))) if args.columns_second is not None else None

if args.delim_first is not None:
    delimiter_first = args.delim_first
if args.delim_second is not None:
    delimiter_second = args.delim_second
if args.join_first is not None:
    join_first = args.join_first
if args.join_second is not None:
    join_second = args.join_second

# print(first_file)
# print(second_file)
#
# print(delimiter_first)
# print(delimiter_second)
#
# print(join_first)
# print(join_second)
#
# print(columns_first)
# print(columns_second)


def LoadDict(file_path: str, delim: str, key_column: int, col_select, suppress_warnings=True):
    result = dict()
    with open(file_path, 'r') as file:
        for line in file:
            line = line.rstrip()

            tokens = line.split(delim)

            if len(tokens) <= key_column:
                if not suppress_warnings:
                    print("Warning: line (len: {}) is shorter than key column index ({}). \n{}".format(key_column, len(tokens), line))
                continue
            # print(col_select)
            if col_select is not None and len(tokens) <= max(col_select):
                if not suppress_warnings:
                    print("Warning: line (len: {}) is shorter than max index in column select ({}). \n{}".format(max(col_select), len(tokens), line))
                continue

            key = tokens[key_column]

            if not key:
                continue

            if key in result:
                print("key {} not unique".format(key))
                exit(9)

            if col_select is not None:
                tokens = [tokens[i] for i in col_select]
                # print(tokens)
            result[key] = tokens

    return result


first_dict = LoadDict(first_file, delimiter_first, join_first, columns_first)

second_dict = LoadDict(second_file, delimiter_second, join_second, columns_second)

# for key, tokens in first_dict.items():
#     print(key, tokens)
#
# exit(8)

# print(second_dict.keys())
# print(first_dict.keys())
# print(set(first_dict.keys()).intersection(set(second_dict.keys())))

for key, tokens in first_dict.items():
    if key not in second_dict:
        # print("{}".format(
        #     delimiter_output.join(tokens)
        # ))
        continue

    tokens_second = second_dict[key]

    print("{}{}{}".format(
        delimiter_output.join(tokens),
        delimiter_output,
        delimiter_output.join(tokens_second)
    ))


