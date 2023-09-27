from openpyxl import Workbook
import sys
class WorkbookLogSample:
    def __init__(self, output_path):
        self.workbook = Workbook()
        self.active_sheet = self.workbook.active
        self.output_path = output_path
        self.row = 2
        self.write_header()
    
    def write_header(self):
        self.active_sheet["A1"] = "Sample"
        self.active_sheet["B1"] = "Taxon"
        self.active_sheet["C1"] = "Rank"
        self.active_sheet["D1"] = "Dataset"
        self.active_sheet["E1"] = "Tool"
        self.active_sheet["F1"] = "Type"
        self.active_sheet["G1"] = "GOLD"
        self.active_sheet["H1"] = "PRED"
    
    def add_row(self, name, taxon, rank, dataset, tool, binary_pred, abundance_gold, abundance_pred):
        self.active_sheet[f"A{self.row}"] = name
        self.active_sheet[f"B{self.row}"] = taxon
        self.active_sheet[f"C{self.row}"] = rank
        self.active_sheet[f"D{self.row}"] = dataset
        self.active_sheet[f"E{self.row}"] = tool
        self.active_sheet[f"F{self.row}"] = binary_pred
        self.active_sheet[f"G{self.row}"] = abundance_gold
        self.active_sheet[f"H{self.row}"] = abundance_pred
        self.row += 1
        
    def add_sample(self, name, rank, dataset, tool, gold_dict, pred_dict):
        tp_taxa = set(gold_dict).intersection(pred_dict)
        fp_taxa = set(pred_dict).difference(gold_dict)
        fn_taxa = set(gold_dict).difference(pred_dict)
        
        for t in tp_taxa:
            self.add_row(name, t, rank, dataset, tool, "TP", gold_dict[t], pred_dict[t])
        for t in fp_taxa:
            self.add_row(name, t, rank, dataset, tool, "FP", 0, pred_dict[t])
        for t in fn_taxa:
            self.add_row(name, t, rank, dataset, tool, "FN", gold_dict[t], 0)
        
    def save(self):
        self.workbook.save(self.output_path)
        
    def write_xsv(self, sep='\t'):
        xsv = open("data.csv", "w+")
        
        for row in self.active_sheet:
            row_list = list(row)
            for i in range(len(row_list)):
                if i == len(row_list) - 1:
                    xsv.write(str(row_list[i].value))
                else:
                    xsv.write(str(row_list[i].value) + sep)
                xsv.write('\n')
                
        xsv.close()
    
    def print_xsv(self, sep='\t', output=sys.stdout):
        for row in self.active_sheet:
            row_list = list(row)
            for i in range(len(row_list)):
                if i == len(row_list) - 1:
                    output.write(str(row_list[i].value))
                else:
                    output.write(str(row_list[i].value) + sep)

            output.write('\n')
            