import sys


class BinaryStatistics:
    def __init__(self):
        self.name = ""
        self.tp = 0
        self.tn = 0
        self.fp = 0
        self.fn = 0

    def Validate(self):
        if self.tp == 0:
            raise Exception("Binary Stats are not right: {}\t{}\t{}\t{}".format(
                self.tp, self.tn, self.fp, self.fn))

    def Sensitivity(self):
        return self.tp / (self.tp + self.fn)

    def Precision(self):
        return self.tp / (self.tp + self.fp)

    def Specificity(self):
        return self.tn / (self.tn + self.fp)

    def Accuracy(self):
        return (self.tn + self.tp) / (self.tn + self.tp + self.fn + self.fp)

    def F1(self):
        return (2 * self.tp) / (2 * self.tp + self.fp + self.fn)

