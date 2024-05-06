class BinaryStatistics:
    def __init__(self, name):
        self.name = name
        self.tp = 0
        self.tn = 0
        self.fp = 0
        self.fn = 0

    def Validate(self):
        if self.tp == 0:
            print("WARNING")
            # input()
            # raise Exception("Binary Stats are not right: {}\t{}\t{}\t{}".format(
            #     self.tp, self.tn, self.fp, self.fn))

    def Sensitivity(self):
        if self.tp + self.fn == 0:
            return 0
        return self.tp / (self.tp + self.fn)

    def Precision(self):
        if self.tp + self.fp == 0:
            return 0
        return self.tp / (self.tp + self.fp)

    def Specificity(self):
        if self.tn + self.fp == 0:
            return 0
        return self.tn / (self.tn + self.fp)

    def Accuracy(self):
        if self.tn + self.tp + self.fn + self.fp == 0:
            return 0
        return (self.tn + self.tp) / (self.tn + self.tp + self.fn + self.fp)

    def F1(self):
        if self.tp + self.fp + self.fn == 0:
            return 0
        return (2 * self.tp) / (2 * self.tp + self.fp + self.fn)

