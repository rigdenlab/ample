#!/usr/bin/env python

import locale
import sys
locale.setlocale(locale.LC_NUMERIC, "")

class Table:

    def format_num(self, num):
        """Format a number according to given places.
        Adds commas, etc. Will truncate floats into ints!"""

        try:
            if "." in num:
                inum = float(num)
                return locale.format("%.2f", (0, inum), True)
            else:
                inum = int(num)
                return locale.format("%.*f", (0, inum), True)

        except (ValueError, TypeError):
            return num.encode('utf-8') if isinstance(num,unicode) else str(num)


    def get_max_width(self, table, index):
        """Get the maximum width of the given column index"""
        return max([len(self.format_num(row[index])) for row in table])


    def pprint_table(self, table):
        """Prints out a table of data, padded for alignment
        @param table: The table to print. A list of lists.
        Each row must have the same number of columns. """
        col_paddings = []

        out = ""
        for i in range(len(table[0])):
            col_paddings.append(self.get_max_width(table, i))

        for row in table:
            # left col
            out += str(row[0]).ljust(col_paddings[0] + 1)
            # rest of the cols
            for i in range(1, len(row)):
                col = self.format_num(row[i]).rjust(col_paddings[i] + 2)
                out +=  col
            out += "\n"

        return out


if __name__ == "__main__":
    T=Table()
    T.bumppath = '/home/jmht/ample-dev1/examples/toxd-example/ROSETTA_MR_3/MRBUMP/cluster_run1'
    T.cluster = True
    table = T.maketable()
    out = sys.stdout
    T.pprint_table(out, table)
    print('done')

