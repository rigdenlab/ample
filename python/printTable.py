#!/usr/bin/env python

import locale
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
           return str(num)
   
   
   def get_max_width(self, table, index):
       """Get the maximum width of the given column index"""
       return max([len(self.format_num(row[index])) for row in table])
   
   
   def pprint_table(self, out, table):
       """Prints out a table of data, padded for alignment
       @param out: Output stream (file-like object)
       @param table: The table to print. A list of lists.
       Each row must have the same number of columns. """
       col_paddings = []
   
       for i in range(len(table[0])):
           col_paddings.append(self.get_max_width(table, i))
   
       for row in table:
           # left col
           print >> out, row[0].ljust(col_paddings[0] + 1),
           # rest of the cols
           for i in range(1, len(row)):
               col = self.format_num(row[i]).rjust(col_paddings[i] + 2)
               print >> out, col,
           print >> out

if __name__ == "__main__":
    table = [["", "Ensemble No", "Truncation level", "Side chain type", "MR program", "Final Rfree", "Shelxe CC"],
        ["spam", 1, 2, "All", "Phaser", 0.25, 23.51],
        ["eggs", 2, 6, "Reliable", "Molrep", 0.45, 4.56],
        ["lumberjacks", 3, 10, "Poly", "Phaser", 0.53, 12.12]]
    import sys
    out = sys.stdout
    T=Table()
    T.pprint_table(out, table)
