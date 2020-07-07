from random import shuffle

#filename = 'out-soc-Pokec.txt'
# filename = 'out-soc-Slashdot0811.txt'
filename = 'soc-pokec-relationships.txt'
filename = 'out-soc-LiveJournal1.txt'
with open(filename, 'r') as f:
  with open('rand-' + filename, 'w') as out:
    lines = f.readlines()
    shuffle(lines)
    out.write("".join(lines))
