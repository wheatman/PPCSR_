# filename = 'soc-Slashdot0811.txt'
filename = 'soc-LiveJournal1.txt'

with open(filename, 'r') as f:
  with open('out-' + filename, 'w') as out:
    for line in f:
      if line[0] != '#':
        parts = line.split()
        src = int(parts[0]) + 1
        dest = int(parts[1]) + 1
        out.write(str(src) + '\t' + str(dest) + '\n')
