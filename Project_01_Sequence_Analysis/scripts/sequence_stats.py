file = open("data/sample.fasta")
header = file.readline()
sequence = file.readline().strip()
length = len(sequence)
gc = sequence.count("G") + sequence.count("C")
gc_content = (gc / length) * 100

print("DNA Sequence Statistics")
print("------------------------")
print("Sequence length:", length)
print("GC content:", round(gc_content,2), "%")

file.close()
