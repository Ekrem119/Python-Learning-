import matplotlib.pyplot as plt

genes = ["Gene A", "Gene B", "Gene C"]

expression = [120, 80, 200]

plt.bar(genes, expression)
plt.title("Gene Expression")
plt.ylabel("Expression Level")
plt.show()