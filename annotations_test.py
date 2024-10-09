import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
from scipy import stats

# Sample data
group1 = np.random.normal(loc=5, scale=2, size=50)  # Group 1 data
group2 = np.random.normal(loc=7, scale=2, size=50)  # Group 2 data
print(group1.shape, group2.shape)

# Perform a t-test
t_stat, p_val = stats.ttest_ind(group1, group2)

# Create a bar plot
means = [np.mean(group1), np.mean(group2)]
errors = [np.std(group1)/np.sqrt(len(group1)), np.std(group2)/np.sqrt(len(group2))]

plt.figure(figsize=(6, 6))
bars = plt.bar(['Group 1', 'Group 2'], means, yerr=errors, capsize=10, color=['lightblue', 'lightgreen'])

# Adding significance annotation
if p_val < 0.05:
    significance = "*"
else:
    significance = "ns"  # not significant

# Add annotation for statistical difference
x1, x2 = 0, 1   # the x locations for the bars
y, h, col = max(means) + max(errors) + 0.1, 0.2, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, color=col)
plt.text((x1 + x2) * .5, y + h, significance, ha='center', va='bottom', color=col)

plt.title('Statistical Comparison of Two Groups')
plt.ylabel('Mean Values')

# Display plot
plt.tight_layout()
plt.show()

